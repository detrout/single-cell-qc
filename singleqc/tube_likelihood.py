 #!/usr/bin/env python3
import argparse
from collections import OrderedDict
from glob import glob
import logging
import numpy
import os
import time
import pandas
import scipy.misc
import scipy.stats
from statsmodels.sandbox.stats.multicomp import multipletests

from singleqc import configure_logging, read_concentrations

logger = logging.getLogger('tube_likelihood')

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)
    configure_logging(args)

    concentrations = read_concentrations(args.concentrations)
    sep = args.sep

    data = []
    pool = read_rsem_quantifications(args.pool, 'pool', args.quantification, concentrations)
    data.append(pool)
    single = read_rsem_quantifications(args.single, 'single', args.quantification, concentrations)
    data.append(single)
    data.extend(read_combined_quantifications(args.combined_pool, 'pool', args.quantification, concentrations, sep))
    data.extend(read_combined_quantifications(args.combined_single, 'single', args.quantification, concentrations, sep))
    data = pandas.concat(data)

    if len(data) == 0:
        parser.error('some libraries need to be specified')

    data = log_likelihood(data)
    data = data.sort_values(by="run_LR", ascending=True)
    data = chi(data)

    if args.output:
        data.to_csv(args.output, sep='\t', index=False)
    else:
        print(data.to_string(index=False))


def make_parser():
    parser = argparse.ArgumentParser()

    group = parser.add_argument_group('combined quantification file')
    group.add_argument('--combined-pool', action='append', default=[],
                        help='file with merged pool-split quantifications to read')
    group.add_argument('--combined-single', action='append', default=[],
                        help='file with merge single cell quantifications to read')

    group = parser.add_argument_group('raw RSEM files')
    group.add_argument('-p', '--pool', action='append', default=[],
                        help='pool-split RSEM quantification files')
    group.add_argument('-s', '--single', action='append', default=[],
                    help='single-cell RSEM quantification files')
    group.add_argument('-q', '--quantification', default='FPKM',
                        help='Which RSEM quantification column to use')

    parser.add_argument('-c', '--concentrations', required=True,
                        help='name of file with concentrations for spike ins')
    parser.add_argument('-o', '--output', help='output name')
    parser.add_argument('--sep', help='quantification file seperator', choices=['\t', ','], default='\t')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-d', '--debug', action='store_true')

    return parser


def read_combined_quantifications(filenames, tube_type, quantification_name, concentrations, sep):
    """Read several combined quantification files.

    this is a gene_id vs library_id tables, if there is a column named "gene_name" it
    will be ignored.
    """
    data = []
    for filename in filenames:
        data.append(read_combined_quantification(filename, tube_type, quantification_name, concentrations, sep))

    return data


def read_combined_quantification(filename, tube_type, quantification_name, concentrations, sep):
    """Read a combined quantification files gene_id vs library_id

    this is a gene_id vs library_id tables, if there is a column named "gene_name" it
    will be ignored.
    """
    quantifications = pandas.read_csv(filename, sep=sep, header=0)
    quantifications = quantifications.set_index('gene_id')

    data = []
    for column in quantifications.columns:
        if column == 'gene_name':
            logger.info('Ignoring gene_name column')
        else:
            spikes = make_spike_success_table(
                quantifications[column].to_frame(quantification_name), concentrations, quantification_name, column, tube_type)
            data.append(spikes)

    return pandas.concat(data)


def read_rsem_quantifications(patterns, tube_type, quantification_name, concentrations):
    """Read a specific quantification type column out of RSEM quantification files.
    """
    if patterns is None or len(patterns) == 0:
        df = pandas.DataFrame(columns=[quantification_name])
        df.index.name = 'gene_id'
        return df

    data = []
    for pattern in patterns:
        filenames = glob(pattern)
        for filename in filenames:
            rsem = pandas.read_csv(filename, sep='\t', header=0, usecols=['gene_id', quantification_name])
            rsem = rsem.set_index('gene_id')
            rsem.columns = [quantification_name]
            spikes = make_spike_success_table(rsem, concentrations, quantification_name, filename, tube_type)
            data.append(spikes)

    return pandas.concat(data)


def make_spike_success_table(library_data, concentrations, quantification_name, run_name, tube_type):
    spikes = concentrations.merge(library_data, how='inner', left_index=True, right_index=True)
    success = spikes[quantification_name] > 0
    spikes = pandas.DataFrame.assign(spikes,
                                     run=run_name,
                                     tube_type=tube_type,
                                     success=success,
    )
    return spikes

def log_likelihood(data, data_runs=None):
    results = []
    if data_runs is None:
        data_runs = data.run.unique()

    likelihoods = compute_log_likelihoods(data)
    for i, run_name in enumerate(data_runs):
        results.append(optimize_by_run(data, likelihoods, run_name))

    if len(data.tube_type.unique()) > 1:
        results.append(optimize_by_tube_type(data, likelihoods))
        
    return pandas.DataFrame(results)

def prob(row, p, K, K_factorial, Threshold):
    pr_k = row.concentration ** K * 2.71 ** -row.concentration/K_factorial
    one_minus_p_k = (1 - p) ** K

    if row.success:
        pr_c_g_k = 1 - one_minus_p_k
    else:
        pr_c_g_k = one_minus_p_k

    return (pr_c_g_k * pr_k).sum().clip(min=Threshold)

def compute_log_likelihoods(data):
    """Return log likelihoods for [0.0, 1.0, step=.01]
    """
    K = numpy.arange(0, 61)
    K_factorial = scipy.misc.factorial(K)
    prob_range = numpy.arange(0.01, 1.01, .01)
    Threshold = .000000001

    vr = []
    for p in prob_range:
        vr.append(data.apply(prob, axis=1, args=(p, K, K_factorial, Threshold)))

    vrmatrix = pandas.DataFrame(dict(zip(prob_range, vr)))
    return numpy.log(vrmatrix)

def optimize_by_run(data, likelihoods, run_name):
    member = (data['run'] == run_name)
    notmember = (data['run'] != run_name)
    return optimize_by(data, likelihoods, run_name, member, notmember)

def optimize_by_tube_type(data, likelihoods):
    member = (data['tube_type'] == 'pool')
    notmember = (data['tube_type'] == 'single')
    return optimize_by(data, likelihoods, 'pool_v_single', member, notmember)

def optimize_by(data, likelihoods, name, member, notmember):
    limit = 100
    vr  = []
    vrs = []
    vrp = []
    
    vr = likelihoods.sum(axis=0)
    vr.name = 'vr'
    vrs = likelihoods[notmember].sum(axis=0)
    vrs.name = 'vrs'
    vrp = likelihoods[member].sum(axis=0)
    vrp.name = 'vrp'

    df = pandas.concat([vr,vrs,vrp], axis=1)

    like_tot = df.vr.max()
    like_non_run = df.vrs.max()
    like_run = df.vrp.max()

    run_LR = -2 * (like_tot - (like_non_run + like_run))

    return pandas.Series(OrderedDict([
        ('run_name', name),
        ('run_LR', run_LR),
        ('like_non_run', like_non_run), #vrs
        ('like_run', like_run), #vrp
        ('like_tot', like_tot), #vr
        ('psmc_non_run', df.vrs.idxmax()),
        ('psmc_run', df.vrp.idxmax()),
        ('psmc_tot', df.vr.idxmax()),
    ]))

def chi(data):
    vchi = 1 - scipy.stats.chi2.cdf(data.run_LR, 1)
    reject, vadj, sidak, bonf = multipletests(vchi, method='bonferroni')
    data['vchi'] = vchi
    data['vadj'] = vadj
    # 
    return data.sort_values(by=['vadj','vchi', ], kind='mergesort')
    
if __name__ == "__main__":
    main()
