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

logger = logging.getLogger('tube_likelihood')

def main(cmdline=None):
    parser = argparse.ArgumentParser()
    #parser.add_argument('filename', nargs=1, help='dump file to read')
    parser.add_argument('-p', '--pool', action='append',
                        help='pool-split RSEM quantification files')
    parser.add_argument('-s', '--single', action='append',
                    help='single-cell RSEM quantification files')
    parser.add_argument('-c', '--concentrations', required=True,
                        help='name of file with concentrations for spike ins')
    parser.add_argument('-o', '--output', help='output name')
    parser.add_argument('-q', '--quantification', default='FPKM',
                        help='Which RSEM quantification column to use')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-d', '--debug', action='store_true')
    
    args = parser.parse_args(cmdline)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)

    concentrations = read_concentrations(args.concentrations)
    pool = read_quantifications(args.pool, 'pool', args.quantification, concentrations)
    single = read_quantifications(args.single, 'single', args.quantification, concentrations)
    data = pandas.concat([pool, single])

    if len(data) == 0:
        parser.error('some libraries need to be specified')

    data = log_likelihood(data)
    data = data.sort_values(by="run_LR", ascending=True)
    data = chi(data)

    if args.output:
        data.to_csv(args.output, sep='\t', index=False)
    else:
        print(data.to_string())

def read_concentrations(filename):
    c = pandas.read_csv(filename, sep='\t', header=0)
    return c

def read_quantifications(patterns, tube_type, quantification, concentrations):
    if patterns is None or len(patterns) == 0:
        return pandas.DataFrame(columns=['gene_id', quantification])

    data = []
    for pattern in patterns:
        filenames = glob(pattern)
        for filename in filenames:
            rsem = pandas.read_csv(filename, sep='\t', header=0, usecols=['gene_id', quantification])
            spikes = concentrations.merge(rsem, how='inner')
            success = spikes.apply(lambda x: 1 if x[quantification] > 0 else 0, axis=1)
            spikes = pandas.DataFrame.assign(spikes,
                                             run=filename,
                                             tube_type=tube_type,
                                             success=success,
            )

            data.append(spikes)

    return pandas.concat(data)

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
    return data.sort_values(by=['vadj'], kind='mergesort')
    
if __name__ == "__main__":
    main()
