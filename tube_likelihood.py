 #!/usr/bin/env python3
import argparse
from collections import OrderedDict
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
    parser.add_argument('filename', nargs=1, help='dump file to read')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-o', '--output', help='output name')
    
    args = parser.parse_args(cmdline)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)

    #data = read_quantification(args.filename[0], 'FPKM')
    data = pandas.read_csv(args.filename[0], sep='\t', header=0)

    data = log_likelihood_by_run(data)
    data = data.sort_values(by="run_LR", ascending=True)
    data = chi(data)

    if args.output:
        data.to_csv(args.output, sep='\t', index=False)
    else:
        print(data.to_string())

def log_likelihood_by_run(data, data_runs=None):
    results = []
    if data_runs is None:
        data_runs = data.run.unique()

    likelihoods = compute_log_likelihoods(data)
    for i, run_name in enumerate(data_runs):
        results.append(optimize_by_run(data, likelihoods, run_name))
        
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
    limit = 100
    vr  = []
    vrs = []
    vrp = []
    
    is_run = (data['run'] == run_name)
    is_not_run = (data['run'] != run_name)

    vr = likelihoods.sum(axis=0)
    vr.name = 'vr'
    vrs = likelihoods[is_not_run].sum(axis=0)
    vrs.name = 'vrs'
    vrp = likelihoods[is_run].sum(axis=0)
    vrp.name = 'vrp'

    df = pandas.concat([vr,vrs,vrp], axis=1)

    if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
        filename = run_name + '.vr.df.csv'
        logger.debug('dumping dataframe to %s', filename)
        df.to_csv(filename, index=False)
    
    like_tot = df.vr.max()
    like_non_run = df.vrs.max()
    like_run = df.vrp.max()

    run_LR = -2 * (like_tot - (like_non_run + like_run))

    return pandas.Series(OrderedDict([
        ('run_name', run_name),
        ('run_LR', run_LR),
        ('like_non_run', like_non_run), #vrs
        ('like_run', like_run), #vrp
        ('like_tot', like_tot), #vr
        ('psmc_non_run', df.vrs.idxmax()),
        ('psmc_run', df.vrp.idxmax()),
        ('psmc_tot', df.vr.idxmax()),
    ]))

def optimize_by_pool(data):
    limit = 100
    vr  = []
    single_vrs = []
    pool_vrp = []
    
    is_poolsplit = data['run'].apply(lambda x: x.startswith('p'))
    is_single = data['run'].apply(lambda x: x.startswith('s'))

    prob_range = numpy.arange(0.01, 1.01, .01)
    for p in prob_range:
        result = numpy.log(data.apply(prob, axis=1, args=(p,)))

        vr.append(result.sum())
        pool_vrp.append(result[is_poolsplit].sum())
        single_vrs.append(result[is_single].sum())

    df  = pandas.DataFrame({'p': prob_range,
                            'vr': vr,
                            'single_vrs': single_vrs,
                            'pool_vrp': pool_vrp,
    })

    like_tot_idx = df.vr.idxmax()
    like_single_idx = df.single_vrs.idxmax()
    like_pool_idx = df.pool_vrp.idxmax()
    
    like_tot = df.vr.iloc[like_tot_idx]
    like_single = df.single_vrs.iloc[like_single_idx]
    like_pool = df.pool_vrp.iloc[like_pool_idx]
    
    pool_LR = -2 * (like_tot - (like_pool + like_single))

    return pandas.Series(OrderedDict([
        ('pool_LR', pool_LR),
        ('like_single', like_single),
        ('like_pool', like_pool),
        ('like_tot', like_tot),
        ('p_single', df.p.iloc[like_single_idx]),
        ('p_pool', df.p.iloc[like_pool_idx]),
        ('p_tot', df.p.iloc[like_tot_idx]),
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
