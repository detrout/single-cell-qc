 #!/usr/bin/env python3
"""Prototype port of original processing pipeline to a single process
"""
import argparse
import logging
import numpy
import os
import time
import pandas
import scipy.misc
import scipy.stats
from math import factorial, log
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from statsmodels.sandbox.stats.multicomp import multipletests

logger = logging.getLogger('tube_likelihood_r')

r("""

prob <- function (p, success, concentration)
{
    r1 = 0

    for (k in 0:60)
    {
        if (1 == success) {pr_c_g_k = 1 - (1 - p)^k}
        else              {pr_c_g_k = (1 - p)^k}

        pr_k = concentration ^ k * 2.71 ^ -concentration/factorial (k)
        r1   = r1 + pr_c_g_k * pr_k
    }

    r1 = max (r1, .000000001)
}

optimize <- function(data, run_name)
{
    vr  = 0
    vrs = 0
    vrp = 0
    n   = length (data$run)

    #table (data$success)

    #summary (lm (data$success ~ data$concentration))

    for (i in 1:100)
    {
        p  = i/100
        r  = 0
        rs = 0
        rp = 0

        for (j in 1:n)
        {
            result = log (prob (p, data$success[j], data$concentration [j]))
            r      = r + result;
            if (data$run[j] == run_name)
               rp = rp + result
            else
               rs = rs + result
        }

        vr[i]  = r
        vrs[i] = rs
        vrp[i] = rp
    }

    inc = 1:100
    p   = inc/100
    df  = data.frame (p, vr, vrs, vrp)

    dfp <- df [order (-vrp),]
    dfs <- df [order (-vrs),]
    dfb <- df [order (-vr ),]

#    dfs$vrs[1]
#    dfp$vrp[1]
#    dfb$vr[1]
#
#    -2 * (dfb$vr[1] - (dfp$vrp[1] + dfs$vrs[1]))
#
#    (dfp$p[1] + dfs$p[1])/2
#    dfs$p[1]
#    dfp$p[1]
#    dfb$p[1]
#
    rst = c (
           run_name,
           dfs$vrs[1],
           dfp$vrp[1],
           dfb$vr[1] ,
           dfs$p[1]  ,
           dfp$p[1]  ,
           dfb$p[1])

#    str (rst)
#    class (rst)

    return(rst)
}

chi <- function (data) 
{
    print(data)
    attach(data)

    vchi = 1 - pchisq (LR,1)
    vadj = p.adjust (vchi, method = "bonferroni")
    df   = data.frame(run_name, vchi, vadj, LR, p_none_run, p_run, p_tot)

    dfp  = df [order (vadj),]

    return(dfp)
#    results = 0
#    for (i in 1:length (vadj))
#       {
#       rst <-  c(as.character (dfp$cell_type[i]),
#                 as.character (dfp$run[i]      ),
#                 as.character (dfp$vchi[i]     ),
#                 as.character (dfp$vadj[i]     ),
#                 as.character (dfp$LR[i]       ))
#       results[i] = rst
#
#       write (rst, file = "chi.tbl",
#                       ncolumns = 5, append = TRUE, sep = "\t")
#       }
#    return(results)
}
""")

def main(cmdline=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', nargs=1, help='dump file to read')
    parser.add_argument('--cell-type')
    parser.add_argument('--run-name')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-d', '--debug', action='store_true')
    
    args = parser.parse_args(cmdline)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)

    data = pandas.read_csv(args.filename[0], sep='\t', header=0)
    data = log_likelihood(data)

    data = data.sort_values(by="LR", ascending=True)
    chi = chi_r(data)
    chi.to_csv('chi.csv')

def log_likelihood(data, data_runs=None):
    results = []
    if data_runs is None:
        data_runs = data.run.unique()
    for i, run_name in enumerate(data_runs):
        logger.info(
            'Processing run {}, {} of {}'.format(run_name, i+1, len(data_runs)))
        start = time.monotonic()
        results.append(optimize_r(data, run_name))
        finish = time.monotonic()
        logger.info("run {} took {} seconds".format(run_name, finish-start))
        start = finish
        
    results = pandas.DataFrame(results)
    #results.to_csv('log_likelihood.csv')
    return results
    
def optimize_r(data, run_name):
    columns = ['run_name', 'LR', 'like_non_run', 'like_run', 'like_tot', 'p_none_run',
               'p_run', 'p_tot']

    results = r['optimize'](data, run_name)
    print(results)
    assert run_name == str(results[0])

    like_non_run = float(results[1])
    like_run = float(results[2])
    like_tot = float(results[3])
    p_none_run = float(results[4])
    p_run = float(results[5])
    p_tot = float(results[6])
    LR = -2 * (like_tot - (like_non_run + like_run))
    return pandas.Series({
        'run_name': run_name,
        'LR': LR,
        'like_non_run': like_non_run,
        'like_run': like_run,
        'like_tot': like_tot,
        'p_none_run': p_none_run,
        'p_run': p_run,
        'p_tot': p_tot,
    }, index=columns)

def score_likelihood(data):
    for row in data.iterrows():
        cell_type = row.cell_type
        
    print(optimize(data, cell_type, run_name))

def chi_r(data):
    data_r = pandas2ri.py2ri(data)
    results = pandas2ri.ri2py(r['chi'](data_r))
    chi = results[['run_name', 'vchi', 'vadj', 'LR']]
    return chi
    
if __name__ == "__main__":
    main()
