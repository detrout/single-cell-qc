import pandas
import unittest
import numpy
import scipy.special

from six import StringIO

from rpy2.robjects import r
from rpy2.robjects import pandas2ri
pandas2ri.activate()

from singleqc import tube_likelihood, read_concentrations

from tests import (
    get_concentration_filename,
    get_mm_dump_filename,
    get_rsem_filename,
    make_combined
)

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

    rst = c (
           run_name,
           dfs$vrs[1],
           dfp$vrp[1],
           dfb$vr[1] ,
           dfs$p[1]  ,
           dfp$p[1]  ,
           dfb$p[1])

    return(rst)
}

optimize_pool <- function(data)
{
    vr  = 0
    vrs = 0
    vrp = 0
    n   = length (data$run)

    tube_type = grepl("s", data$run)

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
            if (FALSE == tube_type[j])
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

    LR <- -2 * (dfb$vr[1] - (dfp$vrp[1] + dfs$vrs[1]))

    rst = c (
           LR,
           dfs$vrs[1],
           dfp$vrp[1],
           dfb$vr[1] ,
           dfs$p[1]  ,
           dfp$p[1]  ,
           dfb$p[1])

    return(rst)
}

chi <- function (data) 
{
    attach(data)

    vchi = 1 - pchisq (run_LR,1)
    vadj = p.adjust (vchi, method = "bonferroni")
    df   = data.frame(vchi, vadj, run_LR)

    dfp  = df [order (vadj,vchi),]

    return(dfp)
}
""")


class TubeLikelihood(unittest.TestCase):
    def test_read_rsem_quantifications(self):
        concentrations = read_concentrations(get_concentration_filename())
        rsem_file = get_rsem_filename()

        single = tube_likelihood.read_rsem_quantifications(
            [rsem_file], 'single', 'FPKM', concentrations)
        self.assertEqual(len(single.columns), 5)
        self.assertIn('FPKM', single.columns)
        self.assertEqual(set(single.run), set((rsem_file,)))

    def test_read_combined_quantifications(self):
        concentrations = read_concentrations(get_concentration_filename())
        pool_stream = StringIO()
        make_combined(pool_stream, scale=0.05)
        pool_stream.seek(0)

        pool = tube_likelihood.read_combined_quantification(pool_stream, 'pool', 'FPKM', concentrations, sep='\t')
        # columns should be "quantification", 'concentration', 'run', 'success', 'tube_type'
        self.assertEqual(len(pool.columns), 5)
        self.assertIn('FPKM', pool.columns)
        self.assertEqual(set(pool['run']), set([str(x) for x in range(10)]))

    def test_prob(self):
        concentrations = [
            0.11, 7.055, 28.219, 56.438, 112.875, 1806
        ]
        succeeded = [0, 1]
        K = numpy.arange(0, 61)
        K_factorial = scipy.special.factorial(K)
        Threshold = .000000001

        for p in [0.05, 0.25, 0.5, 0.75, 1.0]:
            for succ in succeeded:
                for c in concentrations:
                    r_answer = r['prob'](p, succ, c)[0]
                    row = pandas.Series({
                        'success': succ,
                        'concentration': c})
                    py_answer = tube_likelihood.prob(row, p, K, K_factorial, Threshold)
                    self.assertAlmostEqual(r_answer, py_answer)

    def test_optimize_pool(self):
        data = pandas.read_csv(get_mm_dump_filename(), sep='\t')
        r_result = r['optimize_pool'](data)
        likelihoods = tube_likelihood.compute_log_likelihoods(data)
        py_result = tube_likelihood.optimize_by_tube_type(data, likelihoods)
        py_result = py_result[['run_LR',
                               'like_non_run', 'like_run', 'like_tot',
                               'psmc_non_run', 'psmc_run', 'psmc_tot']]
        for x, y in zip(r_result, py_result):
            self.assertAlmostEqual(x, y)


    def test_optimize_run(self):
        data = pandas.read_csv(get_mm_dump_filename(), sep='\t', header=0)
        run_name = 'p15290'
        r_result = r['optimize'](data, run_name)
        # chump off run_name from r result
        r_result = [ float(x) for x in r_result[1:] ]
        likelihoods = tube_likelihood.compute_log_likelihoods(data)
        py_result = tube_likelihood.optimize_by_run(data, likelihoods, run_name)
        py_result = py_result[['like_non_run', 'like_run', 'like_tot',
                               'psmc_non_run', 'psmc_run', 'psmc_tot']]
        # the orders of there values needs to match (see above)
        for x, y in zip(r_result, py_result):
            self.assertAlmostEqual(
                x, y,
                msg="failed x={}[{}] y={}[{}]".format(x, type(x), y, type(y))
            )

    def test_chi(self):
        LR = [3.637978807091713e-12,
              0.0007126937198336236,
              0.025544047941366443,
              0.11297192756683216,
              0.1361518794656149,
              0.1361518794656149,
              0.15117851703325869,
              0.1625158952883794,
              0.17410559771815315,
              0.18334502939978847,
              0.35617773549529375,
              0.4202275662973989,
              0.5116423363615468,
              0.5334788473010121,
              0.7624761984661745,
              0.9247238012248999,
              1.0858579376545094,
              1.1433223595813615,
              1.1545867489003285,
              1.340766585446545,
              1.4091111523812287,
              1.4244564379769145,
              1.4244564379769145,
              1.9246341756479524,
              2.251595334138983,
              2.7844221524319437,
              2.960457152054005,
              2.9678872130680247,
              3.5743446097731066,
              4.316940278207767,
              5.57561769991662,
              9.628733330253453,]
        data = pandas.DataFrame({'run_LR': LR})
        r_result = pandas2ri.ri2py(r['chi'](data))
        py_result = tube_likelihood.chi(data)

        for ((r_i, r_row), (py_i, py_row)) in zip(r_result.iterrows(), py_result.iterrows()):
            for colname in ['vadj', 'vchi', 'run_LR']:
                self.assertAlmostEqual(
                    r_row[colname],
                    py_row[colname],
                    msg='failed row {} col {}'.format(r_i, colname))

if __name__ == '__main__':
    unittest.main()
