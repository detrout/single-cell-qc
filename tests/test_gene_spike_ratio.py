import pandas
import numpy
import unittest
import tempfile

from singleqc import gene_spike_ratio
from six import StringIO

from tests import get_concentration_filename, get_rsem_filename, make_combined

class TestGeneSpikeRatio(unittest.TestCase):
    def setUp(self):
        self.rsem_filename = get_rsem_filename()
        
    def test_read_rsem(self):
        rsem = gene_spike_ratio.compute_ratio_from_rsem(
            rsem=[self.rsem_filename, self.rsem_filename],
            library_ids=['1', '2'],
            pool_ids=['1'],
            single_ids=['2'],
            quantification_name='FPKM'
        )
        self.assertEqual(len(rsem), 2)
        for index in ['gene_sum', 'spike_sum', 'ratio']:
            self.assertEqual(rsem[0][index], rsem[1][index])
        
        self.assertEqual(rsem[0].name, '1')
        self.assertEqual(rsem[0].tube_type, 'pool')
        
        self.assertEqual(rsem[1].name, '2')
        self.assertEqual(rsem[1].tube_type, 'single')

    def test_read_combined(self):
        pool_file = StringIO()
        make_combined(pool_file, scale=0.05)
        pool_file.seek(0)
        pool = gene_spike_ratio.compute_ratio_from_combined([pool_file], 'pool', '\t')
        pool = pandas.DataFrame(pool)
        # this is the mean of ratio from the provided test file
        MEAN = 0.473070617734
        self.assertGreater( pool['ratio'].mean(), MEAN - pool['ratio'].std())
        self.assertLess( pool['ratio'].mean(), MEAN + pool['ratio'].std())

    def test_main(self):
        spikein = get_concentration_filename()
        cmdline = ['-c', spikein]
        for tube_type in ['--pool', '--single']:
            for name in range(10):
                cmdline.extend(['--rsem', self.rsem_filename,
                                '--rsem-library', str(name),
                                tube_type, str(name)])
        ret = gene_spike_ratio.main(cmdline)
        self.assertEqual(ret, 0)
        
        
if __name__ == '__main__':
    unittest.main()

