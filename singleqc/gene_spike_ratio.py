#+/usr/bin/env python3
from __future__ import print_function

import argparse
import logging
import pandas
from scipy.stats import ttest_ind

from singleqc import configure_logging

logger = logging.getLogger('gene_spike_ratio')

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)
    configure_logging(args)

    if len(args.rsem) != len(args.rsem_library):
        parser.error("every rsem filename must be paired with a library ID")

    data = []
    data.extend(compute_ratio_from_rsem(args.rsem, args.rsem_library, args.pool, args.single, args.quantification_name))
    data.extend(compute_ratio_from_combined(args.combined_pool, 'pool', args.sep))
    data.extend(compute_ratio_from_combined(args.combined_single, 'single', args.sep))
    data = pandas.DataFrame(data)
    print(data)

    pool = data[data['tube_type'] == 'pool']['ratio']
    if len(pool) > 0:
        print('Pool-split: mean {:.3} stdev {:.3}'.format(pool.mean(), pool.std()))

    single = data[data['tube_type'] == 'single']['ratio']
    if len(single) > 0:
        print('Single: mean {:.3} stdev {:.3}'.format(single.mean(), single.std()))

    if len(pool) > 0 and len(single) > 0:
        t = ttest_ind(pool, single, equal_var=False)
        print(t)


def make_parser():
    parser = argparse.ArgumentParser()

    group = parser.add_argument_group('combined quantification file')
    group.add_argument('--combined-pool', action='append', default=[],
                        help='file with merged pool-split quantifications to read')
    group.add_argument('--combined-single', action='append', default=[],
                        help='file with merge single cell quantifications to read')

    group = parser.add_argument_group('raw RSEM files')
    group.add_argument('--rsem', action='append', default=[],
                       help='name of rsem quantification files ')
    group.add_argument('--rsem-library', action='append', default=[],
                       help='library id for rsem quantification file')
    group.add_argument('--pool', action='append', default=[],
                       help='pool-split library names')
    group.add_argument('--single', action='append', default=[],
                       help='single-cell library names (default)')

    parser.add_argument('-c', '--concentrations', required=True,
                        help='name of file with concentrations for spike ins')
    parser.add_argument('-q', '--quantification-name', default='FPKM',
                        help='Which RSEM quantification column to use')
    parser.add_argument('-o', '--output', help='output name')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-s', '--sep', choices=[',','\t'], default='\t')

    return parser


def compute_gene_spike_ratios(library_id, quantification, tube_type):
    '''Compare the sum of genes and sum of spiikes.
    '''
    spikes = ['Spikein' in x for x in quantification.index ]
    notspikes = ['Spikein' not in x for x in quantification.index ]

    spike_sum = quantification[spikes].sum(axis=0)
    gene_sum = quantification[notspikes].sum(axis=0)
    ratio = (spike_sum/(spike_sum+gene_sum))*100

    s = pandas.Series([gene_sum, spike_sum, ratio, tube_type],
                      index=['gene_sum', 'spike_sum', 'ratio', 'tube_type'])
    s.name = library_id
    return s


def compute_ratio_from_combined(filenames, tube_type, sep):
    data = []
    for filename in filenames:
        combined = pandas.read_csv(filename, sep=sep, index_col=0)
        for column in combined.columns:
            if column == 'gene_name':
                logger.warning('Ignoring column named gene_name')
            else:
                data.append(compute_gene_spike_ratios(column, combined[column], tube_type))
    return data


def compute_ratio_from_rsem(rsem, library_ids, pool_ids, single_ids, quantification_name):
    data = []
    for filename, library_id in zip(rsem, library_ids):
        current = pandas.read_csv(filename, sep='\t', index_col=0)
        current = current[quantification_name]
        if library_id in pool_ids:
            tube_type = 'pool'
        elif library_id in single_ids:
            tube_type = 'single'
        else:
            tube_type = 'NOT DEFINED'

        ratios = compute_gene_spike_ratios(library_id, current, tube_type)
        data.append(ratios)

    return data

    
if __name__ == '__main__':
    main()
