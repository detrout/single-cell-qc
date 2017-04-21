#+/usr/bin/env python3
from __future__ import print_function

import argparse
import logging
import pandas
from scipy.stats import ttest_ind

from tube_likelihood import read_concentrations

logger = logging.getLogger('gene_spike_ratio')

def main(cmdline=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='*', help='merged quantification file to read')
    parser.add_argument('--rsem', action='append', default=[],
                        help='name of rsem quantification files ')
    parser.add_argument('--rsem-library', action='append', default=[],
                        help='library id for rsem quantification file')
    parser.add_argument('--pool', action='append', default=[],
                        help='pool-split library names')
    parser.add_argument('--single', action='append', default=[],
                        help='single-cell library names (default)')
    parser.add_argument('-c', '--concentrations', required=True,
                        help='name of file with concentrations for spike ins')
    parser.add_argument('-q', '--rsem-quantification', default='FPKM',
                        help='Which RSEM quantification column to use')
    parser.add_argument('-o', '--output', help='output name')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-s', '--sep', choices=[',','\t'], default='\t')
    
    args = parser.parse_args(cmdline)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)

    if len(args.rsem) != len(args.rsem_library):
        parser.error("every rsem filename must be paired with a library ID")

    logger.debug('Filenames: %s', ','.join(args.filenames))
    logger.debug('RSEM: %s', ','.join(args.rsem))
    logger.debug('Poolids: %s', ','.join(args.pool))

    concentrations = read_concentrations(args.concentrations)
    data = read_quantifications(args.filenames, args.rsem, args.rsem_library, args.rsem_quantification, args.sep)

    spikes = ['Spikein' in x for x in data.index ]
    notspikes = ['Spikein' not in x for x in data.index ]

    spike_sum = data[spikes].sum(axis=0)
    gene_sum = data[notspikes].sum(axis=0)
    ratio = (spike_sum/(spike_sum+gene_sum))*100

    if len(args.pool) > 0:
        pool = ratio[args.pool]
        print('pool mean', pool.mean())

    if len(args.single) > 0:
        single = ratio[args.single]
        print('single mean', single.mean())

    
    if len(args.pool) > 0 and len(args.single) > 0:
        t = ttest_ind(pool, single, equal_var=False)
        
        print(t)

def read_quantifications(filenames, rsem, rsem_ids, quantification, sep):
    data = []
    for filename in filenames:
        data.append(pandas.read_csv(filename, sep=sep, index_col=0))

    for filename, library_id in zip(rsem, rsem_ids):
        current = pandas.read_csv(filename, sep='\t',
                                  index_col=0)
        current = current[[quantification]]
        current.columns = [library_id]
        current.index.name = 'gene_id'
        data.append(current)

    return pandas.concat(data, axis=1)
    
if __name__ == '__main__':
    main()
