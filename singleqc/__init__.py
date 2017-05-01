'''Common code shared between single-cell qc modules
'''
import pandas
import logging

def configure_logging(args):
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)


def read_concentrations(filename, sep='\t'):
    '''Read table of spike names and expected copies per cell

    Do allow fractional copies per cell.
    '''
    c = pandas.read_csv(filename, sep=sep, header=0, index_col=0)
    return c
