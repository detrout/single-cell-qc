import pandas
import numpy
import pkg_resources


def get_rsem_filename():
    '''Return path to an example rsem file'''
    return pkg_resources.resource_filename(__name__, 'rsem.genes.results')


def get_concentration_filename():
    '''Return path to the example concentration file'''
    return pkg_resources.resource_filename('singleqc', 'ENCSR535LMC.tsv')


def get_mm_dump_filename():
    return pkg_resources.resource_filename(__name__, 'dump_mm_purkinje.txt')


def make_combined(stream, scale=0.1):
    '''Fill a StringIO with a simulated combined quantification file.

    it shifts the quantifications around using a normal distribution
    centered on the example RSEM file, with a width of scale.

    use with
    stream = StringIO()
    make_combined(stream, scale=0.1)
    stream.seek(0)
    call read_csv with the stream.
    '''
    rsem_filename = get_rsem_filename()
    rsem = pandas.read_csv(rsem_filename, sep='\t', index_col=0)
    fpkm = rsem['FPKM']
    data = {}
    for i in range(10):
        data[i] = numpy.random.normal(loc=fpkm, scale=scale*fpkm, size=fpkm.shape)

    df = pandas.DataFrame(data, columns=range(10), index=fpkm.index)
    df.to_csv(stream, sep='\t')
