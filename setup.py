from setuptools import setup, find_packages
from singleqc.version import get_git_version

setup(
    name='singleqc',
    version=get_git_version(),
    packages = find_packages(),
    package_data={
        'singleqc': [
            'singleqc/RELEASE-VERSION',
            'singleqc/gspikein.txt'],
    },
    entry_points={
        'console_scripts': [
            'tube_likelihood = singleqc.tube_likelihood:main',
            'gene_spike_ratio = singleqc.gene_spike_ratio:main',
        ],
    },
    install_requires=[
        'pandas>=0.18',
        'statsmodels>=0.6.1',
        'scipy>=0.18',
    ],
    tests_require=[
        'py.test',
        'rpy2',
        'six',
    ]
)
   
