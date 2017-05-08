from __future__ import print_function

import os
from setuptools import setup, find_packages
from singleqc.version import get_git_version

REQUIREMENTS = "requirements.txt"
if os.path.exists(REQUIREMENTS):
    required = open(REQUIREMENTS).readlines()
else:
    required = []
    print("WARNING: Can't find requirements file")

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
    install_requires=required,
    tests_require=[
        'py.test',
        'rpy2',
        'six',
    ]
)
   
