#!/usr/bin/env python3
import os
from distutils.core import setup
from setuptools import find_packages
from genomic_address_service.version import __version__
author = 'James Robertson'

classifiers = """
Development Status :: 4 - Beta
Environment :: Console
License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python
Programming Language :: Python :: 3.10
Programming Language :: Python :: 3.11
Programming Language :: Python :: 3.12
Programming Language :: Python :: Implementation :: CPython
Operating System :: POSIX :: Linux
""".strip().split('\n')


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


exec(open('genomic_address_service/version.py').read())

setup(
    name='genomic_address_service',
    include_package_data=True,
    version=__version__,
    python_requires='>=3.10.0,<3.13.0',
    setup_requires=['pytest-runner'],
    tests_require=['pytest', 'pytest-workflow'],
    packages=find_packages(exclude=['tests']),
    url='https://github.com/phac-nml/genomic_address_service',
    license='GPLv3',
    author='James Robertson',
    author_email='james.robertson@phac-aspc.gc.ca',
    description=(
        'Genomic Address Service: De novo clustering and cluster address assignment'),
    keywords='cgMLST, wgMLST, outbreak, surveillance, clustering, nomenclature',
    classifiers=classifiers,
    package_dir={'gas': 'genomic_address_service'},
    package_data={
        "": ["*.txt"],
    },

    install_requires=[
        'pyarrow>=14.0.0',
        'numba==0.59.1',
        'numpy==1.26.4',
        'tables==3.9.1',
        'six>=1.16.0',
        'pandas==2.0.2 ',
        'pytest==8.3.3',
        'scipy==1.14.1',
        'psutil==6.1.0',
        'fastparquet==2023.4.0' #Will drop support of fastparquet in future versions

    ],

    entry_points={
        'console_scripts': [
            'gas=genomic_address_service.main:main',
        ],
    },
)
