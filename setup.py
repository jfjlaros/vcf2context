import sys
from setuptools import setup

if sys.version_info < (2, 6):
    raise Exception('VCF2context requires Python 2.6 or higher.')

# Todo: How does this play with pip freeze requirement files?
requires = ['biopython', 'pyvcf']

# Python 2.6 does not include the argparse module.
try:
    import argparse
except ImportError:
    requires.append('argparse')

import vcf2context as distmeta

setup(
    name='vcf2context',
    version=distmeta.__version__,
    description='VCF to context converter.',
    long_description=distmeta.__doc__,
    author=distmeta.__author__,
    author_email=distmeta.__contact__,
    url=distmeta.__homepage__,
    license='MIT License',
    platforms=['any'],
    packages=['vcf2context'],
    install_requires=requires,
    entry_points = {
        'console_scripts': [
            'vcf2context = vcf2context.vcf2context:main'
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
    ],
    keywords='bioinformatics'
)
