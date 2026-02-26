#!/usr/bin/env python3

import os
import sys
from setuptools import setup, find_packages

# Ensure local package imports work when setup.py is executed in an isolated build env.
sys.path.insert(0, os.path.dirname(__file__))

from gottcha import __version__ as version

# Read README file
readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
if os.path.exists(readme_path):
    with open(readme_path, 'r', encoding='utf-8') as fh:
        long_description = fh.read()
else:
    long_description = "GOTTCHA2: Genomic Origin Through Taxonomic CHAllenge version 2"

setup(
    name="GOTTCHA2",
    version=version,
    author="Po-E Li",
    author_email="po-e@lanl.gov",
    description="GOTTCHA2: Genomic Origin Through Taxonomic CHAllenge v2",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="BSD-3-Clause",
    keywords=['bioinformatics', 'taxonomy', 'profiler', 'metagenomics', 'microbiome'],
    url="https://github.com/poeli/GOTTCHA2",
    project_urls={
        "Bug Reports": "https://github.com/poeli/GOTTCHA2/issues",
        "Source": "https://github.com/poeli/GOTTCHA2",
        "Documentation": "https://github.com/poeli/GOTTCHA2/blob/master/README.md",
    },
    packages=find_packages(exclude=['tests*', 'docs*']),
    package_data={
        'gottcha': ['data/*', 'templates/*'],
    },
    include_package_data=True,
    python_requires=">=3.8",
    install_requires=[
        'numpy>=1.19.0',
        'pandas>=1.2.0',
        'requests>=2.25.0',
        'biom-format>=2.1.7',
        'pysam>=0.22.0',
        'tqdm',
    ],
    entry_points={
        'console_scripts': [
            'gottcha2=gottcha.utils.cmd:gottcha2_command',
        ],
    },
    scripts=[
        'gottcha/gottcha2.py',
        'gottcha/utils/pull_database.py',
        'gottcha/utils/taxonomy.py',
        'gottcha/utils/split_reads.py',
        'gottcha/utils/gottcha_sam_to_bam.py',
        'gottcha/utils/gottcha2_bam.py',
        'gottcha/utils/extract_taxa_reads.py',
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bioinformatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    zip_safe=False,
)
