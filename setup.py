#!/usr/bin/env python3
"""
Set up the package
"""

# Standard inputs
from distutils.util import convert_path
import os

# Third party inputs
from setuptools import setup, find_packages

# Find the version
version = {}
with open(convert_path(os.path.join('allele_tools', 'version.py')), 'r') as version_file:
    exec(version_file.read(), version)

setup(
    name="AlleleFinder",
    version=version['__version__'],
    scripts=[
        # os.path.join('allele_tools', 'allele_finder.py'),
        os.path.join('allele_tools', 'allele_profiler.py'),
        os.path.join('allele_tools', 'allele_updater.py'),
        os.path.join('allele_tools', 'allele_translate_reduce.py'),
        # os.path.join('allele_tools', 'probe_creator.py'),
        os.path.join('allele_tools', 'profile_reduce.py'),
        # os.path.join('allele_tools', 'stec_attributer.py'),
        os.path.join('allele_tools', 'stec.py')
    ],
    packages=find_packages(),
    include_package_data=True,
    author="Adam Koziol",
    author_email="adam.koziol@inspection.gc.ca",
    url="https://github.com/OLC-Bioinformatics/AlleleFinder",
)
