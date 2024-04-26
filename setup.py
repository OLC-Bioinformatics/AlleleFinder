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
with open(convert_path(
        os.path.join('allele_tools', 'version.py')), 'r') as version_file:
    exec(version_file.read(), version)

# Read the contents of your README file
with open('README.md', 'r') as f:
    long_description = f.read()
    
# Open the requirements.txt file
with open('requirements.txt') as f:
    # Read the file and split it into lines
    # Each line should be a separate requirement
    requirements = f.read().splitlines()

setup(
    name="AlleleFinder",
    version=version['__version__'],
    scripts=[
        os.path.join('allele_tools', 'allele_profiler.py'),
        os.path.join('allele_tools', 'allele_updater.py'),
        os.path.join('allele_tools', 'allele_translate_reduce.py'),
        os.path.join('allele_tools', 'profile_reduce.py'),
        os.path.join('allele_tools', 'stec.py')
    ],
    packages=find_packages(),
    include_package_data=True,
    author="Adam Koziol",
    author_email="adam.koziol@inspection.gc.ca",
    url="https://github.com/OLC-Bioinformatics/AlleleFinder",
    # The requirements from requirements.txt
    install_requires=requirements,
    # Classifiers help users find your project by categorizing it.
    classifiers=[
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',
        # Specify the Python versions you support here
        'Programming Language :: Python :: 3.9',
        # Show your project's status - planning, alpha, beta, or stable
        'Development Status :: 5 - Production/Stable',
        # Indicate if your project relates to a particular topic
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
