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
version_file_path = convert_path(
    os.path.join(
        'allele_tools', 'version.py'
    )
)
with open(version_file_path, 'r', encoding='utf-8') as version_file:
    for line in version_file:
        if line.startswith('__version__'):
            DELIMITER = '"' if '"' in line else "'"
            version['__version__'] = line.split(DELIMITER)[1]
            break

# Read the contents of your README file
with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

# Open the requirements.txt file
with open('requirements.txt', 'r', encoding='utf-8') as f:
    # Read the file and split it into lines
    # Each line should be a separate requirement
    requirements = f.read().splitlines()

setup(
    name="AlleleFinder",
    version=version['__version__'],
    author="Adam Koziol",
    author_email="adam.koziol@inspection.gc.ca",
    url="https://github.com/OLC-Bioinformatics/AlleleFinder",
    long_description=long_description,
    # The content type of the long description. Necessary for PyPI
    long_description_content_type='text/markdown',
    packages=find_packages(),
    include_package_data=True,
    scripts=[
        os.path.join('allele_tools', 'allele_updater.py'),
        os.path.join('allele_tools', 'allele_translate_reduce.py'),
        os.path.join('allele_tools', 'profile_reduce.py'),
        os.path.join('allele_tools', 'stec.py'),
        os.path.join('allele_tools', 'stec_combined_subunits.py'),
        os.path.join('allele_tools', 'split_db.py'),
    ],
    entry_points={
        'console_scripts': [
            'allele_updater=allele_tools.allele_updater:cli',
            'allele_translate_reduce=allele_tools.allele_translate_reduce:cli',
            'profile_reduce=allele_tools.profile_reduce:cli',
            'stec=allele_tools.stec:cli',
            'split_db=allele_tools.split_db:cli'
        ]
    },
    install_requires=requirements,
    # Classifiers categorize the project for users.
    classifiers=[
        # Specifies the intended audience of the project
        'Intended Audience :: Science/Research',
        # Defines the license of the project
        'License :: OSI Approved :: MIT License',
        # Specifies the supported Python versions
        'Programming Language :: Python :: 3.9',
        # Indicates the development status of the project
        'Development Status :: 5 - Production/Stable',
        # Specifies the topic related to the project
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
