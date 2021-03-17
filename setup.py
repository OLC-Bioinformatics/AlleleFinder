#!/usr/bin/env python3
from setuptools import setup, find_packages
import os
setup(
    name="AlleleFinder",
    version="0.0.1",
    scripts=[
    	os.path.join('allele_tools', 'allele_finder.py'),
    	os.path.join('allele_tools', 'allele_profiler.py'),
    	os.path.join('allele_tools', 'allele_translate_reduce.py'),
    	os.path.join('allele_tools', 'probe_creator.py'),
    	os.path.join('allele_tools', 'profile_reduce.py'),
    	os.path.join('allele_tools', 'stec_attributer.py'),
    ],
    packages=find_packages(),
    include_package_data=True,
    author="Adam Koziol",
    author_email="adam.koziol@canada.ca",
    url="https://github.com/OLC-Bioinformatics/AlleleFinder",
)
