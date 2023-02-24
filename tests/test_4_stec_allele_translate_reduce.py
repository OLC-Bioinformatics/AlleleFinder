#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/stec.py Only testing the profile_reduce functionality
"""

# Standard imports
from unittest.mock import patch
import argparse
import shutil
import os

# Local imports
from allele_tools.stec import \
    translate_reduce, \
    cli
from .test_1_allele_translate_reduce import \
    setup, \
    clean_outputs

assert vars(setup)['_pytestfixturefunction'].name == 'variables'


def prepare_files(variables):
    if not os.path.isdir(variables.nt_profile_path):
        os.makedirs(variables.nt_profile_path)
    shutil.copyfile(
        src=variables.nt_profile_file,
        dst=os.path.join(variables.nt_profile_path, 'profile.txt')
    )
    if not os.path.isdir(variables.nt_allele_path):
        os.makedirs(variables.nt_allele_path)
    for allele_file in variables.nt_allele_files:
        file_name = os.path.basename(allele_file)
        shutil.copyfile(
            src=allele_file,
            dst=os.path.join(variables.nt_allele_path, file_name)
        )


@patch('argparse.ArgumentParser.parse_args')
def test_allele_translate_reduce_integration(mock_args, variables):
    prepare_files(variables=variables)
    mock_args.return_value = argparse.Namespace(
        allele_path=variables.nt_allele_path,
        profile_file=variables.nt_profile,
        report_path=variables.aa_profile_path,
        translated_path=variables.aa_allele_path
    )
    arguments = cli()
    translate_reduce(args=arguments)


def test_stec_profile_header(variables):
    from .test_1_allele_translate_reduce import test_profile_header
    test_profile_header(variables=variables)


def test_stec_profile_contents(variables):
    from .test_1_allele_translate_reduce import test_profile_contents
    test_profile_contents(variables=variables)


def test_stec_profile_end(variables):
    from .test_1_allele_translate_reduce import test_profile_end
    test_profile_end(variables=variables)


def test_stec_clean_outputs(variables):
    clean_outputs(variables=variables)
