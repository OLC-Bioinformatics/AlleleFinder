#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/stec.py Only testing the
allele_translate_find functionality
"""

# Standard imports
import argparse
from glob import glob
import os
from unittest.mock import patch

# Third-party imports
import pytest

# Local imports
from allele_tools.stec import (
    allele_translate_find,
    cli
)

from .test_2_allele_updater import (
    clean_outputs,
    prepare_files
)


@pytest.fixture(name='variables', scope='module')
def setup():
    """
    Setup fixture for the tests. This fixture is responsible for setting up the
    necessary variables for the tests.

    :return: An instance of the Variables class.
    """
    class Variables:
        def __init__(self):
            # Define the paths and files for the tests
            # Get the absolute path of the directory of this file
            self.test_path = os.path.abspath(os.path.dirname(__file__))
            # Construct the path to the test files
            self.file_path = os.path.join(
                self.test_path,
                'test_files',
                'stec_allele_find'
            )
            # Construct the path to the original files
            self.original_files_path = os.path.join(
                self.file_path,
                'original_files'
            )
            # Get the paths to the original nt allele files
            self.original_nt_allele_files = glob(
                os.path.join(self.original_files_path, 'nt_alleles', 'stx2*'))
            # Get the path to the original nt profile file
            self.original_nt_profile_file = os.path.join(
                self.original_files_path, 'nt_profile', 'profile.txt')
            # Get the paths to the original nt query files
            self.original_nt_query_files = glob(
                os.path.join(self.original_files_path, 'nt_query', '*.fasta'))
            # Get the paths to the original aa allele files
            self.original_aa_allele_files = glob(
                os.path.join(self.original_files_path, 'aa_alleles', 'stx2*'))
            # Get the path to the original aa profile path
            self.original_aa_profile_path = os.path.join(
                self.original_files_path, 'aa_profile')
            # Get the path to the original aa profile file
            self.original_aa_profile_file = os.path.join(
                self.original_aa_profile_path, 'profile.txt')
            # Get the paths to the original aa query files
            self.original_aa_query_files = glob(
                os.path.join(self.original_files_path, 'aa_query', '*.fasta'))
            # Construct the path to the nt allele path
            self.nt_allele_path = os.path.join(self.file_path, 'nt_alleles')
            # Get the path to the nt profile file
            self.nt_profile_file = os.path.join(
                self.file_path,
                'nt_profile',
                'profile.txt'
            )
            # Construct the path to the nt profile path
            self.nt_profile_path = os.path.join(self.file_path, 'nt_profile')
            # Construct the path to the aa allele path
            self.aa_allele_path = os.path.join(self.file_path, 'aa_alleles')
            # Construct the path to the aa profile path
            self.aa_profile_path = os.path.join(self.file_path, 'aa_profile')
            # Get the path to the aa profile file
            self.aa_profile_file = os.path.join(
                self.aa_profile_path,
                'profile.txt'
            )
            # Construct the path to the aa notes path
            self.aa_notes_path = os.path.join(self.file_path, 'aa_notes')
            # Construct the path to the aa reports path
            self.aa_reports_path = os.path.join(self.file_path, 'aa_reports')
            # Get the path to the aa report file
            self.aa_report_file = os.path.join(
                self.aa_reports_path,
                'aa_profiles.tsv'
            )
            # Construct the path to the query path
            self.query_path = os.path.join(self.file_path, 'query')
            # Construct the path to the report path
            self.report_path = os.path.join(self.file_path, 'reports')
            # Get the path to the report file
            self.report_file = os.path.join(self.report_path, 'profiles.tsv')
            # Define the length dictionary
            self.length_dict = {
                'stx1B': 82,
                'stx1A': 313,
                'stx2A': 313,
                'stx2B': 84
            }
            # Define a fake path for testing
            self.fake_path = os.path.join('~', 'completely_fake_path')

    return Variables()


@patch('argparse.ArgumentParser.parse_args')
def test_allele_translate_find_integration_no_files(mock_args, variables):
    """
    Test the allele_translate_find function when no files are provided.

    :param mock_args: The mock arguments.
    :param variables: The test variables.
    """
    # We expect a SystemExit exception when no files are provided
    with pytest.raises(SystemExit):
        # Mock the arguments to simulate no files being provided
        mock_args.return_value = argparse.Namespace(
            nt_profile=variables.nt_profile_file,
            aa_profile=variables.aa_profile_file,
            nt_alleles=variables.nt_allele_path,
            aa_alleles=variables.aa_allele_path,
            query_path=variables.query_path,
            report_path=variables.report_path
        )
        # Call the command line interface function
        arguments = cli()
        # Call the function under test
        allele_translate_find(args=arguments)


@patch('argparse.ArgumentParser.parse_args')
def test_allele_translate_find_integration(mock_args, variables):
    """
    Test the allele_translate_find function with files provided.

    :param mock_args: The mock arguments.
    :param variables: The test variables.
    """
    # Prepare the files for testing
    prepare_files(
        variables=variables,
        links=True
    )
    # Mock the arguments to simulate files being provided
    mock_args.return_value = argparse.Namespace(
        nt_profile=variables.nt_profile_file,
        aa_profile=variables.aa_profile_file,
        nt_alleles=variables.nt_allele_path,
        aa_alleles=variables.aa_allele_path,
        query_path=variables.query_path,
        report_path=variables.report_path
    )
    # Call the command line interface function
    arguments = cli()
    # Call the function under test
    allele_translate_find(args=arguments)
    # Check that the report file was created
    variables.allele_report = os.path.join(
        variables.report_path,
        'stec_report.tsv'
    )
    assert os.path.isfile(variables.allele_report)


def test_report_contents(variables):
    """
    Test the contents of the report.

    :param variables: The test variables.
    """
    # Open the report file and read its contents
    variables.report_contents = open(
        variables.allele_report,
        'r',
        encoding='utf-8'
    ).readlines()
    # Check that the contents of the report are as expected
    assert variables.report_contents[1] == \
        '2013-SEQ-0039\t875\t137\t259156\t6\t2\t68\tNovel nt_seq_type ' \
        '259156 links to aa_seq type 68\n'
    assert variables.report_contents[10] == '2017-SEQ-0617\t0\t0\t1\t0\t0\t1\n'


def test_clean(variables):
    """
    Test the clean_outputs function.

    :param variables: The test variables.
    """
    # Call the function under test
    clean_outputs(variables=variables)
