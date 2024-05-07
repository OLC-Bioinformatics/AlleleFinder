#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/stec.py Only testing the
allele_find functionality
"""

# Standard library imports
import argparse
from glob import glob
import os
from unittest.mock import patch

# Third party imports
import pytest

# Local application imports
from allele_tools.stec import (
    allele_find,
    cli
)
from .test_2_allele_updater import (
    clean_outputs,
    prepare_files
)


@pytest.fixture(name='variables', scope='module')
def setup():
    """
    Pytest fixture for setting up the test environment.
    This fixture is named 'variables' and has module scope.
    """
    class Variables:
        def __init__(self):
            # Extract the connection string prior to running tests
            self.test_path = os.path.abspath(os.path.dirname(__file__))

            # Define file paths for testing
            self.file_path = os.path.join(
                self.test_path, 'test_files', 'stec_allele_find'
            )
            self.original_files_path = os.path.join(
                self.file_path, 'original_files'
            )
            self.original_nt_allele_files = glob(os.path.join(
                self.original_files_path, 'nt_alleles', 'stx2*'
            ))
            self.original_nt_profile_file = os.path.join(
                self.original_files_path, 'nt_profile', 'profile.txt'
            )
            self.original_nt_query_files = glob(os.path.join(
                self.original_files_path, 'nt_query', '*.fasta'
            ))
            self.original_aa_allele_files = glob(os.path.join(
                self.original_files_path, 'aa_alleles', 'stx2*'
            ))
            self.original_aa_profile_path = os.path.join(
                self.original_files_path, 'aa_profile'
            )
            self.original_aa_profile_file = os.path.join(
                self.original_aa_profile_path, 'profile.txt'
            )
            self.original_aa_query_files = glob(os.path.join(
                self.original_files_path, 'aa_query', '*.fasta'
            ))
            self.nt_allele_path = os.path.join(
                self.file_path, 'nt_alleles'
            )
            self.nt_profile_file = os.path.join(
                self.file_path, 'nt_profile', 'profile.txt'
            )
            self.nt_profile_path = os.path.join(
                self.file_path, 'nt_profile'
            )
            self.aa_allele_path = os.path.join(
                self.file_path, 'aa_alleles'
            )
            self.aa_profile_path = os.path.join(
                self.file_path, 'aa_profile'
            )
            self.aa_profile_file = os.path.join(
                self.aa_profile_path, 'profile.txt'
            )
            self.aa_notes_path = os.path.join(
                self.file_path, 'aa_notes'
            )
            self.aa_reports_path = os.path.join(
                self.file_path, 'aa_reports'
            )
            self.aa_report_file = os.path.join(
                self.aa_reports_path, 'aa_profiles.tsv'
            )
            self.query_path = os.path.join(
                self.file_path, 'query'
            )
            self.report_path = os.path.join(
                self.file_path, 'reports'
            )
            self.report_file = os.path.join(
                self.report_path, 'profiles.tsv'
            )

            # Define length dictionary for testing
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
def test_allele_find_integration_no_files(mock_args, variables):
    """
    Test the integration of the allele find functionality when no files
    are provided.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    with pytest.raises(SystemExit):
        # Mock the command line arguments for the allele find functionality
        mock_args.return_value = argparse.Namespace(
            nt_profile=variables.nt_profile_file,
            aa_profile=variables.aa_profile_file,
            nt_alleles=variables.nt_allele_path,
            aa_alleles=variables.aa_allele_path,
            query_path=variables.query_path,
            report_path=variables.report_path
        )

        # Run the command line interface and get the arguments
        arguments = cli()

        # Run the allele find functionality with the arguments
        allele_find(args=arguments)


@patch('argparse.ArgumentParser.parse_args')
def test_allele_find_integration(mock_args, variables):
    """
    Test the integration of the allele find functionality.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    # Prepare the necessary files for the test
    prepare_files(
        variables=variables,
        links=True
    )

    # Mock the command line arguments for the allele find functionality
    mock_args.return_value = argparse.Namespace(
        nt_profile=variables.nt_profile_file,
        aa_profile=variables.aa_profile_file,
        nt_alleles=variables.nt_allele_path,
        aa_alleles=variables.aa_allele_path,
        query_path=variables.query_path,
        report_path=variables.report_path
    )

    # Run the command line interface and get the arguments
    arguments = cli()

    # Run the allele find functionality with the arguments
    allele_find(args=arguments)

    # Check if the allele report file exists
    variables.allele_report = os.path.join(
        variables.report_path,
        'stec_report.tsv'
    )
    assert os.path.isfile(variables.allele_report)


def test_report_contents(variables):
    """
    Test the contents of the report.

    :param variables: An instance of the Variables class.
    """
    # Read the contents of the allele report file
    variables.report_contents = open(
        variables.allele_report,
        'r',
        encoding='utf-8'
    ).readlines()

    # Check the contents of the allele report file
    assert variables.report_contents[1] == \
        '2013-SEQ-0039\t6\t9\t241\t868\t137\t259156\tNovel nt_seq_type ' \
        '241, and novel aa_seq_type 259156\n'
    assert variables.report_contents[16] == \
        '2017-SEQ-0617	0	2	116	0		N/A	stx2B amino acid ' \
        'sequence does not start with M; ' \
        'stx2B amino acid sequence was 13 amino acid residues. ' \
        'Minimum allowed length is 84 amino acid residues\n'


def test_clean(variables):
    """
    Test the clean up of the outputs.

    :param variables: An instance of the Variables class.
    """
    clean_outputs(variables=variables)
