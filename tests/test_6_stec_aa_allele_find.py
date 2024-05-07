#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/stec.py Only testing the
allele_find functionality
"""

# Standard library imports
import argparse
import os
from unittest.mock import patch

# Third party imports
import pytest

# Local application imports
from allele_tools.stec import (
    aa_allele_find,
    cli
)
from .test_2_allele_updater import (
    clean_outputs,
    prepare_files,
    setup
)

# Ensure the setup function is correctly named
assert vars(setup)['_pytestfixturefunction'].name == 'variables'


@patch('argparse.ArgumentParser.parse_args')
def test_allele_find_integration_missing_files(mock_args, variables):
    """
    Test the integration of the allele find functionality when files
    are missing.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    with pytest.raises(SystemExit):
        # Mock the command line arguments for the allele find functionality
        mock_args.return_value = argparse.Namespace(
            aa_alleles=variables.aa_allele_path,
            query_path=variables.query_path,
            report_path=variables.report_path,
            cutoff=90,
        )

        # Run the command line interface and get the arguments
        arguments = cli()

        # Run the allele find functionality with the arguments
        aa_allele_find(args=arguments)


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
        nucleotide=False,
        links=True
    )

    # Mock the command line arguments for the allele find functionality
    mock_args.return_value = argparse.Namespace(
        aa_alleles=variables.aa_allele_path,
        query_path=variables.query_path,
        report_path=variables.report_path,
        cutoff=90,
    )

    # Run the command line interface and get the arguments
    arguments = cli()

    # Run the allele find functionality with the arguments
    aa_allele_find(args=arguments)

    # Check if the allele report file exists
    variables.allele_report = os.path.join(
        variables.report_path,
        'allele_report.tsv'
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
    assert variables.report_contents[3] == 'stx2A_1388\t\t\n'
    assert variables.report_contents[12] == \
           'stx2A_9\tStx2A_868|319aa\tAmino acid allele Stx2A_868|319aa ' \
           'is novel\n'


def test_clean(variables):
    """
    Test the clean up of the outputs.

    :param variables: An instance of the Variables class.
    """
    clean_outputs(variables=variables)
