#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/stec.py Only testing the allele_find functionality
"""

# Standard imports
from unittest.mock import patch
import argparse
import os

# Third-party imports
import pytest

# Local imports
from allele_tools.stec import \
    aa_allele_find, \
    cli
from .test_2_allele_updater import \
    setup, \
    prepare_files, \
    clean_outputs

assert vars(setup)['_pytestfixturefunction'].name == 'variables'


@patch('argparse.ArgumentParser.parse_args')
def test_allele_find_integration_missing_files(mock_args, variables):
    with pytest.raises(SystemExit):
        mock_args.return_value = argparse.Namespace(
            aa_alleles=variables.aa_allele_path,
            query_path=variables.query_path,
            report_path=variables.report_path,
            cutoff=90,
        )
        arguments = cli()
        aa_allele_find(args=arguments)

@patch('argparse.ArgumentParser.parse_args')
def test_allele_find_integration(mock_args, variables):
    prepare_files(
        variables=variables,
        nucleotide=False,
        links=True
    )
    mock_args.return_value = argparse.Namespace(
        aa_alleles=variables.aa_allele_path,
        query_path=variables.query_path,
        report_path=variables.report_path,
        cutoff=90,
    )
    arguments = cli()
    aa_allele_find(args=arguments)
    variables.allele_report = os.path.join(variables.report_path, 'allele_report.tsv')
    assert os.path.isfile(variables.allele_report)


def test_report_contents(variables):
    variables.report_contents = open(variables.allele_report, 'r', encoding='utf-8').readlines()
    assert variables.report_contents[12] == \
           'ECs1205_9	ECs1205_1000000;ECs1205_9	Amino acid matches previous result: ECs1205_9\n'


def test_clean(variables):
    clean_outputs(variables=variables)
