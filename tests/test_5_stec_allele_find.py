#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/stec.py Only testing the allele_find functionality
"""

# Standard imports
from unittest.mock import patch
from glob import glob
import argparse
import os

# Third-party imports
import pytest

# Local imports
from allele_tools.stec import \
    allele_find, \
    cli
from .test_2_allele_updater import \
    clean_outputs, \
    prepare_files


@pytest.fixture(name='variables', scope='module')
def setup():
    class Variables:
        def __init__(self):
            # Extract the connection string prior to running tests
            self.test_path = os.path.abspath(os.path.dirname(__file__))
            self.file_path = os.path.join(self.test_path, 'test_files', 'stec_allele_find')
            self.original_files_path = os.path.join(self.file_path, 'original_files')
            self.original_nt_allele_files = glob(os.path.join(self.original_files_path, 'nt_alleles', 'ECs12*'))
            self.original_nt_profile_file = os.path.join(self.original_files_path, 'nt_profile', 'profile.txt')
            self.original_nt_query_files = glob(os.path.join(self.original_files_path, 'nt_query', '*.fasta'))
            self.original_aa_allele_files = glob(os.path.join(self.original_files_path, 'aa_alleles', 'ECs12*'))
            self.original_aa_profile_path = os.path.join(self.original_files_path, 'aa_profile')
            self.original_aa_profile_file = os.path.join(self.original_aa_profile_path, 'profile.txt')
            self.original_aa_query_files = glob(os.path.join(self.original_files_path, 'aa_query', '*.fasta'))
            self.nt_allele_path = os.path.join(self.file_path, 'nt_alleles')
            self.nt_profile_file = os.path.join(self.file_path, 'nt_profile', 'profile.txt')
            self.nt_profile_path = os.path.join(self.file_path, 'nt_profile')
            self.aa_allele_path = os.path.join(self.file_path, 'aa_alleles')
            self.aa_profile_path = os.path.join(self.file_path, 'aa_profile')
            self.aa_profile_file = os.path.join(self.aa_profile_path, 'profile.txt')
            self.aa_notes_path = os.path.join(self.file_path, 'aa_notes')
            self.aa_reports_path = os.path.join(self.file_path, 'aa_reports')
            self.aa_report_file = os.path.join(self.aa_reports_path, 'aa_profiles.tsv')
            self.query_path = os.path.join(self.file_path, 'query')
            self.report_path = os.path.join(self.file_path, 'reports')
            self.report_file = os.path.join(self.report_path, 'profiles.tsv')
            self.length_dict = {
                'stx1B': 82,
                'stx1A': 313,
                'stx2A': 313,
                'stx2B': 84
            }
            self.fake_path = os.path.join('~', 'completely_fake_path')

    return Variables()


@patch('argparse.ArgumentParser.parse_args')
def test_allele_find_integration_no_files(mock_args, variables):
    with pytest.raises(SystemExit):
        mock_args.return_value = argparse.Namespace(
            nt_profile=variables.nt_profile_file,
            aa_profile=variables.aa_profile_file,
            nt_alleles=variables.nt_allele_path,
            aa_alleles=variables.aa_allele_path,
            query_path=variables.query_path,
            report_path=variables.report_path
        )
        arguments = cli()
        allele_find(args=arguments)


@patch('argparse.ArgumentParser.parse_args')
def test_allele_find_integration(mock_args, variables):
    prepare_files(
        variables=variables,
        links=True
    )
    mock_args.return_value = argparse.Namespace(
        nt_profile=variables.nt_profile_file,
        aa_profile=variables.aa_profile_file,
        nt_alleles=variables.nt_allele_path,
        aa_alleles=variables.aa_allele_path,
        query_path=variables.query_path,
        report_path=variables.report_path
    )
    arguments = cli()
    allele_find(args=arguments)
    variables.allele_report = os.path.join(variables.report_path, 'stec_report.tsv')
    assert os.path.isfile(variables.allele_report)


def test_report_contents(variables):
    variables.report_contents = open(variables.allele_report, 'r', encoding='utf-8').readlines()
    assert variables.report_contents[1] == \
        '2013-SEQ-0039	6	9	241	6	2	68\t\n'
    assert variables.report_contents[16] == \
           '2017-SEQ-0617	0	2	116	0		N/A	ECs1206 amino acid sequence does not start with M; ' \
           'ECs1206 amino acid sequence was 14 amino acid residues. Minimum allowed length is 84 amino acid residues; ' \
           'Novel nt_seq_type 116, and aa_seq_type N/A\n'


def test_clean(variables):
    clean_outputs(variables=variables)
