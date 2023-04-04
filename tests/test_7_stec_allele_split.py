#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/profile_reduce.py
"""

# Standard imports
from unittest.mock import patch
from glob import glob
import argparse
import shutil
import os

# Third-party imports
from Bio import SeqIO
import pytest

# Local imports
from allele_tools.stec import \
    allele_split, \
    cli


@pytest.fixture(name='variables', scope='module')
def setup():
    class Variables:
        def __init__(self):
            # Extract the connection string prior to running tests
            self.test_path = os.path.abspath(os.path.dirname(__file__))
            self.file_path = os.path.join(self.test_path, 'test_files', 'stec_allele_split')
            self.original_file_path = os.path.join(self.file_path, 'original_files')
            self.original_allele_files = glob(os.path.join(self.original_file_path, '*.fasta'))
            self.query_path = os.path.join(self.file_path, 'query')
            self.output_path = os.path.join(self.file_path, 'split_alleles')

    return Variables()


def prepare_files(variables):
    if not os.path.isdir(variables.query_path):
        os.makedirs(variables.query_path)
    for allele_file in variables.original_allele_files:
        shutil.copyfile(
            src=allele_file,
            dst=os.path.join(variables.query_path, os.path.basename(allele_file))
        )


def clean_outputs(variables):
    if os.path.isdir(variables.query_path):
        shutil.rmtree(variables.query_path)
    if os.path.isdir(variables.output_path):
        shutil.rmtree(variables.output_path)


@patch('argparse.ArgumentParser.parse_args')
def test_profile_reduce_missing_tilde_profile_missing_files(mock_args, variables):
    with pytest.raises(SystemExit):
        mock_args.return_value = argparse.Namespace(
            query_path=variables.query_path,
            output_path=variables.output_path
        )
        arguments = cli()
        allele_split(args=arguments)

@patch('argparse.ArgumentParser.parse_args')
def test_profile_reduce_missing_tilde_profile(mock_args, variables):
    prepare_files(variables=variables)
    mock_args.return_value = argparse.Namespace(
        query_path=variables.query_path,
        output_path=variables.output_path
    )
    arguments = cli()
    allele_split(args=arguments)
    variables.output_alleles = sorted(glob(os.path.join(variables.output_path, '*.fasta')))
    assert len(variables.output_alleles) == 1096


def test_allele_contents(variables):
    for record in SeqIO.parse(variables.output_alleles[0], 'fasta'):
        assert record.id == 'ECs1205_1'
        assert str(record.seq) == \
               'MKCILFKWVLCLLLGFSSVSYSREFTIDFSTQQSYVSSLNSIRTEISTPLEHISQGTTSV' \
               'SVINHTPPGSYFAVDIRGLDVYQARFDHLRLIIEQNNLYVAGFVNTATNTFYRFSDFTHI' \
               'SVPGVTTVSMTTDSSYTTLQRVAALERSGMQISRHSLVSSYLALMEFSGNTMTRDASRAV' \
               'LRFVTVTAEALRFRQIQREFRQALSETAPVYTMTPGDVDLTLNWGRISNVLPEYRGEDGV' \
               'RVGRISFNNISAILGTVAVILNCHHQGARSVRAVNEESQPECQITGDRPVIKINNTLWES' \
               'NTAAAFLNRKSQFLYTTGK*'


def test_clean_outputs(variables):
    clean_outputs(variables=variables)
    assert not os.path.isdir(variables.query_path)
    assert not os.path.isdir(variables.output_path)
