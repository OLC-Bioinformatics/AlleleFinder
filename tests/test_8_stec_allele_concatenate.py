#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/stec.py allele_concatenate
"""

# Standard imports
from unittest.mock import patch
import argparse
import shutil
import os

# Third-party imports
import pytest

# Local imports
from allele_tools.allele_profiler import read_profile
from allele_tools.methods import (
    concatenate_alleles,
    load_alleles,
    write_concatenated_sequences
)
from allele_tools.stec import (
    allele_concatenate,
    cli
)


@pytest.fixture(name='variables', scope='module')
def setup():
    class Variables:
        def __init__(self):
            # Extract the connection string prior to running tests
            self.test_path = os.path.abspath(os.path.dirname(__file__))
            self.file_path = os.path.join(self.test_path, 'test_files', 'allele_updater', 'original_files')
            self.nt_allele_path = os.path.join(self.file_path, 'nt_alleles')
            self.aa_allele_path = os.path.join(self.file_path, 'aa_alleles')
            self.nt_profile_file = os.path.join(self.file_path, 'nt_profile', 'profile.txt')
            self.aa_profile_file = os.path.join(self.file_path, 'aa_profile', 'profile.txt')
            self.concatenate_path = os.path.join(self.test_path, 'test_files', 'concatenate_alleles')
            self.linker_length_dict = {
                'stx1': 9,
                'stx2': 12,
            }
            self.allele_order = {
                'stx1': ['stx1B', 'stx1B'],
                'stx2': ['stx2A', 'stx2B']
            }
            self.gene_allele = {
                'stx1': 'stx1A_stx1B',
                'stx2': 'stx1A_stx2B'
            }

    return Variables()


def prepare_files(variables):
    if not os.path.isdir(variables.concatenate_path):
        os.makedirs(variables.concatenate_path)
    assert os.path.isdir(variables.concatenate_path)


def clear_path(variables):
    shutil.rmtree(variables.concatenate_path)
    assert not os.path.isdir(variables.concatenate_path)


def test_read_profile(variables):
    variables.nt_profile_data = read_profile(
        profile_file=variables.nt_profile_file
    )
    assert variables.nt_profile_data['113'] == {
        'ECs1205': '12',
        'ECs1206': '1'
    }


def test_load_alleles(variables):
    variables.gene, variables.nt_alleles = load_alleles(
        allele_path=variables.nt_allele_path,
        allele_order=variables.allele_order
    )
    assert variables.gene == 'stx2'
    assert variables.nt_alleles['ECs1205']['ECs1205_1'][:50] == 'ATGAAGTGTATATTATTTAAATGGGTACTGTGCCTGTTACTGGGTTTTTC'


def test_concatenate_alleles(variables):
    variables.concatenated_nt_seq = concatenate_alleles(
        profile_data=variables.nt_profile_data,
        allele_dict=variables.nt_alleles,
        allele_order=variables.allele_order,
        stx_gene=variables.gene,
        linker_length_dict=variables.linker_length_dict,
        molecule='nt'
    )
    assert str(variables.concatenated_nt_seq[0].seq)[:50] == 'ATGAAGTGTATATTATTTAAATGGGTACTGTGCCTGTTACTGGGTTTTTC'
    assert str(variables.concatenated_nt_seq[0].seq)[-50:] == 'CCTGTGAATCAGGCTCCGGATTTGCTGAAGTGCAGTTTAATAATGACTGA'


def test_write_concatenated_alleles(variables):
    prepare_files(variables=variables)
    write_concatenated_sequences(
        concatenated_sequences=variables.concatenated_nt_seq,
        concatenate_path=variables.concatenate_path,
        file_name=variables.gene_allele[variables.gene],
        molecule='nt'
    )
    assert os.path.isfile(os.path.join(variables.concatenate_path, 'nt', 'ECs1205_ECs1206.fasta'))
    clear_path(variables=variables)


@patch('argparse.ArgumentParser.parse_args')
def test_allele_concatenate_integration(mock_args, variables):
    prepare_files(variables=variables)
    mock_args.return_value = argparse.Namespace(
        nt_profile=variables.nt_profile_file,
        aa_profile=variables.aa_profile_file,
        nt_alleles=variables.nt_allele_path,
        aa_alleles=variables.aa_allele_path,
        concatenate_path=variables.concatenate_path
    )
    arguments = cli()
    allele_concatenate(args=arguments)
    assert os.path.isfile(os.path.join(variables.concatenate_path, 'aa', 'ECs1205_ECs1206.fasta'))
    clear_path(variables=variables)
