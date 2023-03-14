#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/allele_translate_reduce.py
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
from allele_tools.allele_updater import \
    Updater, \
    cli


@pytest.fixture(name='variables', scope='module')
def setup():
    class Variables:
        def __init__(self):
            # Extract the connection string prior to running tests
            self.test_path = os.path.abspath(os.path.dirname(__file__))
            self.file_path = os.path.join(self.test_path, 'test_files', 'allele_updater')
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
                'ECs2973': 90,
                'ECs2974': 316,
                'ECs1205': 320,
                'ECs1206': 88
            }
            self.fake_path = os.path.join('~', 'completely_fake_path')

    return Variables()


def clean_outputs(variables):
    if os.path.isdir(variables.nt_allele_path):
        shutil.rmtree(variables.nt_allele_path)
    if os.path.isdir(variables.nt_profile_path):
        shutil.rmtree(variables.nt_profile_path)
    if os.path.isdir(variables.aa_allele_path):
        shutil.rmtree(variables.aa_allele_path)
    if os.path.isdir(variables.aa_profile_path):
        shutil.rmtree(variables.aa_profile_path)
    if os.path.isdir(variables.query_path):
        shutil.rmtree(variables.query_path)
    if os.path.isdir(variables.report_path):
        shutil.rmtree(variables.report_path)
    if os.path.isdir(variables.aa_notes_path):
        shutil.rmtree(variables.aa_notes_path)
    if os.path.isdir(variables.aa_reports_path):
        shutil.rmtree(variables.aa_reports_path)


def prepare_files(variables, nucleotide=True, links=False):
    if not os.path.isdir(variables.nt_allele_path):
        os.makedirs(name=variables.nt_allele_path)
    for nt_allele in variables.original_nt_allele_files:
        shutil.copyfile(
            src=nt_allele,
            dst=os.path.join(variables.nt_allele_path, os.path.basename(nt_allele))
        )
    if not os.path.isdir(variables.nt_profile_path):
        os.makedirs(name=variables.nt_profile_path)
    shutil.copyfile(
        src=variables.original_nt_profile_file,
        dst=os.path.join(variables.nt_profile_path, os.path.basename(variables.original_nt_profile_file))
    )
    if not os.path.isdir(variables.aa_allele_path):
        os.makedirs(name=variables.aa_allele_path)
    for aa_allele in variables.original_aa_allele_files:
        shutil.copyfile(
            src=aa_allele,
            dst=os.path.join(variables.aa_allele_path, os.path.basename(aa_allele))
        )
    if not os.path.isdir(variables.aa_profile_path):
        os.makedirs(name=variables.aa_profile_path)
    shutil.copyfile(
        src=variables.original_aa_profile_file,
        dst=os.path.join(variables.aa_profile_path, os.path.basename(variables.original_aa_profile_file))
    )
    if links:
        link_file = 'aa_nt_profile_links.tsv'
        shutil.copyfile(
            src=os.path.join(variables.original_aa_profile_path, link_file),
            dst=os.path.join(variables.aa_profile_path, link_file)
        )
    os.makedirs(name=variables.query_path)
    if nucleotide:
        for nt_query in variables.original_nt_query_files:
            shutil.copyfile(
                src=nt_query,
                dst=os.path.join(variables.query_path, os.path.basename(nt_query))
            )
    else:
        for aa_query in variables.original_aa_query_files:
            shutil.copyfile(
                src=aa_query,
                dst=os.path.join(variables.query_path, os.path.basename(aa_query))
            )


@patch('argparse.ArgumentParser.parse_args')
def test_allele_updater_integration(mock_args, variables):
    prepare_files(variables=variables)
    mock_args.return_value = argparse.Namespace(
        path=variables.file_path,
        amino_acid=False,
    )
    cli()
    assert os.path.isfile(variables.report_file)


def test_report_contents(variables):
    variables.report_contents = open(variables.report_file, 'r', encoding='utf-8').readlines()
    assert variables.report_contents[1] == '2013-SEQ-0039	259156	6	8\n'


def test_novel_alleles(variables):
    novel_allele_file = os.path.join(variables.report_path, 'nt_ECs1205_novel_alleles.fasta')
    novel_alleles = open(novel_allele_file, 'r', encoding='utf-8').readlines()
    assert novel_alleles[0] == '>ECs1205_875\n'
    clean_outputs(variables=variables)


def test_tilde_path(variables):
    Updater(
        path=variables.fake_path,
        amino_acid=False
    )
    abs_path = os.path.abspath(os.path.expanduser(os.path.join(variables.fake_path)))
    assert os.path.isdir(abs_path)
    shutil.rmtree(abs_path)


@patch('argparse.ArgumentParser.parse_args')
def test_allele_updater_integration_aa_query(mock_args, variables):
    prepare_files(
        variables=variables,
        nucleotide=False
    )
    mock_args.return_value = argparse.Namespace(
        path=variables.file_path,
        amino_acid=True,
    )
    cli()
    assert os.path.isfile(variables.aa_report_file)


def test_aa_report_contents(variables):
    aa_report_contents = open(variables.aa_report_file, 'r', encoding='utf-8').readlines()
    assert aa_report_contents[2] == 'ECs1205_11	2	11	0\n'
    clean_outputs(variables=variables)


@patch('argparse.ArgumentParser.parse_args')
def test_allele_updater_integration_new_aa_files(mock_args, variables):
    prepare_files(variables=variables)
    shutil.rmtree(variables.aa_allele_path)
    mock_args.return_value = argparse.Namespace(
        path=variables.file_path,
        amino_acid=False,
    )
    cli()
    assert os.path.isfile(variables.report_file)


def test_novel_aa_alleles(variables):
    novel_aa_allele_file = os.path.join(variables.aa_allele_path, 'ECs1206_alleles.fasta')
    records = list(SeqIO.parse(novel_aa_allele_file, 'fasta'))
    assert records[0].seq == \
           'MKKMFMAVLFALASVNAMAADCAKGKIEFSKYNEDDTFTVKVDGKEYWTSRWNLQPLLQSAQLTGMTVTIKSSTCESGSGFAEVQFNND*'
    clean_outputs(variables=variables)
