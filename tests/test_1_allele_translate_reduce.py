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
import pytest

# Local imports
from allele_tools.allele_translate_reduce import \
    Translate, \
    cli


@pytest.fixture(name='variables', scope='module')
def setup():
    class Variables:
        def __init__(self):
            # Extract the connection string prior to running tests
            self.test_path = os.path.abspath(os.path.dirname(__file__))
            self.file_path = os.path.join(self.test_path, 'test_files', 'allele_translate_reduce')
            self.nt_profile_path = os.path.join(self.file_path, 'nt_profile')
            self.nt_profile = os.path.join(self.nt_profile_path, 'profile.txt')
            self.nt_profile_file = os.path.join(self.file_path, 'original_files', 'profile.txt')
            self.nt_allele_files = glob(os.path.join(self.file_path, 'original_files', 'ECs12*'))
            self.nt_allele_path = os.path.join(self.file_path, 'nt_alleles')
            self.aa_allele_path = os.path.join(self.file_path, 'aa_alleles')
            self.aa_profile_path = os.path.join(self.file_path, 'aa_profile')
            self.aa_profile_file = os.path.join(self.aa_profile_path, 'profile.txt')
            self.length_dict = {
                'ECs2973': 90,
                'ECs2974': 316,
                'ECs1205': 316,
                'ECs1206': 88
            }
            self.fake_path = os.path.join('~', 'completely_fake_path')

    return Variables()


def clean_outputs(variables):
    if os.path.isdir(variables.aa_allele_path):
        shutil.rmtree(variables.aa_allele_path)
    if os.path.isdir(variables.aa_profile_path):
        shutil.rmtree(variables.aa_profile_path)
    if os.path.isdir(variables.nt_profile_path):
        shutil.rmtree(variables.nt_profile_path)
    if os.path.isdir(variables.nt_allele_path):
        shutil.rmtree(variables.nt_allele_path)
    allele_files = glob(os.path.join(variables.file_path, '*.fasta'))
    for allele_file in allele_files:
        os.remove(allele_file)
    text_files = glob(os.path.join(variables.file_path, '*.txt'))
    for text_file in text_files:
        os.remove(text_file)
    assert not os.path.isdir(variables.aa_allele_path)
    assert not os.path.isdir(variables.aa_profile_path)


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
            dst=os.path.join(variables.file_path, file_name)
        )


@patch('argparse.ArgumentParser.parse_args')
def test_allele_translate_reduce_integration(mock_args, variables):
    prepare_files(variables=variables)
    mock_args.return_value = argparse.Namespace(
        path=variables.file_path,
        profile=True,
        report_path=variables.aa_profile_path,
        translated_path=variables.aa_allele_path
    )
    variables.arguments = cli()
    assert os.path.isfile(variables.aa_profile_file)


def test_profile_header(variables):
    variables.profiles = open(variables.aa_profile_file, 'r', encoding='utf-8').readlines()
    assert variables.profiles[0] == 'ST	ECs1205	ECs1206\n'


def test_profile_contents(variables):
    assert variables.profiles[1] == '1	0	0\n'


def test_profile_end(variables):
    assert variables.profiles[-1] == '259155	30	13\n'
    clean_outputs(variables=variables)


def test_allele_translate_reduce_length_dict(variables):
    prepare_files(variables=variables)
    translate = Translate(
        path=variables.file_path,
        profile=True,
        report_path=variables.aa_profile_path,
        translated_path=variables.aa_allele_path,
        length_dict=variables.length_dict
    )
    translate.main()
    assert os.path.isfile(os.path.join(variables.aa_profile_path, 'profile.txt'))
    clean_outputs(variables=variables)


def test_allele_translate_reduce_tilde_path(variables):
    with pytest.raises(SystemExit):
        Translate(
            path='~/completely_fake_path',
            profile=True,
            report_path=variables.aa_profile_path,
            translated_path=variables.aa_allele_path,
            length_dict=variables.length_dict
        )


def test_allele_translate_reduce_profile_provided(variables):
    with pytest.raises(SystemExit):
        Translate(
            path='~/completely_fake_path',
            profile=variables.nt_profile_file,
            report_path=variables.aa_profile_path,
            translated_path=variables.aa_allele_path,
            length_dict=variables.length_dict
        )


def test_allele_translate_reduce_tilde_report_path(variables):
    prepare_files(variables=variables)
    translate = Translate(
        path=variables.file_path,
        profile=True,
        report_path=variables.fake_path,
        translated_path=variables.aa_allele_path,
        length_dict=variables.length_dict
    )
    assert os.path.isdir(translate.report_path)
    shutil.rmtree(translate.report_path)
    clean_outputs(variables=variables)


def test_allele_translate_reduce_no_profile(variables):
    prepare_files(variables=variables)
    translate = Translate(
        path=variables.file_path,
        profile=False,
        report_path=variables.aa_profile_path,
        translated_path=variables.fake_path,
        length_dict=variables.length_dict
    )
    assert os.path.isdir(translate.translated_path)
    shutil.rmtree(translate.translated_path)
    clean_outputs(variables=variables)


def test_allele_translate_reduce_missing_alleles(variables):
    if not os.path.isdir(variables.nt_profile_path):
        os.makedirs(variables.nt_profile_path)
    shutil.copyfile(
        src=variables.nt_profile_file,
        dst=os.path.join(variables.file_path, 'nt_profile', 'profile.txt')
    )
    with pytest.raises(SystemExit):
        Translate(
            path=variables.file_path,
            profile=True,
            report_path=variables.aa_profile_path,
            translated_path=variables.aa_allele_path,
            length_dict=variables.length_dict
        )
    clean_outputs(variables=variables)
