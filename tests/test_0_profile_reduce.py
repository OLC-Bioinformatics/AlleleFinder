#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/profile_reduce.py
"""

# Standard imports
from unittest.mock import patch
import argparse
import shutil
import os

# Third-party imports
import pytest

# Local imports
from allele_tools.profile_reduce import \
    cli


@pytest.fixture(name='variables', scope='module')
def setup():
    class Variables:
        def __init__(self):
            # Extract the connection string prior to running tests
            self.test_path = os.path.abspath(os.path.dirname(__file__))
            self.file_path = os.path.join(self.test_path, 'test_files', 'profile_reduce')
            self.profile = os.path.join(self.file_path, 'profiles.list')
            self.names = os.path.join(self.file_path, 'genes.txt')
            self.output_path = os.path.join(self.file_path, 'profile')
            self.output_profiles = os.path.join(self.file_path, 'profile', 'profile.txt')
            self.output_notes = os.path.join(self.file_path, 'profile', 'reducing_notes.txt')

    return Variables()


@patch('argparse.ArgumentParser.parse_args')
def test_profile_reduce_missing_tilde_profile(mock_args, variables):
    mock_args.return_value = argparse.Namespace(
        profile='~/fake_file_profiles.list',
        names=variables.names
    )
    with pytest.raises(SystemExit):
        cli()


@patch('argparse.ArgumentParser.parse_args')
def test_profile_reduce_missing_tilde_names(mock_args, variables):
    mock_args.return_value = argparse.Namespace(
        profile=variables.profile,
        names='~/fake_file_genes.txt'
    )
    with pytest.raises(SystemExit):
        cli()


@patch('argparse.ArgumentParser.parse_args')
def test_profile_reduce_integration(mock_args, variables):
    mock_args.return_value = argparse.Namespace(
        profile=variables.profile,
        names=variables.names
    )
    cli()
    assert os.path.isfile(variables.output_profiles)


def test_output_profile_header(variables):
    variables.profiles = open(variables.output_profiles, 'r', encoding='utf-8').readlines()
    assert variables.profiles[0] == 'ST	ECs1205	ECs1206\n'


def test_output_profile_sequence_type(variables):
    assert variables.profiles[1] == '1	0	0\n'


def test_output_notes_header(variables):
    variables.notes = open(variables.output_notes, 'r', encoding='utf-8').readlines()
    assert variables.notes[0] == 'OriginalSequenceType	ReducedSequenceType	Notes\n'


def test_output_notes_novel(variables):
    assert variables.notes[1] == '1	1\n'


def test_output_notes_duplicate(variables):
    assert variables.notes[2] == '2	0	duplicate\n'


def test_output_notes_length(variables):
    assert len(variables.notes) == 10


def test_output_clean_up(variables):
    shutil.rmtree(variables.output_path)
    assert not os.path.isfile(variables.output_path)
