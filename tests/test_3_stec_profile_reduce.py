#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/stec.py Only testing the profile_reduce functionality
"""
import shutil
# Standard imports
from unittest.mock import patch
import argparse
import os

import pytest

# Local imports
from allele_tools.stec import \
    profile_reduce, \
    cli
from .test_0_profile_reduce import setup

assert vars(setup)['_pytestfixturefunction'].name == 'variables'


def test_missing_tilde(variables):
    from .test_0_profile_reduce import test_profile_reduce_missing_tilde_profile
    test_profile_reduce_missing_tilde_profile(variables=variables)


def test_missing_tilde_names(variables):
    from .test_0_profile_reduce import test_profile_reduce_missing_tilde_names
    test_profile_reduce_missing_tilde_names(variables=variables)


@patch('argparse.ArgumentParser.parse_args')
def test_stec_profile_reduce_integration(mock_args, variables):
    mock_args.return_value = argparse.Namespace(
        profile_file=variables.profile,
        gene_names=variables.names,
        output_folder=variables.output_path
    )
    arguments = cli()
    profile_reduce(args=arguments)
    assert os.path.isfile(variables.output_profiles)


def test_header(variables):
    from .test_0_profile_reduce import test_output_profile_header
    test_output_profile_header(variables=variables)


def test_sequence_type(variables):
    from .test_0_profile_reduce import test_output_profile_sequence_type
    test_output_profile_sequence_type(variables=variables)


def test_notes_header(variables):
    from .test_0_profile_reduce import test_output_notes_header
    test_output_notes_header(variables=variables)


def test_notes_novel(variables):
    from .test_0_profile_reduce import test_output_notes_novel
    test_output_notes_novel(variables=variables)


def test_notes_duplicate(variables):
    from .test_0_profile_reduce import test_output_notes_duplicate
    test_output_notes_duplicate(variables=variables)


def test_notes_length(variables):
    from .test_0_profile_reduce import test_output_notes_length
    test_output_notes_length(variables=variables)


@patch('argparse.ArgumentParser.parse_args')
def test_stec_profile_reduce_integration_no_genes(mock_args, variables):
    with pytest.raises(SystemExit):
        mock_args.return_value = argparse.Namespace(
            profile_file=variables.profile,
            gene_names=variables.names + '.fake',
            output_folder=variables.output_path
        )
        arguments = cli()
        profile_reduce(args=arguments)


@patch('argparse.ArgumentParser.parse_args')
def test_stec_profile_reduce_integration_no_folder(mock_args, variables):
    variables.fake_dir = os.path.join(variables.file_path, 'fake')
    fake_names = os.path.join(variables.fake_dir, variables.names)
    os.makedirs(variables.fake_dir)
    mock_args.return_value = argparse.Namespace(
        profile_file=variables.profile,
        gene_names=fake_names,
        output_folder=variables.output_path
    )
    arguments = cli()
    profile_reduce(args=arguments)
    assert os.path.isdir(variables.fake_dir)


@patch('argparse.ArgumentParser.parse_args')
def test_stec_profile_reduce_integration_no_file_or_folder(mock_args, variables):
    with pytest.raises(SystemExit):
        variables.fake_dir_fake = os.path.join(variables.file_path, 'fake_fake')
        fake_names = os.path.join(variables.fake_dir_fake, 'fake.txt')
        mock_args.return_value = argparse.Namespace(
            profile_file=variables.profile,
            gene_names=fake_names,
            output_folder=variables.output_path
        )
        arguments = cli()
        profile_reduce(args=arguments)


@patch('argparse.ArgumentParser.parse_args')
def test_stec_profile_reduce_integration_empty_file(mock_args, variables):
    with pytest.raises(SystemExit):
        variables.fake_dir_two = os.path.join(variables.file_path, 'fake_path')
        os.makedirs(variables.fake_dir_two)
        fake_names = os.path.join(variables.fake_dir_two, 'fake.txt')
        with open(fake_names, 'w', encoding='utf-8') as fake:
            fake.write('')
        mock_args.return_value = argparse.Namespace(
            profile_file=variables.profile,
            gene_names=fake_names,
            output_folder=variables.output_path
        )
        arguments = cli()
        profile_reduce(args=arguments)


def test_clean_up(variables):
    from .test_0_profile_reduce import test_output_clean_up
    test_output_clean_up(variables=variables)
    shutil.rmtree(variables.fake_dir)
    shutil.rmtree(variables.fake_dir_two)
    os.remove(variables.names + '.fake')
