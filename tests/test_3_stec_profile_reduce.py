#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/stec.py Only testing the
profile_reduce functionality
"""

# Standard library imports
import os
import argparse
import shutil
from unittest.mock import patch

# Third-party imports
import pytest

# Local imports
from allele_tools.stec import (
    cli,
    profile_reduce
)
from .test_0_profile_reduce import setup

# Assert that the setup fixture is correctly defined
assert vars(setup)['_pytestfixturefunction'].name == 'variables'


def test_missing_tilde(variables):
    """
    Test the handling of missing tilde in the profile.

    :param variables: An instance of the Variables class.
    """
    # Import the test function for missing tilde in the profile
    from .test_0_profile_reduce import \
        test_profile_reduce_missing_tilde_profile

    # Run the test function with the provided variables
    test_profile_reduce_missing_tilde_profile(variables=variables)


def test_missing_tilde_names(variables):
    """
    Test the handling of missing tilde in the names.

    :param variables: An instance of the Variables class.
    """
    # Import the test function for missing tilde in the names
    from .test_0_profile_reduce import test_profile_reduce_missing_tilde_names

    # Run the test function with the provided variables
    test_profile_reduce_missing_tilde_names(variables=variables)


@patch('argparse.ArgumentParser.parse_args')
def test_stec_profile_reduce_integration(mock_args, variables):
    """
    Test the integration of the STEC profile reduction functionality.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    # Mock the command line arguments for the profile reduction functionality
    mock_args.return_value = argparse.Namespace(
        profile_file=variables.profile,
        gene_names=variables.names,
        output_folder=variables.output_path
    )

    # Run the command line interface and get the arguments
    arguments = cli()

    # Run the profile reduction functionality with the arguments
    profile_reduce(args=arguments)

    # Assert that the output profiles file is created
    assert os.path.isfile(variables.output_profiles)


def test_header(variables):
    """
    Test the header of the output profile.

    :param variables: An instance of the Variables class.
    """
    # Import the test function for the output profile header
    from .test_0_profile_reduce import test_output_profile_header

    # Run the test function with the provided variables
    test_output_profile_header(variables=variables)


def test_sequence_type(variables):
    """
    Test the sequence type of the output profile.

    :param variables: An instance of the Variables class.
    """
    # Import the test function for the output profile sequence type
    from .test_0_profile_reduce import test_output_profile_sequence_type

    # Run the test function with the provided variables
    test_output_profile_sequence_type(variables=variables)


def test_notes_header(variables):
    """
    Test the notes header of the output profile.

    :param variables: An instance of the Variables class.
    """
    # Import the test function for the output notes header
    from .test_0_profile_reduce import test_output_notes_header

    # Run the test function with the provided variables
    test_output_notes_header(variables=variables)


def test_notes_novel(variables):
    """
    Test the novel notes of the output profile.

    :param variables: An instance of the Variables class.
    """
    # Import the test function for the output novel notes
    from .test_0_profile_reduce import test_output_notes_novel

    # Run the test function with the provided variables
    test_output_notes_novel(variables=variables)


def test_notes_duplicate(variables):
    """
    Test the duplicate notes of the output profile.

    :param variables: An instance of the Variables class.
    """
    # Import the test function for the output duplicate notes
    from .test_0_profile_reduce import test_output_notes_duplicate

    # Run the test function with the provided variables
    test_output_notes_duplicate(variables=variables)


def test_notes_length(variables):
    """
    Test the length notes of the output profile.

    :param variables: An instance of the Variables class.
    """
    # Import the test function for the output length notes
    from .test_0_profile_reduce import test_output_notes_length

    # Run the test function with the provided variables
    test_output_notes_length(variables=variables)


@patch('argparse.ArgumentParser.parse_args')
def test_stec_profile_reduce_integration_no_genes(mock_args, variables):
    """
    Test the integration of the STEC profile reduction functionality
    when no genes are provided.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    with pytest.raises(SystemExit):
        # Mock the command line arguments for the profile reduction
        # functionality
        mock_args.return_value = argparse.Namespace(
            profile_file=variables.profile,
            gene_names=variables.names + '.fake',
            output_folder=variables.output_path
        )
        # Run the command line interface and get the arguments
        arguments = cli()
        # Run the profile reduction functionality with the arguments
        profile_reduce(args=arguments)


@patch('argparse.ArgumentParser.parse_args')
def test_stec_profile_reduce_integration_no_folder(mock_args, variables):
    """
    Test the integration of the STEC profile reduction functionality
    when no folder is provided.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    variables.fake_dir = os.path.join(variables.file_path, 'fake')
    fake_names = os.path.join(variables.fake_dir, variables.names)
    os.makedirs(variables.fake_dir)
    # Mock the command line arguments for the profile reduction functionality
    mock_args.return_value = argparse.Namespace(
        profile_file=variables.profile,
        gene_names=fake_names,
        output_folder=variables.output_path
    )
    # Run the command line interface and get the arguments
    arguments = cli()
    # Run the profile reduction functionality with the arguments
    profile_reduce(args=arguments)
    # Assert that the fake directory is created
    assert os.path.isdir(variables.fake_dir)


@patch('argparse.ArgumentParser.parse_args')
def test_stec_profile_reduce_integration_no_file_or_folder(
        mock_args,
        variables):
    """
    Test the integration of the STEC profile reduction functionality
    when no file or folder is provided.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    with pytest.raises(SystemExit):
        variables.fake_dir_fake = os.path.join(
            variables.file_path,
            'fake_fake'
        )
        fake_names = os.path.join(variables.fake_dir_fake, 'fake.txt')
        # Mock the command line arguments for the profile reduction
        # functionality
        mock_args.return_value = argparse.Namespace(
            profile_file=variables.profile,
            gene_names=fake_names,
            output_folder=variables.output_path
        )
        # Run the command line interface and get the arguments
        arguments = cli()
        # Run the profile reduction functionality with the arguments
        profile_reduce(args=arguments)


@patch('argparse.ArgumentParser.parse_args')
def test_stec_profile_reduce_integration_empty_file(mock_args, variables):
    """
    Test the integration of the STEC profile reduction functionality
    when an empty file is provided.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    with pytest.raises(SystemExit):
        variables.fake_dir_two = os.path.join(variables.file_path, 'fake_path')
        os.makedirs(variables.fake_dir_two)
        fake_names = os.path.join(variables.fake_dir_two, 'fake.txt')
        with open(fake_names, 'w', encoding='utf-8') as fake:
            fake.write('')
        # Mock the command line arguments for the profile reduction
        # functionality
        mock_args.return_value = argparse.Namespace(
            profile_file=variables.profile,
            gene_names=fake_names,
            output_folder=variables.output_path
        )
        # Run the command line interface and get the arguments
        arguments = cli()
        # Run the profile reduction functionality with the arguments
        profile_reduce(args=arguments)


def test_clean_up(variables):
    """
    Test the clean up of the output profile.

    :param variables: An instance of the Variables class.
    """
    from .test_0_profile_reduce import test_output_clean_up
    # Run the test function with the provided variables
    test_output_clean_up(variables=variables)
    # Remove the fake directories and files
    shutil.rmtree(variables.fake_dir)
    shutil.rmtree(variables.fake_dir_two)
    os.remove(variables.names + '.fake')
