#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/allele_translate_reduce.py
"""

# Standard library imports
import argparse
import os
import shutil
from glob import glob
from unittest.mock import patch

# Third-party imports
import pytest
from Bio import SeqIO

# Local imports
from allele_tools.allele_updater import (
    cli,
    Updater
)


@pytest.fixture(name='variables', scope='module')
def setup():
    """
    This fixture is used to setup the test environment.
    It initializes a Variables class with paths and other variables used in
    the tests.
    The scope is set to 'module' so the setup runs once per module.
    """
    class Variables:
        def __init__(self):
            """
            Initialize the Variables class with paths and other variables used
            in the tests.
            """
            # Define the absolute path to the directory containing the
            # test file
            self.test_path = os.path.abspath(os.path.dirname(__file__))

            # Define the path to the directory containing the test files
            self.file_path = os.path.join(
                self.test_path, 'test_files', 'allele_updater'
            )

            # Define the path to the directory containing the original files
            self.original_files_path = os.path.join(
                self.file_path, 'original_files'
            )

            # Define the paths to the original nucleotide allele files
            self.original_nt_allele_files = glob(
                os.path.join(self.original_files_path, 'nt_alleles', 'stx2*')
            )

            # Define the path to the original nucleotide profile file
            self.original_nt_profile_file = os.path.join(
                self.original_files_path, 'nt_profile', 'profile.txt'
            )

            # Define the paths to the original nucleotide query files
            self.original_nt_query_files = glob(
                os.path.join(self.original_files_path, 'nt_query', '*.fasta')
            )

            # Define the paths to the original amino acid allele files
            self.original_aa_allele_files = glob(
                os.path.join(self.original_files_path, 'aa_alleles', 'stx2*')
            )

            # Define the path to the directory containing the original
            # amino acid profile
            self.original_aa_profile_path = os.path.join(
                self.original_files_path, 'aa_profile'
            )

            # Define the path to the original amino acid profile file
            self.original_aa_profile_file = os.path.join(
                self.original_aa_profile_path, 'profile.txt'
            )

            # Define the paths to the original amino acid query files
            self.original_aa_query_files = glob(
                os.path.join(self.original_files_path, 'aa_query', '*.fasta')
            )

            # Define the path to the directory containing the nucleotide
            # alleles
            self.nt_allele_path = os.path.join(
                self.file_path, 'nt_alleles'
            )

            # Define the path to the nucleotide profile file
            self.nt_profile_file = os.path.join(
                self.file_path, 'nt_profile', 'profile.txt'
            )

            # Define the path to the directory containing the
            # nucleotide profile
            self.nt_profile_path = os.path.join(
                self.file_path, 'nt_profile'
            )

            # Define the path to the directory containing the amino
            # acid alleles
            self.aa_allele_path = os.path.join(
                self.file_path, 'aa_alleles'
            )

            # Define the path to the directory containing the amino
            # acid profile
            self.aa_profile_path = os.path.join(
                self.file_path, 'aa_profile'
            )

            # Define the path to the amino acid profile file
            self.aa_profile_file = os.path.join(
                self.aa_profile_path, 'profile.txt'
            )

            # Define the path to the directory containing the amino acid notes
            self.aa_notes_path = os.path.join(
                self.file_path, 'aa_notes'
            )

            # Define the path to the directory containing the query files
            self.query_path = os.path.join(
                self.file_path, 'query'
            )

            # Define the path to the directory containing the report files
            self.report_path = os.path.join(
                self.file_path, 'reports'
            )

            # Define the path to the amino acid report file
            self.aa_report_file = os.path.join(
                self.report_path, 'aa_profiles.tsv'
            )

            # Define the path to the nucleotide report file
            self.report_file = os.path.join(
                self.report_path, 'nt_profiles.tsv'
            )

            # Define a dictionary for the expected length of each allele
            self.length_dict = {
                'stx1B': 82,
                'stx1A': 313,
                'stx2A': 313,
                'stx2B': 84
            }

            # Define a fake path for testing
            self.fake_path = os.path.join('~', 'completely_fake_path')

    return Variables()


def clean_outputs(variables):
    """
    This function cleans up the output directories by removing them if
    they exist.

    :param variables: An instance of the Variables class containing paths to
    various directories that need to be cleaned up.
    """
    # Remove the nucleotide allele directory if it exists
    if os.path.isdir(variables.nt_allele_path):
        shutil.rmtree(variables.nt_allele_path)

    # Remove the nucleotide profile directory if it exists
    if os.path.isdir(variables.nt_profile_path):
        shutil.rmtree(variables.nt_profile_path)

    # Remove the amino acid allele directory if it exists
    if os.path.isdir(variables.aa_allele_path):
        shutil.rmtree(variables.aa_allele_path)

    # Remove the amino acid profile directory if it exists
    if os.path.isdir(variables.aa_profile_path):
        shutil.rmtree(variables.aa_profile_path)

    # Remove the query directory if it exists
    if os.path.isdir(variables.query_path):
        shutil.rmtree(variables.query_path)

    # Remove the report directory if it exists
    if os.path.isdir(variables.report_path):
        shutil.rmtree(variables.report_path)

    # Remove the amino acid notes directory if it exists
    if os.path.isdir(variables.aa_notes_path):
        shutil.rmtree(variables.aa_notes_path)


def prepare_files(variables, nucleotide=True, links=False):
    """
    This function prepares the files for testing by copying them from the
    original locations to the test locations.

    :param variables: An instance of the Variables class containing paths to
    various directories and files.
    :param nucleotide: A boolean indicating whether to prepare
    nucleotide files.
    :param links: A boolean indicating whether to prepare link files.
    """
    # Create the nucleotide allele directory if it doesn't exist
    if not os.path.isdir(variables.nt_allele_path):
        os.makedirs(name=variables.nt_allele_path)

    # Copy the original nucleotide allele files to the test location
    for nt_allele in variables.original_nt_allele_files:
        shutil.copyfile(
            src=nt_allele,
            dst=os.path.join(
                variables.nt_allele_path,
                os.path.basename(nt_allele)
            )
        )

    # Create the nucleotide profile directory if it doesn't exist
    if not os.path.isdir(variables.nt_profile_path):
        os.makedirs(name=variables.nt_profile_path)

    # Copy the original nucleotide profile file to the test location
    shutil.copyfile(
        src=variables.original_nt_profile_file,
        dst=os.path.join(
            variables.nt_profile_path,
            os.path.basename(variables.original_nt_profile_file)
        )
    )

    # Create the amino acid allele directory if it doesn't exist
    if not os.path.isdir(variables.aa_allele_path):
        os.makedirs(name=variables.aa_allele_path)

    # Copy the original amino acid allele files to the test location
    for aa_allele in variables.original_aa_allele_files:
        shutil.copyfile(
            src=aa_allele,
            dst=os.path.join(
                variables.aa_allele_path,
                os.path.basename(aa_allele)
            )
        )

    # Create the amino acid profile directory if it doesn't exist
    if not os.path.isdir(variables.aa_profile_path):
        os.makedirs(name=variables.aa_profile_path)

    # Copy the original amino acid profile file to the test location
    shutil.copyfile(
        src=variables.original_aa_profile_file,
        dst=os.path.join(
            variables.aa_profile_path,
            os.path.basename(variables.original_aa_profile_file)
        )
    )

    # Copy the link file if links is True
    if links:
        link_file = 'aa_nt_profile_links.tsv'
        shutil.copyfile(
            src=os.path.join(variables.original_aa_profile_path, link_file),
            dst=os.path.join(variables.aa_profile_path, link_file)
        )

    # Create the query directory
    os.makedirs(name=variables.query_path)

    # Copy the original query files to the test location based on the
    # nucleotide flag
    if nucleotide:
        for nt_query in variables.original_nt_query_files:
            shutil.copyfile(
                src=nt_query,
                dst=os.path.join(
                    variables.query_path,
                    os.path.basename(nt_query)
                )
            )
    else:
        for aa_query in variables.original_aa_query_files:
            shutil.copyfile(
                src=aa_query,
                dst=os.path.join(
                    variables.query_path,
                    os.path.basename(aa_query)
                )
            )


@patch('argparse.ArgumentParser.parse_args')
def test_allele_updater_integration(mock_args, variables):
    """
    This function tests the integration of the allele updater by preparing the
    files, mocking the command line arguments, running the command line
    interface, and checking that the report file is created.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class containing paths to
    various directories and files.
    """
    # Prepare the files for testing
    prepare_files(variables=variables)

    # Mock the command line arguments
    mock_args.return_value = argparse.Namespace(
        path=variables.file_path,
        amino_acid=False,
    )

    # Run the command line interface
    cli()

    # Check that the report file is created
    assert os.path.isfile(variables.report_file)


def test_report_contents(variables):
    """
    Test the contents of the report file.

    :param variables: An instance of the Variables class.
    """
    # Open the report file and read the lines into report_contents
    variables.report_contents = open(
        variables.report_file,
        'r',
        encoding='utf-8'
    ).readlines()

    # Assert that the second line of the report file is as expected
    assert variables.report_contents[1] == '2013-SEQ-0039\t259156\t6\t8\n'


def test_novel_alleles(variables):
    """
    Test the contents of the novel alleles file.

    :param variables: An instance of the Variables class.
    """
    # Define the path to the novel allele file
    novel_allele_file = os.path.join(
        variables.report_path,
        'nt_stx2A_novel_alleles.fasta'
    )

    # Open the novel allele file and read the lines
    novel_alleles = open(
        novel_allele_file,
        'r',
        encoding='utf-8'
    ).readlines()

    # Assert that the first line of the novel allele file is as expected
    assert novel_alleles[0] == '>Stx2a_875|958bp\n'

    # Clean up the output directories
    clean_outputs(variables=variables)


def test_tilde_path(variables):
    """
    Test the handling of tilde (~) in the path.

    :param variables: An instance of the Variables class.
    """
    # Create an Updater instance with a fake path
    Updater(
        path=variables.fake_path,
        amino_acid=False
    )

    # Expand the tilde in the path and convert it to an absolute path
    abs_path = os.path.abspath(
        os.path.expanduser(
            os.path.join(
                variables.fake_path
            )
        )
    )

    # Assert that the directory at the absolute path exists
    assert os.path.isdir(abs_path)

    # Remove the directory at the absolute path
    shutil.rmtree(abs_path)


@patch('argparse.ArgumentParser.parse_args')
def test_allele_updater_integration_aa_query(mock_args, variables):
    """
    Test the integration of the allele updater with an amino acid query.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    # Prepare the files for testing with a nucleotide query
    prepare_files(
        variables=variables,
        nucleotide=False
    )

    # Mock the command line arguments for an amino acid query
    mock_args.return_value = argparse.Namespace(
        path=variables.file_path,
        amino_acid=True,
    )

    # Run the command line interface
    cli()

    # Assert that the amino acid report file is created
    assert os.path.isfile(variables.aa_report_file)


def test_aa_report_contents(variables):
    """
    Test the contents of the amino acid report file.

    :param variables: An instance of the Variables class.
    """
    # Open the amino acid report file and read the lines
    aa_report_contents = open(
        variables.aa_report_file,
        'r',
        encoding='utf-8'
    ).readlines()

    # Assert that the third line of the amino acid report file is as expected
    assert aa_report_contents[2] == 'stx2A_11\t175\t11\t0\n'

    # Clean up the output directories
    clean_outputs(variables=variables)


@patch('argparse.ArgumentParser.parse_args')
def test_allele_updater_integration_new_aa_files(mock_args, variables):
    """
    Test the integration of the allele updater with new amino acid files.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    # Prepare the files for testing
    prepare_files(variables=variables)

    # Remove the amino acid allele path
    shutil.rmtree(variables.aa_allele_path)

    # Mock the command line arguments for a nucleotide query
    mock_args.return_value = argparse.Namespace(
        path=variables.file_path,
        amino_acid=False,
    )

    # Run the command line interface
    cli()

    # Assert that the report file is created
    assert os.path.isfile(variables.report_file)


def test_novel_aa_alleles(variables):
    """
    Test the contents of the novel amino acid alleles file.

    :param variables: An instance of the Variables class.
    """
    # Define the path to the novel amino acid allele file
    novel_aa_allele_file = os.path.join(
        variables.aa_allele_path,
        'stx2B_alleles.fasta'
    )

    # Parse the fasta file and convert it to a list of records
    records = list(SeqIO.parse(novel_aa_allele_file, 'fasta'))

    # Assert that the sequence of the first record is as expected
    assert records[0].seq == \
           'MKKMFMAVLFALASVNAMAADCAKGKIEFSKYNEDDTFTVKVDGKEYWTSRWN' \
           'LQPLLQSAQLTGMTVTIKSSTCESGSGFAEVQFNND*'

    # Clean up the output directories
    clean_outputs(variables=variables)
