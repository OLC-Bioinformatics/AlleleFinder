#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/stec.py allele_split
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
from allele_tools.stec import (
    allele_split,
    cli
)


@pytest.fixture(name='variables', scope='module')
def setup():
    """
    Setup fixture for the tests. This fixture is responsible for setting up the
    necessary variables for the tests.

    :return: An instance of the Variables class.
    """
    class Variables:
        def __init__(self):
            # Define the paths and files for the tests

            # test_path: the path to the directory containing the test files
            self.test_path = os.path.abspath(os.path.dirname(__file__))

            # file_path: the path to the directory containing the allele files
            self.file_path = os.path.join(
                self.test_path,
                'test_files',
                'stec_allele_split'
            )

            # original_file_path: the path to the directory containing
            # the original allele files
            self.original_file_path = os.path.join(
                self.file_path,
                'original_files'
            )

            # original_allele_files: a list of the original allele files
            self.original_allele_files = glob(
                os.path.join(
                    self.original_file_path,
                    '*.fasta'
                )
            )

            # query_path: the path to the directory where the query files
            # will be stored
            self.query_path = os.path.join(
                self.file_path,
                'query'
            )

            # output_path: the path to the directory where the output
            # files will be stored
            self.output_path = os.path.join(
                self.file_path,
                'split_alleles'
            )

    return Variables()


def prepare_files(variables):
    """
    Prepare the necessary files for the tests. This function is responsible for
    creating the query directory and copying the original allele files to it.

    :param variables: An instance of the Variables class.
    """
    # Create the query directory if it does not exist
    if not os.path.isdir(variables.query_path):
        os.makedirs(variables.query_path)
    # Copy each original allele file to the query directory
    for allele_file in variables.original_allele_files:
        shutil.copyfile(
            src=allele_file,
            dst=os.path.join(
                variables.query_path,
                os.path.basename(allele_file)
            )
        )


def clean_outputs(variables):
    """
    Clean up the outputs after the tests. This function is responsible for
    removing the query and output directories.

    :param variables: An instance of the Variables class.
    """
    # Remove the query directory if it exists
    if os.path.isdir(variables.query_path):
        shutil.rmtree(variables.query_path)
    # Remove the output directory if it exists
    if os.path.isdir(variables.output_path):
        shutil.rmtree(variables.output_path)


@patch('argparse.ArgumentParser.parse_args')
def test_profile_reduce_missing_tilde_profile_missing_files(
        mock_args,
        variables):
    """
    Test the profile reduce functionality when files are missing. This test
    checks if the program exits when the necessary files are missing.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    with pytest.raises(SystemExit):
        # Mock the command line arguments for the profile reduce functionality
        mock_args.return_value = argparse.Namespace(
            query_path=variables.query_path,
            output_path=variables.output_path
        )

        # Run the command line interface and get the arguments
        arguments = cli()

        # Run the profile reduce functionality with the arguments
        allele_split(args=arguments)


@patch('argparse.ArgumentParser.parse_args')
def test_profile_reduce_missing_tilde_profile(mock_args, variables):
    """
    Test the profile reduce functionality. This test checks if the program
    correctly splits the alleles and outputs the correct number of files.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    # Prepare the necessary files for the test
    prepare_files(variables=variables)

    # Mock the command line arguments for the profile reduce functionality
    mock_args.return_value = argparse.Namespace(
        query_path=variables.query_path,
        output_path=variables.output_path
    )

    # Run the command line interface and get the arguments
    arguments = cli()

    # Run the profile reduce functionality with the arguments
    allele_split(args=arguments)

    # Check the number of output alleles
    # output_alleles: a list of the output allele files
    variables.output_alleles = sorted(
        glob(
            os.path.join(
                variables.output_path,
                '*.fasta'
            )
        )
    )
    assert len(variables.output_alleles) == 1096


def test_allele_contents(variables):
    """
    Test the contents of the alleles. This test checks if the first output
    allele has the correct ID and sequence.

    :param variables: An instance of the Variables class.
    """
    # Check the contents of the first output allele
    for record in SeqIO.parse(variables.output_alleles[0], 'fasta'):
        assert record.id == 'stx2A_1'
        assert str(record.seq) == \
               'MKCILFKWVLCLLLGFSSVSYSREFTIDFSTQQSYVSSLNSIRTEISTPLEHISQGTTSV' \
               'SVINHTPPGSYFAVDIRGLDVYQARFDHLRLIIEQNNLYVAGFVNTATNTFYRFSDFTHI' \
               'SVPGVTTVSMTTDSSYTTLQRVAALERSGMQISRHSLVSSYLALMEFSGNTMTRDASRAV' \
               'LRFVTVTAEALRFRQIQREFRQALSETAPVYTMTPGDVDLTLNWGRISNVLPEYRGEDGV' \
               'RVGRISFNNISAILGTVAVILNCHHQGARSVRAVNEESQPECQITGDRPVIKINNTLWES' \
               'NTAAAFLNRKSQFLYTTGK*'


def test_clean_outputs(variables):
    """
    Test the clean up of the outputs. This test checks if the query and output
    directories have been removed.

    :param variables: An instance of the Variables class.
    """
    clean_outputs(variables=variables)

    # Check if the query and output directories have been removed
    assert not os.path.isdir(variables.query_path)
    assert not os.path.isdir(variables.output_path)
