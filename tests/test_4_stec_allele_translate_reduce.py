#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/stec.py Only testing the
profile_reduce functionality
"""

# Standard imports
import argparse
import os
import shutil
import tempfile
from unittest import TestCase
from unittest.mock import patch

# Third-party imports
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

# Local imports
from allele_tools.allele_translate_reduce import Translate
from allele_tools.stec import (
    cli,
    translate_reduce
)
from .test_1_allele_translate_reduce import (
    clean_outputs,
    setup
)

assert vars(setup)['_pytestfixturefunction'].name == 'variables'


def prepare_files(variables):
    """
    Prepare the necessary files for the test.

    :param variables: An instance of the Variables class.
    """
    # Check if the nucleotide profile path exists, if not, create it
    if not os.path.isdir(variables.nt_profile_path):
        os.makedirs(variables.nt_profile_path)

    # Copy the nucleotide profile file to the nucleotide profile path
    shutil.copyfile(
        src=variables.nt_profile_file,
        dst=os.path.join(variables.nt_profile_path, 'profile.txt')
    )

    # Check if the nucleotide allele path exists, if not, create it
    if not os.path.isdir(variables.nt_allele_path):
        os.makedirs(variables.nt_allele_path)

    # Copy each nucleotide allele file to the nucleotide allele path
    for allele_file in variables.nt_allele_files:
        file_name = os.path.basename(allele_file)
        shutil.copyfile(
            src=allele_file,
            dst=os.path.join(variables.nt_allele_path, file_name)
        )


@patch('argparse.ArgumentParser.parse_args')
def test_allele_translate_reduce_integration(mock_args, variables):
    """
    Test the integration of the allele translate reduce functionality.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    # Prepare the necessary files for the test
    prepare_files(variables=variables)

    # Mock the command line arguments for the allele translate
    # reduce functionality
    mock_args.return_value = argparse.Namespace(
        allele_path=variables.nt_allele_path,
        profile_file=variables.nt_profile,
        report_path=variables.aa_profile_path,
        translated_path=variables.aa_allele_path
    )

    # Run the command line interface and get the arguments
    arguments = cli()

    # Run the allele translate reduce functionality with the arguments
    translate_reduce(args=arguments)


def test_stec_profile_header(variables):
    """
    Test the header of the STEC profile.

    :param variables: An instance of the Variables class.
    """
    from .test_1_allele_translate_reduce import test_profile_header
    test_profile_header(variables=variables)


def test_stec_profile_contents(variables):
    """
    Test the contents of the STEC profile.

    :param variables: An instance of the Variables class.
    """
    from .test_1_allele_translate_reduce import test_profile_contents
    test_profile_contents(variables=variables)


def test_stec_profile_end(variables):
    """
    Test the end of the STEC profile.

    :param variables: An instance of the Variables class.
    """
    from .test_1_allele_translate_reduce import test_profile_end
    test_profile_end(variables=variables)


def test_stec_clean_outputs(variables):
    """
    Test the clean up of the STEC outputs.

    :param variables: An instance of the Variables class.
    """
    clean_outputs(variables=variables)


@patch('argparse.ArgumentParser.parse_args')
def test_allele_translate_reduce_integration_no_gene_path(
        mock_args,
        variables):
    """
    Test the integration of the allele translate reduce functionality
    when no gene path is provided.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    with pytest.raises(SystemExit):
        # Mock the command line arguments for the allele translate
        # reduce functionality
        mock_args.return_value = argparse.Namespace(
            allele_path=variables.nt_allele_path,
            profile_file=variables.nt_profile,
            report_path=variables.aa_profile_path,
            translated_path=variables.aa_allele_path
        )

        # Run the command line interface and get the arguments
        arguments = cli()

        # Run the allele translate reduce functionality with the arguments
        translate_reduce(args=arguments)


@patch('argparse.ArgumentParser.parse_args')
def test_allele_translate_reduce_integration_no_gene_file(
        mock_args,
        variables):
    """
    Test the integration of the allele translate reduce functionality
    when no gene file is provided.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    with pytest.raises(SystemExit):
        # Mock the command line arguments for the allele translate
        # reduce functionality
        mock_args.return_value = argparse.Namespace(
            allele_path=variables.nt_allele_path,
            profile_file=variables.nt_profile,
            report_path=variables.aa_profile_path,
            translated_path=variables.aa_allele_path
        )

        # Run the command line interface and get the arguments
        arguments = cli()

        # Run the allele translate reduce functionality with the arguments
        translate_reduce(args=arguments)


class TestFilterProfiles(TestCase):
    def setUp(self):
        # Create a temporary directory for the test files
        self.test_dir = tempfile.TemporaryDirectory()

        # Create a test profile file
        self.profile_file = os.path.join(self.test_dir.name, 'profile.txt')
        with open(self.profile_file, 'w') as f:
            f.write('ST\tgene1\tgene2\n1\t0\t1\n2\t1\t0\n3\t1\t1')

    def tearDown(self):
        # Clean up the temporary directory
        self.test_dir.cleanup()

    def test_filter_profiles(self):
        # Call the function with a set of sequence types to filter
        Translate.filter_profiles({'1', '3'}, self.profile_file)

        # Check the contents of the updated profile file
        with open(self.profile_file, 'r') as f:
            self.assertEqual(f.read(), 'ST\tgene1\tgene2\n2\t1\t0\n')

        # Check the contents of the filtered profile file
        filtered_file = os.path.join(
            os.path.dirname(self.profile_file), 'filtered_profiles.txt')
        with open(filtered_file, 'r') as f:
            self.assertEqual(f.read(), 'ST\tgene1\tgene2\n1\t0\t1\n3\t1\t1')


class TestFindDuplicates(TestCase):
    def test_find_duplicates(self):
        # Create a dictionary of existing sequences
        nt_sequences = {
            'gene1_1': SeqRecord(Seq('ATG')),
            'gene1_2': SeqRecord(Seq('TAA')),
            'gene2_1': SeqRecord(Seq('ATG')),
        }

        # Create a new sequence that matches an existing one
        nt_sequence = SeqRecord(Seq('ATG'))
        allele = 'gene1_3'
        filtered = False
        note = []

        # Call the function
        filtered, note, nt_sequences = Translate.find_duplicates(
            nt_sequences, nt_sequence, allele, filtered, note)

        # Check that the sequence was not added to the dictionary
        self.assertNotIn(allele, nt_sequences)

        # Check that the filtered flag was set to True
        self.assertTrue(filtered)

        # Check that the note was updated correctly
        self.assertEqual(
            note,
            ['Trimmed nt sequence matches previous allele sequence: '
                'gene1_1', 'Trimmed nt sequence matches previous allele '
                'sequence: gene2_1']
        )
