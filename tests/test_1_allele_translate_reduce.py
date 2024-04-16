#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/allele_translate_reduce.py
"""

# Standard imports
import argparse
from glob import glob
import os
import shutil
import unittest
from unittest.mock import (
    Mock,
    patch,
)
# Third-party imports
import pytest

# Local imports
from allele_tools.allele_translate_reduce import (
    cli,
    Translate
)

from allele_tools.methods import (
    create_nt_allele_comprehension
)


@pytest.fixture(name='variables', scope='module')
def setup():
    """
    Pytest fixture for setting up the test environment.
    This fixture is named 'variables' and has module scope.
    """
    class Variables:
        def __init__(self):
            # Extract the path of the tests directory
            self.test_path = os.path.abspath(os.path.dirname(__file__))

            # Define file paths for testing
            self.file_path = os.path.join(
                self.test_path, 'test_files', 'allele_translate_reduce'
            )
            self.nt_profile_path = os.path.join(
                self.file_path, 'nt_profile'
            )
            self.nt_profile = os.path.join(
                self.nt_profile_path, 'profile.txt'
            )
            self.nt_profile_file = os.path.join(
                self.file_path, 'original_files', 'profile.txt'
            )
            self.nt_allele_files = glob(os.path.join(
                self.file_path, 'original_files', 'stx2*'
            ))
            self.nt_allele_path = os.path.join(
                self.file_path, 'nt_alleles'
            )
            self.aa_allele_path = os.path.join(
                self.file_path, 'aa_alleles'
            )
            self.aa_profile_path = os.path.join(
                self.file_path, 'aa_profile'
            )
            self.aa_profile_file = os.path.join(
                self.aa_profile_path, 'profile.txt'
            )

            # Define length dictionary for testing
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
    Function to clean up the output directories and files.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # List of directories to clean
    directories = [
        variables.aa_allele_path,
        variables.aa_profile_path,
        variables.nt_profile_path,
        variables.nt_allele_path
    ]

    # Remove directories if they exist
    for directory in directories:
        if os.path.isdir(directory):
            shutil.rmtree(directory)

    # List of file types to clean
    file_types = ['*.fasta', '*.txt']

    # Remove files of specified types if they exist
    for file_type in file_types:
        files = glob(os.path.join(variables.file_path, file_type))
        for file in files:
            os.remove(file)

    # Assert that the directories were successfully removed
    for directory in directories:
        assert not os.path.isdir(directory)


def prepare_files(variables):
    """
    Function to prepare the necessary files for testing.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Create the nt_profile_path directory if it doesn't exist
    if not os.path.isdir(variables.nt_profile_path):
        os.makedirs(variables.nt_profile_path)

    # Copy the nt_profile_file to the nt_profile_path directory
    shutil.copyfile(
        src=variables.nt_profile_file,
        dst=os.path.join(variables.nt_profile_path, 'profile.txt')
    )

    # Create the nt_allele_path directory if it doesn't exist
    if not os.path.isdir(variables.nt_allele_path):
        os.makedirs(variables.nt_allele_path)

    # Copy each nt_allele_file to the file_path directory
    for allele_file in variables.nt_allele_files:
        file_name = os.path.basename(allele_file)
        shutil.copyfile(
            src=allele_file,
            dst=os.path.join(variables.file_path, file_name)
        )
        assert os.path.isfile(os.path.join(variables.file_path, file_name))


@patch('argparse.ArgumentParser.parse_args')
def test_allele_translate_reduce_integration(mock_args, variables):
    """
    Integration test for the allele_translate_reduce function.

    Parameters:
    mock_args (MagicMock): Mock for the argparse.ArgumentParser.parse_args
    method
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Prepare the necessary files for testing
    prepare_files(variables=variables)

    # Mock the command line arguments
    mock_args.return_value = argparse.Namespace(
        path=variables.file_path,
        profile=True,
        report_path=variables.aa_profile_path,
        translated_path=variables.aa_allele_path
    )

    # Call the command line interface function
    variables.arguments = cli()

    # Assert that the aa_profile_file was created
    assert os.path.isfile(variables.aa_profile_file)


def test_profile_header(variables):
    """
    Test case for validating the header of the profile.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Open and read the profile file
    variables.profiles = open(
        variables.aa_profile_file,
        'r',
        encoding='utf-8').readlines()

    # Assert that the first line (header) of the profile is as expected
    assert variables.profiles[0] == 'ST\tstx2A\tstx2B\n'


def test_profile_contents(variables):
    """
    Test case for validating the contents of the profile.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Assert that the second line (contents) of the profile is as expected
    assert variables.profiles[1] == '1\t0\t0\n'


def test_profile_end(variables):
    """
    Test case for validating the end of the profile.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Assert that the last line (end) of the profile is as expected
    assert variables.profiles[-1] == '259155\t30\t13\n'

    # Clean up the output directories and files
    clean_outputs(variables=variables)


def test_allele_translate_reduce_length_dict(variables):
    """
    Test case for the Translate class with a length dictionary.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Prepare the necessary files for testing
    prepare_files(variables=variables)

    # Create an instance of the Translate class
    translate = Translate(
        path=variables.file_path,
        profile=True,
        report_path=variables.aa_profile_path,
        translated_path=variables.aa_allele_path,
        length_dict=variables.length_dict
    )

    # Call the main method of the Translate class
    translate.main()

    # Assert that the profile file was created
    assert os.path.isfile(os.path.join(
        variables.aa_profile_path,
        'profile.txt'
    ))

    # Clean up the output directories and files
    clean_outputs(variables=variables)


def test_allele_translate_reduce_tilde_path(variables):
    """
    Test case for the Translate class with a non-existent path.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Assert that a SystemExit exception is raised when trying to create an
    # instance of the Translate class with a non-existent path
    with pytest.raises(SystemExit):
        Translate(
            path='~/completely_fake_path',
            profile=True,
            report_path=variables.aa_profile_path,
            translated_path=variables.aa_allele_path,
            length_dict=variables.length_dict
        )


def test_allele_translate_reduce_profile_provided(variables):
    """
    Test case for the Translate class with a provided profile.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Assert that a SystemExit exception is raised when trying to create an
    # instance of the Translate class with a provided profile
    with pytest.raises(SystemExit):
        Translate(
            path='~/completely_fake_path',
            profile=variables.nt_profile_file,
            report_path=variables.aa_profile_path,
            translated_path=variables.aa_allele_path,
            length_dict=variables.length_dict
        )


def test_allele_translate_reduce_tilde_report_path(variables):
    """
    Test case for the Translate class with a tilde in the report path.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Prepare the necessary files for testing
    prepare_files(variables=variables)

    # Create an instance of the Translate class
    translate = Translate(
        path=variables.file_path,
        profile=True,
        report_path=variables.fake_path,
        translated_path=variables.aa_allele_path,
        length_dict=variables.length_dict
    )

    # Assert that the report path directory exists
    assert os.path.isdir(translate.report_path)

    # Remove the report path directory
    shutil.rmtree(translate.report_path)

    # Clean up the output directories and files
    clean_outputs(variables=variables)


def test_allele_translate_reduce_no_profile(variables):
    """
    Test case for the Translate class with no profile.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Prepare the necessary files for testing
    prepare_files(variables=variables)

    # Create an instance of the Translate class
    translate = Translate(
        path=variables.file_path,
        profile=False,
        report_path=variables.aa_profile_path,
        translated_path=variables.fake_path,
        length_dict=variables.length_dict
    )

    # Assert that the translated path directory exists
    assert os.path.isdir(translate.translated_path)

    # Remove the translated path directory
    shutil.rmtree(translate.translated_path)

    # Clean up the output directories and files
    clean_outputs(variables=variables)


def test_allele_translate_reduce_missing_alleles(variables):
    """
    Test case for the Translate class with missing alleles.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Create the nt_profile_path directory if it does not exist
    if not os.path.isdir(variables.nt_profile_path):
        os.makedirs(variables.nt_profile_path)

    # Copy the nt_profile_file to the nt_profile_path directory
    shutil.copyfile(
        src=variables.nt_profile_file,
        dst=os.path.join(variables.file_path, 'nt_profile', 'profile.txt')
    )

    # Assert that a SystemExit exception is raised when trying to create an
    # instance of the Translate class with missing alleles
    with pytest.raises(SystemExit):
        Translate(
            path=variables.file_path,
            profile=True,
            report_path=variables.aa_profile_path,
            translated_path=variables.aa_allele_path,
            length_dict=variables.length_dict
        )

    # Clean up the output directories and files
    clean_outputs(variables=variables)


@patch('builtins.open', new_callable=unittest.mock.mock_open, read_data='data')
def test_aa_profile_key_error(mock_open, variables):
    """
    Test the aa_profile method when a KeyError is raised.

    This test creates a mock self object with a profile_data dictionary that
    contains a gene-allele pair that is not in the allele_links dictionary.
    It then calls the aa_profile method and checks if the method correctly
    handles the KeyError by adding the sequence type to the filtered_list and
    updating the filtered_dict with the filtered geneName_alleleIdentifier.
    """
    # Mock the self object and its attributes
    self_obj = Mock()
    self_obj.profile_data = {'1': {'gene1': 'allele1'}}
    self_obj.allele_links = {'gene1': {'allele2': 'aa_allele1'}}
    self_obj.aa_profile_file = 'aa_profile_file.txt'
    self_obj.gene_names = {'gene1'}
    self_obj.gene_name_file = 'gene_name_file.txt'
    self_obj.profile_file = os.path.join(
        variables.file_path,
        'profile_file.txt'
    )
    self_obj.filtered_list = set()

    # Call the aa_profile method
    Translate.aa_profile(self_obj)

    # Check if the function correctly handled the KeyError
    assert self_obj.profile_data == {}


def test_corresponding_gene_update():
    """
    Test the 'create_nt_allele_comprehension' function to ensure that it
    correctly updates the 'allele_comprehension' dictionary with the new
    gene: allele number for the sample when the gene name is not in the
    corresponding allele.
    """
    # Mock the MetadataObject and its attributes
    runmetadata = Mock()
    runmetadata.samples = [Mock()]
    runmetadata.samples[0].alleles = Mock()
    runmetadata.samples[0].alleles.overlap_dict = {
        'contig1': {
            'range1': {
                ('gene1', 'gene2'): {
                    'allele': ('gene3_allele1', 'gene2_allele2')
                }
            }
        }
    }

    # List of all gene names in the analysis
    gene_names = ['gene1', 'gene2']

    # Call the function
    result = create_nt_allele_comprehension(runmetadata, gene_names)

    # Check the result
    expected_result = {
        'contig1': {
            'range1': {
                'gene2': 'allele2'
            }
        }
    }
    assert result == expected_result


def test_allele_comprehension_update():
    """
    Test 'create_nt_allele_comprehension' function to ensure it correctly
    updates the 'allele_comprehension' dictionary with the negative result
    for each gene when the dictionary is empty.
    """
    # Mock the MetadataObject and its attributes
    runmetadata = Mock()
    runmetadata.samples = [Mock()]
    runmetadata.samples[0].alleles = Mock()
    runmetadata.samples[0].alleles.translated_overlap = {}
    runmetadata.samples[0].alleles.targetsequence = [
        'contig1'
    ]
    runmetadata.samples[0].alleles.overlap_dict = {}

    # List of all gene names in the analysis
    gene_names = ['gene1', 'gene2']

    # Call the function
    result = create_nt_allele_comprehension(
        runmetadata,
        gene_names,
        translated=True
    )

    # Check the result
    expected_result = {
        'contig1': {
            (0, 0): {
                'gene1': '0',
                'gene2': '0'
            }
        }
    }
    assert result == expected_result
