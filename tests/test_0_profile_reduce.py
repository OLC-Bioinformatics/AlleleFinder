#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/profile_reduce.py
"""

# Standard imports
import argparse
import coloredlogs
import logging
import os
import shutil
import tempfile
from typing import Generator
from unittest import (
    mock,
)
from unittest.mock import (
    patch,
)

# Third-party imports
import pytest

# Local imports
from allele_tools.methods import (
    profile_allele_check,
    setup_logging,
    update_profile_file
)
from allele_tools.profile_reduce import (
    cli,
)


@pytest.fixture(name='variables', scope='module')
def setup():
    """
    Sets up the necessary variables for testing.

    Returns:
    Variables: An instance of the Variables class with all the necessary
    attributes set.
    """

    class Variables:
        """
        A class used to represent the necessary variables for testing.
        """

        def __init__(self):
            """
            Constructs all the necessary attributes for the Variables object.
            """

            # Extract the connection string prior to running tests
            self.test_path = os.path.abspath(os.path.dirname(__file__))

            # Set file path
            self.file_path = os.path.join(
                self.test_path,
                'test_files',
                'profile_reduce'
            )

            # Set profile and names files
            self.profile = os.path.join(self.file_path, 'profiles.list')
            self.names = os.path.join(self.file_path, 'genes.txt')

            # Set output paths
            self.output_path = os.path.join(self.file_path, 'profile')
            self.output_profiles = os.path.join(
                self.file_path,
                'profile',
                'profile.txt'
            )
            self.output_notes = os.path.join(
                self.file_path,
                'profile',
                'reducing_notes.txt'
            )
            self.args = argparse.Namespace(
                aa_alleles=None,
                aa_profile=None
            )
    return Variables()


@patch('argparse.ArgumentParser.parse_args')
def test_profile_reduce_missing_tilde_profile(mock_args, variables):
    """
    Test case for missing profile file with tilde (~) in the path.

    Parameters:
    mock_args (MagicMock): Mocked argument parser
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Mock command line arguments
    mock_args.return_value = argparse.Namespace(
        profile='~/fake_file_profiles.list',  # Non-existent file with tilde
        names=variables.names
    )

    # Expect SystemExit when running the command line interface (cli) due to
    # the missing profile file
    with pytest.raises(SystemExit):
        cli()


@patch('argparse.ArgumentParser.parse_args')
def test_profile_reduce_missing_tilde_names(mock_args, variables):
    """
    Test case for missing names file with tilde (~) in the path.

    Parameters:
    mock_args (MagicMock): Mocked argument parser
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Mock command line arguments
    mock_args.return_value = argparse.Namespace(
        profile=variables.profile,
        names='~/fake_file_genes.txt'  # Non-existent file with tilde
    )

    # Expect SystemExit when running the command line interface (cli) due to
    # the missing names file
    with pytest.raises(SystemExit):
        cli()


@patch('argparse.ArgumentParser.parse_args')
def test_profile_reduce_integration(mock_args, variables):
    """
    Integration test case for the profile reduce functionality.

    Parameters:
    mock_args (MagicMock): Mocked argument parser
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Mock command line arguments
    mock_args.return_value = argparse.Namespace(
        profile=variables.profile,
        names=variables.names
    )

    # Run the command line interface (cli)
    cli()

    # Assert that the output profiles file was created
    assert os.path.isfile(variables.output_profiles)


def test_output_profile_header(variables):
    """
    Test case for validating the header of the output profile.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Open and read the output profiles file
    variables.profiles = open(
        variables.output_profiles,
        'r',
        encoding='utf-8'
    ).readlines()

    # Assert that the first line (header) of the profiles is as expected
    assert variables.profiles[0] == 'ST\tECs1205\tECs1206\n'


def test_output_profile_sequence_type(variables):
    """
    Test case for validating the sequence type in the output profile.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Assert that the second line (sequence type) of the profiles is as
    # expected
    assert variables.profiles[1] == '1\t0\t0\n'


def test_output_notes_header(variables):
    """
    Test case for validating the header of the output notes.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Open and read the output notes file
    variables.notes = open(
        variables.output_notes,
        'r',
        encoding='utf-8'
    ).readlines()

    # Assert that the first line (header) of the notes is as expected
    assert variables.notes[0] == \
        'OriginalSequenceType\tReducedSequenceType\tNotes\n'


def test_output_notes_novel(variables):
    """
    Test case for validating the novel notes in the output notes.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Assert that the second line (novel notes) of the notes is as expected
    assert variables.notes[1] == '1\t1\n'


def test_output_notes_duplicate(variables):
    """
    Test case for validating the duplicate notes in the output notes.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Assert that the third line (duplicate notes) of the notes is as expected
    assert variables.notes[2] == '2\t0\tduplicate\n'


def test_output_notes_length(variables):
    """
    Test case for validating the length of the output notes.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Assert that the length of the notes is as expected
    assert len(variables.notes) == 10


def test_output_clean_up(variables):
    """
    Test case for validating the clean up of the output path.

    Parameters:
    variables (Variables): An instance of the Variables class with all the
    necessary attributes set
    """

    # Remove the output path directory and all its contents
    shutil.rmtree(variables.output_path)

    # Assert that the output path no longer exists
    assert not os.path.exists(variables.output_path)


def test_setup_logging():
    """
    Unit test for setup_logging function.
    """
    # Create a mock ArgumentParser object
    mock_args = mock.Mock(spec=argparse.ArgumentParser)
    mock_args.verbosity = 'debug'

    # Call the function with the mock arguments
    setup_logging(mock_args)

    # Check that the logging level was set correctly
    assert logging.getLogger().getEffectiveLevel() == logging.DEBUG

    # Check that the log format was set correctly
    assert coloredlogs.DEFAULT_LOG_FORMAT == '%(asctime)s %(message)s'

    # Check that the level styles were set correctly
    expected_styles = {
        'debug': {'bold': True, 'color': 'green'},
        'info': {'bold': True, 'color': 'blue'},
        'warning': {'bold': True, 'color': 'yellow'},
        'error': {'bold': True, 'color': 'red'},
        'critical': {'bold': True, 'background': 'red'}
    }
    assert coloredlogs.DEFAULT_LEVEL_STYLES == expected_styles


def test_update_profile_file_gene_in_allele():
    """
    Test 'update_profile_file' function to ensure it correctly checks if the
    gene name is in the allele and extracts the allele number.
    """
    # Define test parameters
    profile_file = "test_profile.txt"
    next_seq_type = 1
    allele_dict = {"gene1": "gene1_allele1", "gene2": "gene2_allele2"}
    genes = ["gene1", "gene2"]
    report_path = "."
    molecule = "aa"

    # Create a test profile file
    with open(profile_file, "w") as file:
        file.write("ST\tgene1\tgene2\n")
        file.write("1\tallele1\tallele2\n")

    # Call the function
    result = update_profile_file(
        profile_file, next_seq_type, allele_dict, genes, report_path, molecule
    )

    # Check the result
    assert result is True

    # Check the content of the profile file
    with open(profile_file, "r") as file:
        lines = file.readlines()
        assert lines[-1] == "1\tallele1\tallele2\n"

    # Clean up
    os.remove(profile_file)
    os.remove(f"{molecule}_novel_profiles.txt")


def test_update_profile_file_remove_metadata():
    """
    Test 'update_profile_file' function to ensure it correctly checks if
    there's metadata in the allele name and removes it.
    """
    # Define test parameters
    profile_file = 'test_profile.txt'
    next_seq_type = 1
    allele_dict = {
        'gene1': 'gene1_allele1|metadata',
        'gene2': 'gene2_allele2|metadata'
    }
    genes = ['gene1', 'gene2']
    report_path = '.'
    molecule = 'aa'

    # Create a test profile file
    with open(profile_file, 'w') as file:
        file.write('ST\tgene1\tgene2\n')
        file.write('1\tallele1\tallele2\n')

    # Call the function
    result = update_profile_file(
        profile_file, next_seq_type, allele_dict, genes, report_path, molecule
    )

    # Check the result
    assert result is True

    # Check the content of the profile file
    with open(profile_file, 'r') as file:
        lines = file.readlines()
        assert lines[-1] == '1\tallele1\tallele2\n'

    # Clean up
    os.remove(profile_file)
    os.remove(f'{molecule}_novel_profiles.txt')


def test_update_profile_file_newline_end():
    """
    Test 'update_profile_file' function to ensure it correctly checks if the
    last line in the file ends with a newline and if not, appends one.
    """
    # Define test parameters
    profile_file = 'test_profile.txt'
    next_seq_type = 1
    allele_dict = {
        'gene1': 'gene1_allele1',
        'gene2': 'gene2_allele2'
    }
    genes = ['gene1', 'gene2']
    report_path = '.'
    molecule = 'aa'

    # Create a test profile file without a newline at the end
    with open(profile_file, 'w') as file:
        file.write('ST\tgene1\tgene2\n')
        file.write('1\tallele1\tallele2')

    # Call the function
    result = update_profile_file(
        profile_file, next_seq_type, allele_dict, genes, report_path, molecule
    )

    # Check the result
    assert result is True

    # Check the content of the profile file
    with open(profile_file, 'r') as file:
        lines = file.readlines()
        assert lines[-1] == '1\tallele1\tallele2\n'

    # Clean up
    os.remove(profile_file)
    os.remove(f'{molecule}_novel_profiles.txt')


def test_update_profile_file_empty_last_line():
    """
    Test 'update_profile_file' function to ensure it correctly checks if the
    last line in the file is an empty line and if so, removes it.
    """
    # Define test parameters
    profile_file = 'test_profile.txt'
    next_seq_type = 1
    allele_dict = {
        'gene1': 'gene1_allele1',
        'gene2': 'gene2_allele2'
    }
    genes = ['gene1', 'gene2']
    report_path = '.'
    molecule = 'aa'

    # Create a test profile file with an empty line at the end
    with open(profile_file, 'w') as file:
        file.write('ST\tgene1\tgene2\n')
        file.write('1\tallele1\tallele2\n')
        file.write('\n')

    # Call the function
    result = update_profile_file(
        profile_file, next_seq_type, allele_dict, genes, report_path, molecule
    )

    # Check the result
    assert result is True

    # Check the content of the profile file
    with open(profile_file, 'r') as file:
        lines = file.readlines()
        assert lines[-1] == '1\tallele1\tallele2\n'

    # Clean up
    os.remove(profile_file)
    os.remove(f'{molecule}_novel_profiles.txt')


@pytest.fixture
def setup_temp_dir() -> Generator[str, None, None]:
    """
    Pytest fixture to create a temporary directory for each test.
    The directory is automatically cleaned up after each test.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir


def test_no_fasta_files(setup_temp_dir: str, variables):
    """
    Test profile_allele_check function when no .fasta files exist in the
    amino acid allele folder.
    """
    # Set the amino acid allele folder to the temporary directory
    variables.args.aa_alleles = setup_temp_dir
    # Set the amino acid profile file to a file in the temporary directory
    variables.args.aa_profile = os.path.join(
        setup_temp_dir,
        'test_aa_profile.txt'
    )
    # Create the amino acid profile file
    with open(variables.args.aa_profile, 'w') as file:
        file.write('test')
    result = profile_allele_check(
        args=variables.args,
        errors=[],
    )
    expected_error = (
        'Could not locate sequence files in supplied amino acid '
        f'allele folder: {setup_temp_dir}'
    )
    assert result == [expected_error]


def test_fasta_exists_no_profile(setup_temp_dir, variables):
    """
    Test profile_allele_check function when .fasta file exists but the
    profile file does not.
    """
    # Set the amino acid allele folder to the temporary directory
    variables.args.aa_alleles = setup_temp_dir

    # Create a .fasta file in the amino acid allele folder
    with open(os.path.join(setup_temp_dir, 'test.fasta'), 'w') as file:
        file.write('>test\nATG')

    # Set the amino acid profile file to a file in the temporary directory
    variables.args.aa_profile = os.path.join(
        setup_temp_dir,
        'test_aa_profile.txt'
    )
    result = profile_allele_check(
        args=variables.args,
        errors=[]
    )
    expected_error = (
        f'Could not locate supplied amino acid profile '
        f'file: {os.path.join(setup_temp_dir, "test_aa_profile.txt")}'
    )
    assert result == [expected_error]


def test_fasta_and_profile_exist(setup_temp_dir, variables):
    """
    Test profile_allele_check function when both .fasta file and profile
    file exist.
    """
    # Set the amino acid allele folder to the temporary directory
    variables.args.aa_alleles = setup_temp_dir
    # Set the amino acid profile file to a file in the temporary directory
    variables.args.aa_profile = os.path.join(
        setup_temp_dir,
        'test_aa_profile.txt'
    )
    # Create a .fasta file in the amino acid allele folder
    with open(os.path.join(setup_temp_dir, 'test.fasta'), 'w') as file:
        file.write('>test\nATG')
    # Create the amino acid profile file
    with open(variables.args.aa_profile, 'w') as file:
        file.write('test')
    result = profile_allele_check(
        args=variables.args,
        errors=[])
    assert result == []
