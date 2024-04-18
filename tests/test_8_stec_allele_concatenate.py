#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/stec.py allele_concatenate
"""

# Standard imports
import argparse
import json
import os
import shutil
import tempfile
from unittest.mock import (
    Mock,
    patch
)
# Third-party imports
import pytest

# Local imports
from allele_tools.allele_profiler import (
    append_profiles,
    parseable_blast_outputs,
    profile_alleles,
    read_profile,
    sequence_typer
)
from allele_tools.methods import (
    concatenate_alleles,
    load_alleles,
    write_concatenated_sequences
)
from allele_tools.stec import (
    allele_concatenate,
    cli
)


@pytest.fixture(name='variables', scope='session')
def setup():
    """
    Setup fixture for the tests. This fixture is responsible for setting up the
    necessary variables for the tests.

    :return: An instance of the Variables class.
    """
    class Variables:
        def __init__(self):
            # Define the paths and files for the tests
            # Get the absolute path of the directory of this file
            self.test_path = os.path.abspath(os.path.dirname(__file__))
            # Construct the path to the test files
            self.file_path = os.path.join(
                self.test_path,
                'test_files',
                'allele_updater',
                'original_files'
            )
            # Construct the paths to the nt and aa allele files
            self.nt_allele_path = os.path.join(self.file_path, 'nt_alleles')
            self.aa_allele_path = os.path.join(self.file_path, 'aa_alleles')
            # Construct the paths to the nt and aa profile files
            self.nt_profile_file = os.path.join(
                self.file_path,
                'nt_profile',
                'profile.txt'
            )
            self.aa_profile_file = os.path.join(
                self.file_path,
                'aa_profile',
                'profile.txt'
            )
            # Construct the path to the concatenate directory
            self.concatenate_path = os.path.join(
                self.test_path,
                'test_files',
                'concatenate_alleles'
            )
            # Define the linker lengths for the stx1 and stx2 genes
            self.linker_length_dict = {
                'stx1': 9,
                'stx2': 12,
            }
            # Define the order of the alleles for the stx1 and stx2 genes
            self.allele_order = {
                'stx1': ['stx1B', 'stx1B'],
                'stx2': ['stx2A', 'stx2B']
            }
            # Define the gene-allele mapping for the stx1 and stx2 genes
            self.gene_allele = {
                'stx1': 'stx1A_stx1B',
                'stx2': 'stx2A_stx2B'
            }

    return Variables()


def prepare_files(variables):
    """
    Prepare the necessary files for the tests. This function is responsible for
    creating the concatenate directory if it does not exist.

    :param variables: An instance of the Variables class.
    """
    # Check if the concatenate directory exists
    if not os.path.isdir(variables.concatenate_path):
        # If not, create it
        os.makedirs(variables.concatenate_path)
    # Assert that the concatenate directory now exists
    assert os.path.isdir(variables.concatenate_path)


def clear_path(variables):
    """
    Clean up the outputs after the tests. This function is responsible for
    removing the concatenate directory.

    :param variables: An instance of the Variables class.
    """
    # Remove the concatenate directory
    shutil.rmtree(variables.concatenate_path)
    # Assert that the concatenate directory no longer exists
    assert not os.path.isdir(variables.concatenate_path)


def test_read_profile(variables):
    """
    Test the read_profile function. This test checks if the function correctly
    reads the profile file.

    :param variables: An instance of the Variables class.
    """
    # Read the profile file
    variables.nt_profile_data = read_profile(
        profile_file=variables.nt_profile_file
    )
    # Assert that the profile data is correct
    assert variables.nt_profile_data['113'] == {
        'stx2A': '12',
        'stx2B': '1'
    }


def test_load_alleles(variables):
    """
    Test the load_alleles function. This test checks if the function correctly
    loads the alleles.

    :param variables: An instance of the Variables class.
    """
    # Load the alleles
    variables.gene, variables.nt_alleles = load_alleles(
        allele_path=variables.nt_allele_path,
        allele_order=variables.allele_order
    )
    # Assert that the gene is correct
    assert variables.gene == 'stx2'
    # Assert that the nt alleles are correct
    assert variables.nt_alleles['stx2A']['stx2A_1'][:50] == \
        'ATGAAGTGTATATTATTTAAATGGGTACTGTGCCTGTTACTGGGTTTTTC'


def test_concatenate_alleles(variables):
    """
    Test the concatenate_alleles function. This test checks if the function
    correctly concatenates the alleles.

    :param variables: An instance of the Variables class.
    """
    # Concatenate the alleles
    variables.concatenated_nt_seq = concatenate_alleles(
        profile_data=variables.nt_profile_data,
        allele_dict=variables.nt_alleles,
        allele_order=variables.allele_order,
        stx_gene=variables.gene,
        linker_length_dict=variables.linker_length_dict,
        molecule='nt'
    )
    # Assert that the concatenated sequence is correct
    assert str(variables.concatenated_nt_seq[0].seq)[:50] == \
        'ATGAAGTGTATATTATTTAAATGGGTACTGTGCCTGTTACTGGGTTTTTC'
    assert str(variables.concatenated_nt_seq[0].seq)[-50:] == \
        'CCTGTGAATCAGGCTCCGGATTTGCTGAAGTGCAGTTTAATAATGACTGA'


def test_write_concatenated_alleles(variables):
    """
    Test the write_concatenated_sequences function. This test checks if the
    function correctly writes the concatenated sequences to a file.

    :param variables: An instance of the Variables class.
    """
    # Prepare the files for the test
    prepare_files(variables=variables)
    # Write the concatenated sequences to a file
    write_concatenated_sequences(
        concatenated_sequences=variables.concatenated_nt_seq,
        concatenate_path=variables.concatenate_path,
        file_name=variables.gene_allele[variables.gene],
        molecule='nt'
    )
    # Assert that the file was correctly written
    assert os.path.isfile(os.path.join(
        variables.concatenate_path,
        'nt',
        'stx2A_stx2B.fasta'
    ))
    # Clear the path after the test
    clear_path(variables=variables)


@patch('argparse.ArgumentParser.parse_args')
def test_allele_concatenate_integration(mock_args, variables):
    """
    Test the allele_concatenate function. This test checks if the function
    correctly concatenates the alleles and writes the concatenated sequences
    to a file.

    :param mock_args: A mock object for the command line arguments.
    :param variables: An instance of the Variables class.
    """
    # Prepare the files for the test
    prepare_files(variables=variables)
    # Mock the command line arguments
    mock_args.return_value = argparse.Namespace(
        nt_profile=variables.nt_profile_file,
        aa_profile=variables.aa_profile_file,
        nt_alleles=variables.nt_allele_path,
        aa_alleles=variables.aa_allele_path,
        concatenate_path=variables.concatenate_path
    )
    # Get the command line arguments
    arguments = cli()
    # Run the allele_concatenate function
    allele_concatenate(args=arguments)
    # Assert that the file was correctly written
    assert os.path.isfile(os.path.join(
        variables.concatenate_path,
        'aa',
        'stx2A_stx2B.fasta'
        ))


def test_fieldnames_different_length(variables):
    """
    Test when header and fieldnames have different lengths.
    """
    # Define the path to blast_report directly
    blast_report = os.path.join(
        variables.file_path,
        'blast_report.txt'
    )
    with open(blast_report, 'w') as f:
        # The fields are 'gaps', 'identical', and 'subject_id'
        f.write('0\t100\tsubject1\n10\t90\tsubject2\n')

    # Mock the alleles object and set its blast_report attribute
    alleles = Mock()
    alleles.blast_report = blast_report
    # Mock the sample object and set its alleles attribute
    sample = Mock()
    sample.alleles = alleles
    # Mock the runmetadata object and set its samples attribute
    runmetadata = Mock()
    runmetadata.samples = [sample]
    # Define the fieldnames and extended_fieldnames
    fieldnames = ['gaps', 'identical']
    extended_fieldnames = ['gaps', 'identical', 'subject_id']
    cutoff = 95

    # Create the records dictionary with 'subject1' and 'subject2' as keys
    # and lists of length 100 as values
    records = {'subject1': [0]*100, 'subject2': [0]*100}

    parseable_blast_outputs(
        runmetadata=runmetadata,
        fieldnames=fieldnames,
        extended_fieldnames=extended_fieldnames,
        records=records,
        cutoff=cutoff
    )

    # Check if the function modified the header in the blast_report file
    with open(blast_report, 'r') as f:
        header = f.readline().strip()
    assert header == '\t'.join(extended_fieldnames)
    # Clean up the BLAST report
    os.remove(blast_report)


def test_profile_alleles(variables):
    # Mock the MetadataObject
    runmetadata = Mock()
    # Mock the sample object and set its name and alleles attributes
    sample = Mock()
    sample.name = 'sample1'
    sample.alleles = Mock()

    sample.alleles.blastresults = []
    sample.alleles.targetsequence = {}

    # Define the path to blast_report directly
    blast_report = os.path.join(
        variables.file_path,
        'blast_report.txt'
    )
    # Mock the alleles object and set its blast_report attribute
    sample.alleles.blast_report = blast_report
    with open(blast_report, 'w') as f:
        # The fields are 'gaps', 'identical', 'query_id', and 'subject_id'
        f.write(
            'gaps\tidentical\tquery_id\tsubject_id\n'
            '0\t100\tsubject1\n10\t90\tsubject2\n'
        )

    runmetadata.samples = [sample]
    # Set the profile_dict, profile_set, and records
    profile_dict = {}
    profile_set = []
    records = ['gene3']
    # Call the function with the mock objects and parameters
    profile_dict, profile_set = profile_alleles(
        runmetadata=runmetadata,
        profile_dict=profile_dict,
        profile_set=profile_set,
        records=records,
        novel_alleles=False,
        genome_query=False,
        amino_acid=False,
        allele_path=None,
        report_path=None,
        cutoff=75
    )
    # Check if the function correctly updated profile_set
    assert profile_set == [{'gene3': '0'}]

    # Clean up the BLAST report
    os.remove(blast_report)


def test_profile_alleles_new_hash():
    # Mock the MetadataObject
    runmetadata = Mock()
    # Mock the sample object and set its name and alleles attributes
    sample = Mock()
    sample.name = 'sample1'
    sample.alleles = Mock()
    sample.alleles.blastresults = ['gene1_1', 'gene2_1']
    sample.alleles.targetsequence = {}
    runmetadata.samples = [sample]
    # Set the profile_dict, profile_set, and records
    profile_dict = {}
    profile_set = []
    records = ['gene1', 'gene2']
    # Call the function with the mock objects and parameters
    profile_dict, profile_set = profile_alleles(
        runmetadata=runmetadata,
        profile_dict=profile_dict,
        profile_set=profile_set,
        records=records,
        novel_alleles=False,
        genome_query=False,
        amino_acid=False,
        allele_path=None,
        report_path=None,
        cutoff=75
    )
    # Check if the function correctly updated the profile_dict and profile_set
    expected_hash = hash(
        json.dumps(
            {
                'gene1': '1',
                'gene2': '1'
            },
            sort_keys=True
        )
    )

    assert profile_dict[expected_hash] == ['sample1']
    assert profile_set == [{'gene1': '1', 'gene2': '1'}]


def test_profile_alleles_existing_hash():
    """
    Test the functionality of profile_alleles function when the hash
    already exists.

    This test creates a mock MetadataObject and a temporary file to use as the
    profile_report. It then calls the profile_alleles function with these mock
    objects and checks if the function correctly updates the sample's alleles
    and writes to the report when the hash already exists.

    :return: None
    """
    # Mock the MetadataObject
    runmetadata = Mock()
    # Mock the sample objects and set their name and alleles attributes
    sample1 = Mock()
    sample1.name = 'sample1'
    sample1.alleles = Mock()
    sample1.alleles.blastresults = ['gene1_1', 'gene2_1']
    sample1.alleles.targetsequence = {}
    sample2 = Mock()
    sample2.name = 'sample2'
    sample2.alleles = Mock()
    sample2.alleles.blastresults = ['gene1_1', 'gene2_1']
    sample2.alleles.targetsequence = {}
    runmetadata.samples = [sample1, sample2]
    # Set the profile_dict, profile_set, and records
    profile_dict = {}
    profile_set = []
    records = ['gene1', 'gene2']
    # Call the function with the mock objects and parameters
    profile_dict, profile_set = profile_alleles(
        runmetadata=runmetadata,
        profile_dict=profile_dict,
        profile_set=profile_set,
        records=records,
        novel_alleles=False,
        genome_query=False,
        amino_acid=False,
        allele_path=None,
        report_path=None,
        cutoff=75
    )
    # Check if the function correctly updated the profile_dict and profile_set
    expected_hash = hash(
        json.dumps(
            {
                'gene1': '1',
                'gene2': '1'
            },
            sort_keys=True
        )
    )
    assert profile_dict[expected_hash] == ['sample1', 'sample2']
    assert profile_set == [{'gene1': '1', 'gene2': '1'}]


def test_sequence_typer_no_update():
    """
    Test the sequence_typer function when the update parameter is False.

    This test creates a mock MetadataObject and a temporary file to use as
    the profile_report. It then calls the sequence_typer function with these
    mock objects and checks if the function correctly updates the sample's
    alleles and writes to the report.

    :return: None
    """
    # Mock the MetadataObject
    runmetadata = Mock()
    # Mock the sample object and set its name and alleles attributes
    sample = Mock()
    sample.name = 'sample1'
    sample.alleles = Mock()
    runmetadata.samples = [sample]
    # Set the profile_matches and profile_data
    profile_matches = {'seq_type1': ['sample1']}
    profile_data = {'seq_type1': {'gene1': 'allele1_1', 'gene2': 'allele2_1'}}
    # Create a temporary file to use as the profile_report
    with tempfile.NamedTemporaryFile(delete=False) as temp:
        profile_report = temp.name
    # Call the function with the mock objects and parameters
    runmetadata = sequence_typer(
        profile_report=profile_report,
        data='gene1\tgene2\n',
        runmetadata=runmetadata,
        profile_matches=profile_matches,
        profile_data=profile_data,
        update=False,
        amino_acid=False
    )
    # Check if the function correctly updated the sample's alleles
    assert sample.alleles.nt_st == 'seq_type1'
    assert sample.alleles.nt_profile == profile_data['seq_type1']
    # Check if the function correctly wrote to the report
    with open(profile_report, 'r') as report:
        assert report.read() == \
            'Sample\tgene1\tgene2\nsample1\tseq_type1\t1\t1\n'
    # Delete the temporary file
    os.remove(profile_report)


def test_append_profiles_no_profile_file():
    """
    Test the append_profiles function when the profile_file does not exist.

    This test creates a temporary directory and a non-existent profile_file
    within it. It then calls the append_profiles function with this
    profile_file and checks if the function correctly creates the
    profile_file and writes the data to it.

    :return: None
    """
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Set the profile_file to a non-existent file within the
        # temporary directory
        profile_file = os.path.join(temp_dir, 'non_existent_file.txt')
        # Set the new_profiles and data
        new_profiles = ['profile1', 'profile2']
        data = 'ST\tgene1\tgene2\n'
        # Call the function with the mock objects and parameters
        append_profiles(
            new_profiles=new_profiles,
            profile_file=profile_file,
            data=data,
            novel_profiles=False,
            profile_path=None,
            gene_names=None
        )
        # Check if the function correctly created the profile_file and wrote
        # the data to it
        with open(profile_file, 'r') as file:
            assert file.read() == 'ST\tgene1\tgene2\nprofile1\nprofile2\n'


def test_clear(variables):
    # Clear the path after the test
    clear_path(variables=variables)
