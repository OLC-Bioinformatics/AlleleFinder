#!/usr/bin/env python

"""
Reduce Enterobase wgMLST profiles with a supplied list of genes
"""

# Standard imports
from argparse import ArgumentParser
from csv import DictReader
import logging
import sys
import os

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import make_path, SetupLogging


class ProfileReduce:

    """
    Reduce Enterobase wgMLST profiles
    """

    def main(self):
        """
        Run the required methods in the appropriate order
        """
        self.read_names()
        self.read_profile()

    def read_names(self):
        """
        Read in all the names of the genes of interest
        """
        logging.info('Reading names file')
        with open(self.name_file, 'r', encoding='utf-8') as names:
            self.names = sorted([name.rstrip() for name in names.readlines()])
        logging.debug('Genes used to reduce profile are: %s', '\t'.join(self.names))

    def read_profile(self):
        """
        Load, parse, and reduce the profile information
        """
        logging.info('Reducing profiles')
        with open(self.profile, 'r', encoding='utf-8') as profile:
            with open(self.reduced_profile, 'w', encoding='utf-8') as reduced:
                with open(self.notes_file, 'w', encoding='utf-8') as notes:
                    # Write the header for the notes field
                    notes.write('OriginalSequenceType\tReducedSequenceType\tNotes\n')
                    # Create the header for the reduced profile
                    gene_names = '\t'.join(self.names)
                    reduced_data = f'ST\t{gene_names}\n'
                    logging.info('Loading profile into memory')
                    # Use to DictReader to load the profile file
                    profile_dict = DictReader(profile, dialect='excel-tab')
                    # Initialise a dictionary to store the string of allele numbers: sequence type
                    profile_st = {}
                    logging.info('Parsing profile')
                    # Iterate through all the sequence types in the profile
                    for allele_dict in profile_dict:
                        # Extract the sequence type from the 'ST' column
                        seq_type = allele_dict['ST']
                        # Initialise variables to store the allele numbering information
                        allele_list = []
                        allele_str = ''
                        # Iterate through all the genes of interest
                        for gene in self.names:
                            allele_str += allele_dict[gene]
                            # Add the allele number to the list, and to the string
                            allele_list.append(allele_dict[gene])
                        # Check if the string of allele numbers is already in the dictionary
                        if allele_str not in profile_st:
                            # Update the dictionary with the string of alleles: sequence type
                            profile_st[allele_str] = allele_dict[self.names[0]]
                            alleles = '\t'.join(allele_list)
                            # Add the sequence type allele numbers for this sequence to the string
                            reduced_data += f'{seq_type}\t{alleles.rstrip()}\n'
                            # Updated the notes with the sequence type linking information
                            notes.write(f'{seq_type}\t{seq_type}\n')
                        # Reduced profile already exists
                        else:
                            # Extract the original sequence type with this string of allele numbers
                            original_seq_type = profile_st[allele_str]
                            # Write the sequence type linking information, making note of the fact
                            # that this is a duplicate
                            notes.write(f'{seq_type}\t{original_seq_type}\tduplicate\n')
                    # Write the reduced profile information to file
                    reduced.write(reduced_data)

    def __init__(self, profile, names, output='profile'):
        logging.info('Welcome to profile reducer!')
        if profile.startswith('~'):
            self.profile = os.path.abspath(os.path.expanduser(os.path.join(profile)))
        else:
            self.profile = os.path.abspath(os.path.join(profile))
        try:
            assert os.path.isfile(self.profile)
        except AssertionError as exc:
            logging.error('Cannot locate the specified profile file: %s, (self.profile,)')
            raise SystemExit from exc
        # Create the folder into which the reduced profile and notes are to be placed
        self.report_path = os.path.join(os.path.dirname(self.profile), output)
        make_path(self.report_path)
        self.reduced_profile = os.path.join(self.report_path, 'profile.txt')
        self.notes_file = os.path.join(self.report_path, 'reducing_notes.txt')
        if names.startswith('~'):
            self.name_file = os.path.abspath(os.path.expanduser(os.path.join(names)))
        else:
            self.name_file = os.path.abspath(os.path.join(names))
        try:
            assert os.path.isfile(self.name_file)
        except AssertionError as exc:
            logging.error('Cannot find the supplied file with gene names: %s', self.name_file)
            raise SystemExit from exc
        self.names = []
        self.profile_dict = {}
        self.allele_dict = {}


def cli():
    """
    Collect the arguments, create an object, and run the script
    """
    # Parser for arguments
    parser = ArgumentParser(description='Extract the genes of interest from a profile file')
    parser.add_argument(
        '-p', '--profile',
        metavar='profile',
        required=True,
        help='Name and path of profile file.')
    parser.add_argument(
        '-n', '--names',
        metavar='names',
        required=True,
        help='Name and path to a file containing the gene names (one per line) to be extracted '
             'from the profile')
    # Get the arguments into an object
    arguments = parser.parse_args()
    SetupLogging(debug=True)
    reduce = ProfileReduce(profile=arguments.profile,
                           names=arguments.names)
    reduce.main()
    logging.info('Profile reduction complete!')
    # Prevent the arguments being printed to the console (they are returned in order for the tests to work)
    sys.stderr = open(os.devnull, 'w', encoding='utf-8')
    return arguments


if __name__ == '__main__':
    cli()
