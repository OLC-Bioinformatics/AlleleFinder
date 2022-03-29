#!/usr/bin/env python
from olctools.accessoryFunctions.accessoryFunctions import make_path, SetupLogging
from argparse import ArgumentParser
from csv import DictReader
import logging
import os


class ProfileReduce(object):

    def main(self):
        self.read_names()
        self.read_profile()

    def read_names(self):
        """
        Read in all the names of the genes of interest
        """
        logging.info('Reading names file')
        with open(self.name_file, 'r') as names:
            self.names = [name.rstrip() for name in names.readlines()]

    def read_profile(self):
        """
        Load, parse, and reduce the profile information
        """
        logging.info('Reducing profiles')
        with open(self.profile, 'r') as profile:
            with open(self.reduced_profile, 'w') as reduced:
                with open(self.notes_file, 'w') as notes:
                    # Write the header for the notes field
                    notes.write('OriginalSequenceType\tReducedSequenceType\tNotes\n')
                    # Create the header for the reduced profile
                    reduced_data = 'ST\t{genes}\n'.format(genes='\t'.join(self.names))
                    logging.info('Loading profile into memory')
                    # Use to DictReader to load the profile file
                    profile_dict = DictReader(profile, dialect='excel-tab')
                    # Initialise a dictionary to store the string of allele numbers: sequence type
                    profile_st = dict()
                    logging.info('Parsing profile')
                    # Iterate through all the sequence types in the profile
                    for allele_dict in profile_dict:
                        # Extract the sequence type from the 'ST' column
                        seq_type = allele_dict['ST']
                        # Initialise variables to store the allele numbering information
                        allele_list = list()
                        allele_str = str()
                        # Iterate through all the genes of interest
                        for gene in self.names:
                            # Add the allele number to the list, and to the string
                            allele_list.append(allele_dict[gene])
                            allele_str += allele_dict[gene]
                        # Check if the string of allele numbers is already in the dictionary
                        if allele_str not in profile_st:
                            # Update the dictionary with the string of alleles: sequence type
                            profile_st[allele_str] = allele_dict[gene]
                            alleles = '\t'.join(allele_list)
                            # Add the sequence type allele numbers for this sequence to the string
                            reduced_data += '{st}\t{alleles}\n'.format(st=seq_type,
                                                                       alleles=alleles.rstrip())
                            # Updated the notes with the sequence type linking information
                            notes.write(f'{seq_type}\t{seq_type}\n')
                        # Reduced profile already exists
                        else:
                            # Extract the original sequence type with this string of allele numbers
                            original_seq_type = profile_st[allele_str]
                            # Write the sequence type linking information, making note of the fact that this is a
                            # duplicate
                            notes.write(f'{seq_type}\t{original_seq_type}\tduplicate\n')
                    # Write the reduced profile information to file
                    reduced.write(reduced_data)

    def __init__(self, profile, names):
        logging.info('Welcome to profile reducer!')
        if profile.startswith('~'):
            self.profile = os.path.abspath(os.path.expanduser(os.path.join(profile)))
        else:
            self.profile = os.path.abspath(os.path.join(profile))
        assert os.path.isfile(self.profile), f'Cannot find the supplied profile {self.profile}'
        self.report_path = os.path.join(os.path.dirname(self.profile), 'reports')
        make_path(self.report_path)
        self.reduced_profile = os.path.join(self.report_path, 'profile.txt')
        self.notes_file = os.path.join(self.report_path, 'reducing_notes.txt')
        if names.startswith('~'):
            self.name_file = os.path.abspath(os.path.expanduser(os.path.join(names)))
        else:
            self.name_file = os.path.abspath(os.path.join(names))
        assert os.path.isfile(self.name_file), f'Cannot find the supplied file with gene names: {self.name_file}'
        self.names = list()
        self.profile_dict = dict()
        self.allele_dict = dict()


def cli():
    # Parser for arguments
    parser = ArgumentParser(description='Extract the genes of interest from a profile file')
    parser.add_argument('-p', '--profile',
                        required=True,
                        help='Name and path of profile file.')
    parser.add_argument('-n', '--names',
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


if __name__ == '__main__':
    cli()
