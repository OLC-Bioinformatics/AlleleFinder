#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import combinetargets, make_path, MetadataObject, SetupLogging
from genemethods.geneseekr.blast import BLAST
import allele_finder
from argparse import ArgumentParser
from csv import DictReader
import multiprocessing
from glob import glob
from time import time
import logging
import shutil
import json
import os
__author__ = 'adamkoziol'


class ProfileAlleles(object):

    def main(self):
        self.run_geneseekr()
        self.profile_alleles()
        self.read_profile()
        self.match_profile()
        self.create_profile()
        self.sequence_type()
        self.append_profiles()

    def run_geneseekr(self):
        """
        Run GeneSeekr BLASTn of the allele file against the genomes
        """
        # Run the necessary methods from the GeneSeekr pipeline
        self.blast.blast_db()
        self.blast.run_blast()
        self.blast.parseable_blast_outputs()
        self.blast.parse_results()
        # If the reports are not present, run report creation method - this can take a bit of time every time the
        # script is called, so if it is not required, don't worry about it
        if not os.path.isfile(os.path.join(self.reportpath, 'geneseekr_blastn_detailed.csv')):
            self.blast.create_reports()
        # Create the list of MetadataObject from the blast object metadata
        self.runmetadata.samples = self.blast.metadata

    def profile_alleles(self):
        """
        Create the gene:allele profile from the BLAST outputs from each sample
        """
        logging.info('Determining allele profiles')
        # Iterate through all the samples
        for sample in self.runmetadata.samples:
            # Initialise a dictionary to store the profile information for each samples
            self.profile_dict[sample.name] = dict()
            # Initialise a dictionary to store the gene:allele combinations for each sample
            allele_comprehension = dict()
            # Each gene in the analysis is stored in the list of genes created in allele_finder
            for gene in self.records:
                # Initialise the variable to track whether the current gene is present in the current sample
                present = False
                # Iterate through all the BLAST outputs for the sample
                for allele in sample.geneseekr.blastresults:
                    # If the gene name is present as a substring of the allele e.g. adk in adk_0, then the gene is
                    # present in the BLAST outputs
                    if gene in allele:
                        # Strip off the allele number from the allele e.g. adk_0 yields 0
                        allele_id = allele.split('_')[-1]
                        # Update the dictionary with the new gene: allele number for the sample
                        allele_comprehension.update({gene: allele_id})
                        # Update the gene presence variable
                        present = True
                # If, after iterating through all the BLAST outputs, the gene is not present in the sample, update the
                # gene: allele to reflect this absence
                if not present:
                    # Set missing alleles to 'ND'
                    allele_comprehension.update({gene: 'ND'})
            # In order to hash the dictionary, use JSON, with sorted keys to freeze it
            frozen_allele_comprehension = json.dumps(allele_comprehension, sort_keys=True)
            # Update the dictionary of profiles with the hash of the frozen dictionary: list of samples with that hash
            if hash(frozen_allele_comprehension) not in self.profile_dict:
                self.profile_dict[hash(frozen_allele_comprehension)] = [sample.name]
            else:
                self.profile_dict[hash(frozen_allele_comprehension)].append(sample.name)
            # Add the 'regular' dictionary to the list of all profiles as required
            if allele_comprehension not in self.profile_set:
                self.profile_set.append(allele_comprehension)

    def read_profile(self):
        """
        Load any previously created profiles for this analysis
        """
        # Only load the profile file if it exists
        if os.path.isfile(self.profile_file):
            logging.info('Extracting profiles from profile file')
            # Open an Excel-formatted sequence profile file as a dictionary
            profile = DictReader(open(self.profile_file), dialect='excel-tab')
            # Iterate through the rows
            for row in profile:
                # Populate the profile dictionary with profile number: {gene: allele}.
                allele_comprehension = {gene: allele for gene, allele in row.items() if gene != 'ST'}
                # Extract the sequence type number from the first field name
                st = row[profile.fieldnames[0]]
                # Update the profile data dictionary
                self.profiledata[st] = allele_comprehension

    def match_profile(self):
        """
        Match current profiles to any previously created profiles
        """
        # If the profiledata dictionary was not populated in the read_profiles methods, there is nothing to match
        if self.profiledata:
            logging.info('Matching new profiles against profile file')
            # Extract the sequence type and allele dictionary from the profile file
            for sequence_type, allele_comprehension in self.profiledata.items():
                # Freeze the allele comprehension as above
                frozen_allele_comprehension = json.dumps(allele_comprehension, sort_keys=True)
                try:
                    # Extract the samples that match this profile
                    matches = self.profile_dict[hash(frozen_allele_comprehension)]
                    # Update the dictionary with the matching samples
                    self.profile_matches[sequence_type] = matches
                # The profile will not necessarily match any of the profiles found in the analysis
                except KeyError:
                    pass

    def create_profile(self):
        """
        Create new profiles for novel profiles as required
        """
        # Initialise the sequence type to be 1
        sequence_type = 1
        # If the profiledata dictionary exists, set the sequence type to be the number of entries in the dictionary
        # plus one, as that corresponds to the next sequence type
        if self.profiledata:
            sequence_type = len(self.profiledata) + 1
        # Initialise a list to store the matched samples
        matched = list()
        # Iterate through all the profiles in the analysis
        for allele_comprehension in self.profile_set:
            # Ensure that the allele comprehension (profile) is not already in the profile file
            if allele_comprehension not in [profile_alleles for st, profile_alleles in self.profiledata.items()]:
                # Add the new profile to the list of new profiles
                self.new_profiles.append('{profile_num}\t{alleles}'
                                         .format(profile_num=sequence_type,
                                                 alleles='\t'.join(allele_num for gene, allele_num in
                                                                   sorted(allele_comprehension.items()))))
                # Freeze the comprehension in order to be used as the key in the profile dictionary
                frozen_allele_comprehension = json.dumps(allele_comprehension, sort_keys=True)
                matches = self.profile_dict[hash(frozen_allele_comprehension)]
                # Check to see if this sequence type hasn't already been found in the current analysis
                if matches not in matched:
                    # Update the dictionary with the new sequence type: list of samples
                    self.profile_matches[sequence_type] = matches
                    self.profiledata[sequence_type] = allele_comprehension
                    # Add the matches to the list of matches
                    matched.append(matches)
                # Increment the sequence type number of the next entry
                sequence_type += 1

    def sequence_type(self):
        """
        Perform the final sequence typing, and create the report
        """
        # Open the report
        with open(self.profile_report, 'w') as report:
            # Initialise the header with an extra 'Sample' column plus the comma-separate list of gene names
            data = 'Sample,' + self.data
            for sample in self.runmetadata.samples:
                # Iterate through all the matches to sequence profiles in the analysis
                for sequence_type, sample_names in self.profile_matches.items():
                    # Check if the sample name is in the list of samples names with the current sequence type
                    if sample.name in sample_names:
                        # Add the sample name, sequence type, and all the allele numbers to the report string
                        data += '{sn},{st},{ac}\n'\
                            .format(sn=sample.name,
                                    st=sequence_type,
                                    ac=','.join(allele_num for gene, allele_num in
                                                sorted(self.profiledata[sequence_type].items())))
            # Write the report
            report.write(data)

    def append_profiles(self):
        """
        Add new profiles to the profile file
        """
        # Only try to add new profiles if there are new profiles in the analysis
        if self.new_profiles:
            # Initialise the string to store the new profile
            data = str()
            # If the profile file does not exist, add the string of 'ST', comma-separated gene names to the headers
            if not os.path.isfile(self.profile_file):
                data = self.data
            # Iterate through all the new profiles, and add them to the new profile string
            for profile in self.new_profiles:
                data += '{profile}\n'.format(profile=profile)
            # Open the report with a+ to either create, or append the profile string to it as appropriate
            with open(self.profile_file, 'a+') as profile_file:
                profile_file.write(data)

    def __init__(self, args, records):
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.args = args
        if args.path.startswith('~'):
            self.path = os.path.abspath(os.path.expanduser(os.path.join(args.path)))
        else:
            self.path = os.path.abspath(os.path.join(args.path))
        if args.local_path.startswith('~'):
            self.local_path = os.path.abspath(os.path.expanduser(os.path.join(args.local_path)))
        else:
            self.local_path = os.path.abspath(os.path.join(args.local_path))
        self.sequencepath = os.path.join(self.path, 'geneseekr')
        make_path(self.sequencepath)
        target_files = sorted(glob(os.path.join(self.local_path, '*.tfa')))
        # Create symlinks of the target files in the local path
        for target in target_files:
            try:
                os.symlink(target, os.path.join(self.sequencepath, os.path.basename(target).replace('.tfa', '.fasta')))
            except FileExistsError:
                pass
        self.targetpath = os.path.join(self.sequencepath, 'targets')
        make_path(self.targetpath)
        self.profilepath = os.path.join(self.path, 'profile')
        make_path(self.profilepath)
        self.profile_file = os.path.join(self.profilepath, 'profile.txt')
        self.allelepath = os.path.join(self.path, 'alleles')
        shutil.copyfile(src=os.path.join(self.allelepath, 'combinedtargets.fasta'),
                        dst=os.path.join(self.targetpath, 'combinedtargets.fasta'))
        self.reportpath = os.path.join(self.sequencepath, 'reports')
        self.profile_report = os.path.join(self.reportpath, 'profiles.csv')
        self.cpus = multiprocessing.cpu_count() - 1
        self.starttime = args.starttime
        self.start = self.starttime
        self.runmetadata = MetadataObject()
        self.records = records
        # Create an object for performing BLAST analyses
        self.blast = BLAST(args=self,
                           analysistype='geneseekr',
                           cutoff=100,
                           unique=True)
        self.profile_dict = dict()
        self.profiledata = dict()
        self.profile_set = list()
        self.sequence_profile = dict()
        self.profile_matches = dict()
        self.new_profiles = list()
        # A string of the header to use for formatting the profile file, and the report headers
        self.data = 'ST\t{genes}\n'.format(genes='\t'.join(sorted(self.records)))


def cli():
    # Import the argument parser from allele_finder.py
    parent_parser = allele_finder.cli()
    parser = ArgumentParser(parents=[parent_parser])
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.starttime = time()
    SetupLogging(debug=arguments.verbose)
    # Run the allele-finding pipeline
    finder = allele_finder.AlleleFinder(arguments)
    finder.main()
    logging.info('Allele finding complete')
    # Extract the dictionary of records from the allele finding
    records = finder.records
    # Run the profiling pipeline
    profiler = ProfileAlleles(args=arguments,
                              records=records)
    profiler.main()
    logging.info('Allele Profiling complete')


if __name__ == '__main__':
    cli()
