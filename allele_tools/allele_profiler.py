#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import GenObject, make_path, MetadataObject, run_subprocess, \
    SetupLogging
import allele_finder
from Bio.Blast.Applications import NcbiblastnCommandline, NcbitblastnCommandline
from Bio.Application import ApplicationError
from Bio.Seq import Seq
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
        self.blast()
        self.profile_alleles()
        self.read_profile()
        self.match_profile()
        self.create_profile()
        self.sequence_type()
        self.append_profiles()

    def blast(self):
        """
        BLAST the allele files against the genomes
        """
        self.blast_db()
        self.run_blast()
        self.parseable_blast_outputs()
        self.parse_results()

    def blast_db(self):
        """
        Create the BLAST database of the genomes
        """
        logging.info('Creating BLAST databases from sequence files')
        for query_file in self.query_files:
            # Remove the file extension from the file name
            output = os.path.splitext(query_file)[0]
            cmd = 'makeblastdb -in {fasta} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {output}' \
                .format(fasta=query_file,
                        output=output)
            # Check if database already exists
            if not os.path.isfile('{output}.nhr'.format(output=output)):
                run_subprocess(cmd)

    def run_blast(self):
        """
        BLAST the alleles against the genomes
        """
        logging.info('BLASTing alleles against sequence files')
        for query_file in self.query_files:
            # Create a metadata object to store all the sample-specific information
            sample = MetadataObject()
            sample.alleles = GenObject()
            local_db = os.path.splitext(query_file)[0]
            sample.name = os.path.basename(local_db)
            # Set the name of the BLAST output file
            sample.alleles.blast_report = os.path.join(self.reportpath, '{seq_id}.tsv'.format(seq_id=sample.name))
            # Update the list of metadata objects with this sample
            self.runmetadata.samples.append(sample)
            self.blast_reports.append(sample.alleles.blast_report)
            # Run the appropriate BLAST command: BLASTn for nt; tBLASTn for aa against translated nt
            if self.amino_acid:
                blast = NcbitblastnCommandline(db=local_db,
                                               query=self.target_file,
                                               num_alignments=100000000,
                                               evalue=0.001,
                                               num_threads=self.cpus,
                                               task='tblastn',
                                               outfmt=self.outfmt,
                                               word_size=3,
                                               out=sample.alleles.blast_report)
            else:
                blast = NcbiblastnCommandline(db=local_db,
                                              query=self.target_file,
                                              num_alignments=100000000,
                                              evalue=0.001,
                                              num_threads=self.cpus,
                                              task='blastn',
                                              outfmt=self.outfmt,
                                              out=sample.alleles.blast_report)
            if not os.path.isfile(sample.alleles.blast_report):
                # Run BLAST - supply the record sequence as stdin, so BLAST doesn't look for an input file
                try:
                    blast()
                # BLAST can have issues with genomes that have very large contigs. Retry the analysis using only one
                # thread
                except ApplicationError:
                    os.remove(sample.alleles.blast_report)
                    blast = NcbitblastnCommandline(db=local_db,
                                                   query=self.target_file,
                                                   num_alignments=100000000,
                                                   evalue=0.001,
                                                   num_threads=1,
                                                   task='tblastn',
                                                   outfmt=self.outfmt,
                                                   word_size=3,
                                                   out=sample.alleles.blast_report)
                    blast()

    def parseable_blast_outputs(self):
        """
        Add a header to the BLAST report, so that it is easier to figure out what is in each column
        """
        logging.info('Adding headers to BLAST outputs')
        for sample in self.runmetadata.samples:
            data = list()
            # Load the first line of the report
            with open(sample.alleles.blast_report, 'r') as report:
                header_line = report.readline().strip()
            # Split the header on tabs
            header_list = header_line.split('\t')
            # Check to see if the header has already been added. Skip this step if it has been added.
            if header_list[0] != self.fieldnames[0]:
                with open(sample.alleles.blast_report, 'r') as report:
                    header = [entry for entry in report.readline().split('\t')]
                if len(header) == 15:
                    current_fieldnames = self.fieldnames[:13] + self.fieldnames[14:]
                else:
                    current_fieldnames = self.fieldnames
                blastdict = DictReader(open(sample.alleles.blast_report), fieldnames=current_fieldnames,
                                       dialect='excel-tab')
                # Go through each BLAST result
                for row in blastdict:
                    # Create the query length variable - if the sequences are DNA (blastn), use the query
                    # length as usual; if the sequences are protein (e.g. tblastx), use the query length / 3
                    query_length = float(row['query_length'])
                    # Calculate the percent identity and extract the bitscore from the row
                    # Percent identity is the (length of the alignment - num mismatches) / total query length
                    percentidentity = float('{:0.2f}'.format((float(row['identical']) - float(row['gaps'])) /
                                                             query_length * 100))
                    # Create a percent match entry based on the calculated percent identity match
                    row['percent_match'] = str(percentidentity)
                    # Add the updated row to the list
                    data.append(row)

                # Overwrite the original BLAST outputs to include headers, and the percent match
                with open(sample.alleles.blast_report, 'w') as updated_report:
                    # Add the header
                    updated_report.write('{headers}\n'.format(headers='\t'.join(self.fieldnames)))
                    # Add the results
                    for row in data:
                        for header in self.extended_fieldnames:
                            # Write the value from the row with the header as the key
                            try:
                                updated_report.write('{value}\t'.format(value=row[header]))
                            except KeyError:
                                # noinspection PyTypeChecker
                                updated_report.write('{value}\t'.format(value=''.join(row[None])))
                        # Add a newline for each result
                        updated_report.write('\n')

    def parse_results(self):
        """
        Parse the BLAST results, and populate GenObjects
        """
        logging.info('Parsing BLAST outputs')
        for sample in self.runmetadata.samples:
            # Initialise GenObjects as required
            sample.alleles.blastlist = list()
            sample.alleles.targetsequence = dict()
            # Open the sequence profile file as a dictionary
            blastdict = DictReader(open(sample.alleles.blast_report),
                                   fieldnames=self.extended_fieldnames,
                                   dialect='excel-tab')
            resultdict = dict()
            # Go through each BLAST result
            for row in blastdict:
                # Ignore the headers
                if row['query_id'].startswith(self.fieldnames[0]):
                    pass
                else:
                    # Remove unwanted pipes added to the name
                    target = row['query_id'].lstrip('gb|').rstrip('|') if '|' in row['query_id'] else \
                        row['query_id']
                    row['query_id'] = row['query_id'].lstrip('gb|').rstrip('|') if '|' in row['query_id'] \
                        else row['query_id']
                    # If the percent identity is greater than the cutoff
                    if float(row['percent_match']) == 100:
                        # Append the hit dictionary to the list
                        sample.alleles.blastlist.append(row)
                        # Update the dictionary with the target and percent identity
                        resultdict.update({target: row['percent_match']})
                        # Determine if the orientation of the sequence is reversed compared to the reference
                        if int(row['query_end']) < int(row['query_start']) and not self.amino_acid:
                            # Create a sequence object using Biopython
                            seq = Seq(row['query_sequence'])
                            # Calculate the reverse complement of the sequence
                            querysequence = str(seq.reverse_complement())
                        # If the sequence is not reversed, use the sequence as it is in the output
                        else:
                            querysequence = row['query_sequence']
                        # Add the sequence in the correct orientation to the sample
                        try:
                            sample.alleles.targetsequence[target].append(querysequence)
                        except (AttributeError, KeyError):
                            sample.alleles.targetsequence[target] = list()
                            sample.alleles.targetsequence[target].append(querysequence)
                # Add the percent identity to the object
                sample.alleles.blastresults = resultdict
            # Populate missing results with 'NA' values
            if len(resultdict) == 0:
                sample.alleles.blastresults = 'NA'

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
                for allele in sample.alleles.blastresults:
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
            data = 'Sample\t' + self.data
            for sample in self.runmetadata.samples:
                # Iterate through all the matches to sequence profiles in the analysis
                for sequence_type, sample_names in self.profile_matches.items():
                    # Check if the sample name is in the list of samples names with the current sequence type
                    if sample.name in sample_names:
                        # Add the sample name, sequence type, and all the allele numbers to the report string
                        data += '{sn}\t{st}\t{ac}\n' \
                            .format(sn=sample.name,
                                    st=sequence_type,
                                    ac='\t'.join(allele_num for gene, allele_num in
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

    def __init__(self, path, fasta_path, records, amino_acid):
        if path.startswith('~'):
            self.path = os.path.abspath(os.path.expanduser(os.path.join(path)))
        else:
            self.path = os.path.abspath(os.path.join(path))
        if fasta_path.startswith('~'):
            self.fasta_path = os.path.abspath(os.path.expanduser(os.path.join(fasta_path)))
        else:
            self.fasta_path = os.path.abspath(os.path.join(fasta_path))
        self.working_path = os.path.join(self.path, 'strain_profiles')
        self.sequencepath = os.path.join(self.working_path, 'query')
        make_path(self.sequencepath)
        target_files = [fasta for fasta in sorted(glob(os.path.join(self.fasta_path, '*.fasta')))
                        if os.path.basename(fasta) != 'combinedtargets.fasta']
        self.query_files = list()
        # Create symlinks of the target files in the local path
        for target in target_files:
            try:
                query_file = os.path.join(self.sequencepath, os.path.basename(target).replace('.tfa', '.fasta'))
                self.query_files.append(query_file)
                os.symlink(target, query_file)
            except FileExistsError:
                pass
        self.targetpath = os.path.join(self.working_path, 'targets')
        make_path(self.targetpath)
        self.profilepath = os.path.join(self.working_path, 'sequence_profile')
        make_path(self.profilepath)
        self.profile_file = os.path.join(self.profilepath, 'profile.txt')
        self.target_file = os.path.join(self.targetpath, 'combinedtargets.fasta')
        shutil.copyfile(src=os.path.join(os.path.join(self.path, 'alleles'), 'combinedtargets.fasta'),
                        dst=self.target_file)
        self.reportpath = os.path.join(self.working_path, 'reports')
        make_path(self.reportpath)
        self.strain_profile_path = os.path.join(self.working_path, 'strain_profiles')
        make_path(self.strain_profile_path)
        self.profile_report = os.path.join(self.strain_profile_path, 'profiles.tsv')
        self.cpus = multiprocessing.cpu_count() - 1
        self.starttime = time()
        self.start = self.starttime
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = list()
        self.records = records
        # Create an object for performing BLAST analyses
        if amino_acid:
            self.amino_acid = amino_acid
        else:
            self.amino_acid = None
        if amino_acid:
            self.program = 'tblastn'
        else:
            self.program = 'blastn'
        # Fields used for custom outfmt 6 BLAST output:
        self.fieldnames = ['query_id', 'subject_id', 'identical', 'mismatches', 'gaps',
                           'evalue', 'bit_score', 'query_length', 'subject_length', 'alignment_length',
                           'query_start', 'query_end', 'subject_start', 'subject_end',
                           'query_sequence', 'subject_sequence']
        self.extended_fieldnames = self.fieldnames
        self.extended_fieldnames.insert(14, 'percent_match')
        self.outfmt = "'6 qseqid sseqid nident mismatch gaps evalue bitscore qlen slen length " \
                      "qstart qend sstart send qseq sseq'"
        self.blast_reports = list()
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
    SetupLogging(debug=arguments.verbose)
    # Run the allele-finding pipeline
    finder = allele_finder.AlleleFinder(path=arguments.path,
                                        targetfile=arguments.targetfile,
                                        analysis_type=arguments.blast,
                                        fasta_path=arguments.fasta_path,
                                        genesippr=arguments.genesippr,
                                        metadata_file=arguments.metadatafile,
                                        cutoff=arguments.cutoff,
                                        amino_acid=arguments.amino_acid)
    finder.main()
    # Extract the dictionary of records from the allele finding
    records = finder.records
    logging.info('Allele finding complete')
    # Run the profiling pipeline
    profiler = ProfileAlleles(path=arguments.path,
                              fasta_path=arguments.fasta_path,
                              records=records,
                              amino_acid=arguments.amino_acid)
    profiler.main()
    logging.info('Allele Profiling complete')


if __name__ == '__main__':
    cli()
