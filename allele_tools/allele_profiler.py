#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import combinetargets, GenObject, make_path, MetadataObject, \
    run_subprocess, SetupLogging
from genemethods.geneseekr.geneseekr import GeneSeekr
import allele_finder
from Bio.Blast.Applications import NcbiblastnCommandline, NcbitblastnCommandline
from Bio.Application import ApplicationError
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
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


def allele_prep(allele_path, gene_names, combined_targets, amino_acid):
    """
    Create a 'combinedtargets.fasta' file
    """
    logging.info('Reading allele files')
    custom_targets = os.path.join(allele_path, 'custom.tfa')
    combined_targets_db = glob(os.path.join(allele_path, 'combinedtargets*'))
    records = dict()
    # Clear out any previously created combinedtargets files (and associated BLAST databases)
    clear_alleles(combined_targets_db=combined_targets_db,
                  custom_targets=custom_targets)
    alleles = glob(os.path.join(allele_path, '*.*fa*'))
    # If the dictionary hasn't already been populated by a previous iteration, remove the path and the extension
    # from each of the allele files to determine gene names
    if not gene_names:
        for allele in alleles:
            gene_names.append(os.path.splitext(os.path.basename(allele))[0].replace('_alleles', ''))
    # Populate the header string for the final report
    data = 'ST\t{genes}\n'.format(genes='\t'.join(sorted(gene_names)))
    # Create the combinedtargets file
    if not os.path.isfile(combined_targets):
        if amino_acid:
            combinetargets(targets=alleles,
                           targetpath=allele_path,
                           mol_type='prot')
        else:
            combinetargets(targets=alleles,
                           targetpath=allele_path)
    # Create BLAST databases
    GeneSeekr.makeblastdb(fasta=combined_targets,
                          program='blastn' if not amino_acid else 'blastp')
    # Load all the FASTA records from the combinedtargets file
    for allele_file in alleles:
        for record in SeqIO.parse(allele_file, 'fasta'):
            records[record.id] = str(record.seq)
    return records, gene_names, data


def clear_alleles(combined_targets_db, custom_targets):
    """
    Remove any combinedtargets.fasta or custom.tfa files present in the allele path
    :param combined_targets_db: List of combinedtargets files, including BLAST database files
    :param custom_targets: Name and absolute path to custom.tfa target file
    :return:
    """
    # Iterate through all the combinedtargets files
    for combined_target in combined_targets_db:
        # Remove each file
        os.remove(combined_target)
    # Attempt to remove the targets.tfa file
    try:
        os.remove(custom_targets)
    except FileNotFoundError:
        pass


def read_profile(profile_file):
    """
    Load any previously created profiles for this analysis
    :param profile_file: Name and absolute path to a profile file
    """
    # Initialise the dictionary
    profile_data = dict()
    # Only load the profile file if it exists
    if os.path.isfile(profile_file):
        logging.info('Extracting profiles from profile file')
        # Open an Excel-formatted sequence profile file as a dictionary
        profile = DictReader(open(profile_file), dialect='excel-tab')
        # Iterate through the rows
        for row in profile:
            # Populate the profile dictionary with profile number: {gene: allele}.
            allele_comprehension = {gene: allele for gene, allele in row.items() if gene != 'ST'}
            # Extract the sequence type number from the first field name
            st = row[profile.fieldnames[0]]
            # Update the profile data dictionary
            profile_data[st] = allele_comprehension
    return profile_data


def parseable_blast_outputs(runmetadata, fieldnames, extended_fieldnames, records):
    """
    Add a header to the BLAST report, so that it is easier to figure out what is in each column
    :param runmetadata: Metadata object containing a list of all metadata objects
    :param fieldnames: String of all the field names in the BLAST report
    :param extended_fieldnames: String of the BLAST field names plus the calculated percent identity
    :param records: Dictionary of allele name: allele sequence
    """
    logging.info('Adding headers to BLAST outputs')
    for sample in runmetadata.samples:
        data = list()
        # Load the first line of the report
        with open(sample.alleles.blast_report, 'r') as report:
            header_line = report.readline().strip()
        # Split the header on tabs
        header_list = header_line.split('\t')
        # Check to see if the header has already been added. Skip this step if it has been added.
        if header_list[0] != fieldnames[0]:
            with open(sample.alleles.blast_report, 'r') as report:
                header = [entry for entry in report.readline().split('\t')]
            if len(header) == len(fieldnames):
                current_fieldnames = fieldnames
            else:
                current_fieldnames = extended_fieldnames
            blastdict = DictReader(open(sample.alleles.blast_report), fieldnames=current_fieldnames,
                                   dialect='excel-tab')
            # Go through each BLAST result
            for row in blastdict:
                # Create the target length variable - use query_length or subject length as desired
                # target_length = float(row['query_length']) if not genome_query else float(row['subject_length'])
                # Calculate the percent identity and extract the bitscore from the row
                # Percent identity is the (length of the alignment - num mismatches) / total query length
                percentidentity = float('{:0.2f}'.format((float(row['identical']) - float(row['gaps'])) /
                                                         len(records[row['subject_id']]) * 100))
                # Create a percent match entry based on the calculated percent identity match
                row['percent_match'] = str(percentidentity)
                # Add the updated row to the list
                data.append(row)

            # Overwrite the original BLAST outputs to include headers, and the percent match
            with open(sample.alleles.blast_report, 'w') as updated_report:
                # Add the header
                updated_report.write('{headers}\n'.format(headers='\t'.join(extended_fieldnames)))
                # Add the results
                for row in data:
                    for header in extended_fieldnames:
                        # Write the value from the row with the header as the key
                        try:
                            updated_report.write('{value}\t'.format(value=row[header]))
                        except KeyError:
                            # noinspection PyTypeChecker
                            updated_report.write('{value}\t'.format(value=''.join(row[None])))
                    # Add a newline for each result
                    updated_report.write('\n')


def parse_results(runmetadata, fieldnames, extended_fieldnames, amino_acid, genome_query=False):
    """
    Parse the BLAST results, and populate GenObjects
    :param runmetadata: Metadata object containing list of all metadata objects
    :param fieldnames: String of all the field names in the BLAST report
    :param extended_fieldnames: String of the BLAST field names plus the calculated percent identity
    :param amino_acid: Variable on whether targets are protein
    :param genome_query: Boolean of whether the query is a genome, and the subject is the allele
    :return: Updated runmetadata
    """
    logging.info('Parsing BLAST outputs')
    for sample in runmetadata.samples:
        # Initialise GenObjects as required
        sample.alleles.blastlist = list()
        sample.alleles.targetsequence = dict()
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(open(sample.alleles.blast_report),
                               fieldnames=extended_fieldnames,
                               dialect='excel-tab')
        resultdict = dict()
        # Go through each BLAST result
        for row in blastdict:
            # Ignore the headers
            if row['query_id'].startswith(fieldnames[0]):
                pass
            else:
                target_id = row['query_id'] if not genome_query else row['subject_id']
                target_start = row['subject_start']
                target_end = row['subject_end']
                target_seq = row['query_sequence']
                # Remove unwanted pipes added to the name
                target = target_id.lstrip('gb|').rstrip('|') if '|' in target_id else \
                    target_id
                # If the percent identity is greater than the cutoff
                if float(row['percent_match']) == 100:
                    # Append the hit dictionary to the list
                    sample.alleles.blastlist.append(row)
                    # Update the dictionary with the target and percent identity
                    resultdict.update({target: row['percent_match']})
                    # Determine if the orientation of the sequence is reversed compared to the reference
                    if int(target_end) < int(target_start) and not amino_acid:
                        # Create a sequence object using Biopython
                        seq = Seq(target_seq)
                        # Calculate the reverse complement of the sequence
                        querysequence = str(seq.reverse_complement())
                    # If the sequence is not reversed, use the sequence as it is in the output
                    else:
                        querysequence = target_seq
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
    return runmetadata


def profile_alleles(runmetadata, profile_dict, profile_set, records, amino_acid=False, novel_alleles=False,
                    genome_query=False, allele_path=None, report_path=None):
    """
    Create the gene:allele profile from the BLAST outputs from each sample
    :param runmetadata: Metadata object containing a list of all metadata objects
    :param profile_dict: Dictionary to store gene:allele profile for each sample
    :param profile_set: List of all unique profiles
    :param records: List of all gene names
    :param novel_alleles: Boolean of whether novel alleles should be extracted from BLAST hit if there is no 100% match
    :param genome_query: Boolean of whether the query is a genome, and the subject is the allele
    :param amino_acid: Variable on whether targets are protein
    :param allele_path: Name and absolute path to folder containing allele database files
    :param report_path: Name and absolute path to folder in which reports are to be created
    :return: Updated profile_dict and profile_set
    """
    logging.info('Determining allele profiles')
    # Iterate through all the samples
    for sample in runmetadata.samples:
        # Initialise a dictionary to store the profile information for each samples
        profile_dict[sample.name] = dict()
        # Initialise a dictionary to store the gene:allele combinations for each sample
        allele_comprehension = dict()
        # Each gene in the analysis is stored in the list of genes created in allele_finder
        for gene in records:
            # Initialise the variable to track whether the current gene is present in the current sample
            present = False
            # Iterate through all the BLAST outputs for the sample
            for allele in sample.alleles.blastresults:
                # If the gene name is present as a substring of the allele e.g. adk in adk_1, then the gene is
                # present in the BLAST outputs
                if gene in allele:
                    # Strip off the allele number from the allele e.g. adk_1 yields 1
                    allele_id = allele.split('_')[-1]
                    # Update the dictionary with the new gene: allele number for the sample
                    allele_comprehension.update({gene: allele_id})
                    # Update the gene presence variable
                    present = True
            # If, after iterating through all the BLAST outputs, the gene is not present in the sample, update the
            # gene: allele to reflect this absence
            if not present:
                if novel_alleles:
                    sample, novel_allele, query_sequence = extract_novel_alleles(sample=sample,
                                                                                 gene=gene,
                                                                                 genome_query=genome_query,
                                                                                 amino_acid=amino_acid,
                                                                                 allele_path=allele_path,
                                                                                 report_path=report_path)
                    try:
                        sample.alleles.targetsequence[novel_allele].append(query_sequence)
                    except KeyError:
                        sample.alleles.targetsequence[novel_allele] = [query_sequence]
                    if novel_allele:
                        allele_comprehension.update({gene: novel_allele.split('_')[-1]})
                    else:
                        allele_comprehension.update({gene: '0'})
                else:
                    # Set missing alleles to 'ND'
                    allele_comprehension.update({gene: '0'})
        # In order to hash the dictionary, use JSON, with sorted keys to freeze it
        frozen_allele_comprehension = json.dumps(allele_comprehension, sort_keys=True)
        # Update the dictionary of profiles with the hash of the frozen dictionary: list of samples with that hash
        if hash(frozen_allele_comprehension) not in profile_dict:
            profile_dict[hash(frozen_allele_comprehension)] = [sample.name]
        else:
            profile_dict[hash(frozen_allele_comprehension)].append(sample.name)
        # Add the 'regular' dictionary to the list of all profiles as required
        if allele_comprehension not in profile_set:
            profile_set.append(allele_comprehension)
    return profile_dict, profile_set


def extract_novel_alleles(sample, gene, genome_query, amino_acid, allele_path, report_path, cutoff=75):
    """
    Extract the sequence of novel alleles from samples that do not have a 100% match
    :param sample: Metadata object
    :param gene: Name of current gene
    :param genome_query: Boolean of whether the allele or the genome are the query
    :param amino_acid: Variable indicating whether the current analyses are on DNA or amino acid sequences
    :param allele_path: Name and absolute path to folder containing allele files
    :param report_path: Name and absolute path to folder in which reports are to be created
    :param cutoff: The minimum percent identity cutoff to allow when considering the presence of a sequence in a query
    :return: sample: Updated sample
    :return: novel_allele: Name of novel alleles discovered
    :return: query_sequence: Sequence of novel alleles discovered
    """
    # Open the sequence profile file as a dictionary
    blastdict = DictReader(open(sample.alleles.blast_report),
                           dialect='excel-tab')
    # Initialise a best hit value of 0
    best_hit = 0
    # Initialise strings to store the name and the sequence of novel alleles
    query_sequence = str()
    novel_allele = str()
    # Iterate through all the BLAST hits
    for row in blastdict:
        # Extract the target id with the appropriate key depending on whether genome files are the query or the subject
        target_id = row['query_id'] if not genome_query else row['subject_id']
        # Ensure that the gene name is present in the gene name + allele combination
        if gene in target_id:
            # Create a variable to store the value for percent identity, so it is easier to call
            perc_id = float(row['percent_match'])
            # See if the percent identity for the current match is better than the previous best match, and is above the
            # minimum cutoff threshold
            if perc_id > best_hit and perc_id >= cutoff:
                # Set the start and end variables depending on whether genomes are the query
                target_start = row['query_start'] if not genome_query else row['subject_start']
                target_end = row['query_end'] if not genome_query else row['subject_end']
                target_seq = row['query_sequence']
                # Determine if the orientation of the sequence is reversed compared to the reference
                if int(target_end) < int(target_start) and not amino_acid:
                    # Create a sequence object using Biopython
                    seq = Seq(target_seq)
                    # Calculate the reverse complement of the sequence
                    query_sequence = str(seq.reverse_complement())
                # If the sequence is not reversed, use the sequence as it is in the output
                else:
                    query_sequence = target_seq
                best_hit = perc_id
    # If a query sequence was extracted, use it to update the allele database
    if query_sequence:
        novel_allele = update_allele_database(gene=gene,
                                              query_sequence=query_sequence,
                                              allele_path=allele_path,
                                              report_path=report_path)
    return sample, novel_allele, query_sequence


def update_allele_database(gene, query_sequence, allele_path, report_path):
    """
    Update the allele database with the novel allele extracted above
    :param gene: Name of the current gene being examined
    :param query_sequence: Sequence of the novel allele
    :param allele_path: Name and absolute path to folder containing allele files
    :param report_path: Name and absolute path to folder in which reports are to be created
    :return: novel_allele: Name of the novel allele entered into the database
    """
    # Find the allele database file
    allele_file = glob(os.path.join(allele_path, f'{gene}*.*fa*'))[0]
    # Set the name of the novel allele file in the report path
    new_alleles = os.path.join(report_path, f'{gene}_novel_alleles.fasta')
    # Initialise a variable to store the name of the last allele in the database file
    last_id = str()
    # Create a list to store all the allele records in the database
    records = list()
    # Iterate through all the records in the allele database
    for record in SeqIO.parse(allele_file, 'fasta'):
        # Add the records to the list
        records.append(record)
        # Update the last_id variable
        last_id = record.id
    # Try to separate the gene name from the allele e.g. MutS_1
    try:
        gene_name, allele = last_id.split('_')
    # If there is no allele, set the allele to 1
    except ValueError:
        allele = 1
    # Name the novel allele as the gene name _ allele number + 1
    novel_allele = f'{gene}_{int(allele) + 1}'
    # Create a SeqRecord of the allele using the novel allele name and sequence
    new_record = SeqRecord(seq=Seq(query_sequence),
                           id=novel_allele,
                           name='',
                           description='')
    # Append the SeqRecord to the novel alleles file
    with open(new_alleles, 'a+') as novel:
        SeqIO.write(sequences=new_record,
                    handle=novel,
                    format='fasta')
    # Add the novel allele record to the list of all records
    records.append(new_record)
    # Overwrite the existing allele database file with the updated list of records
    with open(allele_file, 'w') as alleles:
        SeqIO.write(sequences=records,
                    handle=alleles,
                    format='fasta')
    return novel_allele


def match_profile(profile_data, profile_dict, profile_matches):
    """
    Match current profiles to any previously created profiles
    """
    # If the profile_data dictionary was not populated in the read_profiles methods, there is nothing to match
    if profile_data:
        logging.info('Matching new profiles against profile file')
        # Extract the sequence type and allele dictionary from the profile file
        for st, allele_comprehension in profile_data.items():
            # Freeze the allele comprehension as above
            frozen_allele_comprehension = json.dumps(allele_comprehension, sort_keys=True)
            try:
                # Extract the samples that match this profile
                matches = profile_dict[hash(frozen_allele_comprehension)]
                # Update the dictionary with the matching samples
                profile_matches[st] = matches
            # The profile will not necessarily match any of the profiles found in the analysis
            except KeyError:
                pass
    return profile_matches


def create_profile(profile_data, profile_set, new_profiles, profile_dict, profile_matches):
    """
    Create new profiles for novel profiles as required
    """
    # Initialise the sequence type to be 1
    seq_type = 1
    # If the profile_data dictionary exists, set the sequence type to be the last of the entries in the dictionary
    # plus one, as that corresponds to the next sequence type
    if profile_data:
        # seq_type = len(profile_data) + 1
        seq_type = sorted(int(st) for st in profile_data.keys())[-1] + 1
    # Initialise a list to store the matched samples
    matched = list()
    # Iterate through all the profiles in the analysis
    for allele_comprehension in profile_set:
        # Ensure that the allele comprehension (profile) is not already in the profile file
        if allele_comprehension not in [profiled_alleles for st, profiled_alleles in profile_data.items()]:
            # Add the new profile to the list of new profiles
            alleles = '\t'.join(allele_num.split('_')[-1] for gene, allele_num in sorted(allele_comprehension.items()))
            new_profiles.append('{profile_num}\t{alleles}'
                                .format(profile_num=seq_type,
                                        alleles=alleles.rstrip()))
            # Freeze the comprehension in order to be used as the key in the profile dictionary
            frozen_allele_comprehension = json.dumps(allele_comprehension, sort_keys=True)
            matches = profile_dict[hash(frozen_allele_comprehension)]
            # Check to see if this sequence type hasn't already been found in the current analysis
            if matches not in matched:
                # Update the dictionary with the new sequence type: list of samples
                profile_matches[seq_type] = matches
                profile_data[seq_type] = allele_comprehension
                # Add the matches to the list of matches
                matched.append(matches)
            # Increment the sequence type number of the next entry
            seq_type += 1
    return profile_matches, profile_data, new_profiles


def sequence_typer(profile_report, data, runmetadata, profile_matches, profile_data, update=False, amino_acid=False):
    """
    Perform the final sequence typing, and create the report
    """
    # Open the report
    mode = 'w' if not update else 'a+'
    if not update:
        # Initialise the header with an extra 'Sample' column plus the comma-separate list of gene names
        data = 'Sample\t' + data
    else:
        if not os.path.isfile(profile_report):
            # Initialise the header with an extra 'Sample' column plus the comma-separate list of gene names
            data = 'Sample\t' + data
        else:
            data = str()
    with open(profile_report, mode) as report:

        for sample in runmetadata.samples:
            # Iterate through all the matches to sequence profiles in the analysis
            for st, sample_names in profile_matches.items():
                # Check if the sample name is in the list of samples names with the current sequence type
                if sample.name in sample_names:
                    # Add the sample name, sequence type, and all the allele numbers to the report string
                    ac = '\t'.join(allele_num.split('_')[-1] for gene, allele_num in sorted(profile_data[st].items()))
                    data += '{sn}\t{st}\t{ac}\n' \
                        .format(sn=sample.name,
                                st=st,
                                ac=ac.rstrip())
                    # Update the appropriate GenObject based on the current molecule (DNA or amino acid)
                    if not amino_acid:
                        sample.alleles.nt_st = st
                        sample.alleles.nt_profile = profile_data[st]
                    else:
                        sample.alleles.aa_st = st
                        sample.alleles.aa_profile = profile_data[st]
        # Write the report
        report.write(data)
    return runmetadata


def append_profiles(new_profiles, profile_file, data, novel_profiles=False, profile_path=None, gene_names=None):
    """
    Add new profiles to the profile file
    """
    # Only try to add new profiles if there are new profiles in the analysis
    if new_profiles:
        # Initialise the string to store the new profile
        new_data = str()
        # If the profile file does not exist, add the string of 'ST', comma-separated gene names to the headers
        if not os.path.isfile(profile_file):
            new_data = data
        # Iterate through all the new profiles, and add them to the new profile string
        for profile in new_profiles:
            new_data += '{profile}\n'.format(profile=profile)
        # Open the report with a+ to either create, or append the profile string to it as appropriate
        with open(profile_file, 'a+') as profile_file:
            profile_file.write(new_data)
        if novel_profiles:
            novel_profile_file = os.path.join(profile_path, 'novel_profiles.txt')
            if not os.path.isfile(novel_profile_file):
                with open(novel_profile_file, 'w') as novel:
                    novel.write('ST\t{names}\n'.format(names='\t'.join(gene_names)))
            with open(novel_profile_file, 'a+') as novel:
                novel.write(new_data)


class ProfileAlleles(object):

    def main(self):
        self.blast()
        self.profile_dict, self.profile_set = profile_alleles(runmetadata=self.runmetadata,
                                                              profile_dict=self.profile_dict,
                                                              profile_set=self.profile_set,
                                                              records=self.records)
        self.profile_data = read_profile(profile_file=self.profile_file)
        self.profile_matches = match_profile(profile_data=self.profile_data,
                                             profile_dict=self.profile_dict,
                                             profile_matches=self.profile_matches)
        self.profile_matches, self.profile_data, self.new_profiles = \
            create_profile(profile_data=self.profile_data,
                           profile_set=self.profile_set,
                           new_profiles=self.new_profiles,
                           profile_dict=self.profile_dict,
                           profile_matches=self.profile_matches)
        self.runmetadata = sequence_typer(profile_report=self.profile_report,
                                          data=self.data,
                                          runmetadata=self.runmetadata,
                                          profile_matches=self.profile_matches,
                                          profile_data=self.profile_data)
        append_profiles(new_profiles=self.new_profiles,
                        profile_file=self.profile_file,
                        data=self.data,
                        gene_names=self.gene_names)

    def blast(self):
        """
        BLAST the allele files against the genomes
        """
        self.blast_db()
        self.run_blast()
        records, self.gene_names, self.data = \
            allele_prep(allele_path=self.targetpath,
                        gene_names=self.gene_names,
                        combined_targets=self.target_file,
                        amino_acid=self.amino_acid)
        parseable_blast_outputs(runmetadata=self.runmetadata,
                                fieldnames=self.fieldnames,
                                extended_fieldnames=self.extended_fieldnames,
                                records=records)
        self.runmetadata = parse_results(runmetadata=self.runmetadata,
                                         fieldnames=self.fieldnames,
                                         extended_fieldnames=self.extended_fieldnames,
                                         amino_acid=self.amino_acid)

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
                        if os.path.basename(fasta) != 'combinedtargets.fasta'
                        or os.path.basename(fasta) != 'custom.tfa']
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
        self.extended_fieldnames = self.fieldnames.copy()
        self.extended_fieldnames.insert(14, 'percent_match')
        self.outfmt = '6 qseqid sseqid nident mismatch gaps evalue bitscore qlen slen length ' \
                      'qstart qend sstart send qseq sseq'
        self.blast_reports = list()
        self.profile_dict = dict()
        self.profile_data = dict()
        self.profile_set = list()
        self.sequence_profile = dict()
        self.profile_matches = dict()
        self.new_profiles = list()
        # A string of the header to use for formatting the profile file, and the report headers
        genes = '\t'.join(sorted(self.records))
        self.data = 'ST\t{genes}\n'.format(genes=genes.rstrip())
        self.gene_names = list()


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
                                        target_alleles=arguments.no_target_alleles,
                                        allele_hashing=arguments.allele_hashing,
                                        amino_acid=arguments.amino_acid,
                                        one_based=arguments.one_based)
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
