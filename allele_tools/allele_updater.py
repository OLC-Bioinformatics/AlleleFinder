#!/usr/bin/env python
from olctools.accessoryFunctions.accessoryFunctions import GenObject, make_path, MetadataObject, \
    relative_symlink, SetupLogging
from allele_profiler import allele_prep, append_profiles, clear_alleles, create_profile, match_profile, parse_results, \
    parseable_blast_outputs, profile_alleles, read_profile, sequence_typer, update_allele_database
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastpCommandline
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq
from Bio import SeqIO
from argparse import ArgumentParser
import multiprocessing
from glob import glob
import logging
import json
import os


class Updater(object):

    def main(self):
        # Create metadata objects for all files in the query folder
        self.query_prep()
        for sample in self.runmetadata.samples:
            logging.warning('Processing sample {sn}'.format(sn=sample.name))
            if not self.amino_acid:
                records, gene_names, self.data = \
                    allele_prep(allele_path=self.allele_path,
                                gene_names=self.gene_names,
                                combined_targets=self.combined_targets,
                                amino_acid=self.amino_acid)
            else:
                records, gene_names, self.data = \
                    allele_prep(allele_path=self.aa_allele_path,
                                gene_names=self.gene_names,
                                combined_targets=self.combined_targets,
                                amino_acid=self.amino_acid)

            logging.info('Loading profile')
            if not self.amino_acid:
                profile_data = read_profile(profile_file=self.profile_file)
            else:
                profile_data = read_profile(profile_file=self.aa_profile_file)
            self.blast_alleles(runmetadata=sample,
                               amino_acid=self.amino_acid)
            parseable_blast_outputs(runmetadata=sample,
                                    fieldnames=self.fieldnames,
                                    extended_fieldnames=self.extended_fieldnames,
                                    records=records)
            sample = parse_results(runmetadata=sample,
                                   fieldnames=self.fieldnames,
                                   extended_fieldnames=self.extended_fieldnames,
                                   amino_acid=self.amino_acid,
                                   genome_query=True)
            if not self.amino_acid:
                profile_dict, profile_set = profile_alleles(runmetadata=sample,
                                                            profile_dict=dict(),
                                                            profile_set=list(),
                                                            records=self.gene_names,
                                                            novel_alleles=True,
                                                            genome_query=True,
                                                            allele_path=self.allele_path,
                                                            report_path=self.report_path)
            else:
                profile_dict, profile_set = profile_alleles(runmetadata=sample,
                                                            profile_dict=dict(),
                                                            profile_set=list(),
                                                            records=self.gene_names,
                                                            novel_alleles=True,
                                                            genome_query=True,
                                                            allele_path=self.aa_allele_path,
                                                            report_path=self.aa_report_path)
            profile_matches = match_profile(profile_data=profile_data,
                                            profile_dict=profile_dict,
                                            profile_matches=dict())
            profile_matches, profile_data, new_profiles = \
                create_profile(profile_data=profile_data,
                               profile_set=profile_set,
                               new_profiles=list(),
                               profile_dict=profile_dict,
                               profile_matches=profile_matches)
            if not self.amino_acid:
                sample = sequence_typer(profile_report=self.profile_report,
                                        data=self.data,
                                        runmetadata=sample,
                                        profile_matches=profile_matches,
                                        profile_data=profile_data,
                                        update=True)
                append_profiles(new_profiles=new_profiles,
                                profile_file=self.profile_file,
                                data=self.data,
                                novel_profiles=True,
                                profile_path=self.profile_path,
                                gene_names=self.gene_names)
            else:
                sample = sequence_typer(profile_report=self.aa_profile_report,
                                        data=self.data,
                                        runmetadata=sample,
                                        profile_matches=profile_matches,
                                        profile_data=profile_data,
                                        update=True)
                append_profiles(new_profiles=new_profiles,
                                profile_file=self.aa_profile_file,
                                data=self.data,
                                novel_profiles=True,
                                profile_path=self.aa_profile_path,
                                gene_names=self.gene_names)
            if not self.amino_acid:
                # AA
                sample = self.translate(runmetadata=sample)
                self.aa_allele_prep()
                aa_profile_dict, aa_profile_set = self.aa_allele_match(runmetadata=sample,
                                                                       profile_dict=dict(),
                                                                       profile_set=list(),
                                                                       gene_names=gene_names)
                aa_profile_data = read_profile(profile_file=self.aa_profile_file)
                aa_profile_matches = match_profile(profile_data=aa_profile_data,
                                                   profile_dict=aa_profile_dict,
                                                   profile_matches=dict())
                aa_profile_matches, aa_profile_data, aa_new_profiles = \
                    create_profile(profile_data=aa_profile_data,
                                   profile_set=aa_profile_set,
                                   new_profiles=list(),
                                   profile_dict=aa_profile_dict,
                                   profile_matches=aa_profile_matches)

                sample = sequence_typer(profile_report=self.aa_profile_report,
                                        data=self.data,
                                        runmetadata=sample,
                                        profile_matches=aa_profile_matches,
                                        profile_data=aa_profile_data,
                                        update=True,
                                        amino_acid=True)
                make_path(self.aa_profile_path)
                append_profiles(new_profiles=aa_new_profiles,
                                profile_file=self.aa_profile_file,
                                data=self.data,
                                novel_profiles=True,
                                profile_path=self.aa_profile_path,
                                gene_names=self.gene_names)
                self.aa_notes(runmetadata=sample)
                clear_alleles(combined_targets_db=glob(os.path.join(self.allele_path, 'combinedtargets*')),
                              custom_targets=os.path.join(self.allele_path, 'custom.tfa'))

    def query_prep(self):
        """
        Create metadata objects for each sample
        """
        logging.info('Preparing query files')
        # Find all the sequence files in the path
        fastas = sorted(glob(os.path.join(self.query_path, '*.fasta')))
        for fasta in fastas:
            name = os.path.splitext(os.path.basename(fasta))[0]
            if name != 'combinedtargets':
                # Create a metadata object for each sample
                metadata = MetadataObject()
                metadata.samples = list()
                # Populate the metadata object with the required attributes
                metadata.name = name
                metadata.general = GenObject()
                metadata.commands = GenObject()
                metadata.alleles = GenObject()
                metadata.alleles.outputdirectory = os.path.join(self.query_path, metadata.name)
                # Set the name of the BLAST output file
                metadata.alleles.blast_report = os.path.join(metadata.alleles.outputdirectory,
                                                             '{seq_id}.tsv'.format(seq_id=metadata.name))
                try:
                    os.remove(metadata.alleles.blast_report)
                except FileNotFoundError:
                    pass
                make_path(metadata.alleles.outputdirectory)
                metadata.general.bestassemblyfile = relative_symlink(src_file=fasta,
                                                                     output_dir=metadata.alleles.outputdirectory,
                                                                     export_output=True)
                metadata.samples.append(metadata)
                self.runmetadata.samples.append(metadata)

    def blast_alleles(self, runmetadata, amino_acid):
        """
        Run the BLAST analyses on the query
        :param runmetadata: List of metadata objects for each query
        :param amino_acid: Boolean of whether the query sequence is amino acid or nucleotide
        """
        logging.info('Running BLAST analyses')
        for sample in runmetadata.samples:
            if not amino_acid:
                blast = NcbiblastnCommandline(db=os.path.splitext(self.combined_targets)[0],
                                              query=sample.general.bestassemblyfile,
                                              num_alignments=100000000,
                                              evalue=0.001,
                                              num_threads=self.cpus,
                                              task='blastn',
                                              outfmt=self.outfmt,
                                              out=sample.alleles.blast_report)
            else:
                blast = NcbiblastpCommandline(query=sample.general.bestassemblyfile,
                                              db=os.path.splitext(self.combined_targets)[0],
                                              evalue=0.001,
                                              num_alignments=100000000,
                                              num_threads=self.cpus,
                                              outfmt=self.outfmt,
                                              out=sample.alleles.blast_report)
            blast()

    @staticmethod
    def translate(runmetadata):
        """
        Use BioPython to translate DNA to amino acid
        :param runmetadata: List of metadata objects for each query
        :return: Updated list of metadata objects
        """
        logging.info('Translating allele sequences to amino acid')
        for sample in runmetadata.samples:
            # Initialise the dictionary to store the translated sequence
            sample.alleles.nt_alleles_translated = dict()
            for allele, allele_sequence_list in sample.alleles.targetsequence.items():
                for allele_sequence in allele_sequence_list:
                    # Create a sequence object using Biopython
                    seq = Seq(allele_sequence)
                    try:
                        # Translate the sequence
                        aa_seq = str(seq.translate())
                    except TranslationError:
                        allele_seq = allele_sequence.replace('-', '')
                        seq = Seq(allele_seq)
                        aa_seq = str(seq.translate())
                    # Ensure that the allele name exists (isn't an empty string) before adding allele name: translated
                    # sequence to the dictionary
                    if allele:
                        sample.alleles.nt_alleles_translated[allele] = aa_seq
        return runmetadata

    def aa_allele_prep(self):
        """
        Create (first time only) and read the amino acid allele database file
        """
        # Create the amino acid allele database file path as required
        make_path(self.aa_allele_path)
        # Iterate through all the gene in the analysis
        for gene in self.gene_names:
            # Attempt to find the database file
            try:
                allele_file = glob(os.path.join(self.aa_allele_path, f'{gene}*.*fa*'))[0]
            # Create the file if it doesn't exist
            except IndexError:
                allele_file = self.initialise_aa_alleles(gene=gene)
            # Read in and store all the amino acid records in the allele database file
            for record in SeqIO.parse(allele_file, 'fasta'):
                self.aa_allele_dict[record.id] = str(record.seq)

    def initialise_aa_alleles(self, gene):
        """
        Create a gene-specific amino acid allele database file
        :param gene: Name of current gene being analysed
        :return: Name and absolute path to the created database file
        """
        # Find the corresponding gene-specific nucleotide database file
        nt_allele_file = glob(os.path.join(self.allele_path, f'{gene}*.*fa*'))[0]
        # Set the name of the amino acid database file
        aa_allele_file = os.path.join(self.aa_allele_path, f'{gene}_alleles.fasta')
        # Initialise a dictionary to store str(amino acid sequence): allele name to be used in finding duplicate
        # translated alleles
        allele_dict = dict()
        for record in SeqIO.parse(nt_allele_file, 'fasta'):
            # Replace and dashes in the record.id with underscores
            record.id = record.id.replace('-', '_')
            record_id = record.id
            # Translate the sequence to amino acid
            record = record.translate()
            record.id = record_id
            record.description = str()
            # Extract the gene name from the allele number
            gene, allele_id = record.id.split('_')
            # Initialise the gene key in the dictionary as required
            if gene not in self.aa_nt_allele_link_dict:
                self.aa_nt_allele_link_dict[gene] = dict()
            # Check if the translated sequence is not present in the dictionary
            if str(record.seq) not in allele_dict:
                # Update the dictionary with the sequence: allele
                allele_dict[str(record.seq)] = record.id
                # Update the dictionary with gene: allele
                if not self.aa_nt_allele_link_dict[gene]:
                    self.aa_nt_allele_link_dict[gene][record.id] = record.id
                if allele_id == '1':
                    # Write the translated target sequence to file
                    with open(aa_allele_file, 'a+') as aa_targets:
                        SeqIO.write(record, aa_targets, 'fasta')
            else:
                if not self.aa_nt_allele_link_dict[gene]:
                    self.aa_nt_allele_link_dict[gene][record.id] = allele_dict[str(record.seq)]
        return aa_allele_file

    def aa_allele_match(self, runmetadata, profile_dict, profile_set, gene_names):
        """
        Find match the alleles in the current sample to alleles in the database
        :param runmetadata: List of metadata objects for the current sample
        :param profile_dict: Dictionary to store gene:allele profile for each sample
        :param profile_set: List of all unique profiles
        :param gene_names: List of all gene names in the analysis
        :return: profile_dict and profile_set updated with the results from the current sample
        """
        for sample in runmetadata.samples:
            # Initialise a dictionary to store the gene:allele combinations for each sample
            allele_comprehension = dict()
            for allele, allele_seq in sorted(sample.alleles.nt_alleles_translated.items()):
                present = False
                # Strip off the allele number from the allele e.g. adk_1 yields 1
                try:
                    gene, allele_id = allele.split('_')
                except ValueError:
                    gene = str()
                for record_id, record_seq in sorted(self.aa_allele_dict.items()):
                    if allele_seq == record_seq:
                        # Update the dictionary with the new gene: allele number for the sample
                        allele_comprehension.update({gene: record_id.split('_')[-1]})
                        present = True
                # If, after iterating through all the BLAST outputs, the gene is not present in the sample, update the
                # gene: allele to reflect this absence
                if not present:
                    novel_allele = update_allele_database(gene=gene,
                                                          query_sequence=allele_seq,
                                                          allele_path=self.aa_allele_path,
                                                          report_path=self.aa_report_path)
                    if novel_allele:
                        allele_comprehension.update({gene: novel_allele.split('_')[-1]})
                    else:
                        allele_comprehension.update({gene: '0'})
            for gene_name in gene_names:
                if gene_name not in allele_comprehension:
                    allele_comprehension.update({gene_name: '0'})
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

    def aa_notes(self, runmetadata):
        """
        Create (first time only), and update the profile and alleles notes files. These files link amino acid profile(s)
        and allele(s), respectively, to the corresponding nucleotide profile, and allele
        :param runmetadata:
        """
        logging.info('Creating/Updating notes')
        # Read in the profile notes file
        allele_profile_dict = self.read_aa_profile_notes()
        for sample in runmetadata.samples:
            # Check if the nucleotide sequence type has previously been encountered
            try:
                # Extract all the previously encountered amino acid sequence types corresponding to the nucleotide st
                known_allele_profiles = allele_profile_dict[sample.alleles.nt_st]
                # Determine whether the current amino acid st is already in the list of aa sequence types
                if sample.alleles.aa_st not in known_allele_profiles:
                    # If the current st is novel, update the list of profiles, and use the list to update the
                    # notes file
                    allele_profile_dict[sample.alleles.nt_st].append(sample.alleles.aa_st)
                    self.update_aa_profile_notes(aa_profile_dict=allele_profile_dict)
            # If the current nucleotide sequence type is novel, add it to the dictionary
            except KeyError:
                allele_profile_dict[sample.alleles.nt_st] = [sample.alleles.aa_st]
                # Use the dictionary to update the notes file
                self.update_aa_profile_notes(aa_profile_dict=allele_profile_dict)
            # Process the allele file
            for gene, allele in sample.alleles.nt_profile.items():
                # Attempt to read in the previous allele notes file
                allele_notes_dict = self.read_aa_allele_notes(gene=gene)
                # Check if the nucleotide allele is present in the dictionary
                try:
                    # Extract the list of all amino acid alleles corresponding to the nucleotide allele
                    try:
                        known_aa_alleles = allele_notes_dict[int(allele)]
                    # Initialise an empty list if allele_notes_dict is empty
                    except ValueError:
                        known_aa_alleles = list()
                    # Set the aa_allele variable to save typing the long attribute
                    aa_allele = sample.alleles.aa_profile[gene]
                    # Check to see if the amino acid allele has already been encountered
                    if aa_allele not in known_aa_alleles:
                        # Make sure that the allele isn't missing 'ND'
                        if aa_allele != '0':
                            # Add the allele to the list
                            allele_notes_dict[int(allele)].append(aa_allele)
                            # Update the notes file with the
                            self.update_aa_allele_notes(gene=gene,
                                                        allele_notes_dict=allele_notes_dict)
                # If the nucleotide allele is novel, add it to the dictionary
                except KeyError:
                    aa_allele = sample.alleles.aa_profile[gene]
                    if aa_allele != '0':
                        allele_notes_dict[int(allele)] = [aa_allele]
                        self.update_aa_allele_notes(gene=gene,
                                                    allele_notes_dict=allele_notes_dict)

    def read_aa_profile_notes(self):
        """
        Read in all the notes from the profile file (if it exists)
        :return: aa_profile_dict: Dictionary of nucleotide sequence type: amino acid sequence type(s)
        """
        aa_profile_dict = dict()
        try:
            with open(self.aa_profile_notes, 'r') as profile:
                # Ignore the header
                _ = profile.readline()
                # Read in all the lines
                for line in profile:
                    # Split the nucleotide sequence type from the amino acid sequence type(s)
                    nt_profile, aa_profiles = line.rstrip().split('\t')
                    # Create a list of all the amino acid sequence types
                    aa_profile_dict[nt_profile] = aa_profiles.split(';')
        except FileNotFoundError:
            pass
        return aa_profile_dict

    def update_aa_profile_notes(self, aa_profile_dict):
        """
        Update the profile file with a novel entry
        :param aa_profile_dict: Dictionary of nucleotide sequence type: amino acid sequence type(s)
        """
        # Overwrite the previous notes file
        with open(self.aa_profile_notes, 'w') as profile:
            # Create the header for the profile file
            profile.write('nt_profile\taa_profile(s)\n')
            # Iterate through all the sequence type entries in the dictionary
            for nt_profile, aa_profiles in aa_profile_dict.items():
                # Write each nucleotide s
                profile.write('{nt_profile}\t{aa_profiles}\n'
                              .format(nt_profile=nt_profile,
                                      aa_profiles=';'.join([str(profile) for profile in aa_profiles])))

    def read_aa_allele_notes(self, gene):
        """
        Read in the allele notes file. If it doesn't exist, create it
        :param gene: Name of the current gene being analysed
        :return: Dictionary of nucleotide allele: list of corresponding amino acid alleles
        """
        allele_notes_dict = dict()
        # Set the name of the notes file
        notes_file = os.path.join(self.aa_notes_path, f'{gene}_allele_notes.tsv')
        # Attempt to read the notes file
        try:
            with open(notes_file, 'r') as notes:
                # Ignore the header
                _ = notes.readline()
                # Extract the nucleotide allele, ';'-separated list of amino acid alleles from the line
                for line in notes:
                    nt_allele, aa_alleles = line.rstrip().split('\t')
                    # Create a list of all the amino acid alleles
                    allele_notes_dict[int(nt_allele)] = aa_alleles.split(';')
        # Create the file if it doesn't exist
        except FileNotFoundError:
            # Iterate through all the entries in the dictionary of gene name: nt allele: aa allele
            for gene, allele_dict in self.aa_nt_allele_link_dict.items():
                for nt_allele, aa_allele in allele_dict.items():
                    # Split the entries on underscores
                    gene_name, nt_allele_number = nt_allele.split('_')
                    gene_name, aa_allele_number = aa_allele.split('_')
                    # Update the notes dictionary with int(nt allele): list(aa allele)
                    allele_notes_dict[int(nt_allele_number)] = [aa_allele_number]
        return allele_notes_dict

    def update_aa_allele_notes(self, gene, allele_notes_dict):
        """
        Update the amino acid allele notes file with the novel alleles
        :param gene: Current gene being analysed
        :param allele_notes_dict: Dictionary of nucleotide allele: list of corresponding amino acid alleles
        """
        # Set the name and path of the notes file
        notes_file = os.path.join(self.aa_notes_path, f'{gene}_allele_notes.tsv')
        with open(notes_file, 'w') as notes:
            # Create the header for the notes file
            notes.write('nt_allele\taa_alleles\n')
            for nt_allele, aa_alleles in sorted(allele_notes_dict.items()):
                # Write the nucleotide allele and list of corresponding amino acid alleles
                notes.write('{nt_allele}\t{aa_alleles}\n'
                            .format(nt_allele=str(nt_allele),
                                    aa_alleles=';'.join(str(allele) for allele in aa_alleles)))

    def __init__(self, path, amino_acid):
        if path.startswith('~'):
            self.path = os.path.abspath(os.path.expanduser(os.path.join(path)))
        else:
            self.path = os.path.abspath(os.path.join(path))
        self.allele_path = os.path.join(self.path, 'alleles')
        self.aa_allele_path = os.path.join(self.path, 'aa_alleles')
        self.profile_path = os.path.join(self.path, 'profile')
        self.aa_profile_path = os.path.join(self.path, 'aa_profile')
        make_path(self.profile_path)
        self.profile_file = os.path.join(self.profile_path, 'profile.txt')
        self.aa_profile_file = os.path.join(self.aa_profile_path, 'aa_profile.txt')
        self.query_path = os.path.join(self.path, 'query')
        self.report_path = os.path.join(self.path, 'reports')
        self.aa_report_path = os.path.join(self.path, 'aa_reports')
        make_path(self.report_path)
        make_path(self.aa_report_path)
        novel_alleles = glob(os.path.join(self.report_path, '*.fasta'))
        for novel_allele in novel_alleles:
            os.remove(novel_allele)
        self.aa_notes_path = os.path.join(self.path, 'aa_notes')
        make_path(self.aa_notes_path)
        self.aa_profile_notes = os.path.join(self.aa_notes_path, 'aa_profile_notes.tsv')
        self.amino_acid = amino_acid
        if not self.amino_acid:
            self.combined_targets = os.path.join(self.allele_path, 'combinedtargets.fasta')
        else:
            self.combined_targets = os.path.join(self.aa_allele_path, 'combinedtargets.fasta')
        self.gene_names = list()
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = list()
        self.cpus = multiprocessing.cpu_count() - 1
        self.profile_report = os.path.join(self.report_path, 'profiles.tsv')
        self.aa_profile_report = os.path.join(self.aa_report_path, 'aa_profiles.tsv')
        try:
            os.remove(self.profile_report)
        except FileNotFoundError:
            pass
        # Fields used for custom outfmt 6 BLAST output:
        self.fieldnames = ['query_id', 'subject_id', 'identical', 'mismatches', 'gaps',
                           'evalue', 'bit_score', 'query_length', 'subject_length', 'alignment_length',
                           'query_start', 'query_end', 'subject_start', 'subject_end',
                           'query_sequence', 'subject_sequence']
        self.extended_fieldnames = self.fieldnames.copy()
        self.extended_fieldnames.insert(14, 'percent_match')
        self.outfmt = '6 qseqid sseqid nident mismatch gaps evalue bitscore qlen slen length ' \
                      'qstart qend sstart send qseq sseq'
        # A string of the header to use for formatting the profile file, and the report headers
        self.data = str()
        self.aa_allele_dict = dict()
        self.aa_nt_allele_link_dict = dict()


def cli():
    # Parser for arguments
    parser = ArgumentParser(description='Determines profiles of strains against previously calculated allele database '
                                        'and profile. Creates and/or updates both the database of allele definitions '
                                        'and the profile based on novel alleles and/or profiles discovered')
    parser.add_argument('-p', '--path',
                        required=True,
                        help='Specify path. Note that due to code reuse, the query sequence files must be in the '
                             '"query" sub-folder, the alleles must be in the "alleles" sub-folder')
    parser.add_argument('-aa', '--amino_acid',
                        action='store_true',
                        help='The query sequences are protein.')
    # Get the arguments into an object
    arguments = parser.parse_args()
    SetupLogging(debug=True)
    # Run the profiling pipeline
    updater = Updater(path=arguments.path,
                      amino_acid=arguments.amino_acid)
    updater.main()
    logging.info('Allele Updating complete')


if __name__ == '__main__':
    cli()
