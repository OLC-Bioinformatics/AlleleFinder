#!/usr/bin/env python

"""
Determine allele complement and profile of sequences. Update profiles and
databases as necessary
"""

# Standard imports
from argparse import ArgumentParser
from glob import glob
import logging
import json
import multiprocessing
import os


# Third-party imports
from Bio import SeqIO
from olctools.accessoryFunctions.accessoryFunctions import (
    make_path,
    MetadataObject,
    SetupLogging
)

# Local imports
from allele_tools.allele_profiler import (
    allele_prep,
    append_profiles,
    clear_alleles,
    create_profile,
    match_profile,
    parse_results,
    parseable_blast_outputs,
    profile_alleles,
    read_profile,
    sequence_typer
)
from allele_tools.methods import (
    update_allele_database,
    translate,
    query_prep,
    blast_alleles,
    pathfinder
)

__author__ = 'adamkoziol'


class Updater:
    """
    Determine allele complement and profile of sequences. Update profiles
    and databases as necessary
    """

    # pylint: disable=too-many-instance-attributes

    def main(self):
        """
        Run the appropriate methods in the correct order
        """
        # Create metadata objects for all files in the query folder
        self.runmetadata = query_prep(
            query_path=self.query_path,
            runmetadata=self.runmetadata
        )
        for sample in self.runmetadata.samples:
            logging.debug('Processing sample %s', sample.name)
            # Perform necessary prep on the alleles
            if not self.amino_acid:
                records, gene_names, self.data = \
                    allele_prep(
                        allele_path=self.allele_path,
                        gene_names=self.gene_names,
                        combined_targets=self.combined_targets,
                        amino_acid=self.amino_acid
                    )
            else:
                records, gene_names, self.data = \
                    allele_prep(
                        allele_path=self.aa_allele_path,
                        gene_names=self.gene_names,
                        combined_targets=self.combined_targets,
                        amino_acid=self.amino_acid
                    )

            logging.info('Loading profile')
            if not self.amino_acid:
                profile_data = read_profile(profile_file=self.profile_file)
            else:
                profile_data = read_profile(profile_file=self.aa_profile_file)
            # BLAST the query against the allele database
            blast_alleles(
                runmetadata=sample,
                amino_acid=self.amino_acid,
                combined_targets=self.combined_targets,
                cpus=self.cpus,
                outfmt=self.outfmt
            )
            # Add a user-friendly header to the BLAST outputs
            parseable_blast_outputs(
                runmetadata=sample,
                fieldnames=self.fieldnames,
                extended_fieldnames=self.extended_fieldnames,
                records=records
            )
            # Parse the BLAST outputs
            sample = parse_results(
                runmetadata=sample,
                fieldnames=self.fieldnames,
                extended_fieldnames=self.extended_fieldnames,
                amino_acid=self.amino_acid,
                genome_query=True
            )
            # Perform sequence typing of the parsed results
            if not self.amino_acid:
                profile_dict, profile_set = profile_alleles(
                    runmetadata=sample,
                    profile_dict={},
                    profile_set=[],
                    records=self.gene_names,
                    novel_alleles=True,
                    genome_query=True,
                    allele_path=self.allele_path,
                    report_path=self.report_path
                )
            else:
                profile_dict, profile_set = profile_alleles(
                    runmetadata=sample,
                    profile_dict={},
                    profile_set=[],
                    records=self.gene_names,
                    amino_acid=True,
                    novel_alleles=True,
                    genome_query=True,
                    allele_path=self.aa_allele_path,
                    report_path=self.report_path
                )
            # Match the query profile against the profile database
            profile_matches = match_profile(
                profile_data=profile_data,
                profile_dict=profile_dict,
                profile_matches={}
            )
            # Create new profiles as required
            profile_matches, profile_data, new_profiles = \
                create_profile(
                    profile_data=profile_data,
                    profile_set=profile_set,
                    new_profiles=[],
                    profile_dict=profile_dict,
                    profile_matches=profile_matches
                )
            # Perform final sequence typing, and create final report
            if not self.amino_acid:
                sample = sequence_typer(
                    profile_report=self.profile_report,
                    data=self.data,
                    runmetadata=sample,
                    profile_matches=profile_matches,
                    profile_data=profile_data,
                    update=True
                )
                # Write the novel profiles to file
                append_profiles(
                    new_profiles=new_profiles,
                    profile_file=self.profile_file,
                    data=self.data,
                    novel_profiles=True,
                    profile_path=self.profile_path,
                    gene_names=self.gene_names
                )
            else:
                sample = sequence_typer(
                    profile_report=self.aa_profile_report,
                    data=self.data,
                    runmetadata=sample,
                    profile_matches=profile_matches,
                    profile_data=profile_data,
                    update=True
                )
                append_profiles(
                    new_profiles=new_profiles,
                    profile_file=self.aa_profile_file,
                    data=self.data,
                    novel_profiles=True,
                    profile_path=self.aa_profile_path,
                    gene_names=self.gene_names
                )
            if not self.amino_acid:
                # AA
                sample = translate(runmetadata=sample)
                self.aa_allele_prep()
                aa_profile_dict, aa_profile_set = self.aa_allele_match(
                    runmetadata=sample,
                    profile_dict={},
                    profile_set=[],
                    gene_names=gene_names,
                    amino_acid=True
                )
                aa_profile_data = read_profile(
                    profile_file=self.aa_profile_file
                )
                aa_profile_matches = match_profile(
                    profile_data=aa_profile_data,
                    profile_dict=aa_profile_dict,
                    profile_matches={}
                )
                aa_profile_matches, aa_profile_data, aa_new_profiles = \
                    create_profile(
                        profile_data=aa_profile_data,
                        profile_set=aa_profile_set,
                        new_profiles=[],
                        profile_dict=aa_profile_dict,
                        profile_matches=aa_profile_matches
                    )

                sample = sequence_typer(
                    profile_report=self.aa_profile_report,
                    data=self.data,
                    runmetadata=sample,
                    profile_matches=aa_profile_matches,
                    profile_data=aa_profile_data,
                    update=True,
                    amino_acid=True
                    )
                make_path(self.aa_profile_path)
                append_profiles(
                    new_profiles=aa_new_profiles,
                    profile_file=self.aa_profile_file,
                    data=self.data,
                    novel_profiles=True,
                    profile_path=self.aa_profile_path,
                    gene_names=self.gene_names
                )
                self.aa_notes(runmetadata=sample)
                clear_alleles(
                    combined_targets_db=glob(
                        os.path.join(
                            self.allele_path,
                            'combinedtargets*'
                        )
                    ),
                    custom_targets=os.path.join(self.allele_path, 'custom.tfa')
                )

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
                allele_file = glob(
                    os.path.join(
                        self.aa_allele_path,
                        f'{gene}*.*fa*'
                    )
                )[0]
            # Create the file if it doesn't exist
            except IndexError:
                allele_file = self.initialise_aa_alleles(gene=gene)
            # Read in and store all the amino acid records in the allele
            # database file
            for record in SeqIO.parse(allele_file, 'fasta'):
                self.aa_allele_dict[record.id] = str(record.seq)

    def initialise_aa_alleles(self, gene):
        """
        Create a gene-specific amino acid allele database file
        :param gene: Name of current gene being analysed
        :return: Name and absolute path to the created database file
        """
        # Find the corresponding gene-specific nucleotide database file
        nt_allele_file = glob(
            os.path.join(
                self.allele_path,
                f'{gene}*.*fa*'
            )
        )[0]
        # Set the name of the amino acid database file
        aa_allele_file = os.path.join(
            self.aa_allele_path,
            f'{gene}_alleles.fasta'
        )
        # Initialise a dictionary to store str(amino acid sequence): allele
        # name to be used in finding duplicate translated alleles
        allele_dict = {}
        for record in SeqIO.parse(nt_allele_file, 'fasta'):
            # Replace and dashes in the record.id with underscores
            record.id = record.id.replace('-', '_')
            record_id = record.id
            # Translate the sequence to amino acid
            record = record.translate()
            record.id = record_id
            record.description = str()
            # Extract the gene name from the allele number
            gene, allele_id = record.id.rsplit('_', 1)
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
                    with open(
                            aa_allele_file,
                            'a+',
                            encoding='utf-8') as aa_targets:
                        SeqIO.write(record, aa_targets, 'fasta')
            else:
                if not self.aa_nt_allele_link_dict[gene]:
                    self.aa_nt_allele_link_dict[gene][record.id] = \
                        allele_dict[str(record.seq)]
        return aa_allele_file

    def aa_allele_match(
            self,
            runmetadata: list,
            profile_dict: dict,
            profile_set: list,
            gene_names: list,
            amino_acid: str):
        """
        Find match the alleles in the current sample to alleles in the database
        :param runmetadata: List of metadata objects for the current sample
        :param profile_dict: Dictionary to store gene:allele profile for
        each sample
        :param profile_set: List of all unique profiles
        :param gene_names: List of all gene names in the analysis
        :param amino_acid: Variable indicating whether the current analyses
        are on DNA or
        amino acid sequences
        :return: profile_dict and profile_set updated with the results from
        the current sample
        """
        for sample in runmetadata.samples:
            # Initialise a dictionary to store the gene:allele combinations
            # for each sample
            allele_comprehension = {}
            for allele, allele_seq in sorted(
                    sample.alleles.nt_alleles_translated.items()):
                present = False
                # Strip off the allele number from the allele e.g. adk_1
                # yields adk, while EC042_RS26480_2 yields
                # EC042_RS26480
                try:
                    gene, _ = allele.rsplit('_', 1)
                    gene = gene[0].lower() + gene[1:-1] + gene[-1].upper()
                except ValueError:
                    gene = str()
                for record_id, record_seq in sorted(
                        self.aa_allele_dict.items()):
                    if allele_seq == record_seq:
                        # Update the dictionary with the new gene: allele
                        # number for the sample
                        allele_comprehension.update(
                            {
                                gene: record_id.split('_')[-1]
                            }
                        )
                        present = True
                # If, after iterating through all the BLAST outputs, the gene
                # is not present in the sample, update the gene: allele to
                # reflect this absence
                if not present:
                    novel_allele = update_allele_database(
                        gene=gene,
                        query_sequence=allele_seq,
                        allele_path=self.aa_allele_path,
                        report_path=self.report_path,
                        amino_acid=amino_acid
                    )
                    if novel_allele:
                        allele_comprehension.update(
                            {
                                gene: novel_allele.split('_')[-1].split('|')[0]
                            }
                        )
                    else:
                        allele_comprehension.update({gene: '0'})
            profile_dict, profile_set, allele_comprehension = \
                self.profile_dict_set_populate(
                    gene_names=gene_names,
                    allele_comprehension=allele_comprehension,
                    profile_dict=profile_dict,
                    sample=sample,
                    profile_set=profile_set
                )
        return profile_dict, profile_set

    @staticmethod
    def profile_dict_set_populate(
            gene_names: list,
            allele_comprehension: dict,
            profile_dict: dict,
            sample: MetadataObject,
            profile_set: list):
        """
        Update the profile dictionary, profile set, and allele comprehension
        for each gene
        :param gene_names: List of all gene names in the analysis
        :param allele_comprehension: Dictionary of the gene:allele
        combinations for each sample
        :param profile_dict: Dictionary to store gene:allele profile for
        each sample
        :param sample: Metadata objects of current genome
        :param profile_set: List of all unique profiles
        """
        # Iterate through all the genes
        for gene_name in gene_names:
            if gene_name not in allele_comprehension:
                gene_name = gene_name[0].lower() + \
                    gene_name[1:-1] + \
                    gene_name[-1].upper()
                allele_comprehension.update({gene_name: '0'})
            # In order to hash the dictionary, use JSON, with sorted keys to
            # freeze it
            frozen_allele_comprehension = json.dumps(
                allele_comprehension,
                sort_keys=True
            )
            # Update the dictionary of profiles with the hash of the
            # frozen dictionary: list of samples with that hash
            if hash(frozen_allele_comprehension) not in profile_dict:
                profile_dict[hash(frozen_allele_comprehension)] = [sample.name]
            else:
                profile_dict[hash(frozen_allele_comprehension)].append(
                    sample.name
                )
            # Add the 'regular' dictionary to the list of all profiles
            # as required
            if allele_comprehension not in profile_set:
                profile_set.append(allele_comprehension)
        return profile_dict, profile_set, allele_comprehension

    def aa_notes(
            self,
            runmetadata: list):
        """
        Create (first time only), and update the profile and alleles notes
        files. These files link amino acid profile(s) and allele(s),
        respectively, to the corresponding nucleotide profile, and allele
        :param runmetadata: List of metadata objects for the current sample
        """
        logging.info('Creating/Updating notes')
        # Read in the profile notes file
        allele_profile_dict = self.read_aa_profile_notes()
        for sample in runmetadata.samples:
            # Check if the nucleotide sequence type has previously
            # been encountered
            try:
                # Extract all the previously encountered amino acid sequence
                # types corresponding
                # to the nucleotide st
                known_allele_profiles = allele_profile_dict[
                    sample.alleles.nt_st
                ]
                # Determine whether the current amino acid st is already in
                # the list of aa
                # sequence types
                if sample.alleles.aa_st not in known_allele_profiles:
                    # If the current st is novel, update the list of profiles,
                    # and use the list to update the notes file
                    allele_profile_dict[sample.alleles.nt_st].append(
                        sample.alleles.aa_st
                    )
                    self.update_aa_profile_notes(
                        aa_profile_dict=allele_profile_dict
                    )
            # If the current nucleotide sequence type is novel, add it to
            # the dictionary
            except KeyError:
                allele_profile_dict[sample.alleles.nt_st] = \
                    [sample.alleles.aa_st]
                # Use the dictionary to update the notes file
                self.update_aa_profile_notes(
                    aa_profile_dict=allele_profile_dict
                )
            # Process the allele file
            for gene, allele in sample.alleles.nt_profile.items():
                gene = gene[0].lower() + gene[1:-1] + gene[-1].upper()
                # Attempt to read in the previous allele notes file
                allele_notes_dict = self.read_aa_allele_notes(gene=gene)
                # Clean up the allele name
                allele = int(allele.split('|')[0])
                # Check if the nucleotide allele is present in the dictionary
                try:
                    # Extract the list of all amino acid alleles corresponding
                    # to the nucleotide allele
                    try:
                        known_aa_alleles = allele_notes_dict[allele]
                    # Initialise an empty list if allele_notes_dict is empty
                    except ValueError:
                        known_aa_alleles = list()
                    # Set the aa_allele variable to save typing the
                    # long attribute
                    aa_allele = sample.alleles.aa_profile[gene]
                    # Check to see if the amino acid allele has already
                    # been encountered
                    if aa_allele not in known_aa_alleles:
                        # Make sure that the allele isn't missing 'ND'
                        if aa_allele != '0':
                            # Add the allele to the list
                            allele_notes_dict[allele].append(aa_allele)
                            # Update the notes file with the
                            self.update_aa_allele_notes(
                                gene=gene,
                                allele_notes_dict=allele_notes_dict
                            )
                # If the nucleotide allele is novel, add it to the dictionary
                except KeyError:
                    aa_allele = sample.alleles.aa_profile[gene]
                    if aa_allele != '0':
                        allele_notes_dict[allele] = [aa_allele]
                        self.update_aa_allele_notes(
                            gene=gene,
                            allele_notes_dict=allele_notes_dict
                        )

    def read_aa_profile_notes(self) -> dict:
        """
        Read in all the notes from the profile file (if it exists)
        :return: aa_profile_dict: Dictionary of nucleotide seq type:
        amino acid sequence type(s)
        """
        aa_profile_dict = dict()
        try:
            with open(self.aa_profile_notes, 'r', encoding='utf-8') as profile:
                # Ignore the header
                _ = profile.readline()
                # Read in all the lines
                for line in profile:
                    # Split the nucleotide sequence type from the amino
                    # acid sequence type(s)
                    nt_profile, aa_profiles = line.rstrip().split('\t')
                    # Create a list of all the amino acid sequence types
                    aa_profile_dict[nt_profile] = aa_profiles.split(';')
        except FileNotFoundError:
            pass
        return aa_profile_dict

    def update_aa_profile_notes(
            self,
            aa_profile_dict: dict):
        """
        Update the profile file with a novel entry
        :param aa_profile_dict: Dictionary of nucleotide sequence type:
        amino acid sequence type(s)
        """
        # Overwrite the previous notes file
        with open(self.aa_profile_notes, 'w', encoding='utf-8') as profile:
            # Create the header for the profile file
            profile.write('nt_profile\taa_profile(s)\n')
            # Iterate through all the sequence type entries in the dictionary
            for nt_profile, aa_profiles in aa_profile_dict.items():
                aa_profiles = ';'.join(
                    [str(profile) for profile in aa_profiles]
                )
                # Write each nucleotide profile: aa_profile link
                profile.write(f'{nt_profile}\t{aa_profiles}\n')

    def read_aa_allele_notes(
            self,
            gene: str) -> dict:
        """
        Read in the allele notes file. If it doesn't exist, create it
        :param gene: Name of the current gene being analysed
        :return: Dictionary of nucleotide allele: list of corresponding
        amino acid alleles
        """
        allele_notes_dict = {}
        # Set the name of the notes file
        notes_file = os.path.join(
            self.aa_notes_path,
            f'{gene}_allele_notes.tsv'
        )
        # Attempt to read the notes file
        try:
            with open(notes_file, 'r', encoding='utf-8') as notes:
                # Ignore the header
                _ = notes.readline()
                # Extract the nt allele, ';'-separated list of amino acid
                # alleles from the line
                for line in notes:
                    nt_allele, aa_alleles = line.rstrip().split('\t')
                    # Create a list of all the amino acid alleles
                    allele_notes_dict[int(nt_allele)] = aa_alleles.split(';')
        # Create the file if it doesn't exist
        except FileNotFoundError:
            # Iterate through all the entries in the dictionary of gene name:
            # nt allele: aa allele
            for _, allele_dict in self.aa_nt_allele_link_dict.items():
                for nt_allele, aa_allele in allele_dict.items():
                    # Split the entries on underscores
                    _, nt_allele_number = nt_allele.rsplit('_', 1)
                    _, aa_allele_number = aa_allele.rsplit('_', 1)
                    # Update the notes dictionary with int(nt allele):
                    # list(aa allele)
                    allele_notes_dict[int(nt_allele_number)] = \
                        [aa_allele_number]
        return allele_notes_dict

    def update_aa_allele_notes(
            self,
            gene: str,
            allele_notes_dict: dict):
        """
        Update the amino acid allele notes file with the novel alleles
        :param gene: Current gene being analysed
        :param allele_notes_dict: Dictionary of nucleotide allele: list of
        corresponding amino acid alleles
        """
        # Set the name and path of the notes file
        notes_file = os.path.join(
            self.aa_notes_path,
            f'{gene}_allele_notes.tsv'
        )
        with open(notes_file, 'w', encoding='utf-8') as notes:
            # Create the header for the notes file
            notes.write('nt_allele\taa_alleles\n')
            for nt_allele, aa_alleles in sorted(allele_notes_dict.items()):
                # Write the nucleotide allele and list of corresponding amino
                # acid alleles
                aa_alleles = ';'.join(str(allele) for allele in aa_alleles)
                notes.write(f'{str(nt_allele)}\t{aa_alleles}\n')

    def __init__(
            self,
            path: str,
            amino_acid: bool):
        """
        Constructs all the necessary attributes for the Updater object.

        Parameters:
        path (str): Path to the file
        amino_acid (bool): Flag to indicate if the target is an amino acid
        """

        # Set paths
        self.path = pathfinder(path=path)
        self.allele_path = os.path.join(self.path, 'nt_alleles')
        self.aa_allele_path = os.path.join(self.path, 'aa_alleles')
        self.profile_path = os.path.join(self.path, 'nt_profile')
        self.aa_profile_path = os.path.join(self.path, 'aa_profile')

        # Create profile path
        make_path(self.profile_path)

        # Set profile files
        self.profile_file = os.path.join(self.profile_path, 'profile.txt')
        self.aa_profile_file = os.path.join(
            self.aa_profile_path,
            'profile.txt'
        )

        # Set query and report paths
        self.query_path = os.path.join(self.path, 'query')
        self.report_path = os.path.join(self.path, 'reports')

        # Create report path
        make_path(self.report_path)

        # Remove novel alleles
        novel_alleles = glob(os.path.join(self.report_path, '*.fasta'))
        for novel_allele in novel_alleles:
            os.remove(novel_allele)

        # Set notes path and file
        self.aa_notes_path = os.path.join(
            self.path,
            'aa_notes'
        )
        make_path(self.aa_notes_path)
        self.aa_profile_notes = os.path.join(
            self.aa_notes_path,
            'aa_profile_notes.tsv'
        )

        # Set amino acid flag
        self.amino_acid = amino_acid

        # Set combined targets
        if not self.amino_acid:
            self.combined_targets = os.path.join(
                self.allele_path,
                'combinedtargets.fasta'
            )
        else:
            self.combined_targets = os.path.join(
                self.aa_allele_path,
                'combinedtargets.fasta'
            )

        # Initialize other attributes
        self.gene_names = []
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = []
        self.cpus = multiprocessing.cpu_count() - 1
        self.profile_report = os.path.join(
            self.report_path,
            'nt_profiles.tsv'
        )
        self.aa_profile_report = os.path.join(
            self.report_path,
            'aa_profiles.tsv'
        )

        # Remove profile report if it exists
        try:
            os.remove(self.profile_report)
        except FileNotFoundError:
            pass

        # Set field names for BLAST output
        self.fieldnames = [
            'query_id',
            'subject_id',
            'identical',
            'mismatches',
            'gaps',
            'evalue',
            'bit_score',
            'query_length',
            'subject_length',
            'alignment_length',
            'query_start',
            'query_end',
            'subject_start',
            'subject_end',
            'query_sequence',
            'subject_sequence'
        ]
        self.extended_fieldnames = self.fieldnames.copy()
        self.extended_fieldnames.insert(14, 'percent_match')
        self.outfmt = '6 qseqid sseqid nident mismatch gaps evalue bitscore ' \
            'qlen slen length qstart qend sstart send qseq sseq'

        # Initialize data and dictionaries
        self.data = str()
        self.aa_allele_dict = {}
        self.aa_nt_allele_link_dict = {}


def cli():
    """
    Collect the arguments, create an object, and run the script
    """
    # Parser for arguments
    parser = ArgumentParser(
        description='Determines profiles of strains against previously '
        'calculated allele database and profile. Creates and/or updates both '
        'the database of allele definitions and the profile based on novel '
        'alleles and/or profiles discovered'
    )
    parser.add_argument(
        '-p', '--path',
        required=True,
        help='Specify path. Note that due to code reuse, the query sequence '
        'files must be in the "query" sub-folder, nucleotide alleles must be '
        'in the "nt_alleles" sub-folder, the nucleotide profile must be named '
        'profile.txt and be located in the "nt_profile" sub-folder, amino '
        'acid alleles must be in the "aa_allele" sub-folder, and the '
        'amino acid profile must be named profile.txt and be located in the '
        '"aa_profile" sub-folder '
    )
    parser.add_argument(
        '-aa', '--amino_acid',
        action='store_true',
        help='The query sequences are protein.'
    )
    # Get the arguments into an object
    arguments = parser.parse_args()
    SetupLogging(debug=True)
    # Run the profiling pipeline
    updater = Updater(
        path=arguments.path,
        amino_acid=arguments.amino_acid
    )
    updater.main()
    logging.info('Allele Updating complete')


if __name__ == '__main__':
    cli()
