#!/usr/bin/env python

"""
Run allele finding specifically for the co-located subunits of STEC genes
"""

# Standard imports
from argparse import (
    ArgumentParser,
    RawTextHelpFormatter
)
from glob import glob
import logging
import multiprocessing
import sys
import os

# Third party inputs
from olctools.accessoryFunctions.accessoryFunctions import (
    make_path,
    MetadataObject,
)

# Local imports
from allele_tools.allele_profiler import (
    allele_prep,
    parseable_blast_outputs,
    read_profile
)
from allele_tools.allele_translate_reduce import Translate
from allele_tools.profile_reduce import ProfileReduce
from allele_tools.version import __version__
from allele_tools.methods import (
    analyse_aa_alleles,
    blast_alleles,
    blastx_alleles,
    common_allele_find_errors,
    concatenate_alleles,
    create_aa_allele_comprehension,
    create_frozen_allele_comprehension,
    create_gene_names,
    create_nt_allele_comprehension,
    create_stec_report,
    error_print,
    load_alleles,
    match_profile,
    parse_aa_blast,
    parse_colocated_results,
    pathfinder,
    profile_allele_check,
    query_prep,
    report_aa_alleles,
    setup_arguments,
    split_alleles,
    translated_update_nucleotide,
    write_concatenated_sequences
)

__author__ = 'adamkoziol'


class STEC:
    """
    Perform the STEC-specific allele finding
    """

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
            records, gene_names, self.data = \
                allele_prep(
                    allele_path=self.nt_allele_path,
                    gene_names=self.gene_names,
                    combined_targets=self.combined_targets,
                    amino_acid=self.amino_acid
                )
            gene_names = sorted(gene_names)
            logging.info('Loading profile')
            nt_profile_data = read_profile(profile_file=self.nt_profile_file)
            aa_profile_data = read_profile(profile_file=self.aa_profile_file)
            # Extract the sequence type and allele dictionary from the
            # profile file
            blast_alleles(
                runmetadata=sample,
                amino_acid=self.amino_acid,
                combined_targets=self.combined_targets,
                cpus=self.cpus,
                outfmt=self.outfmt
            )
            # Add the header to the BLAST outputs, and filter results based
            # on cutoff
            parseable_blast_outputs(
                runmetadata=sample,
                fieldnames=self.fieldnames,
                extended_fieldnames=self.extended_fieldnames,
                records=records,
                cutoff=90
            )
            # Parse the BLAST results
            sample, notes = parse_colocated_results(
                runmetadata=sample,
                fieldnames=self.fieldnames,
                extended_fieldnames=self.extended_fieldnames,
                amino_acid=self.amino_acid,
                gene_names=gene_names,
                nt_allele_path=self.nt_allele_path,
                aa_allele_path=self.aa_allele_path,
                report_path=self.report_path
            )
            # Create nucleotide allele comprehensions from the BLAST outputs
            nt_allele_comprehension = create_nt_allele_comprehension(
                runmetadata=sample,
                gene_names=gene_names
            )
            # Create an amino acid allele comprehensions from the translated
            # BLAST outputs
            aa_allele_comprehension = create_aa_allele_comprehension(
                runmetadata=sample,
                gene_names=gene_names,
            )
            # Freeze the nucleotide allele comprehension
            nt_frozen_allele_comprehension = \
                create_frozen_allele_comprehension(
                    allele_comprehension=nt_allele_comprehension
                )
            # Freeze the amino acid allele comprehension
            aa_frozen_allele_comprehension = \
                create_frozen_allele_comprehension(
                    allele_comprehension=aa_allele_comprehension
                )
            # Find nucleotide profile matches
            nt_profile_matches, nt_frozen_profiles = match_profile(
                profile_data=nt_profile_data,
                frozen_allele_comprehension=nt_frozen_allele_comprehension,
                report_path=self.report_path,
                profile_file=self.nt_profile_file,
                genes=gene_names,
                allele_comprehension=nt_allele_comprehension,
                molecule='nt'
            )
            # Find amino acid profile matches
            aa_profile_matches, aa_frozen_profiles = match_profile(
                profile_data=aa_profile_data,
                frozen_allele_comprehension=aa_frozen_allele_comprehension,
                report_path=self.report_path,
                profile_file=self.aa_profile_file,
                genes=gene_names,
                allele_comprehension=aa_allele_comprehension,
                molecule='aa'
            )
            # Create the STEC-specific report
            create_stec_report(
                runmetadata=sample,
                nt_profile_matches=nt_profile_matches,
                nt_alleles=nt_allele_comprehension,
                aa_profile_matches=aa_profile_matches,
                aa_alleles=aa_allele_comprehension,
                report_file=self.report_file,
                gene_names=gene_names,
                aa_profile_path=self.aa_profile_path,
                notes=notes
            )

    def __init__(
        self,
        allele_path: str,
        aa_allele_path: str,
        profile_file: str,
        aa_profile_file: str,
        query_path: str,
        report_path: str,
            amino_acid: bool = False):
        """
        Constructs all the necessary attributes for the STEC object.

        Args:
            allele_path (str): Path to the nucleotide allele file.
            aa_allele_path (str): Path to the amino acid allele file.
            profile_file (str): Path to the nucleotide profile file.
            aa_profile_file (str): Path to the amino acid profile file.
            query_path (str): Path to the query file.
            report_path (str): Path to the report directory.
            amino_acid (bool, optional): Flag to indicate if the process is
            for amino acids. Defaults to False.
        """
        # Set paths for allele and profile files
        self.nt_allele_path = pathfinder(path=allele_path)
        self.aa_allele_path = pathfinder(path=aa_allele_path)
        self.nt_profile_file = pathfinder(path=profile_file)
        self.aa_profile_file = pathfinder(path=aa_profile_file)

        # Set paths for profile directories
        self.profile_path = os.path.dirname(self.nt_profile_file)
        self.aa_profile_path = os.path.dirname(self.aa_profile_file)

        # Set paths for query and report files
        self.query_path = pathfinder(path=query_path)
        self.report_path = pathfinder(path=report_path)
        self.aa_report_path = self.report_path

        # Create report directory if it doesn't exist
        make_path(inpath=self.report_path)

        # Remove any existing fasta files in the report directory
        novel_alleles = glob(os.path.join(self.report_path, '*.fasta'))
        for novel_allele in novel_alleles:
            os.remove(novel_allele)

        # Set amino acid flag and combined targets path
        self.amino_acid = amino_acid
        if not self.amino_acid:
            self.combined_targets = os.path.join(
                self.nt_allele_path, 'combinedtargets.fasta')
        else:
            self.combined_targets = os.path.join(
                self.aa_allele_path, 'combinedtargets.fasta')

        # Initialize gene names and metadata
        self.gene_names = []
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = []

        # Set number of CPUs for multiprocessing
        self.cpus = multiprocessing.cpu_count() - 1

        # Set paths for profile report files
        self.profile_report = os.path.join(
            self.report_path, 'nt_profiles.tsv')
        self.aa_profile_report = os.path.join(
            self.aa_report_path, 'aa_profiles.tsv')

        # Remove existing profile report file if it exists
        try:
            os.remove(self.profile_report)
        except FileNotFoundError:
            pass

        # Set path for report file and remove any existing tsv files in the
        # report directory
        self.report_file = os.path.join(self.report_path, 'stec_report.tsv')
        reports = glob(os.path.join(self.report_path, '*.tsv'))
        for report in reports:
            os.remove(report)

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

        # Set format string for BLAST output
        self.outfmt = (
            '6 qseqid sseqid nident mismatch gaps evalue bitscore qlen slen '
            'length qstart qend sstart send qseq sseq'
        )

        # Initialize string for header formatting and dictionaries for alleles
        self.data = str()
        self.aa_allele_dict = {}
        self.aa_nt_allele_link_dict = {}


class STECTranslate:
    """
    Find alleles in genomes against an amino acid database
    """

    def main(self):
        """
        Run the necessary methods in the correct order
        """
        self.runmetadata = query_prep(
            query_path=self.query_path,
            runmetadata=self.runmetadata,
            clear_report=True
        )
        for sample in self.runmetadata.samples:
            logging.debug('Processing sample %s', sample.name)
            notes = {}
            records, gene_names, self.data = \
                allele_prep(
                    allele_path=self.aa_allele_path,
                    gene_names=self.gene_names,
                    combined_targets=self.combined_targets,
                    amino_acid=True
                )
            gene_names = sorted(gene_names)
            logging.info('Loading profile')
            nt_profile_data = read_profile(profile_file=self.nt_profile_file)
            aa_profile_data = read_profile(profile_file=self.aa_profile_file)
            if not os.path.isfile(sample.alleles.blast_report):
                # BLAST the query against the allele database
                blastx_alleles(
                    runmetadata=sample,
                    combined_targets=self.combined_targets,
                    cpus=self.cpus,
                    outfmt=self.outfmt
                )
            # Add headers to the BLAST outputs, and filter based on
            # cutoff value
            parseable_blast_outputs(
                runmetadata=sample,
                fieldnames=self.fieldnames,
                extended_fieldnames=self.extended_fieldnames,
                records=records,
                cutoff=95
            )
            # Parse the amino acid BLAST results
            sample, filtered, notes = parse_aa_blast(
                runmetadata=sample,
                extended_fieldnames=self.extended_fieldnames,
                fieldnames=self.fieldnames,
                gene_names=gene_names,
                notes=notes,
                aa_allele_path=self.aa_allele_path,
                report_path=self.report_path,
                cutoff=95
            )
            # Update the nucleotide database as required
            sample, notes = translated_update_nucleotide(
                runmetadata=sample,
                nt_allele_path=self.nt_allele_path,
                report_path=self.report_path,
                notes=notes,
                gene_names=self.gene_names,
                filtered=filtered
            )
            # Create nucleotide allele comprehensions from the
            # BLAST outputs
            nt_allele_comprehension = create_nt_allele_comprehension(
                runmetadata=sample,
                gene_names=gene_names,
                translated=True
            )
            # Create an amino acid allele comprehensions from the
            # translated BLAST outputs
            aa_allele_comprehension = create_aa_allele_comprehension(
                runmetadata=sample,
                gene_names=gene_names,
            )
            # Freeze the nucleotide allele comprehension
            nt_frozen_allele_comprehension = \
                create_frozen_allele_comprehension(
                    allele_comprehension=nt_allele_comprehension
                )
            # Freeze the amino acid allele comprehension
            aa_frozen_allele_comprehension = \
                create_frozen_allele_comprehension(
                    allele_comprehension=aa_allele_comprehension
                )
            # Find nucleotide profile matches
            nt_profile_matches, _ = match_profile(
                profile_data=nt_profile_data,
                frozen_allele_comprehension=nt_frozen_allele_comprehension,
                report_path=self.report_path,
                profile_file=self.nt_profile_file,
                genes=gene_names,
                allele_comprehension=nt_allele_comprehension,
                molecule='nt'
            )
            # Find amino acid profile matches
            aa_profile_matches, _ = match_profile(
                profile_data=aa_profile_data,
                frozen_allele_comprehension=aa_frozen_allele_comprehension,
                report_path=self.report_path,
                profile_file=self.aa_profile_file,
                genes=gene_names,
                allele_comprehension=aa_allele_comprehension,
                molecule='aa'
            )
            # Create the STEC-specific report
            create_stec_report(
                runmetadata=sample,
                nt_profile_matches=nt_profile_matches,
                nt_alleles=nt_allele_comprehension,
                aa_profile_matches=aa_profile_matches,
                aa_alleles=aa_allele_comprehension,
                report_file=self.report_file,
                gene_names=gene_names,
                aa_profile_path=self.aa_profile_path,
                notes=notes,
                molecule='aa'
            )

    def __init__(
            self,
            allele_path: str,
            aa_allele_path: str,
            profile_file: str,
            aa_profile_file: str,
            query_path: str,
            report_path: str):
        """
        Constructs all the necessary attributes for the STECTranslate object.

        Parameters:
        allele_path (str): Path to the allele file
        aa_allele_path (str): Path to the amino acid allele file
        profile_file (str): Path to the profile file
        aa_profile_file (str): Path to the amino acid profile file
        query_path (str): Path to the query file
        report_path (str): Path to the report file
        """

        # Set paths for allele and profile files
        self.nt_allele_path = pathfinder(path=allele_path)
        self.aa_allele_path = pathfinder(path=aa_allele_path)
        self.nt_profile_file = pathfinder(path=profile_file)
        self.aa_profile_file = pathfinder(path=aa_profile_file)

        # Set paths for profile directories
        self.profile_path = os.path.dirname(self.nt_profile_file)
        self.aa_profile_path = os.path.dirname(self.aa_profile_file)

        # Set paths for query and report files
        self.query_path = pathfinder(path=query_path)
        self.report_path = pathfinder(path=report_path)
        self.aa_report_path = self.report_path

        # Create report directory if it doesn't exist
        make_path(inpath=self.report_path)

        # Remove any existing fasta files in the report directory
        novel_alleles = glob(os.path.join(self.report_path, '*.fasta'))
        for novel_allele in novel_alleles:
            os.remove(novel_allele)

        # Set combined targets path
        self.combined_targets = os.path.join(
            self.aa_allele_path, 'combinedtargets.fasta')

        # Initialize gene names and metadata
        self.gene_names = []
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = []

        # Set number of CPUs for multiprocessing
        self.cpus = multiprocessing.cpu_count() - 1

        # Set paths for profile report files
        self.profile_report = os.path.join(
            self.report_path, 'nt_profiles.tsv')
        self.aa_profile_report = os.path.join(
            self.aa_report_path, 'aa_profiles.tsv')

        # Remove existing profile report file if it exists
        try:
            os.remove(self.profile_report)
        except FileNotFoundError:
            pass

        # Set path for report file and remove any existing tsv files in the
        # report directory
        self.report_file = os.path.join(self.report_path, 'stec_report.tsv')
        reports = glob(os.path.join(self.report_path, '*.tsv'))
        for report in reports:
            os.remove(report)

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

        # Set format string for BLAST output
        self.outfmt = (
            '6 qseqid sseqid nident mismatch gaps evalue bitscore qlen slen '
            'length qstart qend sstart send qseq sseq'
        )

        # Initialize string for header formatting and dictionaries for alleles
        self.data = str()
        self.aa_allele_dict = {}
        self.aa_nt_allele_link_dict = {}


class AASTEC:

    def main(self):
        """
        Run the necessary methods in the correct order
        """
        # Create metadata objects for all files in the query folder
        self.runmetadata = query_prep(
            query_path=self.query_path,
            runmetadata=self.runmetadata,
            clear_report=True
        )

        for sample in self.runmetadata.samples:
            logging.debug('Processing sample %s', sample.name)
            notes = {}
            records, gene_names, self.data = \
                allele_prep(
                    allele_path=self.aa_allele_path,
                    gene_names=self.gene_names,
                    combined_targets=self.combined_targets,
                    amino_acid=self.amino_acid
                )
            gene_names = sorted(gene_names)
            logging.info('Loading profile')
            if not os.path.isfile(sample.alleles.blast_report):
                # BLAST the query against the allele database
                blast_alleles(
                    runmetadata=sample,
                    amino_acid=self.amino_acid,
                    combined_targets=self.combined_targets,
                    cpus=self.cpus,
                    outfmt=self.outfmt
                )
                # Add headers to the BLAST outputs, and filter based on
                # cutoff value
                parseable_blast_outputs(
                    runmetadata=sample,
                    fieldnames=self.fieldnames,
                    extended_fieldnames=self.extended_fieldnames,
                    records=records,
                    cutoff=self.cutoff
                )
            # Parse the amino acid BLAST results
            sample, filtered, notes = parse_aa_blast(
                runmetadata=sample,
                extended_fieldnames=self.extended_fieldnames,
                fieldnames=self.fieldnames,
                gene_names=gene_names,
                notes=notes,
                aa_allele_path=self.aa_allele_path,
                report_path=self.report_path,
                cutoff=self.cutoff
            )
            # Perform content/length checks of the supplied alleles
            sample, notes = analyse_aa_alleles(
                runmetadata=sample,
                gene_names=gene_names,
                notes=notes)
            # Create an amino acid query report
            report_aa_alleles(
                runmetadata=sample,
                report_file=self.report_file,
                notes=notes
            )

    def __init__(
            self,
            aa_allele_path: str,
            query_path: str,
            report_path: str,
            cutoff: float,
            amino_acid: bool = True):
        """
        Constructs all the necessary attributes for the AASTEC object.

        Parameters:
        aa_allele_path (str): Path to the amino acid allele file
        query_path (str): Path to the query file
        report_path (str): Path to the report file
        cutoff (float): Cutoff value for the process
        amino_acid (bool): Flag to indicate if the process involves amino acids
        """

        # Set paths for allele, query, and report files
        self.aa_allele_path = pathfinder(path=aa_allele_path)
        self.query_path = pathfinder(path=query_path)
        self.report_path = pathfinder(path=report_path)
        self.aa_report_path = self.report_path

        # Create report directory if it doesn't exist
        make_path(inpath=self.report_path)

        # Remove any existing fasta files in the report directory
        novel_alleles = glob(os.path.join(self.report_path, '*.fasta'))
        for novel_allele in novel_alleles:
            os.remove(novel_allele)

        # Set cutoff and amino acid flag
        self.cutoff = cutoff
        self.amino_acid = amino_acid

        # Set combined targets path
        self.combined_targets = os.path.join(
            self.aa_allele_path, 'combinedtargets.fasta')

        # Initialize gene names and metadata
        self.gene_names = []
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = []

        # Set number of CPUs for multiprocessing
        self.cpus = multiprocessing.cpu_count() - 1

        # Set path for report file and remove it if it exists
        self.report_file = os.path.join(self.report_path, 'allele_report.tsv')
        try:
            os.remove(self.report_file)
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

        # Set format string for BLAST output
        self.outfmt = (
            '6 qseqid sseqid nident mismatch gaps evalue bitscore qlen slen '
            'length qstart qend sstart send qseq sseq'
        )

        # Initialize string for header formatting
        self.data = str()


class AlleleConcatenate:

    """
    Concatenate stx subunits. Read in profile files. Load alleles. Concatenate
    alleles with appropriate linker
    """

    def main(self):
        """
        Run the necessary methods for AlleleConcatenate
        """
        self.gene, self.nt_alleles = load_alleles(
            allele_path=self.nt_allele_path,
            allele_order=self.allele_order
        )
        self.gene, self.aa_alleles = load_alleles(
            allele_path=self.aa_allele_path,
            allele_order=self.allele_order
        )
        logging.info('Concatenating allele sequences')
        self.concatenated_nt_seq = concatenate_alleles(
            profile_data=self.nt_profile_data,
            allele_dict=self.nt_alleles,
            allele_order=self.allele_order,
            stx_gene=self.gene,
            linker_length_dict=self.linker_length_dict,
            molecule='nt'
        )
        self.concatenated_aa_seq = concatenate_alleles(
            profile_data=self.aa_profile_data,
            allele_dict=self.aa_alleles,
            allele_order=self.allele_order,
            stx_gene=self.gene,
            linker_length_dict=self.linker_length_dict,
            molecule='aa'
        )
        logging.info('Writing concatenated allele sequences to file')
        write_concatenated_sequences(
            concatenated_sequences=self.concatenated_nt_seq,
            concatenate_path=self.concatenate_path,
            file_name=self.gene_allele[self.gene],
            molecule='nt'
        )
        write_concatenated_sequences(
            concatenated_sequences=self.concatenated_aa_seq,
            concatenate_path=self.concatenate_path,
            file_name=self.gene_allele[self.gene],
            molecule='aa'
        )

    def __init__(
            self,
            nt_allele_path: str,
            aa_allele_path: str,
            nt_profile_file: str,
            aa_profile_file: str,
            concatenate_path: str):
        """
        Constructs all the necessary attributes for the AlleleConcatenate
        object.

        Parameters:
        nt_allele_path (str): Path to the nucleotide allele file
        aa_allele_path (str): Path to the amino acid allele file
        nt_profile_file (str): Path to the nucleotide profile file
        aa_profile_file (str): Path to the amino acid profile file
        concatenate_path (str): Path to the concatenate file
        """

        # Set paths for allele and profile files
        self.nt_allele_path = pathfinder(path=nt_allele_path)
        self.aa_allele_path = pathfinder(path=aa_allele_path)
        self.nt_profile_file = pathfinder(path=nt_profile_file)
        self.aa_profile_file = pathfinder(path=aa_profile_file)

        # Set paths for profile directories
        self.nt_profile_path = os.path.dirname(self.nt_profile_file)
        self.aa_profile_path = os.path.dirname(self.aa_profile_file)

        # Set path for concatenate file
        self.concatenate_path = pathfinder(path=concatenate_path)

        # Set linker length dictionary
        self.linker_length_dict = {
            'stx1': 9,
            'stx2': 12,
        }

        # Set the appropriate order for the genes in the report
        self.allele_order = {
            'stx1': ['stx1A', 'stx1B'],
            'stx2': ['stx2A', 'stx2B']
        }

        # Set gene allele dictionary
        self.gene_allele = {
            'stx1': 'stx1A_stx1B',
            'stx2': 'stx2A_stx2B'
        }

        # Read profile data
        self.nt_profile_data = read_profile(profile_file=self.nt_profile_file)
        self.aa_profile_data = read_profile(profile_file=self.aa_profile_file)

        # Initialize gene, alleles, and concatenated sequences
        self.gene = str()
        self.nt_alleles = {}
        self.aa_alleles = {}
        self.concatenated_nt_seq = []
        self.concatenated_aa_seq = []


def profile_reduce(args):
    """
    Reduce the Enterobase profile to only the genes of interest
    :param args: type ArgumentParser arguments
    """
    # Create the gene names file if it doesn't exist or is empty
    genes_path = os.path.dirname(args.gene_names)
    genes_file = os.path.basename(args.gene_names)
    logging.info(genes_path)
    logging.info(genes_file)
    if not os.path.isfile(args.gene_names):
        # Ensure that the path exists
        if not os.path.isdir(genes_path):
            logging.error(
                'Could not locate the supplied path, %s, for the gene file. '
                'Please ensure that it exists, and that it either contains '
                'the gene file, or files with .fasta extensions',
                genes_path
            )
            raise SystemExit
        logging.warning(
            'Could not locate the supplied gene file: %s. Will now attempt to '
            'create it in directory %s',
            args.gene_names, genes_path
        )
        # Attempt to create the file
        create_gene_names(
            path=genes_path,
            name=genes_file
        )
    else:
        # Ensure that the file isn't empty
        if os.stat(args.gene_names).st_size == 0:
            logging.warning(
                'The supplied gene file, %s, is empty. Will now attempt to '
                'populate it in directory %s',
                args.gene_names, genes_path
            )
            # Attempt to repopulate the file
            create_gene_names(
                path=genes_path,
                name=genes_file
            )
    logging.info(
        'Reducing profile file %s to include only the genes found in %s',
        args.profile_file, args.gene_names
        )
    # Create a ProfileReduce object
    profile_reduction = ProfileReduce(
        profile=args.profile_file,
        names=args.gene_names,
        output=args.output_folder
    )
    profile_reduction.main()
    logging.info('Profile reduction complete!')


def translate_reduce(args):
    """
    Translate nucleotide alleles to amino acid, and remove duplicates
    :param args: type ArgumentParser arguments
    """
    log_str = f'Translating and reducing alleles in: {args.allele_path}'
    if args.profile_file:
        log_str += \
            f'. Parsing {args.profile_file} nucleotide profile to create ' \
            f'corresponding, reduced amino acid profile'
    logging.info(log_str)
    length_dict = {
        'stx1B': 82,
        'Stx1B': 82,
        'stx1b': 82,
        'stx1A': 313,
        'Stx1A': 313,
        'stx1a': 313,
        'stx2A': 313,
        'Stx2A': 313,
        'stx2a': 313,
        'stx2B': 84,
        'Stx2B': 84,
        'stx2b': 84,
    }
    allele_translate_reduce = Translate(
        path=args.allele_path,
        profile=args.profile_file,
        report_path=args.report_path,
        translated_path=args.translated_path,
        length_dict=length_dict,
    )
    allele_translate_reduce.main()
    logging.info('Allele translation and reduction complete!')


def allele_find(args):
    """
    Perform allele discovery analyses
    :param args: type ArgumentParser arguments
    """
    log_str = \
        f'Performing STEC allele discovery on sequences in {args.query_path}' \
        f' using nucleotide alleles in {args.nt_alleles}, nucleotide profile' \
        f' in {args.nt_profile}, amino acid alleles in {args.aa_alleles}, ' \
        f'and amino acid profile in {args.aa_profile}'
    logging.info(log_str)
    errors = []

    # Nucleotide allele checks
    if not os.path.isdir(args.nt_alleles):
        errors.append(
            f'Could not find supplied nucleotide allele folder: '
            f'{args.nt_alleles}'
        )
    else:
        if not glob(os.path.join(args.nt_alleles, '*.fasta')):
            errors.append(
                f'Could not locate sequence files in supplied nucleotide '
                f'allele folder: {args.nt_alleles}'
            )
    # Find errors for amino acid and query checks
    errors = common_allele_find_errors(
        args=args,
        errors=errors,
        amino_acid=False
    )
    if not os.path.isfile(args.nt_profile):
        errors.append(
            f'Could not locate supplied nucleotide profile file: '
            f'{args.nt_profile}'
        )
    if errors:
        error_print(errors=errors)
    stec = STEC(
        allele_path=args.nt_alleles,
        aa_allele_path=args.aa_alleles,
        profile_file=args.nt_profile,
        aa_profile_file=args.aa_profile,
        query_path=args.query_path,
        report_path=args.report_path
    )
    stec.main()


def allele_translate_find(args):
    """
    Perform allele discovery analyses. Queries nucleotide sequences against
    amino acid alleles
    :param args: type ArgumentParser arguments
    """
    log_str = \
        f'Performing STEC allele discovery on sequences in {args.query_path}' \
        f' using amino acid alleles in {args.aa_alleles}, and amino acid ' \
        f'profile in {args.aa_profile}'
    logging.info(log_str)
    errors = []

    # Nucleotide allele checks
    if not os.path.isdir(args.nt_alleles):
        os.makedirs(args.nt_alleles)
    else:
        if not glob(os.path.join(args.nt_alleles, '*.fasta')):
            # Create the necessary allele files
            for gene in ['stx1.fasta', 'stx2.fasta']:
                open(os.path.join(args.nt_alleles, gene), 'a').close()
    # Find errors for amino acid and query checks
    errors = common_allele_find_errors(
        args=args,
        errors=errors,
        amino_acid=False
    )
    if not os.path.isfile(args.nt_profile):
        # Create the folder if required
        os.makedirs(os.path.dirname(args.nt_profile), exist_ok=True)
        # Create the empty profile file
        open(args.nt_profile, 'a').close()
    if errors:
        error_print(errors=errors)
    stec = STECTranslate(
        allele_path=args.nt_alleles,
        aa_allele_path=args.aa_alleles,
        profile_file=args.nt_profile,
        aa_profile_file=args.aa_profile,
        query_path=args.query_path,
        report_path=args.report_path
    )
    stec.main()


def aa_allele_find(args):
    """
    Perform allele discovery analyses on amino acid query files
    :param args: type ArgumentParser arguments
    """
    log_str = \
        f'Performing STEC allele discovery on amino acid sequences in ' \
        f'{args.query_path} using amino acid alleles in {args.aa_alleles}'
    logging.info(log_str)
    errors = []
    # Find errors for amino acid and query checks
    errors = common_allele_find_errors(
        args=args,
        errors=errors,
        amino_acid=True
    )
    if errors:
        error_print(errors=errors)
    aa_stec = AASTEC(
        aa_allele_path=args.aa_alleles,
        query_path=args.query_path,
        report_path=args.report_path,
        cutoff=args.cutoff
    )
    aa_stec.main()


def allele_split(args):
    """
    Split files containing multiple alleles into individual files
    :param args: type ArgumentParser arguments
    """
    logging.info(
        'Splitting allele files in %s into individual files', args.query_path
    )
    errors = []
    allele_files = glob(os.path.join(args.query_path, '*.fasta'))
    # Query checks
    if not os.path.isdir(args.query_path):
        errors.append(
            f'Could not find supplied nucleotide allele folder: '
            f'{args.query_path}'
        )
    else:
        if not allele_files:
            errors.append(
                f'Could not locate sequence files in supplied allele folder: '
                f'{args.query_path}'
            )
    if errors:
        error_print(errors=errors)
    split_alleles(
        allele_files=allele_files,
        output_path=args.output_path
    )


def allele_concatenate(args):
    """
    Concatenate subunit files with linkers. Provide linkages between
    nucleotide and amino acid files
    :param args: type ArgumentParser arguments
    """
    logging.info('Concatenating allele subunits')
    errors = []
    # Determine if the profile file, the allele folder, and the alleles
    # all exist
    errors = profile_allele_check(
        args=args,
        errors=errors
    )
    # If there were any errors with the supplied arguments, print them,
    # and exit
    if errors:
        error_print(errors=errors)
    concatenate = AlleleConcatenate(
        nt_allele_path=args.nt_alleles,
        aa_allele_path=args.aa_alleles,
        nt_profile_file=args.nt_profile,
        aa_profile_file=args.aa_profile,
        concatenate_path=args.concatenate_path
    )
    concatenate.main()


def cli():
    """
    Collect the arguments, create an object, and run the script
    """
    # Parser for arguments
    parser = ArgumentParser(
        description='Determines STEC subunit profiles'
        )
    subparsers = parser.add_subparsers(title='Available analyses')
    # Create a parental parser from subparsers can inherit arguments
    parent_parser = ArgumentParser(add_help=False)
    # Add arguments common to all subparsers to the parent parser
    parent_parser.add_argument(
        '-version', '--version',
        action='version',
        version=f'%(prog)s commit {__version__}'
        )
    parent_parser.add_argument(
        '-v', '--verbosity',
        choices=['debug', 'info', 'warning', 'error', 'critical'],
        metavar='verbosity',
        default='info',
        help='Set the logging level. Options are debug, info, warning, error, '
        'and critical. Default is info.'
    )
    profile_reduce_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='profile_reduce',
        description='Reduce full wgMLST profile from Enterobase using genes '
        'of interest',
        formatter_class=RawTextHelpFormatter,
        help='Reduce full wgMLST profile from Enterobase using genes '
        'of interest'
    )
    profile_reduce_subparser.add_argument(
        '-p', '--profile_file',
        metavar='profile_file',
        default='profiles.list',
        help='Specify name and path of profile file. If not provided, the '
        'default "profiles.list" in the current working directory will be used'
        )
    profile_reduce_subparser.add_argument(
        '-g', '--gene_names',
        metavar='gene_names',
        default='genes.txt',
        help='Name and path of text file containing gene names to use to '
        'filter the profile file (one per line). If not provided, the default '
        '"genes.txt" in the current working directory will be used. If the '
        'file does not exist, the program will attempt to create a file using '
        'the .fasta files in the current working directory'
    )
    profile_reduce_subparser.add_argument(
        '-o', '--output_folder',
        metavar='output_folder',
        default='nt_profile',
        help='Name and path of folder into which the reduced profile and notes'
        ' are to be placed. If not provided, the default "nt_profile" folder '
        'in the current working directory will be used'
    )
    profile_reduce_subparser.set_defaults(func=profile_reduce)
    # Create a subparser for allele translation and reduction
    allele_translate_reduce_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='allele_translate_reduce',
        description='Translate allele files in nucleotide format to amino '
        'acid. Remove duplicates. Keep notes',
        formatter_class=RawTextHelpFormatter,
        help='Translate allele files in nucleotide format to amino acid. '
             'Remove duplicates. Keep notes'
    )
    allele_translate_reduce_subparser.add_argument(
        '-a', '--allele_path',
        metavar='allele_path',
        default=os.path.join(os.getcwd(), 'nt_alleles'),
        help='Specify name and path of folder containing allele files. '
        'If not provided, the nt_alleles folder in the current working '
        'directory will be used by default'
    )
    allele_translate_reduce_subparser.add_argument(
        '-p', '--profile_file',
        metavar='profile_file',
        help='Optionally specify name and path of profile file. '
             'Parse the nucleic acid profile, and create the corresponding '
             'reduced amino acid profile'
        )
    allele_translate_reduce_subparser.add_argument(
        '-r', '--report_path',
        metavar='report_path',
        default=os.path.join(os.getcwd(), 'aa_profile'),
        help='Specify the name and path of the folder into which outputs are '
        'to be placed. If not provided, the aa_profile folder in the current '
        'working directory will be used'
    )
    allele_translate_reduce_subparser.add_argument(
        '-t', '--translated_path',
        metavar='translated_path',
        default=os.path.join(os.getcwd(), 'aa_alleles'),
        help='Specify the name and path of the folder into which alleles are '
        'to be placed. If not provided, the aa_alleles folder in the current '
        'working directory will be used'
    )
    allele_translate_reduce_subparser.set_defaults(func=translate_reduce)
    # Create a subparser for allele discovery
    allele_find_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='allele_find',
        description='Analyse sequences to determine allele complement. '
        'Update profiles and databases. Keep notes',
        formatter_class=RawTextHelpFormatter,
        help='Analyse sequences to determine allele complement. Update '
        'profiles and databases. Keep notes'
    )
    allele_find_subparser.add_argument(
        '--nt_profile',
        metavar='nt_profile',
        default=os.path.join(os.getcwd(), 'nt_profile', 'profile.txt'),
        help='Specify name and path of nucleotide profile file. If not '
        'provided, profile.txt in the nt_profile folder in the current '
        'working directory will be used by default'
    )
    allele_find_subparser.add_argument(
        '--aa_profile',
        metavar='aa_profile',
        default=os.path.join(os.getcwd(), 'aa_profile', 'profile.txt'),
        help='Specify name and path of amino acid profile file. If not '
        'provided, profile.txt in the aa_profile folder in the current '
        'working directory will be used by default'
    )
    allele_find_subparser.add_argument(
        '--nt_alleles',
        metavar='nt_alleles',
        default=os.path.join(os.getcwd(), 'nt_alleles'),
        help='Specify name and path of folder containing nucleotide alleles. '
        'If not provided, the nt_allele folder in the current working '
        'directory will be used by default'
    )
    allele_find_subparser.add_argument(
        '--aa_alleles',
        metavar='aa_alleles',
        default=os.path.join(os.getcwd(), 'aa_alleles'),
        help='Specify name and path of folder containing amino acid alleles. '
        'If not provided, the aa_allele folder in the current working '
        'directory will be used by default'
    )
    allele_find_subparser.add_argument(
        '-r', '--report_path',
        metavar='report_path',
        default=os.path.join(os.getcwd(), 'reports'),
        help='Specify name and path of folder into which reports are to be '
        'placed. If not provided, the reports folder in the current working '
        'directory will be used'
    )
    allele_find_subparser.add_argument(
        '-q', '--query_path',
        metavar='query_path',
        default=os.path.join(os.getcwd(), 'query'),
        help='Specify name and path of folder containing query files in FASTA '
        'format. If not provided, the query folder in the current working '
        'directory will be used'
    )
    allele_find_subparser.set_defaults(func=allele_find)
    # Create a subparser for allele discovery using nucleotide inputs against
    # an amino acid database
    allele_translate_find_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='allele_translate_find',
        description='Analyse nucleotide sequences against an amino acid '
        'database to determine allele complement. Update profiles and '
        'databases. Keep notes',
        formatter_class=RawTextHelpFormatter,
        help='Analyse nucleotide sequences against an amino acid database to '
        'determine allele complement. Update profiles and databases. '
        'Keep notes'
    )
    allele_translate_find_subparser.add_argument(
        '--nt_profile',
        metavar='nt_profile',
        default=os.path.join(os.getcwd(), 'nt_profile', 'profile.txt'),
        required=False,
        help='Specify name and path of nucleotide profile file. If not '
        'provided, profile.txt in the nt_profile folder in the current '
        'working directory will be used by default'
    )
    allele_translate_find_subparser.add_argument(
        '--aa_profile',
        metavar='aa_profile',
        default=os.path.join(os.getcwd(), 'aa_profile', 'profile.txt'),
        help='Specify name and path of amino acid profile file. If not '
        'provided, profile.txt in the aa_profile folder in the current '
        'working directory will be used by default'
    )
    allele_translate_find_subparser.add_argument(
        '--nt_alleles',
        metavar='nt_alleles',
        default=os.path.join(os.getcwd(), 'nt_alleles'),
        required=False,
        help='Specify name and path of folder containing nucleotide alleles. '
        'If not provided, the nt_allele folder in the current working '
        'directory will be used by default'
    )
    allele_translate_find_subparser.add_argument(
        '--aa_alleles',
        metavar='aa_alleles',
        default=os.path.join(os.getcwd(), 'aa_alleles'),
        help='Specify name and path of folder containing amino acid alleles. '
        'If not provided, the aa_allele folder in the current working '
        'directory will be used by default'
    )
    allele_translate_find_subparser.add_argument(
        '-r', '--report_path',
        metavar='report_path',
        default=os.path.join(os.getcwd(), 'reports'),
        help='Specify name and path of folder into which reports are to be '
        'placed. If not provided, the reports folder in the current '
        'working directory will be used'
    )
    allele_translate_find_subparser.add_argument(
        '-q', '--query_path',
        metavar='query_path',
        default=os.path.join(os.getcwd(), 'query'),
        help='Specify name and path of folder containing query files in '
        'FASTA format.  If not provided, the query folder in the current '
        'working directory will be used'
    )
    allele_translate_find_subparser.set_defaults(func=allele_translate_find)
    # Create a subparser for allele discovery using amino acid inputs against
    # an amino acid database
    aa_allele_find_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='aa_allele_find',
        description='Analyse amino acid sequences to determine allele '
        'complement. Update profiles and databases. Keep notes',
        formatter_class=RawTextHelpFormatter,
        help='Analyse amino acid sequences to determine allele complement. '
        'Update profiles and databases. Keep notes'
    )
    aa_allele_find_subparser.add_argument(
        '--aa_alleles',
        metavar='aa_alleles',
        default=os.path.join(os.getcwd(), 'aa_alleles'),
        help='Specify name and path of folder containing amino acid alleles. '
        'If not provided, the aa_allele folder in the current working '
        'directory will be used by default'
    )
    aa_allele_find_subparser.add_argument(
        '-r', '--report_path',
        metavar='report_path',
        default=os.path.join(os.getcwd(), 'reports'),
        help='Specify name and path of folder into which reports are to be '
        'placed. If not provided, the reports folder in the current '
        'working directory will be used'
    )
    aa_allele_find_subparser.add_argument(
        '-q', '--query_path',
        metavar='query_path',
        default=os.path.join(os.getcwd(), 'query'),
        help='Specify name and path of folder containing query files in '
        'FASTA format. If not provided, the query folder in the current '
        'working directory will be used'
    )
    aa_allele_find_subparser.add_argument(
        '-c', '--cutoff',
        metavar='cutoff',
        default=90,
        choices=[percent for percent in range(90, 101)],
        help='Specify the percent identity cutoff for matches. Allowed values '
        'are between 90 and 100. Default is 100'
    )
    aa_allele_find_subparser.set_defaults(func=aa_allele_find)
    # Create a subparser for splitting multi-FASTA files of alleles
    allele_split_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='allele_split',
        description='Split combined allele files into individual files',
        formatter_class=RawTextHelpFormatter,
        help='Split combined allele files into individual files'
    )
    allele_split_subparser.add_argument(
        '-q', '--query_path',
        metavar='query_path',
        default=os.path.join(os.getcwd(), 'query'),
        help='Specify name and path of folder containing query files in '
        'FASTA format. If not provided, the query folder in the current '
        'working directory will be used'
    )
    allele_split_subparser.add_argument(
        '-o', '--output_path',
        metavar='output_path',
        default=os.path.join(os.getcwd(), 'split_alleles'),
        help='Specify name and path of folder into which the split allele '
        'files are to be written. If not provided, the split_alleles folder '
        'in the current working directory will be used'
    )
    allele_split_subparser.set_defaults(func=allele_split)
    # Create a subparser for concatenating stx subunits
    allele_concatenate_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='allele_concatenate',
        description='Concatenate stx toxin subunit alleles with linkers',
        formatter_class=RawTextHelpFormatter,
        help='Concatenate stx toxin subunit alleles with linkers'
    )
    allele_concatenate_subparser.add_argument(
        '--nt_profile',
        metavar='nt_profile',
        default=os.path.join(os.getcwd(), 'nt_profile', 'profile.txt'),
        help='Specify name and path of nucleotide profile file. If not '
        'provided, profile.txt in the nt_profile folder in the current '
        'working directory will be used by default'
    )
    allele_concatenate_subparser.add_argument(
        '--aa_profile',
        metavar='aa_profile',
        default=os.path.join(os.getcwd(), 'aa_profile', 'profile.txt'),
        help='Specify name and path of amino acid profile file. If not '
        'provided, profile.txt in the aa_profile folder in the current '
        'working directory will be used by default'
    )
    allele_concatenate_subparser.add_argument(
        '--nt_alleles',
        metavar='nt_alleles',
        default=os.path.join(os.getcwd(), 'nt_alleles'),
        help='Specify name and path of folder containing nucleotide alleles. '
        'If not provided, the nt_allele folder in the current working '
        'directory will be used by default'
    )
    allele_concatenate_subparser.add_argument(
        '--aa_alleles',
        metavar='aa_alleles',
        default=os.path.join(os.getcwd(), 'aa_alleles'),
        help='Specify name and path of folder containing amino acid alleles. '
        'If not provided, the aa_allele folder in the current working '
        'directory will be used by default'
    )
    allele_concatenate_subparser.add_argument(
        '-c', '--concatenate_path',
        metavar='concatenate_path',
        default=os.path.join(os.getcwd(), 'concatenated_alleles'),
        help='Specify name and path of folder into which concatenated subunit '
        'files are to be placed. If not provided, the concatenated_alleles '
        'folder in the current working directory will be used'
    )
    allele_concatenate_subparser.set_defaults(func=allele_concatenate)
    # Get the arguments into an object
    arguments = setup_arguments(parser=parser)
    # Prevent the arguments being printed to the console (they are returned in
    # order for the tests to work)
    sys.stderr = open(os.devnull, 'w', encoding='utf-8')
    return arguments


if __name__ == '__main__':
    cli()
