#!/usr/bin/env python

"""
Run allele finding specifically for the co-located subunits of STEC genes
"""

# Standard imports
from argparse import \
    ArgumentParser, \
    RawTextHelpFormatter
import multiprocessing
from glob import glob
import logging
import sys
import os

# Third party inputs
from olctools.accessoryFunctions.accessoryFunctions import make_path, MetadataObject

# Local imports
from allele_tools.allele_translate_reduce import Translate
from allele_tools.profile_reduce import ProfileReduce
from allele_tools.version import __version__
from allele_tools.methods import \
        setup_arguments, \
        common_allele_find_errors, \
        pathfinder, \
        query_prep, \
        blast_alleles, \
        create_gene_names, \
        parse_colocated_results, \
        create_nt_allele_comprehension, \
        create_aa_allele_comprehension, \
        create_frozen_allele_comprehension, \
        match_profile, \
        create_stec_report, \
        split_alleles, \
        parse_aa_blast, \
        analyse_aa_alleles, \
        report_aa_alleles
from allele_tools.allele_profiler import \
    allele_prep, \
    parseable_blast_outputs, \
    read_profile


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
            # Extract the sequence type and allele dictionary from the profile file
            blast_alleles(
                runmetadata=sample,
                amino_acid=self.amino_acid,
                combined_targets=self.combined_targets,
                cpus=self.cpus,
                outfmt=self.outfmt
            )
            # Add the header to the BLAST outputs, and filter results based on cutoff
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
            # Create an amino acid allele comprehensions from the translated BLAST outputs
            aa_allele_comprehension = create_aa_allele_comprehension(
                runmetadata=sample,
                gene_names=gene_names,
            )
            # Freeze the nucleotide allele comprehension
            nt_frozen_allele_comprehension = create_frozen_allele_comprehension(
                allele_comprehension=nt_allele_comprehension
            )
            # Freeze the amino acid allele comprehension
            aa_frozen_allele_comprehension = create_frozen_allele_comprehension(
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

    def __init__(self, allele_path, aa_allele_path, profile_file, aa_profile_file, query_path, report_path,
                 amino_acid=False):

        if allele_path.startswith('~'):
            self.nt_allele_path = os.path.abspath(os.path.expanduser(os.path.join(allele_path)))
        else:
            self.nt_allele_path = os.path.abspath(os.path.join(allele_path))
        if aa_allele_path.startswith('~'):
            self.aa_allele_path = os.path.abspath(os.path.expanduser(os.path.join(aa_allele_path)))
        else:
            self.aa_allele_path = os.path.abspath(os.path.join(aa_allele_path))
        if profile_file.startswith('~'):
            self.nt_profile_file = os.path.abspath(os.path.expanduser(os.path.join(profile_file)))
        else:
            self.nt_profile_file = os.path.abspath(os.path.join(profile_file))
        if aa_profile_file.startswith('~'):
            self.aa_profile_file = os.path.abspath(
                os.path.expanduser(os.path.join(aa_profile_file))
            )
        else:
            self.aa_profile_file = os.path.abspath(os.path.join(aa_profile_file))
        self.profile_path = os.path.dirname(self.nt_profile_file)
        self.aa_profile_path = os.path.dirname(self.aa_profile_file)
        if query_path.startswith('~'):
            self.query_path = os.path.abspath(os.path.expanduser(os.path.join(query_path)))
        else:
            self.query_path = os.path.abspath(os.path.join(query_path))
        if report_path.startswith('~'):
            self.report_path = os.path.abspath(os.path.expanduser(os.path.join(report_path)))
        else:
            self.report_path = os.path.abspath(os.path.join(report_path))
        self.aa_report_path = self.report_path
        make_path(inpath=self.report_path)
        novel_alleles = glob(os.path.join(self.report_path, '*.fasta'))
        for novel_allele in novel_alleles:
            os.remove(novel_allele)

        self.amino_acid = amino_acid
        if not self.amino_acid:
            self.combined_targets = os.path.join(self.nt_allele_path, 'combinedtargets.fasta')
        else:
            self.combined_targets = os.path.join(self.aa_allele_path, 'combinedtargets.fasta')
        self.gene_names = []
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = []
        self.cpus = multiprocessing.cpu_count() - 1
        self.profile_report = os.path.join(self.report_path, 'nt_profiles.tsv')
        self.aa_profile_report = os.path.join(self.aa_report_path, 'aa_profiles.tsv')
        try:
            os.remove(self.profile_report)
        except FileNotFoundError:
            pass

        self.report_file = os.path.join(self.report_path, 'stec_report.tsv')
        reports = glob(os.path.join(self.report_path, '*.tsv'))
        for report in reports:
            os.remove(report)
        # Fields used for custom outfmt 6 BLAST output:
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
        self.outfmt = '6 qseqid sseqid nident mismatch gaps evalue bitscore qlen slen length ' \
                      'qstart qend sstart send qseq sseq'
        # A string of the header to use for formatting the profile file, and the report headers
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
            notes = []
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
                #
                blast_alleles(
                    runmetadata=sample,
                    amino_acid=self.amino_acid,
                    combined_targets=self.combined_targets,
                    cpus=self.cpus,
                    outfmt=self.outfmt
                )
                # Add headers to the BLAST outputs, and filter based on cutoff value
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

    def __init__(self, aa_allele_path, query_path, report_path, cutoff, amino_acid=True):

        self.aa_allele_path = pathfinder(path=aa_allele_path)
        self.query_path = pathfinder(path=query_path)
        self.report_path = pathfinder(path=report_path)
        self.aa_report_path = self.report_path
        make_path(inpath=self.report_path)
        novel_alleles = glob(os.path.join(self.report_path, '*.fasta'))
        for novel_allele in novel_alleles:
            os.remove(novel_allele)
        self.cutoff = cutoff
        self.amino_acid = amino_acid
        self.combined_targets = os.path.join(self.aa_allele_path, 'combinedtargets.fasta')
        self.gene_names = []
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = []
        self.cpus = multiprocessing.cpu_count() - 1
        self.report_file = os.path.join(self.report_path, 'allele_report.tsv')
        try:
            os.remove(self.report_file)
        except FileNotFoundError:
            pass
        # Fields used for custom outfmt 6 BLAST output:
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
        self.outfmt = '6 qseqid sseqid nident mismatch gaps evalue bitscore qlen slen length ' \
                      'qstart qend sstart send qseq sseq'
        # A string of the header to use for formatting the profile file, and the report headers
        self.data = str()


def profile_reduce(args):
    """
    Reduce the Enterobase profile to only the genes of interest
    :param args: type ArgumentParser arguments
    """
    # Create the gene names file if it doesn't exist or is empty
    genes_path = os.path.dirname(args.gene_names)
    genes_file = os.path.basename(args.gene_names)
    if not os.path.isfile(args.gene_names):
        # Ensure that the path exists
        if not os.path.isdir(genes_path):
            logging.error(
                'Could not locate the supplied path, %s, for the gene file. Please ensure that it'
                ' exists, and that it either contains the gene file, or files with '
                '.fasta extensions',
                genes_path
            )
            raise SystemExit
        logging.warning(
            'Could not locate the supplied gene file: %s. Will now attempt to create it in '
            'directory %s',
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
                'The supplied gene file, %s, is empty. Will now attempt to populate it in '
                'directory %s',
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
        log_str += f'. Parsing {args.profile_file} nucleotide profile to create corresponding,'\
                    ' reduced amino acid profile'
    logging.info(log_str)
    length_dict = {
        'ECs2973': 90,
        'ECs2974': 316,
        'ECs1205': 320,
        'ECs1206': 88
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
    log_str = f'Performing STEC allele discovery on sequences in {args.query_path} using ' \
        f'nucleotide alleles in {args.nt_alleles}, nucleotide profile in {args.nt_profile}, amino' \
        f' acid alleles in {args.aa_alleles}, and amino acid profile in {args.aa_profile}'
    logging.info(log_str)
    errors = []

    # Nucleotide allele checks
    if not os.path.isdir(args.nt_alleles):
        errors.append(f'Could not find supplied nucleotide allele folder: {args.nt_alleles}')
    else:
        if not glob(os.path.join(args.nt_alleles, '*.fasta')):
            errors.append(
                f'Could not locate sequence files in supplied nucleotide allele folder: {args.nt_alleles}'
            )
    # Find errors for amino acid and query checks
    errors = common_allele_find_errors(
        args=args,
        errors=errors,
        amino_acid=False
    )
    if not os.path.isfile(args.nt_profile):
        errors.append(f'Could not locate supplied nucleotide profile file: {args.nt_profile}')
    if errors:
        error_string = '\n'.join(errors)
        was_were = 'was' if len(errors) == 1 else 'were'
        correct = 'error' if len(errors) == 1 else 'errors'
        logging.error(
            'There %s %s %s when attempting to run your command: \n%s', was_were, len(errors), correct, error_string)
        raise SystemExit
    stec = STEC(
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
    log_str = f'Performing STEC allele discovery on amino acid sequences in {args.query_path} using amino' \
              f' acid alleles in {args.aa_alleles}'
    logging.info(log_str)
    errors = []
    # Find errors for amino acid and query checks
    errors = common_allele_find_errors(
        args=args,
        errors=errors,
        amino_acid=True
    )
    if errors:
        error_string = '\n'.join(errors)
        was_were = 'was' if len(errors) == 1 else 'were'
        correct = 'error' if len(errors) == 1 else 'errors'
        logging.error(
            'There %s %s %s when attempting to run your command: \n%s', was_were, len(errors), correct, error_string)
        raise SystemExit
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
    logging.info('Splitting allele files in %s into individual files', args.query_path)
    errors = []
    allele_files = glob(os.path.join(args.query_path, '*.fasta'))
    # Query checks
    if not os.path.isdir(args.query_path):
        errors.append(f'Could not find supplied nucleotide allele folder: {args.query_path}')
    else:
        if not allele_files:
            errors.append(
                f'Could not locate sequence files in supplied allele folder: {args.query_path}'
            )
    if errors:
        error_string = '\n'.join(errors)
        was_were = 'was' if len(errors) == 1 else 'were'
        correct = 'error' if len(errors) == 1 else 'errors'
        logging.error(
            'There %s %s %s when attempting to run your command: \n%s', was_were, len(errors), correct, error_string)
        raise SystemExit
    split_alleles(
        allele_files=allele_files,
        output_path=args.output_path
    )


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
        help='Set the logging level. Options are debug, info, warning, error, and critical. '
             'Default is info.'
    )
    profile_reduce_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='profile_reduce',
        description='Reduce full wgMLST profile from Enterobase using genes of interest',
        formatter_class=RawTextHelpFormatter,
        help='Reduce full wgMLST profile from Enterobase using genes of interest'
    )
    profile_reduce_subparser.add_argument(
        '-p', '--profile_file',
        metavar='profile_file',
        default='profiles.list',
        help='Specify name and path of profile file. If not provided, the default '
             '"profiles.list" in the current working directory will be used'
        )
    profile_reduce_subparser.add_argument(
        '-g', '--gene_names',
        metavar='gene_names',
        default='genes.txt',
        help='Name and path of text file containing gene names to use to filter the profile file '
             '(one per line). If not provided, the default "genes.txt" in the current working '
             'directory will be used. If the file does not exist, the program will attempt to '
             'create a file using the .fasta files in the current working directory'
    )
    profile_reduce_subparser.add_argument(
        '-o', '--output_folder',
        metavar='output_folder',
        default='nt_profile',
        help='Name and path of folder into which the reduced profile and notes are to be placed. '
             'If not provided, the default "nt_profile" folder in the current working directory '
             'will be used'
    )
    profile_reduce_subparser.set_defaults(func=profile_reduce)
    # Create a subparser for allele translation and reduction
    allele_translate_reduce_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='allele_translate_reduce',
        description='Translate allele files in nucleotide format to amino acid. Remove duplicates. Keep notes',
        formatter_class=RawTextHelpFormatter,
        help='Translate allele files in nucleotide format to amino acid. '
             'Remove duplicates. Keep notes'
    )
    allele_translate_reduce_subparser.add_argument(
        '-a', '--allele_path',
        metavar='allele_path',
        default=os.path.join(os.getcwd(), 'nt_alleles'),
        help='Specify name and path of folder containing allele files. If not provided, the '
             'nt_alleles folder in the current working directory will be used by default'
    )
    allele_translate_reduce_subparser.add_argument(
        '-p', '--profile_file',
        metavar='profile_file',
        help='Optionally specify name and path of profile file. '
             'Parse the nucleic acid profile, and create the corresponding reduced amino acid profile'
        )
    allele_translate_reduce_subparser.add_argument(
        '-r', '--report_path',
        metavar='report_path',
        default=os.path.join(os.getcwd(), 'aa_profile'),
        help='Specify the name and path of the folder into which outputs are to be placed. If not '
              'provided, the aa_profile folder in the current working directory will be used'
    )
    allele_translate_reduce_subparser.add_argument(
        '-t', '--translated_path',
        metavar='translated_path',
        default=os.path.join(os.getcwd(), 'aa_alleles'),
        help='Specify the name and path of the folder into which alleles are to be placed. If not '
              'provided, the aa_alleles folder in the current working directory will be used'
    )
    allele_translate_reduce_subparser.set_defaults(func=translate_reduce)
    # Create a subparser for allele discovery
    allele_find_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='allele_find',
        description='Analyse sequences to determine allele complement. Update profiles and '
        'databases. Keep notes',
        formatter_class=RawTextHelpFormatter,
        help='Analyse sequences to determine allele complement. Update profiles and databases. '
        'Keep notes'
    )
    allele_find_subparser.add_argument(
        '--nt_profile',
        metavar='nt_profile',
        default=os.path.join(os.getcwd(), 'nt_profile', 'profile.txt'),
        help='Specify name and path of nucleotide profile file. If not provided, profile.txt in '
             'the nt_profile folder in the current working directory will be used by default'
    )
    allele_find_subparser.add_argument(
        '--aa_profile',
        metavar='aa_profile',
        default=os.path.join(os.getcwd(), 'aa_profile', 'profile.txt'),
        help='Specify name and path of amino acid profile file. If not provided, profile.txt in '
             'the aa_profile folder in the current working directory will be used by default'
    )
    allele_find_subparser.add_argument(
        '--nt_alleles',
        metavar='nt_alleles',
        default=os.path.join(os.getcwd(), 'nt_alleles'),
        help='Specify name and path of folder containing nucleotide alleles. If not provided, the '
             'nt_allele folder in the current working directory will be used by default'
    )
    allele_find_subparser.add_argument(
        '--aa_alleles',
        metavar='aa_alleles',
        default=os.path.join(os.getcwd(), 'aa_alleles'),
        help='Specify name and path of folder containing amino acid alleles. If not provided, the '
             'aa_allele folder in the current working directory will be used by default'
    )
    allele_find_subparser.add_argument(
        '-r', '--report_path',
        metavar='report_path',
        default=os.path.join(os.getcwd(), 'reports'),
        help='Specify name and path of folder into which reports are to be placed. If not '
             'provided, the reports folder in the current working directory will be used'
    )
    allele_find_subparser.add_argument(
        '-q', '--query_path',
        metavar='query_path',
        default=os.path.join(os.getcwd(), 'query'),
        help='Specify name and path of folder containing query files in FASTA format. '
             'If not provided, the query folder in the current working directory will be used'
    )
    allele_find_subparser.set_defaults(func=allele_find)
    # Create a subparser for allele discovery
    aa_allele_find_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='aa_allele_find',
        description='Analyse amino acid sequences to determine allele complement. Update profiles and '
                    'databases. Keep notes',
        formatter_class=RawTextHelpFormatter,
        help='Analyse amino acid sequences to determine allele complement. Update profiles and databases. '
             'Keep notes'
    )
    aa_allele_find_subparser.add_argument(
        '--aa_alleles',
        metavar='aa_alleles',
        default=os.path.join(os.getcwd(), 'aa_alleles'),
        help='Specify name and path of folder containing amino acid alleles. If not provided, the '
             'aa_allele folder in the current working directory will be used by default'
    )
    aa_allele_find_subparser.add_argument(
        '-r', '--report_path',
        metavar='report_path',
        default=os.path.join(os.getcwd(), 'reports'),
        help='Specify name and path of folder into which reports are to be placed. If not '
             'provided, the reports folder in the current working directory will be used'
    )
    aa_allele_find_subparser.add_argument(
        '-q', '--query_path',
        metavar='query_path',
        default=os.path.join(os.getcwd(), 'query'),
        help='Specify name and path of folder containing query files in FASTA format. '
             'If not provided, the query folder in the current working directory will be used'
    )
    aa_allele_find_subparser.add_argument(
        '-c', '--cutoff',
        metavar='cutoff',
        default=90,
        choices=[percent for percent in range(90, 101)],
        help='Specify the percent identity cutoff for matches. Allowed values are between 90 and 100. Default is 100'
    )
    aa_allele_find_subparser.set_defaults(func=aa_allele_find)
    # Create a subparser for allele discovery
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
        help='Specify name and path of folder containing query files in FASTA format. '
             'If not provided, the query folder in the current working directory will be used'
    )
    allele_split_subparser.add_argument(
        '-o', '--output_path',
        metavar='output_path',
        default=os.path.join(os.getcwd(), 'split_alleles'),
        help='Specify name and path of folder into which the split allele files are to be written. '
             'If not provided, the split_alleles folder in the current working directory will be used'
    )
    allele_split_subparser.set_defaults(func=allele_split)
    # Get the arguments into an object
    arguments = setup_arguments(parser=parser)
    # Prevent the arguments being printed to the console (they are returned in order for the tests to work)
    sys.stderr = open(os.devnull, 'w', encoding='utf-8')
    return arguments


if __name__ == '__main__':
    cli()