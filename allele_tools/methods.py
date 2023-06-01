#!/usr/bin/env python

"""
Collection of methods for Allele finding
"""

# Standard imports
from csv import DictReader
from glob import glob
import logging
import json
import math
import os

# Third party inputs
from olctools.accessoryFunctions.accessoryFunctions import GenObject, \
    make_path, \
    MetadataObject, \
    relative_symlink
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastpCommandline
from Bio.Data.CodonTable import TranslationError
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import coloredlogs


def setup_logging(arguments):
    """
    Set the custom colour scheme and message format to used by coloredlogs
    :param arguments: type parsed ArgumentParser object
    """
    # Set up a dictionary of the default colour scheme, and font styles
    coloredlogs.DEFAULT_LEVEL_STYLES = {
        'debug': {
            'bold': True, 'color': 'green'},
        'info': {
            'bold': True, 'color': 'blue'},
        'warning': {
            'bold': True, 'color': 'yellow'},
        'error': {
            'bold': True, 'color': 'red'},
        'critical': {
            'bold': True, 'background': 'red'}

    }
    # Change the default log format to be the time prepended to the appropriately formatted 
    # message string
    coloredlogs.DEFAULT_LOG_FORMAT = '%(asctime)s %(message)s'
    # Set the logging level
    coloredlogs.install(level=arguments.verbosity.upper())


def setup_arguments(parser):
    """
    Finalise setting up the ArgumentParser arguments into an object, and running subparser
    functions, or displaying the help message
    :param parser: type: ArgumentParser object
    :return: parsed ArgumentParser object
    """
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Run the appropriate function for each sub-parser.
    if hasattr(arguments, 'func'):
        # Set up logging
        setup_logging(arguments=arguments)
        arguments.func(arguments)
    # If the 'func' attribute doesn't exist, display the basic help for the appropriate subparser
    # (if any)
    else:
        try:
            # Determine which subparser was called by extracting it from the arguments.
            # Note that this requires the use of the desc keyword when creating subparsers
            command = list(vars(arguments).keys())[0]
            # If the extracted command exists, use the command-specific subparser help
            if command:
                parser.parse_args([command, '-h'])
            # Otherwise, use the basic help
            else:
                parser.parse_args(['-h'])
        # If there were no subparsers specified (the list of keys in the arguments is empty),
        # use the basic help
        except IndexError:
            parser.parse_args(['-h'])
    return arguments


def common_allele_find_errors(args, errors, amino_acid):
    """
    Perform checks for arguments shared between allele finding scripts
    :param args: type ArgumentParser arguments
    :param errors: List of errors with supplied arguments
    :param amino_acid: Boolean of whether the query sequence is amino acid or nucleotide
    :return: Updated list of errors
    """
    # Query folder checks
    if not os.path.isdir(args.query_path):
        errors.append(f'Could not find supplied query folder: {args.query_path}')
    else:
        if not glob(os.path.join(args.query_path, '*.fasta')):
            errors.append(
                f'Could not locate sequence files in supplied query folder: {args.query_path}'
            )
        else:
            if amino_acid:
                errors = detect_protein(
                    query_path=args.query_path,
                    errors=errors
                )
    # Amino acid allele checks
    if not os.path.isdir(args.aa_alleles):
        errors.append(f'Could not find supplied amino acid allele folder: {args.aa_alleles}')
    else:
        if not glob(os.path.join(args.aa_alleles, '*.fasta')):
            errors.append(
                f'Could not locate sequence files in supplied amino acid allele folder: {args.aa_alleles}'
            )
    try:
        if not os.path.isfile(args.aa_profile):
            errors.append(f'Could not locate supplied amino acid profile file: {args.aa_profile}')
    except AttributeError:
        pass
    return errors


def detect_protein(query_path, errors):
    """
    Attempt to determine whether a supplied file contains protein sequence
    :param query_path: String of absolute path to folder containing sequence files
    :param errors: List of errors with supplied arguments
    :return: Updated list of errors
    """
    # Create a list of all the FASTA files in the query path
    seq_files = glob(os.path.join(query_path, '*.fasta'))
    # Iterate through all the files
    for seq_file in seq_files:
        # Initialise a boolean to track whether the sequences appear to be amino acid or nucleotide
        aa = False
        for record in SeqIO.parse(seq_file, format='fasta'):
            # Convert the sequence object to a string
            seq = str(record.seq)
            # Create a set of all the characters in the string of the sequence
            seq_set = set(seq)
            # Since there are only 4(5 with N, 6 with N and -) possible characters in DNA, but 20 in protein, I chose
            # a length of 10 to allow for relatively low complexity protein sequences to pass, but DNA to fail
            if len(seq_set) > 10:
                # Update the boolean - note that only a single sequence in the file needs to be considered protein for
                # the entire file to pass
                aa = True
        # Update the errors if the file appears to be DNA
        if not aa:
            errors.append(f'Query file {seq_file} does not appear to be protein')
    return errors


def pathfinder(path):
    """
    Create absolute path user-supplied path. Allows for tilde expansion from
    :param path: String of path supplied by user. Could be relative, tilde expansion, or absolute
    :return: out_path: String of absolute path provided by user.
    """
    # Determine if the path requires path expansion
    if path.startswith('~'):
        # Create the absolute path of the tilde expanded path
        out_path = os.path.abspath(os.path.expanduser(os.path.join(path)))
    else:
        # Create the absolute path from the path
        out_path = os.path.abspath(os.path.join(path))
    return out_path


def query_prep(query_path, runmetadata, clear_report=True):
    """
    Create metadata objects for each sample
    :param query_path: String of absolute path to folder containing sequence files
    :param runmetadata: List of metadata objects for each query
    :param clear_report: Boolean of whether to clear previous iterations of BLAST reports. Default is True
    :return runmetadata: Metadata object updated with query information
    """
    logging.info('Preparing query files')
    # Find all the sequence files in the path
    fasta_files = sorted(glob(os.path.join(query_path, '*.fasta')))
    for fasta in fasta_files:
        name = os.path.splitext(os.path.basename(fasta))[0]
        if name != 'combinedtargets':
            # Create a metadata object for each sample
            metadata = MetadataObject()
            metadata.samples = []
            # Populate the metadata object with the required attributes
            metadata.name = name
            metadata.general = GenObject()
            metadata.commands = GenObject()
            metadata.alleles = GenObject()
            metadata.alleles.outputdirectory = os.path.join(query_path, metadata.name)
            # Set the name of the BLAST output file
            metadata.alleles.blast_report = os.path.join(
                metadata.alleles.outputdirectory,
                f'{metadata.name}.tsv'
            )
            # As the name and number of alleles can change over multiple iterations of the
            # program, it's best to clear out any old reports before running again
            if clear_report:
                try:
                    os.remove(metadata.alleles.blast_report)
                except FileNotFoundError:
                    pass
            make_path(metadata.alleles.outputdirectory)
            # Create a symlink of the sequence file in its own subdirectory
            metadata.general.bestassemblyfile = relative_symlink(
                src_file=fasta,
                output_dir=metadata.alleles.outputdirectory,
                export_output=True
            )
            metadata.samples.append(metadata)
            runmetadata.samples.append(metadata)
    return runmetadata


def blast_alleles(runmetadata, amino_acid, combined_targets, cpus, outfmt):
    """
    Run the BLAST analyses on the query
    :param runmetadata: List of metadata objects for each query
    :param amino_acid: Boolean of whether the query sequence is amino acid or nucleotide
    :param combined_targets: String of absolute path to file containing all sequences in other files in folder
    :param cpus: Integer of number of threads to use for BLAST analyses
    :param outfmt: String of BLAST fields to include in the report
    """
    logging.info('Running BLAST analyses')
    # Iterate through the samples
    for sample in runmetadata.samples:
        # Run the appropriate BLAST command BLASTN for nucleotide, BLASTP for protein
        if not amino_acid:
            blast = NcbiblastnCommandline(
                db=os.path.splitext(combined_targets)[0],
                query=sample.general.bestassemblyfile,
                num_alignments=100000000,
                evalue=0.001,
                num_threads=cpus,
                task='blastn',
                outfmt=outfmt,
                out=sample.alleles.blast_report
            )
        else:
            blast = NcbiblastpCommandline(
                query=sample.general.bestassemblyfile,
                db=os.path.splitext(combined_targets)[0],
                evalue=0.001,
                num_alignments=100000000,
                num_threads=cpus,
                outfmt=outfmt,
                out=sample.alleles.blast_report
            )
        blast()


def create_gene_names(path=os.getcwd(), name='genes.txt'):
    """
    Create a file with gene names to use in reducing a wgMLST profile by finding any .fasta files
    in a folder and adding them to the file (one per line)
    :param path: type: String of the path in which the genes file is to be created
    :param name: type: String of the name of the gene file
    """
    # Find all the .fasta files in the path
    fasta_files = glob(os.path.join(path, '*.fasta'))
    # Set the name and path of the gene file
    gene_file = os.path.join(path, name)
    # Open the gene file to write
    with open(gene_file, 'w', encoding='utf-8') as genes_file:
        for fasta in fasta_files:
            # Remove the path information and the file extension. Print to file
            genes_file.write(os.path.basename(os.path.splitext(fasta)[0]) + '\n')
            logging.debug(
                'Adding %s to %s',
                os.path.basename(os.path.splitext(fasta)[0]), gene_file
            )
    # Check to see if the file is empty
    if os.stat(gene_file).st_size == 0:
        # Log an error stating that the file could not be properly populated
        logging.error(
            'Created gene file, %s, is empty. Please ensure that directory %s has files with '
            '.fasta extensions', gene_file, path
            )
        raise SystemExit


def create_blast_dict(sample, extended_fieldnames):
    """
    Use DictReader to open and read a BLAST report into a dictionary
    :param sample: Metadata object
    :param extended_fieldnames: List of the BLAST fields used, as well as the additional percent
    match in index 14
    """
    # Open the BLAST report as a dictionary
    blastdict = DictReader(
        open(sample.alleles.blast_report,
             encoding='utf-8'),
        fieldnames=extended_fieldnames,
        dialect='excel-tab'
    )
    return blastdict


def parse_colocated_results(runmetadata, fieldnames, extended_fieldnames, amino_acid, gene_names, nt_allele_path,
                            aa_allele_path, report_path, overlap_range=50, cutoff=90):
    """
    Parse BLAST outputs. Ensure co-location of genes that must be co-located
    :param runmetadata: List of metadata objects
    :param fieldnames: List of the BLAST fields used
    :param extended_fieldnames: List of the BLAST fields used, as well as the additional percent
    match in index 14
    :param amino_acid: Variable on whether targets are protein
    :param gene_names: List of all gene names in the analysis
    :param nt_allele_path: String of the absolute path to the folder containing nucleotide allele files
    :param aa_allele_path: String of the absolute path to the folder containing amino acid allele files
    :param report_path: String of the absolute path to the folder into which reports are to be written
    :param overlap_range: Integer of the maximum distance allowed between two genes in order for
    them to be considered co-located. Default is 50 bp
    :param cutoff: Integer of the minimum percent identity between query and subject sequence.
    Default is 100%
    :return: runmetadata: Updated list of metadata objects
    :return: notes: List of contig:query_range-specific notes
    """
    logging.info('Parsing BLAST outputs')
    notes = {}
    for sample in runmetadata.samples:
        # Initialise GenObjects as required
        sample.alleles.blastlist = []
        sample.alleles.targetsequence = {}
        # Read the BLAST outputs into a dictionary
        blastdict = create_blast_dict(
            sample=sample,
            extended_fieldnames=extended_fieldnames
        )
        # Initialise dictionaries to store parsed outputs
        resultdict = {}
        colocation_dict = {}
        processed_range_dict = {}
        # Go through each BLAST result
        for row in blastdict:
            # Ignore the headers
            if row['query_id'].startswith(fieldnames[0]):
                continue
            target_id = row['subject_id']
            target_start = row['subject_start']
            target_end = row['subject_end']
            target_seq = row['query_sequence']
            high = max([int(row['query_start']), int(row['query_end'])])
            low = min([int(row['query_start']), int(row['query_end'])])
            # Create a list of the properly ordered start and stop points of the match
            query_range = [low, high]
            # Remove unwanted pipes added to the name
            nt_allele = target_id.lstrip('gb|').rstrip('|') if '|' in target_id else \
                target_id
            # If the percent identity is equal to the cutoff
            if float(row['percent_match']) >= cutoff:
                # Append the hit dictionary to the list
                sample.alleles.blastlist.append(row)
                # Determine if the orientation of the sequence is reversed compared to
                # the reference sequence
                if int(target_end) < int(target_start) and not amino_acid:
                    seq = Seq(target_seq)
                    # Calculate the reverse complement of the sequence
                    nt_querysequence = str(seq.reverse_complement())
                # If the sequence is not reversed, use the sequence as it is in the output
                else:
                    nt_querysequence = target_seq
                # Create a variable to avoid unnecessary typing
                contig = row['query_id']
                # Create a boolean to track whether this contig:range combination has already been processed
                processed = False
                # Add the contig key to the dictionary as required
                if contig not in processed_range_dict:
                    processed_range_dict[contig] = set()
                # Check the processed range dictionary to see if the current range is present
                if processed_range_dict[contig]:
                    for previous_range in processed_range_dict[contig]:
                        # Allow a small overlap of five bases in case the range of one query is slightly different
                        overlap = query_range[1] + 5 >= previous_range[0] and \
                                  previous_range[1] + 5 >= query_range[0]
                        # If the range is already present in the dictionary, update the tracking boolean
                        if overlap:
                            processed = True
                # Add the range to the set if it is empty
                else:
                    processed_range_dict[contig].add(tuple(query_range))
                # If range has already been processed, we can skip this iteration of it
                if processed:
                    continue
                # Update the processed ranges dictionary with the current range
                processed_range_dict[contig].add(tuple(query_range))
                # Create a tuple of the query range list to allow it to be used as a dictionary key
                query_range_tuple = tuple(query_range)
                # Add keys to the targetsequence dictionary as required
                if contig not in sample.alleles.targetsequence:
                    sample.alleles.targetsequence[contig] = {}
                if query_range_tuple not in sample.alleles.targetsequence[contig]:
                    sample.alleles.targetsequence[contig][query_range_tuple] = {}
                # Populate the percent match dictionary
                if contig not in resultdict:
                    resultdict[contig] = {}
                if query_range_tuple not in resultdict[contig]:
                    resultdict[contig][query_range_tuple] = {nt_allele: row['percent_match']}
                # Populate the notes dictionary
                if contig not in notes:
                    notes[contig] = {}
                if query_range_tuple not in notes[contig]:
                    notes[contig][query_range_tuple] = []
                # Determine the name of the gene corresponding to the allele e.g. if the allele
                # ECs1206_138, the corresponding gene is ECs1206_
                base_gene = [gene_name for gene_name in gene_names if gene_name in nt_allele][0]
                # Translate the query sequence to protein
                aa_querysequence = translate_sequence(
                    nt_seq=nt_querysequence
                )
                # Find the amino acid allele corresponding to this sequence
                returned = aa_allele_lookup(
                    aa_seq=aa_querysequence,
                    gene=base_gene,
                    aa_allele_path=aa_allele_path,
                    notes=notes[contig][query_range_tuple]
                )
                # Initialise a string to hold the amino acid allele identifier
                aa_allele = str()
                # Create a boolean to track whether the amino acid sequence fails the screen criteria
                filtered = False
                # If a perfect match to a previous allele was found, a string is returned
                if isinstance(returned, str):
                    aa_allele = returned
                # If not perfect match was found to a previous allele, a tuple of the updated filtered boolean, and
                # notes on the sequence are returned
                else:
                    filtered, notes[contig][query_range_tuple] = returned
                # Add unfiltered imperfect nt alleles to the database
                if float(row['percent_match']) < 100 and not filtered:
                    # Find the next allele identifier for the database
                    nt_allele_id = find_next_allele(
                        gene=base_gene,
                        allele_path=nt_allele_path
                    )
                    # Add the base gene name to the allele identifier
                    nt_allele = f'{base_gene}_{nt_allele_id}'
                    # Update the allele database with the new allele
                    notes[contig][query_range_tuple] = update_allele_databases(
                        query_sequence=nt_querysequence,
                        header=nt_allele,
                        filtered=filtered,
                        gene=base_gene,
                        report_path=report_path,
                        allele_path=nt_allele_path,
                        notes=notes[contig][query_range_tuple],
                        molecule='Nucleotide'
                    )
                # Add unfiltered novel aa alleles to the database
                if not aa_allele and not filtered:
                    aa_allele_id = find_next_allele(
                        gene=base_gene,
                        allele_path=aa_allele_path
                    )
                    aa_allele = f'{base_gene}_{aa_allele_id}'
                    notes[contig][query_range_tuple] = update_allele_databases(
                        query_sequence=aa_querysequence,
                        header=aa_allele,
                        filtered=filtered,
                        gene=base_gene,
                        report_path=report_path,
                        allele_path=aa_allele_path,
                        notes=notes[contig][query_range_tuple],
                        molecule='Amino acid'
                    )
                # Populate the targetsequence dictionary with information on the nt and aa alleles
                if base_gene not in sample.alleles.targetsequence[contig][query_range_tuple]:
                    sample.alleles.targetsequence[contig][query_range_tuple][base_gene] = {
                        'nt': {
                            'allele': nt_allele,
                            'sequence': nt_querysequence
                        },
                        'aa': {
                            'allele': aa_allele,
                            'sequence': aa_querysequence
                        }
                    }
                # Populate the co-location dictionary with the required keys as necessary
                if contig not in colocation_dict:
                    colocation_dict[contig] = {}
                # The query_ranges and target keys both correspond to lists of values
                if 'query_ranges' not in colocation_dict[contig]:
                    colocation_dict[contig] = {
                        'query_ranges': [query_range],
                        'target': [nt_allele]
                    }
                # If the keys already exist, append to the lists
                else:
                    colocation_dict[contig]['query_ranges'].append(query_range)
                    colocation_dict[contig]['target'].append(nt_allele)
            # Store the BLAST outputs in the metadata object
            sample.alleles.blastresults = resultdict
        # Populate missing results with 'NA' values
        if len(resultdict) == 0:
            sample.alleles.blastresults = 'NA'
        sample.alleles.overlap_dict = colocation_calculation(
            colocation_dict=colocation_dict,
            gene_names=gene_names,
            overlap_range=overlap_range
        )
    return runmetadata, notes


def translate_sequence(nt_seq):
    """
    Uses BioPython to translate a nucleotide sequence to protein, and trims it to the first stop
    codon
    :param nt_seq: String of the nucleotide sequence
    :return aa_seq: String of the trimmed amino acid sequence
    """
    # Create a sequence object from the nucleotide sequence
    nt_seq_object = Seq(nt_seq)
    # Translate the sequence to protein
    try:
        # Translate the sequence
        aa_seq_object = nt_seq_object.translate()
    # BioPython cannot translate a sequence with gaps (-)
    except TranslationError:
        allele_seq = str(nt_seq_object).replace('-', '')
        seq = Seq(allele_seq)
        aa_seq_object = str(seq.translate())
    # Split the sting on stop codons, keep only the first part of the split
    aa_seq = str(aa_seq_object).split('*', maxsplit=1)[0] + '*'
    return str(aa_seq)


def aa_allele_lookup(aa_seq, gene, aa_allele_path, notes):
    """
    Read in the amino acid allele file. Search for exact matches to the current sequence
    :param aa_seq: String of the amino acid sequence
    :param gene: Sting of the gene name
    :param aa_allele_path: String of the absolute path to the folder containing the amino acid
    allele files
    :param notes: List of notes for the current contig: query_range
    :return record.id: Allele identifier corresponding to the sequence matching the aa_seq
    if no matches:
    :return filtered: Boolean of whether the amino acid sequence passes length thresholds
    :return notes: Populated notes
    """
    # Set the name of the amino acid allele file by joining the allele folder path to the gene name
    aa_allele_file = os.path.join(aa_allele_path, f'{gene}.fasta')
    # Iterate through all the alleles in the file
    for record in SeqIO.parse(aa_allele_file, 'fasta'):
        # If the sequence in the file matches the current sequence, return the allele identifier
        if aa_seq == str(record.seq):
            return record.id
    # If no records match, evaluate whether the aa allele passes necessary length thresholds
    filtered, notes = evaluate_translated_allele(
        aa_seq=aa_seq,
        gene=gene,
        notes=notes
    )
    return filtered, notes


def evaluate_translated_allele(aa_seq, gene, notes, aa=False):
    """
    Evaluate whether an aa sequence passes the necessary length thresholds after trimming of an interior stop codons
    :param aa_seq: String of the amino acid sequence to evaluate
    :param gene: String of the name of the gene (no allele information) being evaluated
    :param notes: List of notes for the current contig: query_range
    :param aa: Boolean of whether the query sequence is amino acid. Triggers filtering if sequence doesn't end with a
    stop codon
    :return filtered: Boolean of whether the amino acid sequence passes length thresholds
    :return notes: Populated notes
    """
    # Dictionary of minimum acceptable lengths for each of the STEC genes
    length_dict = {
        'ECs2973': 90,
        'ECs2974': 316,
        'ECs1205': 316,
        'ECs1206': 88
    }
    filtered = False
    if not aa_seq.endswith('*'):
        notes.append(f'{gene} trimmed sequence did not end with a stop codon')
        if aa:
            filtered = True
    # Remove all sequence after a stop codon (*)
    aa_seq = aa_seq.split('*', maxsplit=1)[0] + '*'
    # Evaluate the translated length of the sequence
    filtered, notes = evaluate_translated_length(
        aa_seq=aa_seq,
        length_dict=length_dict,
        gene=gene,
        notes=notes,
        filtered=filtered
    )
    return filtered, notes


def update_allele_databases(query_sequence, header, filtered, gene, report_path, allele_path, notes, molecule):
    """
    Update the appropriate allele file depending on quality filter status and molecule
    :param query_sequence: SEQIO sequence object of the novel allele
    :param header: String of the allele name (gene_allele ID)
    :param filtered: Boolean of whether the allele has been quality filtered
    :param gene: String of the name of the gene
    :param report_path: String of the absolute path to the folder into which the reports are to be written
    :param allele_path: String of the absolute path to the folder containing the allele database
    :param notes: List of notes on the alleles
    :param molecule: String of the current molecule. Options are Nucleotide and Amino acid
    :return: notes: Updated list of notes
    """
    # Create a SeqRecord of the allele using the novel allele name and sequence
    new_record = SeqRecord(seq=Seq(query_sequence),
                           id=header,
                           name='',
                           description='')
    # Create a string to prepend to allele file names
    molecule_str = 'nt' if molecule == 'Nucleotide' else 'aa'
    # Set the correct files depending on the filtering status
    if not filtered:
        new_alleles = os.path.join(report_path, f'{molecule_str}_{gene}_novel_alleles.fasta')
        allele_file = os.path.join(allele_path, f'{gene}.fasta')
    else:
        new_alleles = os.path.join(report_path, f'{molecule_str}_{gene}_filtered_alleles.fasta')
        allele_file = os.path.join(allele_path, f'{gene}_filtered.txt')
    records = []
    # Iterate through all the records in the allele database
    if os.path.isfile(allele_file):
        for record in SeqIO.parse(allele_file, 'fasta'):
            # Append all the records to the list
            records.append(record)
    # Check to see if the query sequence is novel in the database
    if query_sequence not in [str(seq.seq) for seq in records]:
        # Append the SeqRecord to the novel alleles file
        with open(new_alleles, 'a+', encoding='utf-8') as novel:
            SeqIO.write(sequences=new_record,
                        handle=novel,
                        format='fasta')
        records.append(new_record)
        # Overwrite the existing allele database file with the updated list of records
        with open(allele_file, 'w', encoding='utf-8') as alleles:
            SeqIO.write(sequences=records,
                        handle=alleles,
                        format='fasta')
        remove_combined_db_files(allele_path=allele_path)
        notes.append(f'{molecule} allele {header} is novel')
    # Non-novel sequences will have updated notes with the match
    else:
        for record in records:
            if str(query_sequence) == record.seq:
                # Append the previous finding to the notes
                notes.append(f'{molecule} matches previous result: {record.id}')
    return notes


def colocation_calculation(colocation_dict, gene_names, overlap_range):
    """
    Determine if gene results are co-located on a contig
    :param colocation_dict: Dictionary of contig: {'query_ranges': [query_range],
    'target': [allele_id]}
    :param gene_names: List of all genes in the analysis
    :param overlap_range: Integer of the maximum distance allowed between two separate hits before
    they can no longer be considered co-located on a contig
    :return overlap_dict: Dictionary of contig:full_range:gene_pair: {'overlap': overlap,
    'allele': [allele_identifiers]}
    """
    # Initialise a dictionary to store the details of any co-located sequences
    overlap_dict = {}
    # Iterate over all the contigs with hits
    for contig, info_dict in colocation_dict.items():
        # Update the overlap dictionary with the contig name as required
        if contig not in overlap_dict:
            overlap_dict[contig] = {}
        # Extract the query range and the allele identifiers from info_dict
        query_ranges = info_dict['query_ranges']
        targets = info_dict['target']
        # Iterate over all the query ranges with hits on the current contig
        for query_iterator, query_range in enumerate(query_ranges):
            # Create a variable to track whether the current contig:query range combination has
            # been added to the overlap dictionary
            processed = False
            # Extract the name of the current allele from the list of all alleles in the range
            current_allele = targets[query_iterator]
            # Create a dictionary of tuple of other ranges present in the list of ranges: iterator
            other_ranges = {
                tuple(other): j for j, other in enumerate(query_ranges) if other != query_range
            }
            # Iterate over these other ranges
            for other_range, other_iterator in other_ranges.items():
                # Calculate whether the current range overlaps with this other range
                # e.g. query_range = (100, 500), other_range = (525, 1000), overlap_range = 50
                # Check if 500 + 50 >= 525 and 1000 + 50 >= 100
                overlap = query_range[1] + overlap_range >= other_range[0] and \
                        other_range[1] + overlap_range >= query_range[0]
                # If these ranges overlap, populate the overlap dictionary
                if overlap:
                    overlap_dict, processed = positive_overlap(
                        info_dict=info_dict,
                        other_iterator=other_iterator,
                        query_range=query_range,
                        other_range=other_range,
                        overlap_dict=overlap_dict,
                        current_allele=current_allele,
                        contig=contig,
                        gene_names=gene_names,
                        overlap=overlap
                    )
            # If the current contig: range was not entered into the overlap dictionary, there were
            # either no other hits, or the hits did not overlap
            if not processed:
                # Create a tuple containing only the current range
                tuple_range = tuple(query_range)
                # Update the dictionary as required
                if tuple_range not in overlap_dict[contig]:
                    overlap_dict[contig][tuple_range] = {}
                # Extract the gene name corresponding to the allele identifier
                # e.g. gene = ECs1206 allele = ECs1206_138 will create ECs1206
                gene = [gene_name for gene_name in gene_names if gene_name in current_allele][0]
                # Add the gene name to the dictionary, and update create the overlap and allele keys
                if gene not in overlap_dict[contig][tuple_range]:
                    overlap_dict[contig][tuple_range][gene] = {
                        'overlap': False,
                        'allele': []
                    }
                # Append the current allele identifier to the list of alleles
                overlap_dict[contig][tuple_range][gene]['allele'].append(current_allele)
    return overlap_dict


def positive_overlap(
        info_dict, other_iterator, query_range, other_range, overlap_dict,
        current_allele, contig, gene_names, overlap):
    """
    Determine the combined range of two overlapping ranges, extract gene names corresponding to
    allele names, populate dictionary of range overlaps
    :param info_dict: Dictionary of {'query_ranges': [query_range], 'target': [allele_id]}
    :param other_iterator: Integer of the iterator corresponding to the current other_range from the
    dictionary of other_range: iterator
    :param query_range: Range of hit corresponding to current_allele in info_dict
    e.g. info_dict['query_ranges'][query_iterator] and info_dict['target'][query_iterator]
    :param other_range: Range of hit corresponding to non-current allele in info_dict
    e.g. info_dict['query_ranges'][other_iterator] and info_dict['target'][other_iterator]
    :param overlap_dict: Dictionary of to be populated with overlap information
    e.g. contig:full_range:gene_pair: {'overlap': overlap, 'allele': [allele_identifiers]}
    :param current_allele: String of the name of the allele extracted from info_dict
    e.g. info_dict['target'][query_iterator]
    :param contig: Name of the contig within which the BLAST hits are located
    :param gene_names: List of all gene names in the analysis
    :param overlap: Boolean of whether query_range overlaps with other_range
    """
    # Extract the name of the other allele from info_dict using the iterator of the other range
    other_allele = info_dict['target'][other_iterator]
    full_range = calculate_full_range(
        query_range=query_range,
        other_range=other_range
    )
    # Update the overlap dictionary with the full range as required
    if full_range not in overlap_dict[contig]:
        overlap_dict[contig][full_range] = {}
    # Create a sorted tuple of the allele names
    alleles = tuple(sorted([current_allele, other_allele]))
    # Create a set of all the genes with hits in the current overlap
    genes = set()
    # Iterate over all the alleles
    for allele in alleles:
        # Add the gene from the list of genes if it is present in the allele
        # identifier e.g. gene = ECs1206 allele = ECs1206_138 will add ECs1206
        genes.add([gene for gene in gene_names if gene in allele][0])
    # Create a tuple of the sorted list of genes present in the set
    gene_pair = tuple(sorted(list(genes)))
    # Update the dictionary as required
    if gene_pair not in overlap_dict[contig][full_range]:
        overlap_dict[contig][full_range][gene_pair] = {
            'overlap': overlap,
            'allele': []
        }
    # Append the current allele to the overlap dictionary
    overlap_dict[contig][full_range][gene_pair]['allele'].append(current_allele)
    # Set processed to True to indicate that there was an overlap and that the
    # dictionary was populated
    processed = True
    return overlap_dict, processed


def calculate_full_range(query_range, other_range):
    """
    Determine if two ranges overlap
    :param query_range: Range of hit corresponding to current_allele
    :param other_range: Range of hit corresponding to non-current allele
    :return full_range: Tuple of minimum coordinate from both ranges, maximum coordinate from
    both ranges
    """
    # Determine in the minimum and maximum coordinates of the two ranges
    # e.g. query_range = (100, 500), other_range = (525, 1000)
    # min_range = 100
    min_range = (min(query_range[0], other_range[0]))
    # max_range = 1000
    max_range = (max(query_range[1], other_range[1]))
    # The full range is a tuple of (min_range, max_range)
    # full_range = (100, 1000)
    full_range = tuple(sorted([min_range, max_range]))
    return full_range


def evaluate_translated_length(aa_seq, length_dict, gene, notes, filtered):
    """
    Evaluate whether a translated sequence passes a length filter and starts with a methionine
    residue
    :param aa_seq: String of the amino acid sequence to evaluate
    :param length_dict: Dictionary of minimum acceptable length for each gene in the analysis
    :param gene: String of the name of the gene
    :param notes: List of notes on the alleles
    :param filtered: Boolean of whether the allele has been filtered based on length or content
    :return: filtered: Updated filtering boolean
    :return: notes: Updated list of notes
    """
    # Proper protein sequences must start with a methionine (M)
    if not aa_seq.startswith('M'):
        filtered = True
        notes.append(f'{gene} amino acid sequence does not start with M')
    # The length of the sequence must also be greater than the minimum gene-specific length
    if len(aa_seq) < length_dict[gene]:
        filtered = True
        notes.append(
            f'{gene} amino acid sequence was {len(aa_seq)} amino acid residues. Minimum allowed '
            f'length is {length_dict[gene]} amino acid residues')
    return filtered, notes


def generic_evaluate_translated_length(aa_seq, sequence, gene, notes, filtered, cutoff=0.95):
    """
    Evaluate whether a translated sequence passes a generic length filter and starts with a methionine
    residue
    :param aa_seq: String of the amino acid sequence to evaluate
    :param sequence: String of untrimmed nucleotide sequence
    :param gene: String of the name of the gene
    :param notes: List of notes on the alleles
    :param filtered: Boolean of whether the allele has been filtered based on length or content
    :param cutoff: Float of minimum cutoff value to be used for filtering trimmed sequences. Default is 0.95
    :return: filtered: Updated filtering boolean
    :return: notes: Updated list of notes
    """
    # Proper protein sequences must start with a methionine (M)
    if not aa_seq.startswith('M'):
        filtered = True
        notes.append(f'{gene} amino acid sequence does not start with M')
    # Minimum length of a trimmed amino acid allele permitted is 95% the length of the theoretical length of
    # the translated nucleotide sequence e.g. a 99 bp nt sequence would be 33 amino acid residues, and 95% of
    # that is 31.35 -> 31 (rounded down)
    minimum_length = math.floor(len(sequence) / 3 * cutoff)
    aa_seq_length = len(aa_seq)
    if aa_seq_length < minimum_length:
        filtered = True
        notes.append(
            f'{gene} amino acid sequence was trimmed to {aa_seq_length} residues '
            f'the minimum length allowed is {minimum_length} residues')
    return filtered, notes


def find_next_allele(gene, allele_path, extension='.fasta'):
    """
    Update the allele database with the novel allele extracted above
    :param gene: Name of the current gene being examined
    :param allele_path: Name and absolute path to folder containing allele files
    :param extension: String of the file extension. Default is .fasta
    :return: last_id: Number of the last alleles in the current database
    """
    # Find the allele database file
    allele_file = os.path.join(allele_path, f'{gene}{extension}')
    # Initialise a variable to store the name of the last allele in the database file
    last_id = int()
    records = []
    if os.path.isfile(allele_file):
        # Iterate through all the records in the allele database
        for record in SeqIO.parse(allele_file, 'fasta'):
            # Update the last_id variable
            last_id = int(record.id.split('_')[-1])
            records.append(record)
    else:
        last_id = 0
    # Make it clear that these are novel profiles by starting at 1000000
    if last_id < 1000000:
        last_id = 999999
    return last_id + 1


def remove_combined_db_files(allele_path):
    """
    Remove all the combined gene files used in BLAST analyses
    :param allele_path: String of the absolute path to the folder containing the alleles
    """
    # Find all the files in the directory with the word combined in the name
    combined_files = glob(os.path.join(allele_path, 'combined*'))
    # Remove each of the files
    for file in combined_files:
        os.remove(file)


def create_nt_allele_comprehension(runmetadata, gene_names):
    """
    Create gene: nucleotide allele ID comprehensions for each contig: range combination with hits
    :param runmetadata: List of metadata objects
    :param gene_names: List of all gene names in the analysis
    :return: allele_comprehension: nucleotide allele comprehension. allele_comprehension[contig][full_range] =
    {gene:allele}
    """
    logging.info('Determining nucleotide allele profiles')
    # Initialise an empty allele comprehension dictionary
    allele_comprehension = {}
    # Iterate through all the samples
    for sample in runmetadata.samples:
        # Extract hit information from the overlap dictionary
        for contig, range_dict in sample.alleles.overlap_dict.items():
            # Update the allele comprehension dictionary with the contig key as required
            if contig not in allele_comprehension:
                allele_comprehension[contig] = {}
            # Iterate through each query range with a hit in the current contig
            for query_range, gene_dict in range_dict.items():
                # Update the dictionary with the query range key as required
                if query_range not in allele_comprehension[contig]:
                    allele_comprehension[contig][query_range] = {}
                # Iterate over each gene with a hit within this contig:query_range combination
                for gene_pair, info_dict in gene_dict.items():
                    # If the gene_pair variable is a string (instead of a tuple), there is only a single gene present
                    if isinstance(gene_pair, str):
                        # Extract the gene_allele ID from the dictionary
                        corresponding_allele = info_dict['allele'][0]
                        # Remove the gene name (and an underscore) from the corresponding_allele variable
                        # (leaving only the allele ID)
                        allele_number = corresponding_allele.replace(f'{gene_pair}_', '')
                        # Update the allele comprehension dictionary with the gene name: alleleID
                        allele_comprehension[contig][query_range].update(
                            {gene_pair: allele_number}
                        )
                        # Determine which gene(s) are missing from this contig:query_range
                        missing = [other_gene for other_gene in gene_names if other_gene != gene_pair]
                        # Populate the allele comprehension dictionary with the missing genes
                        for missing_gene in missing:
                            allele_comprehension[contig][query_range].update({missing_gene: '0'})
                    # A tuple of gene names indicates that multiple co-located genes are present in this
                    # contig:query_range combination
                    else:
                        # Iterate over each gene in the gene_pair tuple
                        for i, gene_name in enumerate(gene_pair):
                            # Extract the gene_alleleID from the dictionary for this gene
                            corresponding_allele = info_dict['allele'][i]
                            # Remove the gene name information from the corresponding_allele variable
                            allele_number = corresponding_allele.replace(f'{gene_name}_', '')
                            # Update the dictionary with the new gene: allele number for the sample
                            allele_comprehension[contig][query_range].update(
                                {gene_name: allele_number}
                            )
        # If the allele_comprehension dictionary exists, it doesn't need to be further populated
        if allele_comprehension:
            continue
        # Otherwise iterate through the targetsequence dictionary
        for contig, range_dict in sample.alleles.targetsequence.items():
            # Update the allele comprehension dictionary as required
            if contig not in allele_comprehension:
                allele_comprehension[contig] = {}
            # Iterate through all genes in the analysis
            for gene in gene_names:
                # Set an 'empty' range as (0, 0)
                full_range = (0, 0)
                # Add the range to the dictionary as required
                if full_range not in allele_comprehension[contig]:
                    allele_comprehension[contig][full_range] = {}
                # Update the dictionary with the negative result
                allele_comprehension[contig][full_range].update(
                                {gene: '0'}
                            )
    return allele_comprehension


def create_aa_allele_comprehension(runmetadata, gene_names):
    """
    Create gene: amino acid allele ID comprehensions for each contig: range combination with hits
    :param runmetadata: List of metadata objects
    :param gene_names: List of all gene names in the analysis
    :return: allele_comprehension: amino acid allele comprehension. allele_comprehension[contig][full_range] =
    {gene:allele}
    """
    logging.info('Determining amino acid allele profiles')
    # Initialise a dictionary to store contig:query_range gene:alleleID results
    allele_comprehension = {}
    # Iterate through all the samples
    for sample in runmetadata.samples:
        # Iterate through all the contigs in the targetsequence dictionary
        for contig, range_dict in sample.alleles.targetsequence.items():
            # Update the dictionary as required
            if contig not in allele_comprehension:
                allele_comprehension[contig] = {}
            # If the current contig is not in the overlap dictionary, populate the allele comprehension with
            # negative values
            if contig not in sample.alleles.overlap_dict:
                # Iterate over every gene in the analysis
                for gene in gene_names:
                    # Set the 'empty' value to (0, 0)
                    full_range = (0, 0)
                    # Update the dictionary with the negative values
                    if full_range not in allele_comprehension[contig]:
                        allele_comprehension[contig][full_range] = {}
                    allele_comprehension[contig][full_range].update(
                                    {gene: '0'}
                                )
                # The dictionary has been populated, so continue
                continue
            # Extract all the ranges with hits in the overlap dictionary
            for full_range in sample.alleles.overlap_dict[contig]:
                # Extract the query ranges with hits on the current contig in the targetsequence dictionary
                for query_range, gene_dict in range_dict.items():
                    # Determine if these two ranges have an overlap
                    overlap = query_range[1] >= full_range[0] and full_range[1] >= query_range[0]
                    # If they do not overlap, they are not the same
                    if not overlap:
                        continue
                    # Update the dictionary as required
                    if full_range not in allele_comprehension[contig]:
                        allele_comprehension[contig][full_range] = {}
                    # Create a list to store all genes in the analysis that do not have hits in the range
                    missing_genes = []
                    # Iterate over all genes in the analysis
                    for gene in gene_names:
                        # If the gene is not in the gene dictionary contained in the targetsequence dictionary,
                        # add it to the list of missing genes
                        if gene not in gene_dict:
                            missing_genes.append(gene)
                        # Otherwise, update the dictionary
                        else:
                            # Extract the name of the gene:alleleID from gene_dict in the targetsequence dictionary
                            full_allele = gene_dict[gene]['aa']['allele']
                            # Remove the gene name information from the full_allele variable
                            allele_number = full_allele.replace(f'{gene}_', '')
                            # Add the gene:alleleID to the dictionary
                            allele_comprehension[contig][full_range].update(
                                {gene: allele_number}
                            )
                    # Add any missing genes to the dictionary with negative values
                    for gene in missing_genes:
                        # Ensure that the gene isn't already present in the dictionary
                        if gene not in allele_comprehension[contig][full_range]:
                            allele_comprehension[contig][full_range].update(
                                    {gene: '0'}
                                )
    return allele_comprehension


def create_frozen_allele_comprehension(allele_comprehension):
    """
    Freeze allele comprehension dictionaries
    :param allele_comprehension: Dictionary of contig:full_range: {gene:allele}
    :return: frozen_allele_comprehension: Dictionary of contig:query_range: json.dumps({gene:allele}, sort_keys=True)
    """
    # Initialise a dictionary to store the frozen allele comprehensions
    frozen_allele_comprehension = {}
    # Iterate over all the contigs with hits
    for contig, query_dict in allele_comprehension.items():
        # Update the frozen allele dictionary as required
        if contig not in frozen_allele_comprehension:
            frozen_allele_comprehension[contig] = {}
        # Iterate over all the ranges and allele comprehensions on the current contig with hits
        for query_range, allele_dict in query_dict.items():
            # Freeze the allele comprehension
            frozen_allele_dict = json.dumps(allele_dict, sort_keys=True)
            # Update the dictionary with the range and the frozen allele string
            if query_range not in frozen_allele_comprehension[contig]:
                frozen_allele_comprehension[contig][query_range] = frozen_allele_dict
    return frozen_allele_comprehension


def extract_novel_alleles(sample, gene, genome_query, amino_acid, allele_path, report_path,
                          cutoff=75):
    """
    Extract the sequence of novel alleles from samples that do not have a 100% match
    :param sample: Metadata object
    :param gene: Name of current gene
    :param genome_query: Boolean of whether the allele or the genome are the query
    :param amino_acid: Variable indicating whether the current analyses are on DNA or
    amino acid sequences
    :param allele_path: Name and absolute path to folder containing allele files
    :param report_path: Name and absolute path to folder in which reports are to be created
    :param cutoff: The minimum percent identity cutoff to allow when considering the presence of a
    sequence in a query
    :return: sample: Updated sample
    :return: novel_allele: Name of novel alleles discovered
    :return: query_sequence: Sequence of novel alleles discovered
    """
    # Open the sequence profile file as a dictionary
    blastdict = DictReader(open(sample.alleles.blast_report, encoding='utf-8'),
                           dialect='excel-tab')
    # Initialise the best hit value of 0
    best_hit = 0
    # Initialise strings to store the name and the sequence of novel alleles
    query_sequence = str()
    novel_allele = str()
    # Iterate through all the BLAST hits
    for row in blastdict:
        # Extract the target id with the appropriate key depending on whether genome files
        # are the query or the subject
        target_id = row['query_id'] if not genome_query else row['subject_id']
        # Ensure that the gene name is present in the gene name + allele combination
        if gene in target_id:
            # Create a variable to store the value for percent identity, so it is easier to call
            perc_id = float(row['percent_match'])
            # See if the percent identity for the current match is better than the previous best
            # match, and is above the
            # minimum cutoff threshold
            if perc_id > best_hit and perc_id >= cutoff:
                # Set the start and end variables depending on whether genomes are the query
                target_start = row['query_start'] if not genome_query else row['subject_start']
                target_end = row['query_end'] if not genome_query else row['subject_end']
                target_seq = row['query_sequence']
                # Determine if the orientation of the sequence is reversed compared to the reference
                if int(target_end) < int(target_start) and not amino_acid:
                    # Create a sequence object using BioPython
                    seq = Seq(target_seq)
                    # Calculate the reverse complement of the sequence
                    query_sequence = str(seq.reverse_complement())
                # If the sequence is not reversed, use the sequence as it is in the output
                else:
                    query_sequence = target_seq
                best_hit = perc_id

    # If a query sequence was extracted, use it to update the allele database
    if query_sequence:
        novel_allele = update_allele_database(
            gene=gene,
            query_sequence=query_sequence,
            allele_path=allele_path,
            report_path=report_path,
            amino_acid=amino_acid,
        )
    return sample, novel_allele, query_sequence


def update_allele_database(gene, query_sequence, allele_path, report_path, amino_acid):
    """
    Update the allele database with the novel allele extracted above
    :param gene: Name of the current gene being examined
    :param query_sequence: Sequence of the novel allele
    :param allele_path: Name and absolute path to folder containing allele files
    :param report_path: Name and absolute path to folder in which reports are to be created
    :param amino_acid: Variable indicating whether the current analyses are on DNA or
    amino acid sequences
    :return: novel_allele: Name of the novel allele entered into the database
    """
    # Find the allele database file
    allele_file = glob(os.path.join(allele_path, f'{gene}*.*fa*'))[0]
    # Set the appropriate molecule type based on the current analysis
    molecule = 'nt' if not amino_acid else 'aa'
    # Set the name of the novel allele file in the report path
    new_alleles = os.path.join(report_path, f'{molecule}_{gene}_novel_alleles.fasta')
    # Initialise a variable to store the name of the last allele in the database file
    last_id = str()
    # Create a list to store all the allele records in the database
    records = []
    # Iterate through all the records in the allele database
    for record in SeqIO.parse(allele_file, 'fasta'):
        # Add the records to the list
        records.append(record)
        # Update the last_id variable
        last_id = record.id
    # Try to separate the gene name from the allele e.g. MutS_1
    try:
        _, allele = last_id.rsplit('_', 1)
    # If there is no allele, set the allele to 1
    except ValueError:
        allele = 1
    # Typecase the variable to an integer
    allele = int(allele)
    # If the sequence type corresponds to an Enterobase number, use our local numbering scheme instead
    if allele < 1000000:
        allele = 999999
    # Name the novel allele as the gene name _ allele number + 1
    novel_allele = f'{gene}_{int(allele) + 1}'
    # Create a SeqRecord of the allele using the novel allele name and sequence
    new_record = SeqRecord(
        seq=Seq(query_sequence),
        id=novel_allele,
        name='',
        description=''
    )
    # Append the SeqRecord to the novel alleles file
    with open(new_alleles, 'a+', encoding='utf-8') as novel:
        SeqIO.write(
            sequences=new_record,
            handle=novel,
            format='fasta'
        )
    # Add the novel allele record to the list of all records
    records.append(new_record)
    # Overwrite the existing allele database file with the updated list of records
    with open(allele_file, 'w', encoding='utf-8') as alleles:
        SeqIO.write(
            sequences=records,
            handle=alleles,
            format='fasta'
        )
    return novel_allele


def translate(runmetadata):
    """
    Use BioPython to translate DNA to amino acid
    :param runmetadata: List of metadata objects for each query
    :return: Updated list of metadata objects
    """
    logging.info('Translating allele sequences to amino acid')
    for sample in runmetadata.samples:
        # Initialise the dictionary to store the translated sequence
        sample.alleles.nt_alleles_translated = {}
        for allele, allele_sequence_list in sample.alleles.targetsequence.items():
            for allele_sequence in allele_sequence_list:
                # Create a sequence object using BioPython
                seq = Seq(allele_sequence)
                try:
                    # Translate the sequence
                    aa_seq = str(seq.translate())
                # BioPython cannot translate sequences with gaps (-)
                except TranslationError:
                    # Remove all - from the sequence
                    allele_seq = allele_sequence.replace('-', '')
                    seq = Seq(allele_seq)
                    aa_seq = str(seq.translate())
                # Ensure that the allele name exists (isn't an empty string) before adding
                # allele name: translated sequence to the dictionary
                if allele:
                    sample.alleles.nt_alleles_translated[allele] = aa_seq
    return runmetadata


def match_profile(profile_data, frozen_allele_comprehension, report_path, profile_file, genes,
                  allele_comprehension, molecule):
    """
    Match current profiles to any previously created profiles
    :param profile_data: Dictionary of seq_type: {gene name: allele ID}
    :param frozen_allele_comprehension: Dictionary of json.dumps({gene name: allele ID}, sort_keys=True): seq_type
    :param report_path: Name and absolute path to folder in which reports are to be created
    :param profile_file: Name and path of file containing reduced profiles
    :param genes: List of all genes in the analysis
    :param allele_comprehension: Dictionary of contig:full_range: {gene:allele}
    :param molecule: String of the current molecule being processed. Options are "aa" and "nt"
    :return: profile_matches: Dictionary of contig:query_range:seq_type_match
    """
    # If the profile_data dictionary was not populated in the read_profiles methods,
    # there is nothing to match
    if not profile_data:
        return
    logging.info('Matching new %s profiles against profile file', molecule)
    profile_matches = {}
    # Extract all the profiles from the profile file (as a frozen string)
    frozen_profiles = freeze_profiles(
        profile_data=profile_data
    )
    # Iterate over all the contigs with hits
    for contig, query_dict in frozen_allele_comprehension.items():
        # Iterate over the query ranges with hits on the current contig
        for query_range, frozen_allele_dict in query_dict.items():
            try:
                # Extract the samples that match this profile
                seq_type_match = frozen_profiles[frozen_allele_dict]
                # Update the dictionary with the matching samples
                if contig not in profile_matches:
                    profile_matches[contig] = {}
                if query_range not in profile_matches[contig]:
                    profile_matches[contig][query_range] = seq_type_match
            # The profile will not necessarily match any of the profiles found in the analysis
            except KeyError:
                if contig not in profile_matches:
                    profile_matches[contig] = {}
                # Update the profile file with this novel profile
                profile_matches[contig][query_range] = update_profiles(
                    profile_file=profile_file,
                    report_path=report_path,
                    genes=genes,
                    allele_dict=allele_comprehension[contig][query_range],
                    molecule=molecule
                )
    return profile_matches, frozen_profiles


def freeze_profiles(profile_data):
    """
    Freeze profiles, so that the frozen {gene:allele} dictionary can be used as the key and the corresponding sequence
    type as the value
    :param profile_data: Dictionary of all profiles in seq_type: {gene name: allele ID} format
    :return: frozen_profiles: Dictionary of json.dumps({gene name: allele ID}, sort_keys=True): seq_type
    """
    # Initialise a dictionary to store the frozen profiles information
    frozen_profiles = {}
    # Iterate over all the sequence type: {gene name: allele ID} pairs in the dictionary
    for seq_type, allele_comprehension in profile_data.items():
        # Freeze the allele comprehension
        frozen_allele_comprehension = json.dumps(allele_comprehension, sort_keys=True)
        # Populate the dictionary with frozen_allele_comprehension: seq_type
        frozen_profiles[frozen_allele_comprehension] = seq_type
    return frozen_profiles


def update_profiles(profile_file, report_path, genes, allele_dict, molecule):
    """
    Run methods to add novel profiles to the profile file. Determine the sequence type to use, and update the file
    :param profile_file: Name and path of file containing reduced profiles
    :param report_path: Name and absolute path to folder in which reports are to be created
    :param genes: List of all genes in the analysis
    :param allele_dict: Dictionary of a single allele comprehension. Extracted from
    allele_comprehension[contig][query_range] = {gene: allele}
    :param molecule: String of the current molecule being processed. Options are "aa" and "nt"
    :return: next_seq_type: Integer of the sequence type assigned to the novel profile
    """
    # Extract the sequence type to use for the novel profile
    next_seq_type = return_next_seq_type(
        profile_file=profile_file
    )
    # Update the profile file
    update_profile_file(
        profile_file=profile_file,
        next_seq_type=next_seq_type,
        allele_dict=allele_dict,
        genes=genes,
        report_path=report_path,
        molecule=molecule
    )
    return next_seq_type


def return_next_seq_type(profile_file):
    """
    Parse the profile file, and return the value for the next sequence type to be used. Local profiles will start at
    1000000 in order to be distinct from Enterobase profiles
    :param profile_file: Name and path of file containing reduced profiles
    :return: last_seq_type + 1: Integer of the sequence type to be assigned to the novel profile
    """
    # Open the profile file
    with open(profile_file, 'r', encoding='utf-8') as profile:
        # Create a list of all the lines in the file
        lines = profile.read().splitlines()
        # Extract the last value from the list of lines
        last_line = lines[-1]
        # Split the line on tabs, and set the last_seq_type variable to the first entry e.g. 22\t3\t2\n yields a
        # sequence type of 22
        last_seq_type = last_line.split('\t')[0]
    # Typecase the variable to an integer
    int_last_seq_type = int(last_seq_type)
    # If the sequence type corresponds to an Enterobase number, use our local numbering scheme instead
    if int_last_seq_type < 1000000:
        int_last_seq_type = 999999
    # Return the last sequence type + 1 to give the next sequence type
    return int_last_seq_type + 1


def update_profile_file(profile_file, next_seq_type, allele_dict, genes, report_path, molecule):
    """
    Update the profile file with novel profile. Additionally, either create or update the novel_profiles.txt file
    with the same profile
    :param profile_file: Name and path of file containing reduced profiles
    :param next_seq_type: Integer of the sequence type to be assigned to the novel profile
    :param allele_dict: Dictionary of a single allele comprehension. Extracted from
    allele_comprehension[contig][query_range] = {gene: allele}
    :param genes: List of all genes in the analysis
    :param report_path: Name and absolute path to folder in which reports are to be created
    :param molecule: String of the current molecule being processed. Options are "aa" and "nt"
    """

    # Initialise a string to store the profile information with the novel sequence type
    seq_type_str = f'{next_seq_type}'
    # Initialise a header to store 'ST  gene1   gene2.......geneX\n'
    header = 'ST'
    # Iterate over all the genes in the analysis
    for gene in genes:
        # Extract the allele ID for each gene in the analysis
        allele = allele_dict[gene]
        # Update the header with the gene
        header += f'\t{gene}'
        # Update the profile string with the allele ID
        seq_type_str += f'\t{allele}'
    # Open the profile file (to update) and write the novel profile
    with open(profile_file, 'a+', encoding='utf-8') as profile:
        profile.write(seq_type_str + '\n')
    # Set the name of the file containing novel profiles using the molecule variable ('aa' or 'nt')
    novel_profile_file = os.path.join(report_path, f'{molecule}_novel_profiles.txt')
    # Check to see if the novel profile file exists
    if not os.path.isfile(novel_profile_file):
        # If it does not exist, create it, and write the header line before the novel profile
        with open(novel_profile_file, 'w', encoding='utf-8') as novel_profile:
            novel_profile.write(header + '\n')
            novel_profile.write(seq_type_str + '\n')
    # Otherwise, update the existing file with the novel profile
    else:
        with open(novel_profile_file, 'a+', encoding='utf-8') as novel_profile:
            novel_profile.write(seq_type_str + '\n')


def create_stec_report(runmetadata, nt_profile_matches, nt_alleles, aa_profile_matches,
                       aa_alleles, report_file, gene_names, aa_profile_path, notes):
    """
    Create a STEC-specific report including the allele matches for each gene and sequence type for both nucleotide and
    amino acid sequence information
    :param runmetadata: List of metadata objects for each query
    :param nt_profile_matches: Dictionary of contig:query_range:nucleotide seq_type_match
    :param nt_alleles: Dictionary of nucleotide contig:full_range:nucleotide seq_type_match
    :param aa_profile_matches: Dictionary of contig:query_range:amino acid seq_type_match
    :param aa_alleles: Dictionary of amino acid contig:full_range: {gene:allele}
    :param report_file: String of the name and path of the report file
    :param gene_names: List of all gene names in the analysis
    :param aa_profile_path: String of the absolute path of the folder in which the amino acid profile file is located
    :param notes: List of notes on the alleles
    """
    logging.info('Creating report')
    # Set the appropriate order for the genes in the report (stx1 genes are not in numerical order)
    gene_order = {
        'stx1': ['ECs2974', 'ECs2973'],
        'stx2': ['ECs1205', 'ECs1206']
    }
    # Create a list to store the ordered genes
    ordered_genes = []
    # Set the ordered genes according to the genes used in the current analysis (stx1 or stx2)
    for _, gene_list in gene_order.items():
        # If the sorted list matches the list of genes in the analysis, use the unsorted list as the gene order
        if sorted(gene_list) == gene_names:
            ordered_genes = gene_list
    # Create a header for the report. Includes which alleles are present and the sequence type for both the nucleotide
    # and amino acid sequences of the query and notes
    header = f'Sample\tnt_{ordered_genes[0]}\tnt_{ordered_genes[1]}\tnt_seq_type\t' \
             f'aa_{ordered_genes[0]}\taa_{ordered_genes[1]}\taa_seq_type\tnotes\n'
    # Create a string to store the query information
    data = str()
    # Iterate over the samples
    for sample in runmetadata.samples:
        # Iterate over all the contigs that had hits
        for contig, range_dict in nt_profile_matches.items():
            # Iterate over the ranges: nucleotide profiles that had hits on this contig
            for query_range, nt_profile in range_dict.items():
                # Update the data string with the sample name
                data += f'{sample.name}\t'
                # Extract the corresponding amino acid profile from aa_profile_matches
                aa_profile = aa_profile_matches[contig][query_range]
                # Extract the allele dictionaries ({gene name: allele ID}) using the contig and query range
                nt_allele_dict = nt_alleles[contig][query_range]
                aa_allele_dict = aa_alleles[contig][query_range]
                # Iterate over the genes in the analysis to extract their corresponding nucleotide alleles
                for gene in ordered_genes:
                    # Update the string with the nucleotide allele ID
                    data += f'{nt_allele_dict[gene]}\t'
                # Update the string with the nucleotide sequence type
                data += f'{nt_profile}\t'
                # Iterate over the genes in the analysis to extract their corresponding amino acid alleles
                for gene in ordered_genes:
                    # Update the string with the amino acid allele ID
                    data += f'{aa_allele_dict[gene]}\t'
                # Update the string with the amino acid sequence type
                data += f'{aa_profile}\t'
                # Create a list to store sample:contig:query_range-specific notes
                note_list = []
                # Determine if there are already notes for this contig in the notes dictionary
                if contig in notes:
                    # Determine if there are already notes for this contig: range in the notes dictionary
                    if query_range in notes[contig]:
                        # Update the profile linking file. Use notes in the notes dictionary
                        note_list = update_profile_link_file(
                            nt_seq_type=nt_profile,
                            aa_seq_type=aa_profile,
                            aa_profile_path=aa_profile_path,
                            note=notes[contig][query_range]
                        )
                    # If there are no notes for the contig:range, create notes from scratch
                    else:
                        # Update the profile linking file. Use notes in the notes_list list
                        note_list = update_profile_link_file(
                            nt_seq_type=nt_profile,
                            aa_seq_type=aa_profile,
                            aa_profile_path=aa_profile_path,
                            note=note_list
                        )
                # If there are no notes for the contig, create notes from scratch
                else:
                    # Update the profile linking file. Use notes in the notes_list list
                    note_list = update_profile_link_file(
                        nt_seq_type=nt_profile,
                        aa_seq_type=aa_profile,
                        aa_profile_path=aa_profile_path,
                        note=note_list
                    )
                # Join all the notes from the list with semicolons
                note_str = '; '.join(note_list)
                # Update the data string with the notes
                data += f'{note_str}'
                # Add a newline to the data string
                data += '\n'
        # If there were no hits for the sample, add negative values to the data string
        if not data:
            data = f'{sample.name}\t0\t0\t1\t0\t0\t1\n'
    # If the report file does not already exist, write the header and data strings
    if not os.path.isfile(report_file):
        with open(report_file, 'w', encoding='utf-8') as report:
            report.write(header)
            report.write(data)
    # If the report already exists, write only the data string
    else:
        with open(report_file, 'a+', encoding='utf-8') as report:
            report.write(data)


def update_profile_link_file(nt_seq_type, aa_seq_type, note, aa_profile_path):
    """
    Update the file linking amino acid sequence type to the (multiple) corresponding nucleotide sequence type(s)
    :param nt_seq_type: String of the nucleotide sequence type
    :param aa_seq_type: String of the amino acid sequence type
    :param note: List of notes on the alleles
    :param aa_profile_path: String of the absolute path of the folder in which the amino acid profile file is located
    :return: note: Update list of notes
    """
    # Set the name of the link file
    link_file = os.path.join(aa_profile_path, 'aa_nt_profile_links.tsv')
    # Create a dictionary of nucleotide sequence type matches
    nt_match = {}
    # Initialise a boolean to track whether the amino acid sequence type is already present in the profile link file
    aa_match = False
    # Initialise a dictionary to store the aa_seq_type: nt_seq_type(s)
    links = {}
    # Initialise a list of all amino acid sequence types present the file
    records = []
    # Open the profile link file to read in the contents
    with open(link_file, 'r', encoding='utf-8') as profile_link:
        for line in profile_link:
            # Split the amino acid sequence type from the nucleotide sequence type(s)
            # e.g 1 1 or 116 116;125;187;39973;92286;1000005
            aa_seq, nt_seq = line.split('\t')
            # Check to see if the extracted amino acid sequence type matches the sequence type of the sample
            if aa_seq_type == aa_seq:
                # Set the match boolean to True (there is a match)
                aa_match = True
                # Check to see if the nucleotide sequence type of sample is in the semicolon-separated list of
                # nucleotide sequence types corresponding to the amino acid sequence type
                if nt_seq_type in nt_seq.rstrip().split(';'):
                    # Update the nucleotide sequence type match dictionary
                    nt_match[aa_seq] = nt_seq.rstrip()
            # Update the link dictionary
            links[aa_seq] = nt_seq.rstrip()
            # Add the amino acid sequence type to the list
            records.append(aa_seq)
    # Check if the amino acid of the sample matched a previous sequence type
    if aa_match:
        # Check if there was a match of the sample's nucleotide sequence type
        if not nt_match:
            # Append the nucleotide sequence type to the string of nucleotide sequence type matches
            links[aa_seq_type] += f';{nt_seq_type}'
            # Update the note
            note.append(f'Novel nt_seq_type {nt_seq_type} links to aa_seq type {aa_seq_type}')
    # If no match, this is a novel amino acid sequence type
    else:
        # Update the link dictionary novel amino acid sequence type: novel nucleotide sequence type
        links[aa_seq_type] = nt_seq_type
        # Add the novel sequence type to the list
        records.append(aa_seq_type)
        # Update the notes
        note.append(f'Novel nt_seq_type {nt_seq_type}, and aa_seq_type {aa_seq_type}')
    # Overwrite the profile link file with the updated links
    with open(link_file, 'w', encoding='utf-8') as profile_link:
        for record in records:
            profile_link.write(f'{record}\t{links[record]}\n')
    return note


def split_alleles(allele_files, output_path):
    """
    Split FASTA files into individual sequences
    :param allele_files: List of absolute path to FASTA-formatted allele sequence files
    :param output_path: String of the absolute path into which the individual sequence files are to be written
    """
    # Create the output path if it doesn't already exist
    make_path(inpath=output_path)
    # Ensure that the path could be created
    if not os.path.isdir(output_path):
        logging.error('Could not create desired output folder: %s', output_path)
        raise SystemExit
    # Iterate over all the allele files
    for allele_file in allele_files:
        # Use SeqIO to load all the records in the file
        for record in SeqIO.parse(allele_file, 'fasta'):
            # Set the name of the file to be the FASTA header
            output_file = os.path.join(output_path, f'{record.id}.fasta')
            # Write the record to the new file
            with open(output_file, 'w', encoding='utf-8') as output:
                SeqIO.write(record, output, 'fasta')


def parse_aa_blast(runmetadata, extended_fieldnames, fieldnames, gene_names, notes, aa_allele_path,
                   report_path, cutoff):
    """
    Parse amino acid BLAST results
    :param runmetadata: List of metadata objects
    :param extended_fieldnames: List of the BLAST fields used, as well as the additional percent
    match in index 14
    :param fieldnames: List of the BLAST fields used
    :param gene_names: List of all gene names in the analysis
    :param notes: List of sample-specific notes
    :param aa_allele_path: String of the absolute path to the folder containing amino acid allele files
    :param report_path: String of the absolute path to the folder into which reports are to be written
    them to be considered co-located. Default is 50 bp
    :param cutoff: Integer of the minimum percent identity between query and subject sequence.
    Default is 90
    :return: runmetadata: Updated list of metadata objects
    :return: filtered: Boolean of whether the sample fails quality/length checks
    :return: notes: Updated list of sample-specific notes
    """
    logging.info('Parsing BLAST outputs')
    # Initialise a boolean to track if the sequence fails checks
    filtered = False
    for sample in runmetadata.samples:
        # Initialise GenObjects as required
        sample.alleles.blastlist = []
        sample.alleles.targetsequence = {}
        # Read the BLAST outputs into a dictionary
        blastdict = create_blast_dict(
            sample=sample,
            extended_fieldnames=extended_fieldnames
        )
        # Initialise a boolean to track whether this contig:query_range has already been processed
        processed = False
        # Go through each BLAST result
        for row in blastdict:
            # Ignore the headers
            if row['query_id'].startswith(fieldnames[0]):
                continue
            # Create variables to reduce extra typing and for extra clarity
            target_id = row['subject_id']
            percent_id = float(row['percent_match'])
            # If the match is perfect
            if percent_id == 100:
                # Add the name of the matching allele to the list
                sample.alleles.blastlist.append(target_id)
                # Update the processed boolean to indicate that this region has been processed
                processed = True
            # If the match is imperfect, but greater than the cutoff
            elif cutoff < percent_id < 100 and not processed:
                # Determine which gene is being processed by finding the match of the genes against the allele
                gene = [gene for gene in gene_names if gene in target_id][0]
                query_seq = row['query_sequence']
                # Evaluate the sequence for length, as well as required start/stop codons
                filtered, notes = evaluate_translated_allele(
                    aa_seq=query_seq,
                    gene=gene,
                    notes=notes,
                    aa=True
                )
                # Find the next available allele identifier in the database
                aa_allele_id = find_next_allele(
                    gene=gene,
                    allele_path=aa_allele_path
                )
                # Set the name of the allele as gene_alleleID
                aa_allele = f'{gene}_{aa_allele_id}'
                # Update the allele database with the novel allele
                notes = update_allele_databases(
                    query_sequence=query_seq,
                    header=aa_allele,
                    filtered=filtered,
                    gene=gene,
                    report_path=report_path,
                    allele_path=aa_allele_path,
                    notes=notes,
                    molecule='Amino acid'
                )
                # If the allele passes the necessary checks, update the list of results
                if not filtered:
                    sample.alleles.blastlist.append(aa_allele)
                # Update the processed boolean
                processed = True
    return runmetadata, filtered, notes


def analyse_aa_alleles(runmetadata, gene_names, notes):
    """
    Analyse supplied amino acid alleles to ensure that they pass quality checks
    :param runmetadata: List of metadata objects
    :param gene_names: List of all gene names in the analysis
    :param notes: List of sample-specific notes
    :return: runmetadata: Updated list of metadata objects
    :return: notes: Updated list of sample-specific notes
    """
    # Iterate through all the samples
    for sample in runmetadata.samples:
        # Iterate over all the records in the file (should only be one, as these files must be split with the
        # split_alleles function)
        for record in SeqIO.parse(sample.general.bestassemblyfile, 'fasta'):
            # Determine to which gene the allele corresponds
            gene = [gene for gene in gene_names if gene in record.id][0]
            # Perform content/length checks
            filtered, notes = evaluate_translated_allele(
                aa_seq=record.seq,
                gene=gene,
                notes=notes,
                aa=True
            )
    return runmetadata, notes


def report_aa_alleles(runmetadata, report_file, notes):
    """
    Create an amino acid query-specific report with sample name, allele match, and notes
    :param runmetadata: List of metadata objects
    :param report_file: String of absolute path to the report file
    :param notes: List of sample-specific notes
    """
    # Initialise the header string
    header = 'Sample\tMatch\tNotes\n'
    # Create an empty string to store the sample-specific results
    data = str()
    # Iterate over all the samples
    for sample in runmetadata.samples:
        # Extract the list of hits from the metadata object, and join with semicolons
        matches = ';'.join(sample.alleles.blastlist)
        # Join the list of notes with
        note = ';'.join(notes)
        # Populate the data string with the matches and notes
        data += f'{sample.name}\t{matches}\t{note}\n'
    # If the report doesn't already exist write the header and data string
    if not os.path.isfile(report_file):
        with open(report_file, 'w', encoding='utf-8') as report:
            report.write(header)
            report.write(data)
    # Otherwise write only the data string
    else:
        with open(report_file, 'a+', encoding='utf-8') as report:
            report.write(data)
