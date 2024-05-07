#!/usr/bin/env python3

"""
Perform allele profiling
"""

# Standard imports
from csv import DictReader
from glob import glob
import logging
import json
import os

# Third party imports
from Bio.Seq import Seq
from Bio import SeqIO

from olctools.accessoryFunctions.accessoryFunctions import MetadataObject


# Local imports
from allele_tools.methods import \
     extract_novel_alleles
from olctools.accessoryFunctions.accessoryFunctions import combinetargets
from genemethods.geneseekr.geneseekr import GeneSeekr


__author__ = 'adamkoziol'


def allele_prep(allele_path, gene_names, combined_targets, amino_acid):
    """
    Create a 'combinedtargets.fasta' file
    :param allele_path: Name and absolute path to folder containing allele
    database files
    :param gene_names: List of all gene names in the analysis
    :param combined_targets: String of absolute path to the
    combinedtargets.fasta file
    :param amino_acid: Boolean of whether the query sequences are amino acid
    or nucleotide
    :return: records: Dictionary of allele_header:allele_sequence
    :return: gene_names: List of all gene names extracted from allele headers
    :return: data: String of header including all gene names. To be used in
    creating final report
    """
    logging.info('Reading allele files')
    custom_targets = os.path.join(allele_path, 'custom.tfa')
    combined_targets_db = glob(os.path.join(allele_path, 'combinedtargets*'))
    records = {}
    # Clear out any previously created combinedtargets files (and associated
    # BLAST databases)
    clear_alleles(combined_targets_db=combined_targets_db,
                  custom_targets=custom_targets)
    alleles = glob(os.path.join(allele_path, '*.*fa*'))
    # If the dictionary hasn't already been populated by a previous iteration,
    # remove the path and
    # the extension from each of the allele files to determine gene names
    if not gene_names:
        for allele in alleles:
            gene_names.append(
                os.path.splitext(
                    os.path.basename(
                        allele
                    )
                )[0].replace('_alleles', '')
            )
    # Populate the header string for the final report
    genes = '\t'.join(sorted(gene_names))
    data = f'ST\t{genes}\n'
    # Create the combinedtargets file
    if not os.path.isfile(combined_targets):
        if amino_acid:
            combinetargets(
                targets=alleles,
                targetpath=allele_path,
                mol_type='prot'
            )
        else:
            combinetargets(
                targets=alleles,
                targetpath=allele_path
            )
    # Create BLAST databases
    GeneSeekr.makeblastdb(
        fasta=combined_targets,
        program='blastn' if not amino_acid else 'blastp'
    )
    # Load all the FASTA records from the combinedtargets file
    for allele_file in alleles:
        for record in SeqIO.parse(allele_file, 'fasta'):
            records[record.id] = str(record.seq)
    return records, gene_names, data


def clear_alleles(combined_targets_db, custom_targets):
    """
    Remove any combinedtargets.fasta or custom.tfa files present in the
    allele path
    :param combined_targets_db: List of combinedtargets files, including BLAST
    database files
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
    profile_data = {}
    # Only load the profile file if it exists
    if os.path.isfile(profile_file):
        logging.info('Extracting profiles from profile file')
        # Open an Excel-formatted sequence profile file as a dictionary
        profile = DictReader(
            open(
                profile_file,
                encoding='utf-8'
            ), dialect='excel-tab'
        )
        # Iterate through the rows
        for row in profile:
            # Populate the profile dictionary with profile number:
            # {gene: allele}.
            allele_comprehension = {
                gene: allele for gene, allele in row.items() if gene != 'ST'
            }
            # Extract the sequence type number from the first field name
            seq_type = row[profile.fieldnames[0]]
            # Update the profile data dictionary
            profile_data[seq_type] = allele_comprehension
    return profile_data


def parseable_blast_outputs(
        runmetadata: MetadataObject,
        fieldnames: str,
        extended_fieldnames: str,
        records: dict,
        cutoff: int = 90):
    """
    Add a header to the BLAST report, so that it is easier to figure out what
    is in each column
    :param runmetadata: Metadata object containing a list of all
    metadata objects
    :param fieldnames: String of all the field names in the BLAST report
    :param extended_fieldnames: String of the BLAST field names plus the
    calculated percent identity
    :param records: Dictionary of allele name: allele sequence
    :param cutoff: Integer of the minimum percent identity between query and
    subject sequence. Default is 90
    """
    logging.info('Adding headers to BLAST outputs')
    for sample in runmetadata.samples:
        data = []
        # Load the first line of the report
        with open(
                sample.alleles.blast_report,
                'r',
                encoding='utf-8') as report:
            header_line = report.readline().strip()
        # Split the header on tabs
        header_list = header_line.split('\t')
        # Check to see if the header has already been added. Skip this step
        # if it has been added.
        if header_list[0] != fieldnames[0]:
            with open(
                sample.alleles.blast_report,
                'r',
                    encoding='utf-8') as report:
                header = list(report.readline().split('\t'))
            if len(header) == len(fieldnames):
                current_fieldnames = fieldnames
            else:
                current_fieldnames = extended_fieldnames
            blastdict = DictReader(
                open(
                    sample.alleles.blast_report, encoding='utf-8'
                ),
                fieldnames=current_fieldnames,
                dialect='excel-tab'
            )
            # Go through each BLAST result
            for row in blastdict:
                # Calculate the percent identity and extract the bit score
                # from the row Percent identity is the (length of the
                # alignment - num mismatches) / total query length
                percent_identity = float(
                    '{:0.2f}'.format(
                        (float(
                            row['identical']) - float(row['gaps'])
                         ) / len(records[row['subject_id']]) * 100
                    )
                )
                # Filter the results based on the cutoff value
                if percent_identity < cutoff:
                    continue
                # Create a percent match entry based on the calculated
                # percent identity match
                row['percent_match'] = str(percent_identity)
                # Add the updated row to the list
                data.append(row)

            # Overwrite the original BLAST outputs to include headers,
            # and the percent match
            with open(
                sample.alleles.blast_report,
                'w',
                    encoding='utf-8') as updated_report:
                # Add the header
                headers = '\t'.join(extended_fieldnames)
                updated_report.write(f'{headers}\n')
                # Add the results
                for row in data:
                    for header in extended_fieldnames:
                        # Write the value from the row with the header as
                        # the key
                        updated_report.write(
                            '{value}\t'.format(value=row[header])
                        )
                    # Add a newline for each result
                    updated_report.write('\n')


def parse_results(
        runmetadata: MetadataObject,
        fieldnames: str,
        extended_fieldnames: str,
        amino_acid: bool,
        genome_query: bool = False):
    """
    Parse the BLAST results, and populate GenObjects
    :param runmetadata: Metadata object containing list of all metadata objects
    :param fieldnames: String of all the field names in the BLAST report
    :param extended_fieldnames: String of the BLAST field names plus the
    calculated percent identity
    :param amino_acid: Variable on whether targets are protein
    :param genome_query: Boolean of whether the query is a genome, and the
    subject is the allele
    :return: Updated runmetadata
    """
    logging.info('Parsing BLAST outputs')
    for sample in runmetadata.samples:
        # Initialise GenObjects as required
        sample.alleles.blastlist = []
        sample.alleles.targetsequence = {}
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(
            open(
                sample.alleles.blast_report, encoding='utf-8'
            ),
            fieldnames=extended_fieldnames,
            dialect='excel-tab'
        )
        resultdict = {}
        # Go through each BLAST result
        for row in blastdict:
            # Ignore the headers
            if row['query_id'].startswith(fieldnames[0]):
                pass
            else:
                target_id = row['query_id'] if not genome_query else \
                    row['subject_id']
                target_start = row['subject_start']
                target_end = row['subject_end']
                target_seq = row['query_sequence']
                # Remove unwanted pipes added to the name
                target = target_id.lstrip('gb|').rstrip('|') if '|' \
                    in target_id else target_id
                # If the percent identity is equal to the cutoff
                if float(row['percent_match']) == 100:
                    # Append the hit dictionary to the list
                    sample.alleles.blastlist.append(row)
                    # Update the dictionary with the target and
                    # percent identity
                    resultdict.update({target: row['percent_match']})
                    # Determine if the orientation of the sequence is
                    # reversed compared to
                    # the reference sequence
                    if int(target_end) < int(target_start) and not amino_acid:
                        # Create a sequence object using BioPython
                        seq = Seq(target_seq)
                        # Calculate the reverse complement of the sequence
                        querysequence = str(seq.reverse_complement())
                    # If the sequence is not reversed, use the sequence as it
                    # is in the output
                    else:
                        querysequence = target_seq
                    # Add the sequence in the correct orientation to the sample
                    try:
                        sample.alleles.targetsequence[target].append(
                            querysequence
                        )
                    except (AttributeError, KeyError):
                        sample.alleles.targetsequence[target] = []
                        sample.alleles.targetsequence[target].append(
                            querysequence
                        )
            # Add the percent identity to the object
            sample.alleles.blastresults = resultdict
        # Populate missing results with 'NA' values
        if len(resultdict) == 0:
            sample.alleles.blastresults = 'NA'
    return runmetadata


def profile_alleles(
        runmetadata: MetadataObject,
        profile_dict: dict,
        profile_set: list,
        records: list,
        amino_acid: bool = False,
        novel_alleles: bool = False,
        genome_query: bool = False,
        allele_path: str = None,
        report_path: str = None,
        cutoff: int = 75):
    """
    Create the gene:allele profile from the BLAST outputs from each sample
    :param runmetadata: Metadata object containing a list of all
    metadata objects
    :param profile_dict: Dictionary to store gene:allele profile for
    each sample
    :param profile_set: List of all unique profiles
    :param records: List of all gene names
    :param novel_alleles: Boolean of whether novel alleles should be extracted
    from BLAST hit if there is no 100% match
    :param genome_query: Boolean of whether the query is a genome, and the
    subject is the allele
    :param amino_acid: Variable on whether targets are protein
    :param allele_path: Name and absolute path to folder containing allele
    database files
    :param report_path: Name and absolute path to folder in which reports are
    to be created
    :param cutoff: Integer of the minimum percent identity between query and
    subject sequence. Default is 75
    :return: Updated profile_dict and profile_set
    """
    logging.info('Determining allele profiles')
    # Iterate through all the samples
    for sample in runmetadata.samples:
        # Initialise a dictionary to store the profile information for
        # each sample
        profile_dict[sample.name] = {}
        # Initialise a dictionary to store the gene:allele combinations for
        # each sample
        allele_comprehension = {}
        # Each gene in the analysis is stored in the list of genes created
        # in allele_finder
        for gene in records:
            gene = gene[0].lower() + gene[1:-1] + gene[-1].upper()
            # Create a variable to track whether the current gene is present
            # in the current sample
            present = False
            # Iterate through all the BLAST outputs for the sample
            for allele in sample.alleles.blastresults:
                # If the gene name is present as a substring of the allele
                # e.g. adk in adk_1, then the gene is present in the
                # BLAST outputs
                if gene.lower() in allele.lower():
                    # Strip off the allele number from the allele e.g.
                    # adk_1 yields 1
                    allele_id = allele.split('_')[-1]
                    # Update the dictionary with the new gene: allele number
                    # for the sample
                    allele_comprehension.update({gene: allele_id})
                    # Update the gene presence variable
                    present = True
            # If, after iterating through all the BLAST outputs, the gene is
            # not present in the sample, update the gene: allele to reflect
            # this absence
            if not present:
                if novel_alleles:
                    sample, novel_allele, query_sequence = \
                        extract_novel_alleles(
                            sample=sample,
                            gene=gene,
                            genome_query=genome_query,
                            amino_acid=amino_acid,
                            allele_path=allele_path,
                            report_path=report_path,
                            cutoff=cutoff,
                        )
                    novel_allele = novel_allele.split('|')[0]
                    try:
                        sample.alleles.targetsequence[novel_allele].append(
                            query_sequence
                        )
                    except KeyError:
                        sample.alleles.targetsequence[novel_allele] = \
                            [query_sequence]
                    if novel_allele:
                        allele_comprehension.update(
                            {
                                gene: novel_allele.split('_')[-1]
                                }
                        )
                    else:
                        allele_comprehension.update({gene: '0'})
                else:
                    # Set missing alleles to 'ND'
                    allele_comprehension.update({gene: '0'})
        # In order to hash the dictionary, use JSON, with sorted keys
        # to freeze it
        frozen_allele_comprehension = json.dumps(
            allele_comprehension,
            sort_keys=True
        )
        # Update the dictionary of profiles with the hash of the frozen
        # dictionary: list of samples with that hash
        if hash(frozen_allele_comprehension) not in profile_dict:
            profile_dict[hash(frozen_allele_comprehension)] = [sample.name]
        else:
            profile_dict[hash(frozen_allele_comprehension)].append(sample.name)
        # Add the 'regular' dictionary to the list of all profiles as required
        if allele_comprehension not in profile_set:
            profile_set.append(allele_comprehension)
    return profile_dict, profile_set


def match_profile(
        profile_data: dict,
        profile_dict: dict,
        profile_matches: dict):
    """
    Match current profiles to any previously created profiles
    :param profile_data: Dictionary of seq_type: {gene name:allele ID}
    :param profile_dict: Dictionary of gene:allele profile for each sample
    :param profile_matches: Dictionary of seq_type: matching profile
    :return: profile_matches: Updated dictionary of seq_type: matching profiles
    """
    # If the profile_data dictionary was not populated in the read_profiles
    # methods, there is nothing to match
    if profile_data:
        logging.info('Matching new profiles against profile file')
        # Extract the sequence type and allele dictionary from the profile file
        for seq_type, allele_comprehension in profile_data.items():
            # Freeze the allele comprehension as above
            frozen_allele_comprehension = json.dumps(
                allele_comprehension,
                sort_keys=True
            )
            try:
                # Extract the samples that match this profile
                matches = profile_dict[hash(frozen_allele_comprehension)]
                # Update the dictionary with the matching samples
                profile_matches[seq_type] = matches
            # The profile will not necessarily match any of the profiles
            # found in the analysis
            except KeyError:
                pass
    return profile_matches


def create_profile(
        profile_data: dict,
        profile_set: list,
        new_profiles: list,
        profile_dict: dict,
        profile_matches: dict):
    """
    Create new profiles for novel profiles as required
    :param profile_data: Dictionary of seq_type: {gene name:allele ID}
    :param profile_set: List of all unique profiles
    :param new_profiles: List of novel profiles
    :param profile_dict: Dictionary of gene:allele profile for each sample
    :param profile_matches: Dictionary of seq_type: matching profile
    :return: profile_matches: Updated dictionary of seq_type: matching profiles
    :return: profile_data: Updated dictionary of seq_type: {gene name:
    allele ID}
    :return: new_profiles: Updated list of novel profiles
    """
    # Initialise the sequence type to be 1
    seq_type = 1
    # If the profile_data dictionary exists, set the sequence type to be the
    # last of the entries in the dictionary plus one, as that corresponds to
    # the next sequence type
    if profile_data:
        # seq_type = len(profile_data) + 1
        seq_type = sorted(int(st) for st in profile_data.keys())[-1] + 1
    # Initialise a list to store the matched samples
    matched = []
    # Iterate through all the profiles in the analysis
    for allele_comprehension in profile_set:
        # Ensure that the allele comprehension (profile) is not already in
        # the profile file
        if allele_comprehension not in [
                profiled_alleles for _,
                profiled_alleles in profile_data.items()]:
            # Add the new profile to the list of new profiles
            alleles = '\t'.join(
                allele_num.split('_')[-1] for gene, allele_num in sorted(
                    allele_comprehension.items()
                )
            )
            new_profiles.append(
                f'{seq_type}\t{alleles.rstrip()}'
            )
            # Freeze the comprehension in order to be used as the key in the
            # profile dictionary
            frozen_allele_comprehension = json.dumps(
                allele_comprehension,
                sort_keys=True
            )
            matches = profile_dict[hash(frozen_allele_comprehension)]
            # Check to see if this sequence type hasn't already been found in
            # the current analysis
            if matches not in matched:
                # Update the dictionary with the new sequence type: list
                # of samples
                profile_matches[seq_type] = matches
                profile_data[seq_type] = allele_comprehension
                # Add the matches to the list of matches
                matched.append(matches)
            # Increment the sequence type number of the next entry
            seq_type += 1
    return profile_matches, profile_data, new_profiles


def sequence_typer(
            profile_report: str,
            data: str,
            runmetadata: MetadataObject,
            profile_matches: dict,
            profile_data: dict,
            update: bool = False,
            amino_acid: bool = False):
    """
    Perform the final sequence typing, and create the report
    :param profile_report: String of absolute path to report file
    :param data: String of header including all gene names. To be used in
    creating final report
    :param runmetadata: Metadata object containing a list of all
    metadata objects
    :param profile_matches: Dictionary of seq_type: matching profile
    :param profile_data: Dictionary of seq_type: {gene name:allele ID}
    :param update: Boolean of whether the report is to be created or updated.
    Default is False (created)
    :param amino_acid: Boolean of whether the query sequences are amino acid.
    Default is False
    """
    # Open the report
    mode = 'w' if not update else 'a+'
    if not update:
        # Initialise the header with an extra 'Sample' column plus the
        # comma-separated list
        # of gene names
        data = 'Sample\t' + data
    else:
        if not os.path.isfile(profile_report):
            # Initialise the header with an extra 'Sample' column plus the
            # comma-separated list of gene names
            data = 'Sample\t' + data
        else:
            data = str()
    with open(profile_report, mode, encoding='utf-8') as report:
        for sample in runmetadata.samples:
            # Iterate through all the matches to sequence profiles in
            # the analysis
            for seq_type, sample_names in profile_matches.items():
                # Check if the sample name is in the list of samples names
                # with the current sequence type
                if sample.name in sample_names:
                    # Add the sample name, sequence type, and all the allele
                    # numbers to the
                    # report string
                    complement = '\t'.join(
                        allele_num.split('_')[-1] for _, allele_num in sorted(
                            profile_data[seq_type].items())
                    )
                    data += \
                        f'{sample.name}\t{seq_type}\t{complement.rstrip()}\n'
                    # Update the appropriate GenObject based on the current
                    # molecule (DNA or amino acid)
                    if not amino_acid:
                        sample.alleles.nt_st = seq_type
                        sample.alleles.nt_profile = profile_data[seq_type]
                    else:
                        sample.alleles.aa_st = seq_type
                        sample.alleles.aa_profile = profile_data[seq_type]
        # Write the report
        report.write(data)
    return runmetadata


def append_profiles(
        new_profiles: list,
        profile_file: str,
        data: str,
        novel_profiles: bool = False,
        profile_path: str = None,
        gene_names: str = None):
    """
    Add new profiles to the profile file
    :param new_profiles: List of all novel profiles in this analysis
    :param profile_file: Name and absolute path to a profile file
    :param data: String of header including all gene names. To be used in
    creating final report
    :param novel_profiles: Boolean of whether the novel_profiles.txt file is
    to be populated. Default is False
    :param profile_path: String of absolute path of folder in which profiles
    are located
    :param gene_names: List of all genes in the analysis
    """
    # Only try to add new profiles if there are new profiles in the analysis
    if new_profiles:
        # Initialise the string to store the new profile
        new_data = str()
        # If the profile file does not exist, add the string of 'ST',
        # comma-separated gene names to the headers
        if not os.path.isfile(profile_file):
            new_data = data
        # Iterate through all the new profiles, and add them to the new
        # profile string
        for profile in new_profiles:
            new_data += f'{profile}\n'
        # Open the report with a+ to either create, or append the profile
        # string to it
        with open(profile_file, 'a+', encoding='utf-8') as profile:
            profile.write(new_data)
        if novel_profiles:
            novel_profile_file = os.path.join(
                profile_path,
                'novel_profiles.txt'
            )
            if not os.path.isfile(novel_profile_file):
                with open(novel_profile_file, 'w', encoding='utf-8') as novel:
                    novel.write(
                        'ST\t{names}\n'.format(names='\t'.join(gene_names))
                    )
            with open(novel_profile_file, 'a+', encoding='utf-8') as novel:
                novel.write(new_data)
