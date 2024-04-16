#!/usr/bin/env python

"""
Translate nucleotide alleles to amino acid, and remove duplicates
"""

# Standard imports
from argparse import ArgumentParser
from glob import glob
import logging
import os
import shutil
import sys
from typing import TextIO

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    make_path,
    SetupLogging
)
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Local imports
from allele_tools.allele_profiler import read_profile
from allele_tools.profile_reduce import ProfileReduce
from allele_tools.methods import (
    evaluate_translated_length,
    generic_evaluate_translated_length,
    remove_combined_db_files,
    pathfinder
)


class Translate:

    """
    Translate and reduce alleles
    """

    def main(self):
        """
        Run the appropriate methods in the correct order
        """
        # Read in all the allele files
        self.load_alleles()
        # Parse the alleles
        self.parse_alleles(length_dict=self.length_dict)
        # If a profile file has been provided, run the profiling methods
        if self.profile_file:
            # Extract seq_type:allele comprehensions from the nucleotide
            # profile file
            self.profile_data = read_profile(profile_file=self.profile_file)
            # Create the amino acid profiles
            self.aa_profile()
            # Remove duplicate profiles
            reduce = ProfileReduce(
                profile=self.aa_profile_file,
                names=self.gene_name_file
            )
            reduce.main()
            # Read in the profile data again now that the profile file has
            # been updated
            self.aa_profile_data = read_profile(
                profile_file=self.aa_profile_file
            )
            # Find linkages between nucleotide and amino acid profiles
            self.profile_link()
            # Write profile links to file
            self.link_file()
            # Copy and rename the reduced profile file
            self.copy_profile()

    def load_alleles(self):
        """
        Use SeqIO to read in all the gene sequences
        """
        for allele_file in self.sequence_files:
            gene_name = os.path.splitext(os.path.basename(allele_file))[0]
            # Add the gene name to the set of names
            self.gene_names.add(gene_name)
            self.allele_dict[gene_name] = SeqIO.to_dict(
                SeqIO.parse(allele_file, 'fasta')
            )

    def parse_alleles(
            self,
            length_dict: dict = None):
        """
        Parse the allele files to translate the amino acid sequence using
        BioPython. Write the amino acid sequence to file. Store the allele
        name in the notes. Find duplicates, and link the nucleotide allele
        name to the amino acid allele name
        :param length_dict: Dictionary of gene name: minimum acceptable
        length of translated sequence
        """
        logging.info('Translating and parsing alleles')
        for gene_name, allele_dict in self.allele_dict.items():
            # Initialise the dictionary to store the links between the nt
            # alleles and the aa alleles with the name of the gene
            self.allele_links[gene_name] = {}
            #  Initialise a set of allele that do not pass the required
            # filters, and must be removed from the database
            remove_list = set()
            logging.info('Processing %s', gene_name)
            # Initialise a dictionary to store the translated allele
            # sequence: allele name
            seq_allele = {}
            # Set the name and path of the allele and notes files
            allele_file = os.path.join(
                self.translated_path,
                f'{gene_name}.fasta'
            )
            notes_file = os.path.join(
                self.notes_path,
                f'{gene_name}_notes.txt'
            )
            # Open the file to store the translated alleles
            with open(allele_file, 'w', encoding='utf-8') as aa_alleles:
                # Open the notes file
                with open(notes_file, 'w', encoding='utf-8') as notes:
                    # Create the header for the notes file
                    notes.write('nt_allele\taa_allele\tnote\n')
                    # Iterate through all the alleles in the dictionary
                    for allele, details in allele_dict.items():
                        # Create a list to store notes
                        note = []
                        # Create a boolean to track whether the sequence
                        # is filtered
                        filtered = False
                        # Create a string to store the untrimmed nucleotide
                        # sequence
                        original_nt_sequence = str(details.seq)
                        # Calculate the translated sequence
                        translated_allele = details.seq.translate()
                        # Remove all sequence after a stop codon (*)
                        split_allele = translated_allele.split('*')
                        # Create a string to hold the trimmed sequence
                        trimmed_seq = str()
                        # If there are multiple stop codons in the sequence,
                        # trim to the first one
                        if len(split_allele) > 2:
                            if split_allele[1]:
                                for trimmed in split_allele[1:]:
                                    trimmed_seq += trimmed
                            else:
                                for trimmed in split_allele[2:]:
                                    trimmed_seq += trimmed
                            note.append(f'Trimmed {trimmed_seq} from end')
                        elif len(split_allele) == 1:
                            pass
                        else:
                            if split_allele[-1] and not \
                                    str(translated_allele).endswith('*'):
                                trimmed_seq += split_allele[-1]
                                note.append(f'Trimmed {trimmed_seq} from end')
                        # Create a dictionary to store the allele nt sequences
                        # to check to see if any are duplicates following
                        # trimming
                        nt_sequences = {}
                        filtered, note, nt_sequences, translated_allele = \
                            Translate.trim_alleles(
                                note=note,
                                allele=allele,
                                sequence=details,
                                gene_name=gene_name,
                                nt_allele_path=self.path,
                                trim_length=len(trimmed_seq),
                                length_dict=length_dict,
                                filtered=filtered,
                                nt_sequences=nt_sequences,
                                original_nt_sequence=original_nt_sequence
                            )
                        # If the allele has not been filtered through the
                        # trim_alleles
                        if not filtered:
                            # Determine if this amino acid allele is new
                            if str(translated_allele) not in seq_allele:
                                # Add the string of the amino acid sequence to
                                # the dictionary
                                seq_allele[str(translated_allele)] = allele
                                # Create a SeqRecord of the translated allele
                                seq_record = SeqRecord(
                                    seq=translated_allele,
                                    id=allele,
                                    name=str(),
                                    description=str()
                                )
                                # Write the SeqRecord to file
                                SeqIO.write(
                                    sequences=seq_record,
                                    handle=aa_alleles,
                                    format='fasta'
                                )
                                # Update the notes with the allele naming
                                # information
                                notes.write(
                                    f'{allele}\t{allele}\t{";".join(note)}\n'
                                )
                                # Populate the linking dictionary with the nt
                                # allele: aa allele
                                self.allele_links[
                                    gene_name][allele.split('_')[-1]] = \
                                    allele.split('_')[-1]
                            # Amino acid allele already exists
                            else:
                                # Extract the allele name corresponding to the
                                # translated sequence
                                aa_allele = seq_allele[str(translated_allele)]
                                # Update the notes, including that this allele
                                # is a duplicate, and a pointer to the original
                                notes.write(
                                    f'{allele}\t{aa_allele}\tDuplicate'
                                )
                                if not note:
                                    notes.write('\n')
                                else:
                                    notes.write(f'; {";".join(note)}\n')
                                # Populate the linking dictionary with the nt
                                # allele: aa allele
                                self.allele_links[
                                    gene_name][allele.split('_')[-1]] = \
                                    aa_allele.split('_')[-1]
                                self.allele_links[gene_name]['0'] = '0'
                        # Filtered alleles must be removed from the database
                        else:
                            Translate.write_filtered_allele_notes(
                                notes=notes,
                                allele=allele,
                                note=note,
                            )
                            remove_list.add(allele)
                            # Create a SeqRecord of the translated allele
                            seq_record = SeqRecord(
                                seq=translated_allele,
                                id=allele,
                                name=str(),
                                description=str()
                            )
                            # Write the filtered alleles to the filtered
                            # alleles file
                            Translate.create_or_update_filtered_files(
                                gene=gene_name,
                                allele_path=self.translated_path,
                                records_filter=[seq_record]
                            )
            # Remove the filtered sequences from the database
            Translate.remove_filtered_alleles_from_db(
                gene_name=gene_name,
                allele_list=remove_list,
                nt_allele_path=self.path,
            )
            self.load_alleles()

    @staticmethod
    def trim_alleles(
            note: list,
            allele: str,
            sequence: SeqIO,
            gene_name: str,
            nt_allele_path: str,
            trim_length: int,
            length_dict: dict,
            filtered: bool,
            nt_sequences: dict,
            original_nt_sequence: str):
        """
        Trim the alleles based on location of stop codons and whether the
        sequence is a multiple of three nucleotides. Evaluate the trimmed
        sequence based on length and contents.
        :param note: List of sequence-specific notes
        :param allele: String of the allele identifier
        :param sequence: SeqIO sequence object of nucleotide sequence
        :param gene_name: String of the gene name to which the
        allele corresponds
        :param nt_allele_path: String of the absolute path to the folder in
        which the nucleotide alleles are located
        :param trim_length: Integer of the number of nucleotides to be trimmed
        (due to internal stop codons)
        :param length_dict: Dictionary of minimum acceptable length for each
        gene in the analysis
        :param filtered: Boolean to track whether the sequence fails the
        quality/length checks
        :param nt_sequences: Dictionary of allele: sequence
        :param original_nt_sequence: String of the untrimmed nucleotide
        sequence of the allele
        :return: filtered: Updated boolean of whether the sequence fails
        quality/length checks
        :return: note: Update list of sequence-specific notes
        :return: nt_sequences: Updated dictionary of allele: sequence
        :return: translated_allele: SeqIO sequence object of trimmed,
        translated allele
        """
        # Determine the length of sequence to trim from the end of the sequence
        # Multiply the number of amino acid residues to trim by three to get
        # the number of nucleotides. Add the modulo three of the sequence to
        # yield a final length divisible by three
        nt_to_trim = 3 * trim_length + len(sequence) % 3
        # Check if trimming is required
        if nt_to_trim:
            # Slice the string by the number of required bases at the end
            sequence = sequence[:-nt_to_trim]
            # Update the nucleotide sequence in the database with the trimmed
            # version
            Translate.update_allele_db(
                nt_allele_path=nt_allele_path,
                gene=gene_name,
                allele=allele,
                nt_sequence=sequence
            )
        # Translate the sequence to amino acid
        translated_allele = sequence.seq.translate()
        # Perform content and length checks of the protein sequence
        if length_dict:
            filtered, note = evaluate_translated_length(
                aa_seq=str(translated_allele),
                length_dict=length_dict,
                gene=gene_name,
                notes=note,
                filtered=filtered
            )
        else:
            filtered, note = generic_evaluate_translated_length(
                aa_seq=str(translated_allele),
                sequence=original_nt_sequence,
                gene=gene_name,
                notes=note,
                filtered=filtered,
                cutoff=0.95
            )
        # Search the amino acid database for matches
        filtered, note, nt_sequences = Translate.find_duplicates(
            nt_sequences=nt_sequences,
            nt_sequence=sequence,
            allele=allele,
            filtered=filtered,
            note=note
        )
        return filtered, note, nt_sequences, translated_allele

    @staticmethod
    def update_allele_db(
            nt_allele_path: str,
            gene: str,
            allele: str,
            nt_sequence: SeqIO):
        """
        Update nucleotide allele files with newly trimmed sequences
        :param nt_allele_path: String of the absolute path to the folder
        containing the nucleotide allele database
        :param gene: String of the name of the current gene
        :param allele: String of the allele header (geneName_alleleID)
        :param nt_sequence: SeqIO object of the nucleotide sequence
        """
        # Set the name of the allele database file to update by joining the
        # allele path and the gene name
        nt_allele_file = os.path.join(nt_allele_path, f'{gene}.fasta')
        # Create a list to store all database records (including modified ones)
        updated_records = []
        # Iterate over all the sequences in the nucleotide allele file
        for record in SeqIO.parse(handle=nt_allele_file, format='fasta'):
            # Check if the header of the database sequence matches the header
            # of the query sequence
            if record.id == allele:
                # Update the record to be the query allele sequence object
                record = nt_sequence
            # Append the record to the list of records
            updated_records.append(record)
        # Overwrite the allele file with the records
        SeqIO.write(updated_records, handle=nt_allele_file, format='fasta')

    @staticmethod
    def find_duplicates(
            nt_sequences: dict,
            nt_sequence: SeqIO,
            allele: str,
            filtered: bool,
            note: list):
        """
        Match a query amino acid sequence against the protein allele
        database file
        :param nt_sequences: Dictionary of allele: sequence
        :param nt_sequence: SeqIO sequence object of the trimmed
        nucleotide sequence
        :param allele: String of the allele header (geneName_alleleID)
        :param filtered: Boolean of whether the sequence passes
        content/length filters
        :param note: List of sequence-specific notes
        :return: filtered: Updated boolean of sequence-specific
        content/length filtering results
        :return: note: Updated list of sequence-specific notes
        :return: nt_sequences: Updated list of allele: sequence
        """
        # Initialise a list to store any matches
        matches = []
        # Iterate over allele header, sequence in the dictionary
        for prev_allele, nt_seq in nt_sequences.items():
            # Check if the current nucleotide sequence matches a previous entry
            if nt_sequence.seq == nt_seq.seq:  # Compare sequence strings
                # Append matches to the list
                matches.append(prev_allele)
        # Check if the list has been populated with matches
        if not matches:
            # If no matches, update the dictionary with the novel sequence
            nt_sequences[allele] = nt_sequence
            return filtered, note, nt_sequences
        # If there are matches
        for match in matches:
            # Set the filtering boolean to True (this is not a novel sequence,
            # so do not update the database)
            filtered = True
            # Update the note
            note.append(
                f'Trimmed nt sequence matches previous '
                f'allele sequence: {match}'
            )
        return filtered, note, nt_sequences

    @staticmethod
    def write_filtered_allele_notes(
            notes: TextIO,
            allele: str,
            note: list):
        """
        Write notes for filtered sequences to file
        :param notes: File handle for notes file
        :param allele: String of the allele header (geneName_alleleID)
        :param note: List of sequence-specific notes
        """
        # Write the allele \t ND (because the sequence is filtered,there is
        # no amino acid allele) \t
        notes.write(f'{allele}\tND\tFiltered: {"; ".join(note)}\n')

    @staticmethod
    def remove_filtered_alleles_from_db(
            gene_name: str,
            allele_list: list,
            nt_allele_path: str):
        """
        Remove filtered sequences from the allele database file
        :param gene_name: String of the gene name currently being analysed
        :param allele_list: List of alleles to be removed from the database
        :param nt_allele_path: String of the absolute path to the folder
        containing the nucleotide allele database
        """
        # Remove the combinedtargets.fasta and BLAST database files from
        # the folder
        remove_combined_db_files(
            allele_path=nt_allele_path
        )
        # Set the name and path of the allele file
        nt_allele_file = os.path.join(nt_allele_path, f'{gene_name}.fasta')
        # Create a list of all the records in the database file using SeqIO
        records = SeqIO.parse(nt_allele_file, 'fasta')
        # Initialise a list to store unfiltered records
        records_keep = []
        # Initialise a list to store filtered records
        records_filter = []
        # Iterate over all the records in the database
        for record in records:
            # Check if the allele header matches any of the headers of alleles
            # to be filtered
            if record.id not in allele_list:
                # If the record is not to be filtered, add it to the keep list
                records_keep.append(record)
            # If the record header matches any of the alleles to be filtered,
            # add it to the filtered list
            else:
                records_filter.append(record)
        # Overwrite the nucleotide allele database file with all the
        # unfiltered records
        with open(nt_allele_file, 'w', encoding='utf-8') as allele_file:
            for record in records_keep:
                SeqIO.write(
                    sequences=record,
                    handle=allele_file,
                    format='fasta'
                )
        # Write the filtered records to the filtered records file
        if records_filter:
            Translate.create_or_update_filtered_files(
                gene=gene_name,
                allele_path=nt_allele_path,
                records_filter=records_filter
            )

    @staticmethod
    def create_or_update_filtered_files(
            gene: str,
            allele_path: str,
            records_filter: list):
        """
        Write the filtered alleles to the filtered alleles file
        :param gene: String of the gene name currently being analysed
        :param allele_path: String of the absolute path to the folder
        containing the allele database
        :param records_filter: List of SeqIO sequence objects for alleles to
        be added to the filtered alleles file
        """
        # Set the name and path of the filtered alleles file
        filtered_allele_file = os.path.join(
            allele_path,
            f'{gene}_filtered.txt'
        )
        # Append the filtered alleles to the file
        with open(
                filtered_allele_file,
                'a+',
                encoding='utf-8') as filtered_alleles:
            # Iterate and write all the records in the list of filtered records
            for record in records_filter:
                SeqIO.write(
                    sequences=record,
                    handle=filtered_alleles,
                    format='fasta'
                )

    def aa_profile(self):
        """
        Create the amino acid profile
        """
        # Create a list to store profiles containing alleles that have been
        # quality filtered
        filtered_list = set()
        # Initialise a dictionary to store sequenceType:geneName_alleleID
        filtered_dict = {}
        # Initialise a string to hold the profile
        profile_str = ''
        # Iterate through all the sequence types in the profile_data dictionary
        for seq_type in sorted(int(st) for st in self.profile_data):
            # Create a string to store the allele information
            allele_string = str()
            # Iterate through all the gene: allele entries in the nested
            # dictionary
            for gene_name, allele in self.profile_data[str(seq_type)].items():
                try:
                    # Extract the linked amino acid allele from the dictionary
                    allele_string += self.allele_links[
                        gene_name][str(allele)] + '\t'
                # If the gene name + allele ID is not in the allele link
                # dictionary, add the sequence type to the set
                except KeyError:
                    filtered_list.add(seq_type)
                    # Initialise the sequence type key in the filtered
                    # dictionary as required
                    if seq_type not in filtered_dict:
                        filtered_dict[seq_type] = []
                    # Update the dictionary with the filtered
                    # geneName_alleleIdentifier
                    filtered_dict[seq_type].append(f'{gene_name}_{allele}')
            if allele_string:
                # Add the sequence type to the profile string
                profile_str += f'{str(seq_type)}\t{allele_string}'
                # Remove trailing whitespace and add a newline for proper
                # formatting
                profile_str = profile_str.rstrip()
                profile_str += '\n'
        # Create the amino acid profile file
        with open(self.aa_profile_file, 'w', encoding='utf-8') as aa_profile:
            # Create a string of tab-delimited gene names to be used in
            # the header
            names = '\t'.join(sorted(list(self.gene_names)))
            # Write the header string
            aa_profile.write(f'ST\t{names.rstrip()}\n')
            # Write the tab-delimited profile string
            aa_profile.write(profile_str)
        # Create the names file
        with open(self.gene_name_file, 'w', encoding='utf-8') as gene_file:
            # Write the names to file
            gene_file.write('\n'.join(sorted(list(self.gene_names))))
        if filtered_list:
            filtered_list = sorted(list(filtered_list))
            # Remove the filtered profiles from the profile files
            Translate.filter_profiles(
                filtered_list=filtered_list,
                profile_file=self.profile_file,
            )
            # Write the notes to file
            Translate.filter_notes(
                filtered_list=filtered_list,
                filtered_dict=filtered_dict,
                profile_path=os.path.dirname(self.profile_file)
            )
            self.profile_data = read_profile(profile_file=self.profile_file)

    @staticmethod
    def filter_profiles(
            filtered_list: set,
            profile_file: str):
        """
        Remove filtered profiles from the profile file. Write them to the
        filtered profile file
        :param filtered_list: Set of filtered sequence types
        :param profile_file: String of the absolute path to the profile file
        """
        # Initialise a list to store the filtered sequence types
        filtered_rows = []
        # Initialise a list to store the unfiltered sequence types
        keep_rows = []
        # Read in the profile files
        with open(profile_file, 'r', encoding='utf-8') as profiles:
            # Iterate over all the profiles in the profile files
            for row in profiles:
                # Extract the sequence type from the row
                seq_type = row.split('\t')[0]
                # Check if it is the header, which starts with 'ST'
                if seq_type == 'ST':
                    # Filter the header
                    filtered_rows.append(row)
                # Check to see if the sequence type of the row is present in
                # the list of filtered sequence types
                if str(seq_type) in [
                        str(filtered) for filtered in filtered_list]:
                    # Add the row to the list of filtered rows
                    filtered_rows.append(row)
                # Unfiltered rows are added to the list of unfiltered rows
                else:
                    keep_rows.append(row)
        # Extract the absolute path of the folder in which the profile file
        # is located
        profile_path = os.path.dirname(profile_file)
        # Overwrite the profile file with the unfiltered profiles
        with open(os.path.join(
                profile_path,
                'profile.txt'
                ), 'w', encoding='utf-8') as updated_profile:
            updated_profile.write(''.join(keep_rows))
        # Overwrite the filtered profile file with the filtered profiles
        with open(os.path.join(
                profile_path,
                'filtered_profiles.txt'
                ), 'w', encoding='utf-8') as filtered_profiles:
            filtered_profiles.write(''.join(filtered_rows))

    @staticmethod
    def filter_notes(
            filtered_list: set,
            filtered_dict: dict,
            profile_path: str):
        """
        Write the notes regarding profile filtering to file
        :param filtered_list: Set of filtered sequence types
        :param filtered_dict: Dictionary of sequenceType:geneName_alleleID
        :param profile_path: String of the absolute path to the folder
        containing the profile file
        """
        # Set the absolute path of the file containing the profile filtering
        # notes
        filtered_notes = os.path.join(profile_path, 'filtering_notes.txt')
        # Write the notes to file
        with open(filtered_notes, 'w', encoding='utf-8') as notes:
            # Create the header
            notes.write('SequenceType\tFilteringNote\n')
            # Create a new line for each filtered profile
            for seq_type in filtered_list:
                # Create a sting of the list of filtered alleles present in
                # this profile
                gene_allele = ';'.join(filtered_dict[seq_type])
                # Write the sequence type and all the missing alleles to the
                # note
                notes.write(f'{seq_type}\t{gene_allele} missing\n')

    def profile_link(self):
        """
        Link nucleotide and amino profiles
        """
        # Initialise a dictionary to store to number of matches between a
        # query profile and the profiles in the database
        match_score = {}
        # Iterate over all the profiles in the profile file
        for seq_type, gene_dict in self.profile_data.items():
            # Initialise the seq_type key in the dictionary
            self.profile_matches[seq_type] = set()
            match_score[seq_type] = {}
            # Iterate over all the gene name, allele ID combinations in the
            # profile dictionary
            for gene, allele in gene_dict.items():
                # Iterate over all the profiles in the amino acid profile file
                for aa_st, aa_gene_dict in self.aa_profile_data.items():
                    # Use the gene name to extract the amino acid allele ID
                    # from the profile file.
                    # Also extract the amino acid allele ID from the linking
                    # dictionary with the gene name and
                    # nucleotide allele ID. Check if they match
                    if aa_gene_dict[gene] == self.allele_links[gene][allele]:
                        # Initialise the amino acid sequence type in
                        # the dictionary
                        if aa_st not in match_score[seq_type]:
                            match_score[seq_type][aa_st] = 0
                        # Increment the number of matches to the profile
                        match_score[seq_type][aa_st] += 1
        # Iterate over all the matches to the profiles in the profile file
        for seq_type, aa_st_dict in match_score.items():
            # Iterate over the amino acid sequence type matches
            for aa_st, matches, in aa_st_dict.items():
                # Check if the number of matches observed is equal to the
                # required number of matches (one for each gene)
                if matches == len(self.gene_names):
                    # Update the dictionary of matches with the amino acid
                    # sequence type: nucleotide sequence type
                    self.profile_matches[aa_st].add(seq_type)

    def link_file(self):
        """
        Write linking details between nucleotide and amino acid profiles
        to file
        """
        # Open the profile link file
        with open(self.aa_nt_profile_link_file, 'w', encoding='utf-8') as link:
            # Write the header information
            link.write('aa_seq_type\tnt_seq_types\n')
            # Iterate over all the profile matches in the dictionary
            for seq_type, match_set in self.profile_matches.items():
                # Check the length of the match
                if len(match_set) == 1:
                    # With a single match, convert the set to a string
                    nt_st = ''.join(match_set)
                    link.write(f'{seq_type}\t{nt_st}\n')
                # Multiple profile matches
                else:
                    # Semicolon-join the set of matches. Since they are
                    # strings, first typecase to int for proper sorting before
                    # converting back to strings
                    nt_st = ';'.join(
                        str(nt_st) for nt_st in sorted(
                            int(nt_st) for nt_st in match_set
                        )
                    )
                    link.write(f'{seq_type}\t{nt_st}\n')

    def copy_profile(self):
        """
        Copy the reduced profile file from the profile folder to the base
        aa_profile folder
        """
        # Use shutil to copy and rename the profile file to the root of
        # the report path
        shutil.copyfile(
            src=os.path.join(self.report_path, 'profile', 'profile.txt'),
            dst=os.path.join(self.report_path, 'profile.txt')
        )

    def __init__(
            self,
            path: str,
            profile: str,
            report_path: str = os.path.join(os.getcwd(), 'aa_profile'),
            translated_path: str = os.path.join(os.getcwd(), 'aa_alleles'),
            length_dict: dict = None):
        """
        Constructs all the necessary attributes for the Translate object.

        Parameters:
        path (str): Path to the file
        profile (str): Profile file. If not provided, it will look for
        'nt_profile/profile.txt' in the path
        report_path (str): Path to the report. Default is 'aa_profile' in
        the current working directory
        translated_path (str): Path to the translated file. Default is
        'aa_alleles' in the current working directory
        length_dict (dict): Dictionary of lengths. Default is None
        """

        logging.info('Welcome to the allele translator!')
        self.path = pathfinder(path=path)

        # Set profile file
        if profile:
            if not os.path.isfile(profile):
                self.profile_file = os.path.join(
                    self.path,
                    'nt_profile',
                    'profile.txt'
                )
            else:
                self.profile_file = profile
            try:
                assert os.path.isfile(self.profile_file)
            except AssertionError as exc:
                logging.error(
                    'Cannot locate the required profile file: %s. '
                    'Please ensure that the file name and path of your file '
                    'is correct', self.profile_file
                )
                raise SystemExit from exc
        else:
            self.profile_file = None

        # Set sequence files
        self.sequence_files = glob(os.path.join(self.path, '*.fasta'))
        try:
            assert self.sequence_files
        except AssertionError as exc:
            logging.error(
                'Could not locate alleles in provided allele path: %s',
                self.path
            )
            raise SystemExit from exc

        # Set report path
        self.report_path = pathfinder(path=report_path)
        make_path(inpath=self.report_path)

        # Set translated path
        self.translated_path = pathfinder(path=translated_path)
        make_path(inpath=self.translated_path)

        # Set notes path
        self.notes_path = os.path.join(self.translated_path, 'notes')
        make_path(inpath=self.notes_path)

        # Initialize other attributes
        self.length_dict = length_dict
        self.allele_dict = {}
        self.profile_data = {}
        self.allele_links = {}
        self.aa_profile_file = os.path.join(
            self.report_path,
            'aa_full_profile.txt'
        )
        self.gene_names = set()
        self.gene_name_file = os.path.join(
            self.report_path, 'gene_names.txt'
        )
        self.aa_profile_data = {}
        self.profile_matches = {}
        self.aa_nt_profile_link_file = os.path.join(
            self.report_path,
            'aa_nt_profile_links.tsv'
        )


def cli():
    """
    Collect the arguments, create an object, and run the script
    """
    # Parser for arguments
    parser = ArgumentParser(
        description='Translate allele files in nucleotide format to '
        'amino acid. Remove duplicates. Keep notes.'
    )
    parser.add_argument(
        '-p', '--path',
        required=True,
        help='Specify path containing allele files.'
    )
    parser.add_argument(
        '--profile',
        action='store_true',
        help='Optionally parse the nucleic acid profile, and create the '
        'corresponding reduced amino acid profile. The profile must be named '
        'profile.txt and be located in nt_profile folder in the path'
    )
    parser.add_argument(
        '--report_path',
        default=os.path.join(os.getcwd(), 'aa_profile'),
        help='Optionally provide the name (and path, if desired) of the '
        'folder into which the amino acid profile and related files are to be '
        'written. Default folder is "aa_profile" in your current '
        'working directory'
    )
    parser.add_argument(
        '--translated_path',
        default=os.path.join(os.getcwd(), 'aa_alleles'),
        help='Optionally provide the name (and path, if desired) of the '
        'folder into which the amino acid alleles and notes are to be written.'
        ' Default folder is "aa_alleles" in your current working directory'
    )
    # Get the arguments into an object
    arguments = parser.parse_args()
    SetupLogging(debug=True)
    translate = Translate(
        path=arguments.path,
        profile=arguments.profile,
        report_path=arguments.report_path,
        translated_path=arguments.translated_path
    )
    translate.main()
    logging.info('Allele translation complete!')
    # Prevent the arguments being printed to the console (they are returned in
    # order for the tests to work)
    sys.stderr = open(os.devnull, 'w', encoding='utf-8')
    return arguments


if __name__ == '__main__':
    cli()
