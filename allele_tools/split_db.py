#!/usr/bin/env python3
"""
This script processes a FASTA file and outputs each record into separate files
based on the FASTA header. Records with headers starting with 'Stx1' or 'Stx2'
(case insensitive) are written to separate files in 'stx1' and 'stx2'
subdirectories respectively.

The script uses the BioPython SeqIO module to parse the FASTA file and argparse
to handle command line arguments.

The script defines a FastaProcessor class with methods to process the file and
write records to separate files. The main function parses the command line
arguments and creates an instance of FastaProcessor to process the file.

Usage:
    python script_name.py filename

Where:
    filename is the name of the FASTA file to process.
"""

__author__ = 'adamkoziol'

# Standard imports
import argparse
import os

# Third party imports
from Bio import SeqIO


class FastaProcessor:
    'Class to process a fasta file and output records into separate files'
    def __init__(
        self,
        *,  # Enforce keyword-only arguments
        filename: str,
        split_stx: bool
    ) -> None:
        """
        Initialize FastaProcessor with a filename.

        :param filename: The name of the fasta file to process.
        :param split_stx: If True, split records into 'stx1' and 'stx2'
        """
        self.filename = filename
        self.split_stx = split_stx

    def process(self) -> None:
        """
        Reads the file and writes each record to a separate file.
        Records with headers starting with 'stx1' or 'stx2' are written to
        separate folders.
        """
        # Open the fasta file in read mode
        with open(self.filename, 'r', encoding='utf-8') as file:
            # Iterate over each record in the fasta file
            for record in SeqIO.parse(file, 'fasta'):
                # Get the description of the record, which serves as the header
                header = record.description
                # If split_stx is False, write the record to a file in the
                # base_dir folder with the header as the filename
                if not self.split_stx:
                    # Write the record to a file in the "split" directory
                    self._write_record(
                        folder='split',
                        header=header,
                        record=record
                    )
                    continue
                # Check if the header starts with 'stx1' (case insensitive)
                if header.lower().startswith('stx1'):
                    # If it does, write the record to a file in the 'stx1'
                    # subfolder. The header is passed to the _write_record
                    # method to be used as the filename
                    self._write_record(
                        folder='stx1',
                        header=header,
                        record=record
                    )
                # Check if the header starts with 'stx2' (case insensitive)
                elif header.lower().startswith('stx2'):
                    # If it does, write the record to a file in the 'stx2'
                    # subfolder. The header is passed to the _write_record
                    # method to be used as the filename
                    self._write_record(
                        folder='stx2',
                        header=header,
                        record=record
                    )

    def _write_record(
        self,
        *,  # Enforce keyword arguments
        folder: str,
        header: str,
        record: SeqIO.SeqRecord
    ) -> None:
        """
        Writes a record to a file in specified folder.

        :param folder: The name of the subfolder to write the record to.
        :param header: The header of the record.
        :param record: The record to write.
        """
        # Get the directory of the input file
        base_dir = os.path.dirname(self.filename)

        # Combine the base directory with the subfolder name
        full_folder_path = os.path.join(base_dir, folder)

        # Create the subfolder if it doesn't exist
        os.makedirs(full_folder_path, exist_ok=True)

        # Remove any character including and after the |
        # This is done by splitting the header at the | character and
        # taking the first part.The split method returns a list of parts,
        # and index 0 is the part before the |
        header = header.split('|')[0]

        # Combine the folder path with the sanitized header and the .fasta
        # extension to get the output file path
        output_file_path = os.path.join(full_folder_path, f'{header}.fasta')

        # Open the output file in write mode
        with open(output_file_path, 'w', encoding='utf-8') as output_file:
            # Write the record to the output file in fasta format
            SeqIO.write(record, output_file, 'fasta')


def cli() -> None:
    'Main function to process command line arguments and run the processor'
    parser = argparse.ArgumentParser(description='Process a fasta file.')
    parser.add_argument(
        '--filename',
        type=str,
        help='The name of the fasta file to process.'
    )
    parser.add_argument(
        '--stx',
        action='store_true',
        help='Split the fasta file into stx1 and stx2 subfolders.'
    )
    args = parser.parse_args()

    processor = FastaProcessor(
        filename=args.filename,
        split_stx=args.stx
    )
    processor.process()


if __name__ == '__main__':
    cli()
