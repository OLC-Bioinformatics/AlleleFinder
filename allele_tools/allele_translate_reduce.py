#!/usr/bin/env python
from olctools.accessoryFunctions.accessoryFunctions import make_path, SetupLogging
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from argparse import ArgumentParser
from glob import glob
import logging
import os


class Translate(object):

    def main(self):
        self.load_alleles()
        self.parse_alleles()

    def load_alleles(self):
        """
        Use SeqIO to read in all the gene sequences
        """
        for allele_file in self.sequence_files:
            gene_name = os.path.splitext(os.path.basename(allele_file))[0]
            self.allele_dict[gene_name] = SeqIO.to_dict(SeqIO.parse(allele_file, 'fasta'))

    def parse_alleles(self):
        """
        Parse the allele files to translate the amino acid sequence using BioPython. Write the amino acid sequence to
        file. Store the allele name in the notes. Find duplicates, and link the nucleotide allele name to the amino
        acid allele name
        """
        logging.info('Translating and parsing alleles')
        for gene_name, allele_dict in self.allele_dict.items():
            logging.info(f'Processing {gene_name}')
            # Initialise a dictionary to store the translated allele sequence: allele name
            seq_allele = dict()
            # Open the file to store the translated alleles
            with open(os.path.join(self.translated_path, '{gn}_aa.fasta'.format(gn=gene_name)), 'w') as aa_alleles:
                # Open the notes file
                with open(os.path.join(self.notes_path, '{gn}_notes.txt'.format(gn=gene_name)), 'w') as notes:
                    # Create the header for the notes file
                    notes.write('nt_allele\taa_allele\tnote\n')
                    # Iterate through all the alleles in the dictionary
                    for allele, details in allele_dict.items():
                        # Calculate the translated sequence
                        translated_allele = details.seq.translate()
                        # Determine if this amino acid allele is new
                        if str(translated_allele) not in seq_allele:
                            # Add the string of the amino acid sequence to
                            seq_allele[str(translated_allele)] = allele
                            # Create a SeqRecord of the translated allele
                            seq_record = SeqRecord(seq=translated_allele,
                                                   id=allele,
                                                   name=str(),
                                                   description=str())
                            # Write the SeqRecord to file
                            SeqIO.write(sequences=seq_record,
                                        handle=aa_alleles,
                                        format='fasta')
                            # Update the notes with the allele naming information
                            notes.write('{nt}\t{aa}\n'.format(nt=allele,
                                                              aa=allele))
                        # Amino acid allele already exists
                        else:
                            # Extract the allele name corresponding to the translated sequence
                            aa_allele = seq_allele[str(translated_allele)]
                            # Update the notes, including that this allele is a duplicate, and a pointer to the original
                            notes.write('{nt}\t{aa}\tduplicate\n'.format(nt=allele,
                                                                         aa=aa_allele))

    def __init__(self, path):
        logging.info('Welcome to the allele translator!')
        if path.startswith('~'):
            self.path = os.path.abspath(os.path.expanduser(os.path.join(path)))
        else:
            self.path = os.path.abspath(os.path.join(path))
        self.sequence_files = glob(os.path.join(self.path, '*.fasta'))
        self.translated_path = os.path.join(self.path, 'translated')
        self.notes_path = os.path.join(self.path, 'notes')
        make_path(inpath=self.translated_path)
        make_path(inpath=self.notes_path)
        self.allele_dict = dict()


def cli():
    # Parser for arguments
    parser = ArgumentParser(description='Translate allele files in nucleotide format to amino acid. '
                                        'Remove duplicates. Keep notes.')
    parser.add_argument('-p', '--path',
                        required=True,
                        help='Specify path containing allele files.')
    # Get the arguments into an object
    arguments = parser.parse_args()
    SetupLogging(debug=True)
    translate = Translate(path=arguments.path)
    translate.main()
    logging.info('Allele translation complete!')


if __name__ == '__main__':
    cli()
