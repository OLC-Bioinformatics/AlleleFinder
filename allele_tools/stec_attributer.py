#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import make_path, SetupLogging
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from argparse import ArgumentParser
from glob import glob
import hashlib
import logging
import os
__author__ = 'adamkoziol'


class Attribute(object):

    def main(self):
        """
        Run the methods in the appropriate order
        """
        self.read_targets()
        self.read_alleles()
        self.attribute_targets()
        self.attribute_alleles()
        self.write_alleles()
        self.align_alleles()

    def read_targets(self):
        """
        Use SeqIO to read in all the target records into a dictionary
        """
        logging.info('Parsing target files')
        self.targetrecords = SeqIO.to_dict(SeqIO.parse(self.targetfile, 'fasta'))

    def read_alleles(self):
        """
        Use SeqIO to read in the raw sequence of all the alleles into a set
        """
        logging.info('Parsing allele files')
        # Iterate through all the allele files
        for allele_file in self.allelefiles:
            for record in SeqIO.parse(allele_file, 'fasta'):
                self.alleleset.add(str(record.seq))

    def attribute_targets(self):
        """
        Find the alleles that correspond to the target sequence
        """
        logging.info('Adding target alleles')
        for allele, record in self.targetrecords.items():
            # Add the target record to the list of the attributed alleles
            self.attributedalleles.append(record)
            # Add the raw sequence of the target to the set of processed allele sequences
            self.complete.add(str(record.seq))

    def attribute_alleles(self):
        """
        Find and name alleles without a corresponding target sequence. Create SeqRecords for these alleles
        """
        logging.info('Adding non-target allele sequences')
        for allele in self.alleleset:
            # Ensure that the allele is not already in the set of processed alleles
            if allele not in self.complete:
                # Use MD5 hashing of the UTF-8 encoded string of the allele sequence.
                # Extract only the first 10 characters
                hash_str = hashlib.md5(str(allele).encode('utf-8')).hexdigest()[:10]
                # Ensure that the allele in not in the set of alleles detected in the incorrect genera
                header = '{gene}_{hash_str}'.format(gene=str(self.gene),
                                                    hash_str=hash_str)
                # Create a SeqRecord for each allele
                fasta = SeqRecord(Seq(allele),
                                  # Without this, the header will be improperly formatted
                                  description=str(),
                                  # Use the gene name_iteration as the allele name
                                  id=header)
                self.attributedalleles.append(fasta)

    def write_alleles(self):
        """
        Use SeqIO to write the list of alleles to file
        """
        logging.info('Writing alleles to {fasta}'
                     .format(fasta=self.unaligned_alleles))
        with open(self.unaligned_alleles, 'w') as fasta:
            SeqIO.write(self.attributedalleles, fasta, 'fasta')

    def align_alleles(self):
        """

        """
        logging.info('Aligning extracted alleles with MUSCLE')
        cline = MuscleCommandline(input=self.unaligned_alleles,
                                  out=self.aligned_alleles)
        cline()

    def __init__(self, allelepath, targetpath, reportpath, gene):
        logging.info('Welcome to the CFIA Allele Attributer')
        if allelepath.startswith('~'):
            self.allelepath = os.path.abspath(os.path.expanduser(os.path.join(allelepath)))
        else:
            self.allelepath = os.path.abspath(os.path.join(allelepath))
        self.allelefiles = glob(os.path.join(self.allelepath, '*.tfa'))
        self.alleleset = set()
        if targetpath.startswith('~'):
            self.targetpath = os.path.expanduser(os.path.abspath(os.path.join(targetpath)))
        else:
            self.targetpath = os.path.abspath(os.path.join(targetpath))
        self.targetfile = sorted(glob(os.path.join(self.targetpath, '*.fasta')))[0]
        self.targetrecords = dict()
        if reportpath.startswith('~'):
            self.reportpath = os.path.expanduser(os.path.abspath(os.path.join(reportpath, 'reports')))
        else:
            self.reportpath = os.path.abspath(os.path.join(reportpath, 'reports'))
        make_path(self.reportpath)
        self.gene = gene
        self.unaligned_alleles = os.path.join(self.reportpath, '{gene}_unaligned_alleles.fasta'.format(gene=self.gene))
        self.aligned_alleles = os.path.join(self.reportpath, '{gene}_aligned_alleles.fasta'.format(gene=self.gene))
        self.attributedalleles = list()
        self.complete = set()


def cli():
    # Parser for arguments
    parser = ArgumentParser(description='Finds the target sequences in allele files. Useful if you have an allele '
                                        'database, and want to attribute subtypes e.g. STEC subtyping to your newly '
                                        'expanded alleles')
    parser.add_argument('-a', '--allelepath',
                        required=True,
                        help='Name and path of folder containing generated allele files')
    parser.add_argument('-t', '--targetpath',
                        required=True,
                        help='Name and path of folder containing sequencing target sequences')
    parser.add_argument('-r', '--reportpath',
                        required=True,
                        help='Name and path of folder in which reports are to be created')
    parser.add_argument('-g', '--gene',
                        required=True,
                        choices=['stx1A', 'stx1B', 'stx2A', 'stx2B'],
                        help='Name of gene being profiled')
    SetupLogging()
    arguments = parser.parse_args()
    # Run the pipeline
    attributer = Attribute(allelepath=arguments.allelepath,
                           targetpath=arguments.targetpath,
                           reportpath=arguments.reportpath,
                           gene=arguments.gene)
    attributer.main()
    logging.info('Allele Attribution complete!')


if __name__ == '__main__':
    cli()
