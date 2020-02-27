#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import combinetargets, make_path, SetupLogging
from genemethods.geneseekr.geneseekr import GeneSeekr
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIWWW
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from argparse import ArgumentParser
from threading import Thread
from csv import DictReader
from queue import Queue
import multiprocessing
from glob import glob
from time import time
import logging
import os

__author__ = 'adamkoziol'


class AlleleFinder(object):

    def main(self):
        """
        Run the required methods in the appropriate order
        """
        self.record_extraction()
        self.parameters()
        self.prep_targets()
        self.populate_local_dict()
        self.local_blast()
        self.parse()
        self.create_allele_file()

    def record_extraction(self):
        """
        Parse the input FASTA file, and create a dictionary of header: sequence for each entry
        """
        for record in SeqIO.parse(self.file, 'fasta'):
            # Replace and dashes in the record.id with underscores
            record.id = record.id.replace('-', '_')
            self.records[record.id] = str(record.seq)

    def parameters(self):
        """
        Optimise BLAST input parameters for sequences below 50 bp
        """
        logging.info('Optimising parameters for targets')
        for record, record_seq in self.records.items():
            # Initialise the set of alleles with the input sequence
            self.alleleset[record] = set()
            self.illegal_alleleset[record] = set()
            self.alleleset[record].add(Seq(record_seq))
            if len(record_seq) < 70:
                self.expect[record] = 0.001
                self.word_size[record] = 4
                self.filter_low_complexity[record] = 'F'
            else:
                self.expect[record] = 0.00001
                self.word_size[record] = 11
                self.filter_low_complexity[record] = 'T'
            logging.debug((record, self.expect[record], self.word_size[record], self.filter_low_complexity[record]))

    def remote_blast(self):
        """
        Create threads for each BLAST database: output file combination. Run multi-threaded BLAST
        """
        logging.info('Running remote nr BLAST')
        # Create and start threads for each fasta file in the list
        for i in range(self.cpus):
            # Send the threads to makeblastdb
            threads = Thread(target=self.run_remote_blast, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for record, record_seq in self.records.items():
            result_file = os.path.join(self.reportpath, '{record}_nr.xml'.format(record=record))
            # Update the dictionary of BLAST output files with the database, and output file
            self.blast_outputs[record] = [{'nr': result_file}]
            self.queue.put((result_file, record, record_seq))
        self.queue.join()

    def run_remote_blast(self):
        """
        Run NCBIWWW BLAST against the appropriate remote database and save results to a local file
        """
        while True:
            result_file, record, record_seq = self.queue.get()
            if not os.path.isfile(result_file):
                # Create the BLAST request with the appropriate database
                result = NCBIWWW.qblast(program='blastn',
                                        database='nr',
                                        format_type='XML',
                                        hitlist_size=1000000,
                                        expect=self.expect[record],
                                        word_size=self.word_size[record],
                                        filter=self.filter_low_complexity[record],
                                        sequence=Seq(record_seq))
                # Write the results to file
                with open(result_file, 'w') as results:
                    results.write(result.read())
            self.queue.task_done()

    def prep_targets(self):
        """

        """
        target_files = sorted(glob(os.path.join(self.path, '*.tfa')))
        logging.debug('Target files: {target_files}'.format(target_files=','.join(target_files)))
        combined_targets = os.path.join(self.path, 'combinedtargets.fasta')
        if not os.path.isfile(combined_targets):
            combinetargets(targets=target_files,
                           targetpath=self.path)
        #
        for record in SeqIO.parse(combined_targets, 'fasta'):
            # Replace and dashes in the record.id with underscores
            record.id = record.id.replace('-', '_')
            self.genus_alleles['custom'] = list()
            self.gene_dict[record.id] = 'custom'
            self.mismatches[record.id] = int(0.25 * len(record.seq))

        GeneSeekr.makeblastdb(fasta=combined_targets)

    def populate_local_dict(self):
        """
        Populate self.local_dict with the
        """
        seq_files = sorted(glob(os.path.join(self.local_path, '*.tfa')))
        logging.debug('Query files: {seq_files}'.format(seq_files=','.join(seq_files)))
        combined_sequences = os.path.join(self.local_path, 'combinedtargets.fasta')
        if not os.path.isfile(combined_sequences):
            combinetargets(targets=seq_files,
                           targetpath=self.local_path)
        # Populate self.local_dict with the 'custom' genus and the name and path of the combined targets file
        self.local_dict['custom'] = combined_sequences
        # Create the BLAST database from the combined file
        GeneSeekr.makeblastdb(fasta=combined_sequences)

    def create_local_db(self):
        """
        Calls BLAST database creation method for local assembly and refseq databases
        """
        for db_type, output_file in sorted(self.local_dict.items()):
            # Set the extension to look for with glob
            if 'assemblies' in output_file:
                glob_command = os.path.join(self.local_path, '{db}*.fasta'.format(db=db_type.split('_')[0]))
            else:
                # glob_extension = 'fna'
                glob_command = os.path.join(self.local_path, '*.fna')
            # Get a list of files for the appropriate database
            db_files = glob(glob_command)
            if db_files:
                # Create the database as necessary
                self.run_create_local_db(db_type=db_type,
                                         db_files=db_files,
                                         output_file=output_file)

    @staticmethod
    def run_create_local_db(db_type, db_files, output_file):
        """
        Create a combined file from all .fasta files in the local database path. Use the combined file to make a
        BLAST database
        :param db_type: string of the database type: either assemblies or refseq
        :param db_files: list of FASTA files to use in creating the BLAST database
        :param output_file: string of path and file name of the combined FASTA file
        """
        logging.info('Reading {db_type} FASTA files'.format(db_type=db_type))
        if not os.path.isfile(output_file):
            records = list()
            # A set to ensure that all record IDs are unique
            ids = set()
            # Use SeqIO to create a properly-formatted multifasta file
            for db_file in sorted(db_files):
                iterator = SeqIO.parse(db_file, "fasta")
                for record in iterator:
                    # Extract the sequence record from each entry in the multifasta
                    # Replace and dashes in the record.id with underscores
                    record.id = record.id.replace('-', '_')
                    # Remove and dashes or 'N's from the sequence data - makeblastdb can't handle sequences
                    # with gaps
                    # noinspection PyProtectedMember
                    record.seq._data = record.seq._data.replace('-', '').replace('N', '')
                    # Clear the name and description attributes of the record
                    record.name = ''
                    record.description = ''
                    # Only add the record to the list of records if the id is unique (i.e. not in the set of ids)
                    if record.id not in ids:
                        # Add the id to the list of ids
                        ids.add(record.id)
                        # Add the record to the list of records
                        records.append(record)
            logging.info('Merging {db_type} FASTA files'.format(db_type=db_type))
            # Write each record to the combined file
            SeqIO.write(records, output_file, 'fasta')
        # Create the BLAST database from the combined file
        GeneSeekr.makeblastdb(fasta=output_file)

    def local_blast(self):
        """
        Calls BLAST method for local assembly and refseq databases
        """
        for db_type, db_file in sorted(self.local_dict.items()):
            if os.path.isfile(db_file):
                for record, record_seq in self.records.items():
                    # Set the name of the BLAST output file
                    local_output = os.path.join(self.reportpath, '{record}_{db_type}.csv'.format(record=record,
                                                                                                 db_type=db_type))
                    # Run BLAST with the appropriate parameters
                    self.run_local_blast(db_type=db_type,
                                         local_db=db_file,
                                         local_output=local_output,
                                         record=record,
                                         record_seq=record_seq)
                    # Add the local database to the dictionary
                    try:
                        self.blast_outputs[record].append({db_type: local_output})
                    except KeyError:
                        self.blast_outputs[record] = list()
                        self.blast_outputs[record].append({db_type: local_output})

    def run_local_blast(self, db_type, local_db, local_output, record, record_seq):
        """
        BLAST the probe file against the local database
        :param db_type: string of the database type: either assemblies or refseq
        :param local_db: name and path of the BLAST database file
        :param local_output: BLAST output file
        :param record: Name of the record
        :param record_seq: Sequence of the record
        """
        logging.info('Running {record} local {db_type} BLAST'.format(record=record,
                                                                     db_type=db_type))
        # BLAST command line call.
        blastn = NcbiblastnCommandline(db=os.path.splitext(local_db)[0],
                                       num_alignments=100000000,
                                       evalue=self.expect[record],
                                       num_threads=self.cpus,
                                       task='blastn',
                                       outfmt=self.outfmt,
                                       perc_identity=75,
                                       word_size=self.word_size[record],
                                       out=local_output
                                       )
        if not os.path.isfile(local_output):
            # Run BLAST - supply the record sequence as stdin, so BLAST doesn't look for an input file
            blastn(stdin=record_seq)

    def load_genera(self):
        """
        Parse the combined metadata file, and extract seqid: calculated genus key: value pairs
        """
        with open(self.metadata_file, 'r') as metadata:
            # Skip the header
            next(metadata)
            # Iterate through all the entries in the metadata report
            for line in metadata:
                # Split on commas
                data = line.split(',')
                try:
                    # Extract the genus from the ReferenceGenome column e.g.
                    # Escherichia_coli_O26_H11_11368_DNA_260853213 becomes Escherichia
                    genus = data[6].split('_')[0]
                    strain = data[0]
                    # Ensure that the genus is one of Escherichia, Listeria, or Salmonella, and that the strain is
                    # one of our locally-assembled strain with a SEQID similar to: 2013-SEQ-0001
                    if genus in self.genera and 'SEQ' in strain:
                        self.strain_genera[strain] = genus
                # Ignore strains without results
                except IndexError:
                    pass

    def parse(self):
        """
        Call the report parsing method for all the BLAST output files
        """
        logging.info('Parsing outputs')
        # Call parse_report for every file
        for database, result_files in sorted(self.blast_outputs.items()):
            for dictionary in result_files:
                for db, result_file in dictionary.items():
                    if os.path.isfile(result_file):
                        self.parse_report(database, db, result_file)

    def parse_report(self, gene_name, db, result_file):
        """
        Extract matches from the BLAST results
        :param gene_name: Gene name
        :param db: database name
        :param result_file: name and path of the result file to be parsed
        """
        # Read in the BLAST results
        try:
            # Open the sequence profile file as a dictionary
            logging.info('Parsing {gene_name} {db} report'.format(gene_name=gene_name,
                                                                  db=db))
            blastdict = DictReader(open(result_file), fieldnames=self.fieldnames, dialect='excel-tab')
            # Go through each BLAST result
            for row in blastdict:
                if '-' not in row['subject_sequence']:
                    deviation = int(0.05 * len(self.records[gene_name]))
                    # Only retrieve sequences that are as long as the query sequence, and do not have gaps
                    if len(self.records[gene_name]) - deviation <= len(row['subject_sequence']) <= \
                            len(self.records[gene_name]) + deviation:
                        # Try to use the gene-specific number of mismatches
                        try:
                            if int(row['positives']) >= len(self.records[gene_name]) - self.mismatches[gene_name]:
                                # Create a Seq object to add to the set
                                self.alleleset[gene_name].add(Seq(row['subject_sequence'].replace('-', '')))
                        # Do not allow for more than five mismatches
                        except KeyError:
                            if int(row['positives']) >= len(self.records[gene_name]) - 5:
                                # Create a Seq object to add to the set
                                self.alleleset[gene_name].add(Seq(row['subject_sequence'].replace('-', '')))
        except FileNotFoundError:
            pass

    def create_allele_file(self):
        """
        Create and populate the allele file with all alleles in the set
        """
        for record, record_seq in sorted(self.records.items()):
            logging.info('Creating {record} allele FASTA file'.format(record=record))
            # Set the name of the allele file
            allelefile = '{allele_path}_alleles.tfa'.format(allele_path=os.path.join(self.allelepath, record))
            with open(allelefile, 'w') as alleles:
                fastalist = list()
                # Iterate through the set of alleles
                count = 0
                for alleleseq in self.alleleset[record]:
                    # Ensure that the allele in not in the set of alleles detected in the incorrect genera
                    header = '{record}_{count}'.format(record=str(record),
                                                       count=count)
                    # Create a SeqRecord for each allele
                    fasta = SeqRecord(alleleseq,
                                      # Without this, the header will be improperly formatted
                                      description='',
                                      # Use the gene name_iteration as the allele name
                                      id=header)
                    fastalist.append(fasta)
                    self.all_alleles.append(fasta)
                    # Add the fasta record to the appropriate genus
                    self.genus_alleles[self.gene_dict[record]].append(fasta)
                    count += 1
                # Write the alleles to file
                SeqIO.write(fastalist, alleles, 'fasta')
        # Write all the alleles to the combined allele file
        with open(os.path.join(self.allelepath, 'combinedtargets.fasta'), 'w') as combined:
            SeqIO.write(self.all_alleles, combined, 'fasta')
        # Create the genus-specific allele file
        for genus, allele_list in self.genus_alleles.items():
            allele_file = os.path.join(self.allelepath, '{genus}.tfa'.format(genus=genus))
            with open(allele_file, 'w') as genus_file:
                SeqIO.write(allele_list, genus_file, 'fasta')

    def __init__(self, args):
        logging.info('Welcome to the CFIA Allele Finder (CAlF)')
        self.path = os.path.join(args.path)
        self.start = args.starttime
        self.file = os.path.join(self.path, args.filename)
        assert os.path.isfile(self.file), 'Cannot find the supplied FASTA file: {fn}'.format(fn=self.file)
        self.reportpath = os.path.join(self.path, 'reports')
        self.allelepath = os.path.join(self.path, 'alleles')
        self.local_path = os.path.join(args.local_path)
        self.records = dict()
        self.record_parameters = dict()
        self.expect = dict()
        self.word_size = dict()
        self.filter_low_complexity = dict()
        self.blast_outputs = dict()
        self.alleleset = dict()
        self.illegal_alleleset = dict()
        self.strain_genera = dict()
        self.all_alleles = list()
        self.devnull = open(os.devnull, 'wb')
        make_path(self.reportpath)
        make_path(self.allelepath)
        self.cpus = multiprocessing.cpu_count()
        self.queue = Queue()
        # Fields used for custom outfmt 6 BLAST output:
        self.fieldnames = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps',
                           'evalue', 'bit_score', 'subject_length', 'alignment_length',
                           'query_start', 'query_end', 'subject_start', 'subject_end',
                           'query_sequence', 'subject_sequence']
        self.outfmt = "'6 qseqid sseqid positive mismatch gaps evalue bitscore slen length qstart qend sstart send " \
                      "qseq sseq'"
        self.genesippr = args.genesippr
        if self.genesippr:
            self.metadata_file = os.path.join(self.path, args.metadata)
            self.local_dict = {
                '20_assemblies': os.path.join(self.local_path, 'combined_assemblies.tfa'),
                '2013_assemblies': os.path.join(self.local_path, 'combined_2013_assemblies.tfa'),
                '2014_assemblies': os.path.join(self.local_path, 'combined_2014_assemblies.tfa'),
                '2015_assemblies': os.path.join(self.local_path, 'combined_2015_assemblies.tfa'),
                '2016_assemblies': os.path.join(self.local_path, 'combined_2016_assemblies.tfa'),
                '2017_assemblies': os.path.join(self.local_path, 'combined_2017_assemblies.tfa'),
                '2018_assemblies': os.path.join(self.local_path, 'combined_2018_assemblies.tfa'),
                'refseq': os.path.join(self.local_path, 'combined_refseq.tfa')
            }
            self.genera = ['Bacillus', 'Campylobacter, ''Escherichia', 'Listeria', 'Salmonella', 'Staphylococcus',
                           'Vibrio']
            self.gene_dict = {
                'IGS': 'Listeria',
                'hlyA': 'Listeria',
                'inlJ': 'Listeria',
                'invA': 'Salmonella',
                'stn': 'Salmonella',
                'VT1': 'Escherichia',
                'VT2': 'Escherichia',
                'VT2f': 'Escherichia',
                'uidA': 'Escherichia',
                'eae': 'Escherichia',
                'hylA': 'Escherichia',
                'aggR': 'Escherichia',
                'tlh': 'Vibrio'
            }
            self.mismatches = {
                'IGS': 5,
                'hlyA': 7,
                'inlJ': 7,
                'invA': 7,
                'stn': 11,
                'VT1': 7,
                'VT2': 7,
                'VT2f': 7,
                'uidA': 7,
                'eae': 7,
                'hylA': 7,
                'aggR': 7,
                'tlh': 5
            }
            self.genus_alleles = {
                'Bacillus': list(),
                'Campylobacter': list(),
                'Escherichia': list(),
                'Listeria': list(),
                'Salmonella': list(),
                'Staphylococcus': list(),
                'Vibrio': list()
            }
            self.record_extraction()
            self.parameters()
            self.remote_blast()
            self.create_local_db()
            self.local_blast()
            self.load_genera()
            self.parse()
            self.create_allele_file()
        else:
            # TODO self.strain_genera needs to be populated properly
            self.metadata_file = str()
            self.local_dict = dict()
            self.genera = str()
            self.gene_dict = dict()
            self.mismatches = dict()
            self.genus_alleles = dict()


def cli():
    # Parser for arguments
    parser = ArgumentParser(add_help=False)
    parser.add_argument('-p', '--path',
                        required=True,
                        help='Specify path.')
    parser.add_argument('-f', '--filename',
                        required=True,
                        help='Name of file containing probe sequence to search. The file can be a multi-FASTA. The '
                             'header for each sequence must be unique, as it will be used as the name of the gene.'
                             'This file must be located in the supplied path folder.')
    parser.add_argument('-l', '--local_path',
                        required=True,
                        help='Path to folder containing local files to BLAST. These files should have a .tfa extension')
    parser.add_argument('-g', '--genesippr',
                        default=False,
                        action='store_true',
                        help='Enable mode to specifically create alleles for the defined set of genes used in '
                             'the GeneSippr analysis')
    parser.add_argument('-m', '--metadata',
                        help='Name of combined metadata file used to parse the genus of each local assembly. This '
                             'file must be located in the supplied path folder. NOTE: This is only required if '
                             'performing the "GeneSippr-specific" analysis')
    parser.add_argument('-v', '--verbose',
                        default=False,
                        action='store_true',
                        help='Enable verbose mode')
    arg_parser = ArgumentParser(parents=[parser])
    # Get the arguments into an object
    arguments = arg_parser.parse_args()
    arguments.starttime = time()
    SetupLogging(debug=arguments.verbose)
    # Run the pipeline
    pipeline = AlleleFinder(arguments)
    pipeline.main()
    logging.info('Allele finding complete')
    return parser


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    cli()