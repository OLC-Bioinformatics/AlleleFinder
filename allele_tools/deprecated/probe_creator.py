#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import GenObject, make_path, MetadataObject, SetupLogging
import allele_finder
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from argparse import ArgumentParser
from threading import Thread
from queue import Queue
import multiprocessing
import logging
import shutil
import numpy
import os
__author__ = 'adamkoziol'


class Probes(object):

    def main(self):
        """

        """
        self.record_extraction()
        if self.blast:
            self.parse()
        self.create_allele_file()
        self.allelealigner()
        self.probefinder()
        self.probe_location()
        self.probes()

    def record_extraction(self):
        """
        Parse the input FASTA file, and create a dictionary of header: sequence for each entry
        """
        for record in SeqIO.parse(self.file, 'fasta'):
            metadata = MetadataObject()
            metadata.name = record.id
            metadata.records = record.seq.upper()
            # self.records[record.id] = record.seq.upper()
            self.samples.append(metadata)

    def remote_blast(self):
        """
        Create threads for each BLAST database: output file combination. Run multithreaded BLAST
        """
        logging.info('Running remote nr BLAST')
        for sample in self.samples:
            # Initialise the set of alleles with the input sequence
            try:
                sample.alleleset.add(str(sample.records))
            except AttributeError:
                sample.alleleset = set()
                sample.alleleset.add(str(sample.records))
            sample.blast_outputs = os.path.join(self.path, '{sn}_nr.xml'.format(sn=sample.name))
            if not os.path.isfile(sample.blast_outputs):
                # Create the BLAST request with the appropriate database
                result = NCBIWWW.qblast(program='blastn',
                                        database='nr',
                                        format_type='XML',
                                        expect=0.00001,
                                        word_size=28,
                                        filter='F',
                                        hitlist_size=1000000,
                                        sequence=str(sample.records))
                # Write the results to file
                with open(sample.blast_outputs, 'w') as results:
                    results.write(result.read())

    def parse(self):
        """
        Call the report parsing method for all the BLAST output files
        """
        logging.info('Parsing outputs')
        # Call parse_report for every file
        for sample in self.samples:
            if os.path.isfile(sample.blast_outputs):
                # Read in the BLAST results
                try:
                    with open(sample.blast_outputs, 'r') as result_handle:
                        logging.info('Parsing {sn} nr report'.format(sn=sample.name))
                        blast_record = NCBIXML.read(result_handle)
                        # Iterate through all the alignments
                        for alignment in blast_record.alignments:
                            # Iterate through each HSP per alignment
                            for hsp in alignment.hsps:
                                # Only retrieve sequences that are as long as the query sequence, and do not have gaps
                                if len(hsp.sbjct) == len(sample.records) and '-' not in hsp.sbjct:
                                    # Do not allow for more than five mismatches
                                    if hsp.identities >= len(sample.records) * (self.cutoff / 100):
                                        # Create a Seq object to add to the set
                                        sample.alleleset.add(hsp.sbjct)
                except FileNotFoundError:
                    pass

    def create_allele_file(self):
        """
        Create and populate the allele file with all alleles in the set
        """
        for sample in self.samples:
            logging.info('Creating {sn} allele FASTA file'.format(sn=sample.name))
            # Set the name of the allele file
            sample.allelefile = '{sn}_alleles.tfa'.format(sn=os.path.join(self.probepath, sample.name))
            try:
                sample.allelefiles.append(sample.allelefile)
            except AttributeError:
                sample.allelefiles = list()
                sample.allelefiles.append(sample.allelefile)
            with open(sample.allelefile, 'w') as alleles:
                fastalist = list()
                # Iterate through the set of alleles
                count = 0
                for alleleseq in sample.alleleset:
                    header = '{sn}_{count}'.format(sn=str(sample.name),
                                                   count=count)
                    # Create a SeqRecord for each allele
                    fasta = SeqRecord(Seq(alleleseq),
                                      # Without this, the header will be improperly formatted
                                      description='',
                                      # Use the gene name_iteration as the allele name
                                      id=header)
                    fastalist.append(fasta)
                    count += 1
                # Write the alleles to file
                SeqIO.write(fastalist, alleles, 'fasta')

    def allelealigner(self):
        """
        Perform a multiple sequence alignment of the allele sequences
        """
        logging.info('Aligning alleles')
        # Create the threads for the analysis
        for i in range(self.cpus):
            threads = Thread(target=self.alignthreads, args=())
            threads.setDaemon(True)
            threads.start()
        for sample in self.samples:
            sample.alignpath = os.path.join(self.path, 'alignedalleles')
            make_path(sample.alignpath)
            # Create a list to store objects
            sample.alignedalleles = list()
            for outputfile in sample.allelefiles:
                aligned = os.path.join(sample.alignpath, os.path.basename(outputfile))
                sample.alignedalleles.append(aligned)
                # Create the command line call
                clustalomega = ClustalOmegaCommandline(infile=outputfile,
                                                       outfile=aligned,
                                                       threads=4,
                                                       auto=True)
                sample.clustalomega = str(clustalomega)
                self.queue.put((sample, clustalomega, outputfile, aligned))
        self.queue.join()

    def alignthreads(self):
        while True:
            sample, clustalomega, outputfile, aligned = self.queue.get()
            if not os.path.isfile(aligned):
                # Perform the alignments
                # noinspection PyBroadException
                try:
                    clustalomega()
                # Files with a single sequence cannot be aligned. Copy the original file over to the aligned folder
                except:
                    shutil.copyfile(outputfile, aligned)
            self.queue.task_done()

    def probefinder(self):
        """
        Find the longest probe sequences
        """
        logging.info('Finding and filtering probe sequences')
        for sample in self.samples:
            # A list to store the metadata object for each alignment
            sample.gene = list()
            for align in sample.alignedalleles:
                # Create an object to store all the information for each alignment file
                metadata = GenObject()
                metadata.name = os.path.splitext(os.path.basename(align))[0]
                metadata.alignmentfile = align
                # Create an alignment object from the alignment file
                try:
                    metadata.alignment = AlignIO.read(align, 'fasta')
                except ValueError:
                    # If a ValueError: Sequences must all be the same length is raised, pad the shorter sequences
                    # to be the length of the longest sequence
                    # https://stackoverflow.com/q/32833230
                    records = SeqIO.parse(align, 'fasta')
                    # Make a copy, otherwise our generator is exhausted after calculating maxlen
                    records = list(records)
                    # Calculate the length of the longest sequence
                    maxlen = max(len(record.seq) for record in records)
                    # Pad sequences so that they all have the same length
                    for record in records:
                        if len(record.seq) != maxlen:
                            sequence = str(record.seq).ljust(maxlen, '.')
                            record.seq = Seq(sequence)
                    assert all(len(record.seq) == maxlen for record in records)
                    # Write to file and do alignment
                    metadata.alignmentfile = '{}_padded.tfa'.format(os.path.splitext(align)[0])
                    with open(metadata.alignmentfile, 'w') as padded:
                        SeqIO.write(records, padded, 'fasta')
                    # Align the padded sequences
                    metadata.alignment = AlignIO.read(metadata.alignmentfile, 'fasta')

                metadata.summaryalign = AlignInfo.SummaryInfo(metadata.alignment)
                # The dumb consensus is a very simple consensus sequence calculated from the alignment. Default
                # parameters of threshold=.7, and ambiguous='X' are used
                consensus = metadata.summaryalign.dumb_consensus()
                metadata.consensus = str(consensus)
                # The position-specific scoring matrix (PSSM) stores the frequency of each based observed at each
                # location along the entire consensus sequence
                metadata.pssm = metadata.summaryalign.pos_specific_score_matrix(consensus)
                metadata.identity = list()
                # Find the prevalence of each base for every location along the sequence
                for line in metadata.pssm:
                    try:
                        bases = [line['A'], line['C'], line['G'], line['T'], line['-']]
                        # Calculate the frequency of the most common base - don't count gaps
                        metadata.identity.append(float('{:.2f}'.format(max(bases[:4]) / sum(bases) * 100)))
                    except KeyError:
                        bases = [line['A'], line['C'], line['G'], line['T']]
                        # Calculate the frequency of the most common base - don't count gaps
                        metadata.identity.append(float('{:.2f}'.format(max(bases) / sum(bases) * 100)))
                # List to store metadata objects
                metadata.windows = list()
                # Variable to store whether a suitable probe has been found for the current organism + gene pair.
                # As the probe sizes are evaluated in descending size, as soon as a probe has been discovered, the
                # search for more probes can stop, and subsequent probes will be smaller than the one(s) already found
                passing = False
                # Create sliding windows of size self.max - self.min from the list of identities for each column
                # of the alignment
                for i in reversed(range(self.min, self.max + 1)):
                    if not passing:
                        windowdata = MetadataObject()
                        windowdata.size = i
                        windowdata.max = 0
                        windowdata.sliding = list()
                        # Create a counter to store the starting location of the window in the sequence
                        n = 0
                        # Create sliding windows from the range of sizes for the list of identities
                        windows = self.window(metadata.identity, i)
                        # Go through each window from the collection of sliding windows to determine which window(s)
                        # has (have) the best results
                        for window in windows:
                            # Create another object to store all the data for the window
                            slidingdata = MetadataObject()
                            # Only consider the window if every position has a percent identity greater than the cutoff
                            if min(window) > self.cutoff:
                                # Populate the object with the necessary variables
                                slidingdata.location = '{}:{}'.format(n, n + i)
                                slidingdata.min = min(window)
                                slidingdata.mean = float('{:.2f}'.format(numpy.mean(window)))
                                slidingdata.sequence = str(consensus[n:n+i])
                                # Create attributes for evaluating windows. A greater/less windowdata.max/windowdata.min
                                #  means a better/less overall percent identity, respectively
                                windowdata.max = slidingdata.mean if slidingdata.mean >= windowdata.max \
                                    else windowdata.max
                                windowdata.min = slidingdata.mean if slidingdata.mean <= windowdata.max \
                                    else windowdata.min
                                # Add the object to the list of objects
                                windowdata.sliding.append(slidingdata)
                                passing = True
                            n += 1
                        # All the object to the list of objects
                        metadata.windows.append(windowdata)
                # All the object to the list of objects
                sample.gene.append(metadata)

    def probe_location(self):
        """
        Find the 'best' probes for each gene by evaluating the percent identity of the probe to the best recorded
        percent identity for that organism + gene pair. Extract the location of the probe, so that all the alleles
        in that location can be recovered
        """
        logging.info('Determining optimal probe sequences')
        for sample in self.samples:
            # Make a folder to store the probes
            for gene in sample.gene:
                for window in gene.windows:
                    # Variable to record whether a probe has already been identified from this gene
                    passed = False
                    for sliding in window.sliding:
                        # Only consider the sequence if the sliding object has data, if the probe in question
                        # has a mean identity equal to the highest observed identity for that probe size, and
                        # if the mean identity is greater or equal than the lowest observed identity
                        if sliding.datastore and sliding.mean == window.max and sliding.mean >= window.min \
                                and not passed:
                            sample.location = sliding.location
                            passed = True

    def probes(self):
        """
        Extract the probe sequence from all the alleles
        """
        for sample in self.samples:
            # Determine the starting and ending locations for the slice of the probe sequence
            start, stop = sample.location.split(':')
            sample.probeallele = set()
            # Open the allele file
            sample.probeoutputfile = os.path.join(self.probepath, '{}_probe_alleles.tfa'.format(sample.name))
            for allele in sample.alleleset:
                sample.probeallele.add(allele[int(start):int(stop)])
            alleles = list()
            with open(sample.probeoutputfile, 'w') as allelefile:
                for count, probe_allele in enumerate(sample.probeallele):
                    # Use 1-based indexing if requested
                    if self.one_based:
                        count += 1
                    header = '{sn}_{count}'.format(sn=str(sample.name),
                                                   count=count)
                    # Create a SeqRecord for each allele
                    fasta = SeqRecord(Seq(probe_allele),
                                      # Without this, the header will be improperly formatted
                                      description='',
                                      # Use the gene name_iteration as the allele name
                                      id=header)
                    alleles.append(fasta)
                    # Write the alleles to file
                SeqIO.write(alleles, allelefile, 'fasta')

    @staticmethod
    def window(iterable, size):
        """
        https://coderwall.com/p/zvuvmg/sliding-window-in-python
        :param iterable: string from which sliding windows are to be created
        :param size: size of sliding window to create
        """
        i = iter(iterable)
        win = []
        for e in range(0, size):
            win.append(next(i))
        yield win
        for e in i:
            win = win[1:] + [e]
            yield win

    def __init__(self, path, targetfile, min_length, max_length, cutoff, perc_gc, blast, one_based):
        # Determine the path in which the sequence files are located. Allow for ~ expansion
        if path.startswith('~'):
            self.path = os.path.abspath(os.path.expanduser(os.path.join(path)))
        else:
            self.path = os.path.abspath(os.path.join(path))
        self.file = os.path.join(self.path, targetfile)
        assert os.path.isfile(self.file), 'Cannot find the supplied FASTA file: {fn}'.format(fn=self.file)
        self.reportpath = os.path.join(self.path, 'reports')
        self.probepath = os.path.join(self.path, 'probes')
        make_path(self.reportpath)
        self.min = min_length
        self.max = max_length
        self.cutoff = cutoff
        self.perc_gc = perc_gc
        self.blast = blast
        self.one_based = one_based
        self.cpus = multiprocessing.cpu_count()
        self.queue = Queue()
        self.samples = list()


def cli():
    # Import the argument parser from allele_finder.py
    parent_parser = allele_finder.cli()
    parser = ArgumentParser(parents=[parent_parser])
    parser.add_argument('-min', '--min',
                        default=20,
                        type=int,
                        help='Minimum size of probe to create')
    parser.add_argument('-max', '--max',
                        default=50,
                        type=int,
                        help='Maximum size of probe to create')
    parser.add_argument('-c', '--cutoff',
                        default=70,
                        help='Cutoff percent identity of a nucleotide location to use')
    parser.add_argument('-gc', '--percentgc',
                        default=50,
                        type=int,
                        help='Desired percent GC of the probe')
    parser.add_argument('-r', '--runblast',
                        action='store_true',
                        help='Run BLAST analyses on the supplied target file. If not enabled, then the program assumes '
                             'that the supplied file includes all the desired alleles to use to create the probe')
    parser.add_argument('-aa', '--amino_acid',
                        choices=['targets_nt', 'targets_aa'],
                        help='Find the amino acid sequence of alleles. The target alleles supplied can either be '
                             'nucleotide or amino acid. Default is nucleotide')
    parser.add_argument('-o', '--one_based',
                        action='store_true',
                        help='Use 1-based indexing rather than the default 0-based')
    # Get the arguments into an object
    arguments = parser.parse_args()
    SetupLogging(debug=arguments.verbose)
    if arguments.runblast:
        # Run the allele-finding pipeline
        finder = allele_finder.AlleleFinder(path=arguments.path,
                                            targetfile=arguments.targetfile,
                                            analysis_type=arguments.blast,
                                            fasta_path=arguments.fasta_path,
                                            genesippr=arguments.genesippr,
                                            metadata_file=arguments.metadatafile,
                                            cutoff=arguments.cutoff,
                                            amino_acid=arguments.amino_acid,
                                            one_based=arguments.one_based)
        finder.main()
    # Run the pipeline
    probes = Probes(path=arguments.path,
                    targetfile=arguments.targetfile,
                    min_length=arguments.min,
                    max_length=arguments.max,
                    cutoff=arguments.cutoff,
                    perc_gc=arguments.percentgc,
                    blast=arguments.runblast,
                    one_based=arguments.one_based)
    probes.main()
    logging.info('Probe finding complete')


if __name__ == '__main__':
    cli()
