#!/usr/bin/env python
from olctools.accessoryFunctions.accessoryFunctions import make_path, SetupLogging
from profile_reduce import ProfileReduce
from allele_profiler import read_profile
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
        if self.profile_file:
            self.profile_data = read_profile(profile_file=self.profile_file)
            self.aa_profile()
            reduce = ProfileReduce(profile=self.aa_profile_file,
                                   names=self.gene_name_file)
            reduce.main()
            self.aa_profile_data = read_profile(profile_file=self.aa_profile_file)
            self.profile_link()
            self.link_file()

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
            # Initialise the dictionary to store the links between the nt alleles and the aa alleles with the name of
            # the gene
            self.allele_links[gene_name] = dict()
            logging.info(f'Processing {gene_name}')
            # Initialise a dictionary to store the translated allele sequence: allele name
            seq_allele = dict()
            # Open the file to store the translated alleles
            with open(os.path.join(self.translated_path, '{gn}.fasta'.format(gn=gene_name)), 'w') as aa_alleles:
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
                            # Add the string of the amino acid sequence to the dictionary
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
                            # Populate the linking dictionary with the nt allele: aa allele
                            self.allele_links[gene_name][allele.split('_')[-1]] = allele.split('_')[-1]
                        # Amino acid allele already exists
                        else:
                            # Extract the allele name corresponding to the translated sequence
                            aa_allele = seq_allele[str(translated_allele)]
                            # Update the notes, including that this allele is a duplicate, and a pointer to the original
                            notes.write('{nt}\t{aa}\tduplicate\n'.format(nt=allele,
                                                                         aa=aa_allele))
                            # Populate the linking dictionary with the nt allele: aa allele
                            self.allele_links[gene_name][allele.split('_')[-1]] = aa_allele.split('_')[-1]
                            if self.one_based:
                                self.allele_links[gene_name]['0'] = '0'
                            else:
                                self.allele_links[gene_name]['0'] = '0'

    def aa_profile(self):
        """

        """
        profile_str = str()
        for seq_type in sorted(int(st) for st in self.profile_data.keys()):
            profile_str += str(seq_type) + '\t'
            for gene_name, allele in self.profile_data[str(seq_type)].items():
                self.gene_names.add(gene_name)
                profile_str += self.allele_links[gene_name][str(allele)] + '\t'
            profile_str = profile_str.rstrip()
            profile_str += '\n'
        with open(self.aa_profile_file, 'w') as aa_profile:
            names = '\t'.join(sorted(list(self.gene_names)))
            aa_profile.write('ST\t{names}\n'.format(names=names.rstrip()))
            aa_profile.write(profile_str)

        with open(self.gene_name_file, 'w') as gene_file:
            gene_file.write('\n'.join(sorted(list(self.gene_names))))

    def profile_link(self):
        """

        """
        match_score = dict()
        for st in self.profile_data:
            self.profile_matches[st] = set()
            match_score[st] = dict()
            for gene, allele in self.profile_data[st].items():

                for aa_st in self.aa_profile_data:
                    # print('\t' + aa_st, self.aa_profile_data[aa_st][gene])
                    # if allele == self.allele_links[gene][self.aa_profile_data[aa_st][gene]]:
                    if self.aa_profile_data[aa_st][gene] == self.allele_links[gene][allele]:
                        if aa_st not in match_score[st]:
                            match_score[st][aa_st] = 0
                        match_score[st][aa_st] += 1
        for st, aa_st_dict in match_score.items():
            for aa_st, matches, in aa_st_dict.items():
                if matches == len(self.gene_names):
                    self.profile_matches[st].add(aa_st)

    def link_file(self):
        """

        :return:
        """
        with open(self.aa_nt_profile_link_file, 'w') as link:
            for st, match_set in self.profile_matches.items():
                if len(match_set) == 1:
                    link.write('{nt_st}\t{aa_st}\n'.format(nt_st=st,
                                                           aa_st=''.join(match_set)))
                else:
                    link.write('{nt_st}\t{aa_st}\n'.format(
                        nt_st=st,
                        aa_st=';'.join(str(aa_st) for aa_st in sorted(int(aa_st) for aa_st in match_set))))

    def __init__(self, path, profile, one_based):
        logging.info('Welcome to the allele translator!')
        if path.startswith('~'):
            self.path = os.path.abspath(os.path.expanduser(os.path.join(path)))
        else:
            self.path = os.path.abspath(os.path.join(path))
        if profile:
            self.profile_file = os.path.join(self.path, 'profile', 'profile.txt')
            assert os.path.isfile(self.profile_file), 'Cannot locate the required profile file: {profile}. Please ' \
                                                      'ensure that the file name and path of your file is correct'\
                .format(profile=self.profile_file)
        else:
            self.profile_file = None
        self.one_based = one_based
        self.sequence_files = glob(os.path.join(self.path, '*.fasta'))
        self.translated_path = os.path.join(self.path, 'aa_alleles')
        self.notes_path = os.path.join(self.path, 'notes')
        make_path(inpath=self.translated_path)
        make_path(inpath=self.notes_path)
        self.allele_dict = dict()
        self.profile_data = dict()
        self.allele_links = dict()
        self.aa_profile_path = os.path.join(self.path, 'aa_profile')
        make_path(self.aa_profile_path)
        self.aa_profile_file = os.path.join(self.aa_profile_path, 'aa_profile.txt')
        self.gene_names = set()
        self.gene_name_file = os.path.join(self.aa_profile_path, 'gene_names.txt')
        self.aa_profile_data = dict()
        self.profile_matches = dict()
        self.aa_nt_profile_link_file = os.path.join(self.aa_profile_path, 'reports', 'aa_nt_profile_links.tsv')


def cli():
    # Parser for arguments
    parser = ArgumentParser(description='Translate allele files in nucleotide format to amino acid. '
                                        'Remove duplicates. Keep notes.')
    parser.add_argument('-p', '--path',
                        required=True,
                        help='Specify path containing allele files.')
    parser.add_argument('--profile',
                        action='store_true',
                        help='Optionally parse the nucleic acid profile, and create the corresponding reduced amino '
                             'acid profile')
    parser.add_argument('-o', '--one_based',
                        action='store_true',
                        help='Use 1-based indexing rather than the default 0-based')
    # Get the arguments into an object
    arguments = parser.parse_args()
    SetupLogging(debug=True)
    translate = Translate(path=arguments.path,
                          profile=arguments.profile,
                          one_based=arguments.one_based)
    translate.main()
    logging.info('Allele translation complete!')


if __name__ == '__main__':
    cli()
