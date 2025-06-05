#!/usr/bin/env python3

"""
Find duplicate alleles in the provided allele files.
Report which alleles were removed.
"""

# Standard imports
from argparse import ArgumentParser
from collections import defaultdict
import os

# Third-party imports
from Bio import SeqIO


class FindDuplicateAlleles:
    """
    Class to find duplicate alleles in the provided allele files.
    """

    def __init__(self, args):
        self.allele_file = args.allele_file
        self.allele_path = os.path.dirname(self.allele_file)

    def main(self):
        records, full_records = self._extract_records(
            allele_file=self.allele_file
        )
        duplicates = self._find_duplicates(seq_to_ids=records)

        print("Duplicate alleles found:")
        for group in duplicates:
            print(f" - {', '.join(group)}")

        self._write_duplicates(
            duplicates=duplicates,
            output_file=os.path.join(self.allele_path, "duplicate_alleles.txt")
        )
        self._write_unique_alleles(
            seq_to_ids=records,
            id_to_record=full_records,
            output_file=os.path.join(self.allele_path, "unique_alleles.fasta"),
        )

    @staticmethod
    def _extract_records(
        *,  # Enforce keyword arguments
        allele_file: str
    ) -> tuple[defaultdict[str, list[str]], dict[str, SeqIO.SeqRecord]]:
        """
        Extract records and group IDs by sequence.

        :param allele_file: Path to the allele FASTA file.
        :return: A tuple of (seq_to_ids, id_to_record).
        """
        seq_to_ids = defaultdict(list)
        id_to_record = {}
        for record in SeqIO.parse(allele_file, "fasta"):
            seq_to_ids[str(record.seq)].append(record.id)
            id_to_record[record.id] = record
        return seq_to_ids, id_to_record

    @staticmethod
    def _find_duplicates(
        *,  # Enforce keyword arguments
        seq_to_ids: defaultdict[str, list[str]]
    ) -> list[list[str]]:
        """
        Return a list of lists of IDs for duplicate sequences.
        """
        return [ids for ids in seq_to_ids.values() if len(ids) > 1]

    @staticmethod
    def _write_duplicates(
        *,  # Enforce keyword arguments
        duplicates: list[list[str]], output_file: str
    ):
        """
        Write duplicate allele groups to a file.
        """
        with open(output_file, "w", encoding="utf-8") as f:
            f.write("Duplicate_Alleles\n")
            for group in duplicates:
                f.write(";".join(group) + "\n")

    @staticmethod
    def _write_unique_alleles(
        *,  # Enforce keyword arguments
        seq_to_ids: defaultdict[str, list[str]],
        id_to_record: dict[str, SeqIO.SeqRecord],
        output_file: str
    ):
        """
        Write unique alleles (first occurrence of each sequence) to a FASTA
        file.
        """
        unique_records = [id_to_record[ids[0]] for ids in seq_to_ids.values()]
        with open(output_file, "w", encoding="utf-8") as f:
            SeqIO.write(unique_records, f, "fasta")


def main():
    parser = ArgumentParser(
        description="Find duplicate alleles in FASTA files."
    )
    parser.add_argument(
        "-a", "--allele_file",
        metavar="allele_file",
        required=True,
        help="Path to the allele FASTA file."
    )
    args = parser.parse_args()

    finder = FindDuplicateAlleles(args)
    finder.main()


if __name__ == "__main__":
    main()
