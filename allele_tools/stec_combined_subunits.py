#!/usr/bin/env python3

"""
Processes STEC alleles with both A and B subunits.
"""

# Standard imports
from argparse import ArgumentParser
from collections import defaultdict
from csv import DictReader
from glob import glob
import hashlib
import logging
import multiprocessing
import os
import re
import shutil
import subprocess
import time
from typing import (
    Dict,
    Optional,
    TypedDict
)

# Third-party imports
from Bio.Application import ApplicationError
from Bio.Blast.Applications import (
    NcbiblastnCommandline,
    NcbiblastxCommandline,
    NcbitblastxCommandline
)
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Local imports
from allele_tools.methods import (
    error_print,
    pathfinder,
    setup_logging
)
from allele_tools.version import __version__

__authors__ = "adamkoziol"


class AlleleHitDict(TypedDict, total=False):
    """
    Dictionary to store the best BLAST hit row for each subunit (A and B).
    Keys are "A" and "B", values are the row dictionaries or None.
    """

    A: Optional[dict]
    B: Optional[dict]


class AlleleStatsDict(TypedDict):
    """
    Dictionary to store cumulative statistics for a given allele's subunits.
    - identical: Total number of identical matches across both subunits.
    - aligned: Total alignment length across both subunits.
    - mismatches: Total mismatches across both subunits.
    - gaps: Total gaps across both subunits.
    - query_seqs: Dictionary with keys "A" and "B" for the query sequences.
    """

    identical: int
    aligned: int
    mismatches: int
    gaps: int
    query_seqs: Dict[str, str]  # keys: "A" and "B"


class CombinedSubunits:
    """
    Class to handle combined subunits for STEC alleles.
    """

    # Use class variables for the registry and counter
    _novel_registry = {}
    _novel_counter = 1

    def __init__(self, args):
        self.args = args
        self.allele_path = pathfinder(path=args.allele_path)
        self.report_path = pathfinder(path=args.report_path)
        self.query_path = pathfinder(path=args.query_path)

        # Set BLAST mode and molecule type
        if args.blast_mode == 'blastn':
            self.blast_mode = args.blast_mode
        elif args.blast_mode == 'tblastx':
            self.blast_mode = args.blast_mode
        elif args.blast_mode == 'blastx':
            self.blast_mode = args.blast_mode
        else:
            self.blast_mode = ['blastn', 'tblastx']

        if self.blast_mode == 'blastn':
            self.molecule = 'nt'
        elif self.blast_mode in ['blastx', 'tblastx']:
            self.molecule = 'aa'
        else:
            self.molecule = ['nt', 'aa']

        # Store the preliminary boolean
        self.preliminary = args.preliminary

        # Set the number of alignments
        self.num_alignments = args.num_alignments

        # Set the percent identity cutoff
        self.cutoff = args.cutoff

        # Create a list to store error messages
        errors = []

        # Nucleotide allele checks
        if not os.path.isdir(self.allele_path):
            errors.append(
                "Could not find supplied nucleotide allele folder: "
                f"{self.allele_path}"
            )
        else:
            # Initialize query files
            self.query_files, errors = self._locate_fasta_files(
                directory=self.query_path,
                errors=errors
            )

        # Query folder checks
        if not os.path.isdir(self.query_path):
            errors.append(
                f"Could not find supplied query folder: {self.query_path}"
            )
        else:
            allele_files, errors = self._locate_fasta_files(
                directory=self.allele_path,
                errors=errors
            )
            # Only take the first allele file found
            self.allele_file = allele_files[0]
            logging.debug("Found allele file: %s", self.allele_file)

        # Try to create the report folder (if it doesn't already exist)
        try:
            os.makedirs(self.report_path, exist_ok=True)
        except (OSError, RuntimeError):
            error_str = f"Error creating report folder: {self.report_path}"
            errors.append(error_str)

        # If there are any errors, print them
        if errors:
            error_print(errors=errors)

        # Set number of CPUs for multiprocessing
        self.cpus = multiprocessing.cpu_count() - 1

        # Set field names for BLAST output
        self.fieldnames = [
            "query_id",
            "subject_id",
            "identical",
            "mismatches",
            "gaps",
            "evalue",
            "bit_score",
            "query_length",
            "subject_length",
            "alignment_length",
            "query_start",
            "query_end",
            "subject_start",
            "subject_end",
            "query_sequence",
            "subject_sequence",
        ]
        self.extended_fieldnames = self.fieldnames.copy()
        self.extended_fieldnames.insert(14, "percent_match")

        # Set format string for BLAST output
        self.outfmt = (
            "6 qseqid sseqid nident mismatch gaps evalue bitscore qlen slen "
            "length qstart qend sstart send qseq sseq"
        )

        # Initialise variables to store the split aa database information
        self.split_aa_dbs = {}
        self.split_aa_dbs_paths = []

        # Find split AA databases (if specified)
        if getattr(args, "split_aa_db_dir", None):
            self.split_aa_dbs = self._find_split_aa_dbs(
                split_aa_db_dir=args.split_aa_db_dir
            )

            # Create the necessary blast databases
            for _, allele_file in self.split_aa_dbs.items():

                # Get the directory of the allele file
                allele_path = os.path.dirname(allele_file)

                # Make the BLAST database for the allele file
                blast_db_path = self.make_blast_db(
                    allele_file=allele_file,
                    allele_path=allele_path,
                    molecule="aa",
                )
                # Store the path to the split AA database
                self.split_aa_dbs_paths.append(blast_db_path)
        else:
            self.split_aa_dbs = {}

        # Initialise the BLAST database
        if self.blast_mode == 'blastx':
            self.blast_db_path = self.make_blast_db(
                allele_file=self.allele_file,
                allele_path=self.allele_path,
                blastx=True,
                molecule='aa'
            )
        else:
            self.blast_db_path = self.make_blast_db(
                allele_file=self.allele_file,
                allele_path=self.allele_path,
                molecule='nt'
            )

        # Initialize metadata dictionary
        self.metadata = {}

    @staticmethod
    def _locate_fasta_files(
        *,  # Enforce keyword arguments
        directory: str,
        errors: list[str],
    ) -> tuple[list[str], list[str]]:
        """
        Locate FASTA files in a directory.

        :param directory: Path to the directory to search
        :param errors: List to store error messages
        :return: A tuple containing a list of found FASTA file paths and a
        list of error messages
        """
        fasta_extensions = (
            "*.fa", "*.fasta", "*.fna", "*.ffn", "*.faa", "*.tfa"
        )

        # Initialize a list to store the FASTA files
        fasta_files = []

        # Iterate over the FASTA extensions and glob for files
        for ext in fasta_extensions:
            fasta_files.extend(glob(os.path.join(directory, ext)))

        # Check if any FASTA files were found
        if not fasta_files:
            error_string = (
                f"No FASTA files found in {directory} with extensions: "
                f"{', '.join(fasta_extensions)}"
            )
            errors.append(error_string)

        return sorted(fasta_files), errors

    @staticmethod
    def _find_split_aa_dbs(
        *,  # Enforce keyword arguments
        split_aa_db_dir: str
    ) -> dict:
        """
        Locate the split AA allele databases for Stx1A, Stx1B, Stx2A, Stx2B.

        :param split_aa_db_dir: Path to the directory containing the split AA
        databases
        :return: A dictionary mapping subunit names to their FASTA file paths
        """
        # Initialize a dictionary to store the database paths
        dbs = {}

        # Use pathfinder to return the absolute path to the directory
        dir_path = pathfinder(path=split_aa_db_dir)

        # Look for Stx1A, Stx1B, Stx2A, and Stx2B databases
        for stx in ["Stx1", "Stx2"]:
            for subunit in ["A", "B"]:
                pattern = os.path.join(dir_path, f"{stx}{subunit}_aa_*.fasta")
                matches = glob(pattern)
                if not matches:
                    raise FileNotFoundError(
                        f"No database found for {stx}{subunit} in {dir_path}"
                    )
                dbs[f"{stx.lower()}{subunit}"] = matches[0]
        return dbs

    @staticmethod
    def make_blast_db(
        *,  # Enforce keyword arguments
        allele_file: str,
        allele_path: str,
        blastx: bool = False,
        molecule: str = 'nt'
    ) -> str:
        """
        Process the allele file with makeblastdb (if necessary)

        :param allele_file: Path to the allele file
        :param allele_path: Path to the allele directory
        :param blastx: Whether to use blastx (default is False)
        :param molecule: Molecule type (default is 'nt')
        :return: Path to the created BLAST database
        """
        if blastx:
            allele_basename = os.path.join(
                allele_path, 'blastx_split_allele_database'
            )
            CombinedSubunits.split_fasta(
                input_fasta=allele_file,
                output_fasta=allele_basename + ".txt",
                linker_char="X",
                linker_min=4
            )
            allele_file = allele_basename + ".txt"
        else:
            # Split off the extension of the allele file
            allele_basename = os.path.join(
                allele_path,
                os.path.splitext(os.path.basename(allele_file))[0]
            )

        # Use makeblastdb to create a BLAST database (if necessary)
        if molecule == 'nt':
            if not os.path.isfile(allele_basename + '.ndb'):
                try:
                    subprocess.run(
                        [
                            "makeblastdb",
                            "-in", allele_file,
                            "-dbtype", "nucl",
                            "-out", allele_basename,
                            "-title", "STEC Allele Database"
                        ],
                        check=True
                    )
                except subprocess.CalledProcessError as exc:
                    logging.error(
                        "Error creating BLAST database: %s", exc
                    )
                    raise SystemExit from exc
        elif molecule == 'aa':
            if not os.path.isfile(allele_basename + '.pdb'):
                try:
                    subprocess.run(
                        [
                            "makeblastdb",
                            "-in", allele_file,
                            "-dbtype", "prot",
                            "-out", allele_basename,
                            "-title", "STEC Allele Database"
                        ],
                        check=True
                    )
                except subprocess.CalledProcessError as exc:
                    logging.error(
                        "Error creating BLAST database: %s", exc
                    )
                    raise SystemExit from exc

        return allele_basename

    @staticmethod
    def split_fasta(
        *,  # Enforce keyword arguments
        input_fasta: str,
        output_fasta: str,
        linker_char: str = "X",
        linker_min: int = 4
    ):
        """
        Splits concatenated subunit FASTA records into separate subunits.
        Each subunit is written as a new FASTA entry, omitting the linker.
        """
        with open(input_fasta, "r") as infile, open(
            output_fasta, "w"
        ) as outfile:
            for record in SeqIO.parse(infile, "fasta"):
                seq = str(record.seq)
                # Find the linker region
                match = re.search(
                    f"{linker_char}{{{linker_min},}}", seq.upper()
                )
                if not match:
                    # No linker found, write as a single subunit
                    SeqIO.write(record, outfile, "fasta")
                    continue
                start, end = match.span()

                # Write subunit A
                record_a = record[:start]
                record_a.id = f"{record.id}_A"
                record_a.description = ""
                record_a.seq = record.seq[:start]
                SeqIO.write(record_a, outfile, "fasta")

                # Write subunit B
                record_b = record[end:]
                record_b.id = f"{record.id}_B"
                record_b.description = ""
                record_b.seq = record.seq[end:]
                SeqIO.write(record_b, outfile, "fasta")

    def process(self):
        """
        Process the combined subunits.
        """
        # Implement processing logic here
        logging.info(
            "Processing combined subunits with the following parameters:"
            "Allele Path: %s, Report Path: %s, Query Path: %s",
            self.allele_path,
            self.report_path,
            self.query_path
        )

        # Perform pre-BLAST preparation of query files
        self.metadata = self._prep_query_files(
            query_files=self.query_files,
            query_path=self.query_path
        )
        for query_name, query_info in self.metadata.items():
            logging.debug(
                "Prepared query file '%s' at path: %s",
                query_name, query_info['path']
            )

        # Process the dual analysis types (e.g. a list of both molecules)
        # Check to see if the molecule is a list or a string
        if isinstance(self.molecule, list):
            for iterator, molecule in enumerate(self.molecule):
                blast_mode = self.blast_mode[iterator]
                self.run(
                    blast_mode=blast_mode,
                    molecule=molecule
                )

            # Print the report
            self._combined_report(
                metadata=self.metadata,
                report_path=self.report_path
            )
        else:
            self.run(
                blast_mode=str(self.blast_mode),
                molecule=self.molecule
                )

    def run(
        self,
        *,  # Enforce keyword arguments
        blast_mode: str,
        molecule: str
    ):
        """
        Run the combined subunits analysis pipeline.
        """
        # Run BLAST
        self._run_blast(
            allele_file=self.blast_db_path,
            blast_mode=blast_mode,
            cpus=self.cpus,
            metadata=self.metadata,
            molecule=molecule,
            outfmt=self.outfmt,
        )

        # Set the blastx variable
        if blast_mode == 'blastx':
            blastx = True
        else:
            blastx = False

        # Format the BLAST output
        self._parseable_blast_outputs(
            blastx=blastx,
            cutoff=self.cutoff,
            fieldnames=self.fieldnames,
            extended_fieldnames=self.extended_fieldnames,
            metadata=self.metadata,
            molecule=molecule,
        )

        # Find the best BLAST hit
        self.metadata = self._find_best_blast_hits(
            blastx=blastx,
            cutoff=self.cutoff,
            metadata=self.metadata,
            molecule=molecule
        )
        for name, info in self.metadata.items():
            novel = info.get(f"{molecule}_novel_hits", [])
            logging.debug(
                "Best BLAST hits for %s: %s. Novel hits: %s",
                name, info[f"{molecule}_best_hits"],
                novel
            )

        self._write_novel_alleles_to_db(
            db_fasta_file=self.allele_file,
            report_path=self.report_path
        )

        # For nucleotide analyses, proceed to find the amino acid alleles for
        # each novel hit
        if molecule == 'nt' and not isinstance(self.molecule, list):

            if not self.preliminary:

                self.metadata = self._export_novel_alleles(
                    blastx=blastx,
                    cutoff=self.cutoff,
                    metadata=self.metadata,
                    molecule="nt",
                    report_path=self.report_path
                )
                self.metadata = self._find_aa_hits_novel_alleles(
                    cpus=self.cpus,
                    metadata=self.metadata,
                    num_alignments=self.num_alignments,
                    outfmt=self.outfmt,
                    split_aa_dbs=self.split_aa_dbs
                )
                self._add_headers_novel_aa_reports(
                    blastx=False,
                    cutoff=self.cutoff,
                    extended_fieldnames=self.extended_fieldnames,
                    fieldnames=self.fieldnames,
                    metadata=self.metadata,
                )

                self.metadata = self._extract_novel_aa_best_hits(
                    metadata=self.metadata,
                )
                self.metadata = self._assign_aa_allele_profiles_and_operons(
                    metadata=self.metadata,
                )
                self._export_aa_novel_operons(
                    metadata=self.metadata,
                    report_path=self.report_path
                )
                self._export_aa_novel_alleles(
                    cutoff=self.cutoff,
                    metadata=self.metadata,
                    molecule='aa',
                    report_path=self.report_path
                )

            # Print the report
            self._nt_report(
                metadata=self.metadata,
                report_path=self.report_path
            )
        elif molecule == 'aa' and not isinstance(self.molecule, list):
            # Export novel alleles
            self._export_novel_alleles(
                blastx=blastx,
                cutoff=self.cutoff,
                metadata=self.metadata,
                molecule=molecule,
                report_path=self.report_path
            )

            # Handle other molecule types
            self._aa_report(
                metadata=self.metadata,
                report_path=self.report_path
            )

    @staticmethod
    def _prep_query_files(
        *,  # Enforce keyword arguments
        query_files: list[str],
        query_path: str
    ):
        """
        Prepare query files for BLAST. Create a subfolder for each query file
        and symlink the query file into that subfolder.

        :param query_files: List of query file paths
        :param query_path: Path to the query directory
        :return: A dictionary containing metadata about the query files
        """
        # Create a dictionary to store file-specific information
        metadata = {}

        for query_file in query_files:
            logging.debug("Preparing query file: %s", query_file)

            # Get the base name of the query file
            query_basename = os.path.splitext(os.path.basename(query_file))[0]

            # Split off any |length
            query_basename = query_basename.split("|")[0].replace(".", "_")

            # Create a subfolder for the query file
            query_subfolder = os.path.join(query_path, query_basename)
            os.makedirs(query_subfolder, exist_ok=True)

            # Create a symlink to the query file in the subfolder
            query_symlink = os.path.join(
                query_subfolder, os.path.basename(query_file)
            )

            # Store metadata about the query file
            metadata[query_basename] = {
                'path': query_subfolder,
                'file': query_symlink
            }

            # Only create the symlink if it doesn't already exist
            if not os.path.exists(query_symlink):
                os.symlink(query_file, query_symlink)

        return metadata

    @staticmethod
    def _run_blast(
        *,  # Enforce keyword arguments
        allele_file: str,
        blast_mode: str,
        cpus: int,
        metadata: dict[str, dict[str, str]],
        molecule: str,
        outfmt: str,
    ) -> dict[str, dict[str, str]]:
        """
        Run BLAST for the given query files against the specified database.

        :param allele_file: Path to the allele file
        :param blast_mode: BLAST mode to use (e.g., "blastn" or "tblastx")
        :param cpus: Number of CPUs to use for BLAST
        :param metadata: Metadata about the query files
        :param molecule: Molecule type (e.g., "nt" or "aa")
        :param outfmt: Output format for BLAST
        :return: Updated metadata with BLAST results
        """
        for file, info in metadata.items():
            blast_output = os.path.join(
                info['path'], f"{file}_{molecule}_blast.tsv")

            blast_cmd = CombinedSubunits._run_blast_cmd(
                allele_file=allele_file,
                blast_mode=blast_mode,
                blast_output=blast_output,
                cpus=cpus,
                outfmt=outfmt,
                query=info['file']
            )

            # Update the metadata with the BLAST output path and command
            info[f'{molecule}_blast_output'] = blast_output
            info[f'{molecule}_blast_cmd'] = blast_cmd

        return metadata

    @staticmethod
    def _run_blast_cmd(
        *,  # Enforce keyword arguments
        allele_file: str,
        blast_mode: str,
        blast_output: str,
        cpus: int,
        outfmt: str,
        query: str,
        num_alignments: int = 10
    ) -> str:
        """
        Run a BLAST command and handle errors.

        :param allele_file: Path to the allele file
        :param blast_mode: BLAST mode to use (e.g., "blastn" or "tblastx")
        :param blast_output: Path to the BLAST output file
        :param cpus: Number of CPUs to use for BLAST
        :param outfmt: Output format for BLAST
        :param query: Query sequence or file
        :param num_alignments: Number of alignments to return
        """
        if blast_mode == 'blastn':
            blast_cmd = NcbiblastnCommandline(
                db=allele_file,
                evalue=1e-5,
                num_alignments=num_alignments,
                num_threads=cpus,
                outfmt=outfmt,
                out=blast_output,
                query=query,
                task="blastn"
            )
        elif blast_mode == 'tblastx':
            blast_cmd = NcbitblastxCommandline(
                db=allele_file,
                evalue=1e-5,
                num_alignments=num_alignments,
                num_threads=cpus,
                outfmt=outfmt,
                out=blast_output,
                query=query,
            )
        elif blast_mode == 'blastx':
            blast_cmd = NcbiblastxCommandline(
                db=allele_file,
                evalue=1e-5,
                num_alignments=num_alignments,
                num_threads=cpus,
                outfmt=outfmt,
                out=blast_output,
                query=query,
                seg='no'
            )

        # Skip BLAST if output already exists
        if os.path.isfile(blast_output):
            return str(blast_cmd)

        logging.info(
            "Running BLAST for query file %s with database %s. BLAST "
            "outputs saved to %s", query, allele_file, blast_output
        )
        try:
            blast_cmd()
        except ApplicationError as exc:
            logging.error("BLAST application error for %s: %s", query, exc)

        return str(blast_cmd)

    @staticmethod
    def _parseable_blast_outputs(
        *,  # Enforce keyword arguments
        blastx: bool,
        cutoff: float,
        fieldnames: list[str],
        extended_fieldnames: list[str],
        metadata: dict[str, dict[str, str]],
        molecule: str
    ):
        """
        Add a header to the BLAST report, so that it is easier to figure out
        what is in each column
        :param blastx: Boolean indicating if this is a BLASTx report
        :param fieldnames: String of all the field names in the BLAST report
        :param extended_fieldnames: String of the BLAST field names plus the
        calculated percent identity
        :param molecule: Molecule type (e.g., "nt" or "aa")
        :param cutoff: Float of the minimum percent identity between query
        and subject sequence.
        """

        logging.info("Adding headers to BLAST outputs")
        for _, info in metadata.items():

            # Set a variable for the blast output file
            blast_output = info[f"{molecule}_blast_output"]

            # Add the headers
            CombinedSubunits._add_headers(
                blast_output=blast_output,
                blastx=blastx,
                fieldnames=fieldnames,
                extended_fieldnames=extended_fieldnames,
                cutoff=cutoff,
                molecule=molecule
            )

    @staticmethod
    def _add_headers(
        *,  # Enforce keyword arguments
        blast_output: str,
        blastx: bool,
        cutoff: float,
        extended_fieldnames: list[str],
        fieldnames: list[str],
        molecule: str = "nt",
    ):
        """
        Add headers to the BLAST output file.

        :param blast_output: Path to the BLAST output file
        :param blastx: Boolean indicating if this is a BLASTx report
        :param cutoff: Float of the minimum percent identity between query
        and subject sequence.
        :param extended_fieldnames: List of extended field names for the
        BLAST report
        :param fieldnames: List of field names for the BLAST report
        :param molecule: Molecule type (e.g., "nt" or "aa")
        """
        # Read all rows, skipping if already has header
        with open(blast_output, "r", encoding="utf-8") as report:
            lines = report.readlines()
        if not lines:
            return  # skip empty files

        # Extract header information
        header_list = lines[0].strip().split("\t")

        # If header already present, skip
        if header_list[0] == fieldnames[0]:
            return

        # Use correct fieldnames
        current_fieldnames = (
            fieldnames if len(header_list) == len(fieldnames) else
            extended_fieldnames
        )

        # Initialise a list to store the updated rows
        data = []

        # Parse rows
        for row in DictReader(
            lines, fieldnames=current_fieldnames, dialect="excel-tab"
        ):
            try:
                # Calculate the length of the subject sequence
                if molecule == 'nt' or blastx:
                    subject_length = int(row["subject_length"])
                else:
                    subject_length = int(row["subject_length"]) / 3
                # Calculate percent identity by subtracting gaps from
                # identical matches and dividing by the length of the
                # subject sequence
                percent_identity = (
                    (float(row["identical"]) - float(row["gaps"]))
                    / subject_length
                    * 100
                )
            except (KeyError, ValueError, ZeroDivisionError):
                continue

            # Filter out low percent identity matches
            if percent_identity < cutoff:
                continue

            # Add percent identity to row
            row["percent_match"] = f"{percent_identity:.2f}"
            data.append(row)

        # Write header and filtered rows
        with open(blast_output, "w", encoding="utf-8") as updated_report:
            updated_report.write("\t".join(extended_fieldnames) + "\n")
            for row in data:
                updated_report.write(
                    "\t".join(row.get(h, "") for h in extended_fieldnames)
                    + "\n"
                )

    @staticmethod
    def _find_best_blast_hits(
        *,
        blastx: bool,
        cutoff: float,
        metadata: dict,
        molecule: str
    ):
        """
        Find the best BLAST hits for each query sequence based on percent
        identity. Look for stx1 and stx2 separately
        :param blastx: Boolean indicating if this is a BLASTx report
        :param cutoff: Float of the minimum percent identity between query
        :param metadata: Dictionary of strain metadata
        :param molecule: Molecule type (e.g., "nt" or "aa")
        """

        logging.info("Finding best BLAST hits")
        for _, info in metadata.items():
            # Get the BLAST output file
            blast_output = info[f"{molecule}_blast_output"]

            # Extract the best hits from the BLAST output
            if blastx:
                best_hits, query_sequences = \
                    CombinedSubunits._extract_best_hits_subunits(
                        blast_output=blast_output,
                        cutoff=cutoff
                    )
                info[f"{molecule}_best_hits"] = best_hits
                info[f"{molecule}_query_sequences"] = query_sequences
            else:
                best_hits, novel_hits, query_sequences, query_ids = \
                    CombinedSubunits._extract_best_hits(
                        blast_output=blast_output,
                        cutoff=cutoff
                    )
                info[f"{molecule}_best_hits"] = best_hits
                info[f"{molecule}_novel_hits"] = novel_hits
                info[f"{molecule}_query_sequences"] = query_sequences
                info[f"{molecule}_query_ids"] = query_ids

        return metadata

    @staticmethod
    def _extract_best_hits(
        *,  # Enforce keyword arguments
        blast_output: str,
        cutoff: float
    ) -> tuple[dict, dict, dict, dict]:
        """
        Extract the best BLAST hits from the BLAST output file.

        :param blast_output: Path to the BLAST output file
        :param cutoff: Minimum percent identity for a hit to be considered
        Returns:
            - best_hits: closest allele match in the database (by name)
            - novel_hits: novel allele name (if <100% identity), else same as
            best_hits
            - query_sequences: dict of query sequences for novel alleles
            - query_ids: dict of query ids
        """
        best_hits = {"stx1": [], "stx2": []}
        novel_hits = {"stx1": [], "stx2": []}
        all_hits = {"stx1": [], "stx2": []}
        query_sequences = {"stx1": {}, "stx2": {}}
        query_ids = {"stx1": {}, "stx2": {}}

        with open(blast_output, "r", encoding="utf-8") as report:
            reader = DictReader(report, dialect="excel-tab")
            for row in reader:
                percent_identity = float(row.get("percent_match", 0))
                subject_id = row.get("subject_id", "").lower().split("|")[0]
                query_sequence = row.get("query_sequence", "")
                query_id = row.get("query_id", "")

                for stx in ("stx1", "stx2"):
                    if stx not in subject_id:
                        continue
                    all_hits[stx].append(
                        (
                            subject_id,
                            percent_identity,
                            query_sequence,
                            query_id
                        )
                    )

        for stx in ("stx1", "stx2"):
            if not all_hits[stx]:
                continue

            # Extract the maximum percent identity
            max_perc = max(hit[1] for hit in all_hits[stx])

            # Iterate over all hits to find the best one
            for subject_id, percent_identity, query_sequence, query_id in \
                    all_hits[stx]:

                # Check if this hit is the best one
                if percent_identity == max_perc:

                    # best_hits always points to the closest database match
                    best_hits[stx].append((subject_id, percent_identity))
                    query_ids[stx][subject_id] = query_id

                    # Novel_hits points to the novel allele name if <100%
                    # identity, else same as best_hits
                    if (
                        percent_identity < 100
                        and percent_identity >= cutoff
                        and query_sequence
                    ):
                        novel_name = CombinedSubunits._register_novel_allele(
                            sequence=query_sequence, stx_type=stx
                        )

                        # Add the novel hit to the dictionary with a 100%
                        # identity
                        novel_hits[stx].append((novel_name, 100))
                        query_sequences[stx][subject_id] = query_sequence
                        query_ids[stx][novel_name] = query_id
                    else:
                        novel_hits[stx].append((subject_id, 100))

        return best_hits, novel_hits, query_sequences, query_ids

    @staticmethod
    def _novel_allele_registry():
        """
        Singleton-style static registry for novel alleles.
        Returns a dict: {hash: {"name": ..., "sequence": ...}}
        """
        return CombinedSubunits._novel_registry

    @staticmethod
    def _register_novel_allele(
        *,  # Enforce keyword arguments
        sequence: str,
        stx_type: str,
    ) -> str:
        """
        Register or retrieve a unique name for a novel allele sequence.
        Returns the unique allele name (existing or new).

        :param sequence: The nucleotide or protein sequence
        :param stx_type: The type of Shiga toxin (e.g., "stx1" or "stx2")
        """
        registry = CombinedSubunits._novel_allele_registry()
        seq_hash = hashlib.sha256(sequence.encode()).hexdigest()
        if seq_hash in registry:
            return registry[seq_hash]["name"]
        else:
            name = f"novel_{stx_type}_{CombinedSubunits._novel_counter}"
            registry[seq_hash] = {"name": name, "sequence": sequence}
            CombinedSubunits._novel_counter += 1
            return name

    @staticmethod
    def _extract_best_hits_gene_subunits(
        *,  # Enforce keyword arguments
        blast_output: str,
        gene_subunit: str,
    ) -> tuple[dict, dict, dict]:
        """
        Extract the best BLAST hits for each gene subunit
        (Stx1A, Stx1B, Stx2A, Stx2B)
        from the BLAST output file, using subject_length as a tie breaker.

        :param blast_output: Path to the BLAST output file
        :return: A tuple containing the best hits, query sequences, and query
        IDs keyed by gene subunit (e.g., 'Stx1A', 'Stx1B', etc.)
        """
        # Initialize dictionaries for best hits, max identity, max subject
        # length, query sequences, and query IDs
        best_hits = {gene_subunit: []}
        max_identity = {gene_subunit: 0.0}
        max_length = {gene_subunit: 0}
        query_sequences = {gene_subunit: {}}
        query_ids = {gene_subunit: {}}

        with open(blast_output, "r", encoding="utf-8") as report:
            reader = DictReader(report, dialect="excel-tab")
            for row in reader:
                # Extract relevant fields from the BLAST report
                percent_identity = float(row.get("percent_match", 0))
                subject_id = row.get("subject_id", "").split("|")[0]
                query_sequence = row.get("query_sequence", "")
                query_id = row.get("query_id", "")
                try:
                    subject_length = int(row.get("subject_length", 0))
                except (TypeError, ValueError):
                    subject_length = 0

                # Find the best hit for the gene subunit, using subject_length
                # as tie breaker
                if percent_identity > max_identity[gene_subunit]:
                    max_identity[gene_subunit] = percent_identity
                    max_length[gene_subunit] = subject_length
                    best_hits[gene_subunit] = [(subject_id, percent_identity)]
                    query_sequences[gene_subunit][subject_id] = query_sequence
                # Check if percent_identity is equal to max_identity
                elif percent_identity == max_identity[gene_subunit]:
                    # Use subject_length as tie breaker
                    if subject_length > max_length[gene_subunit]:
                        max_length[gene_subunit] = subject_length
                        best_hits[gene_subunit] = [
                            (subject_id, percent_identity)
                        ]
                        query_sequences[gene_subunit] = {
                            subject_id: query_sequence
                        }
                    elif subject_length == max_length[gene_subunit]:
                        best_hits[gene_subunit].append(
                            (subject_id, percent_identity)
                        )
                        query_sequences[gene_subunit][subject_id] = \
                            query_sequence

                # Store query_id for each subject_id
                query_ids[gene_subunit][subject_id] = query_id

        return best_hits, query_sequences, query_ids

    @staticmethod
    def _extract_best_hits_subunits(
        *,  # Enforce keyword arguments
        blast_output: str,
        cutoff: float,
    ) -> tuple[dict, dict]:
        """
        Find the best combined A+B subunit hit for each allele.
        Prefer perfect matches; otherwise, find the closest combined match.
        Returns best_hits and query_sequences as before.

        :param blast_output: Path to the BLAST output file
        :param cutoff: Percent identity cutoff for considering hits
        :return: A tuple containing the best hits and query sequences
        """
        # Use the new TypedDicts for type safety and linter friendliness
        allele_hits: defaultdict[str, AlleleHitDict] = defaultdict(
            lambda: {"A": None, "B": None}
        )
        allele_stats: defaultdict[str, AlleleStatsDict] = defaultdict(
            lambda: {
                "identical": 0,
                "aligned": 0,
                "mismatches": 0,
                "gaps": 0,
                "query_seqs": {"A": "", "B": ""}
            }
        )

        with open(blast_output, "r", encoding="utf-8") as report:
            reader = DictReader(report, dialect="excel-tab")
            for row in reader:
                subject_id = row.get("subject_id", "")
                query_sequence = row.get("query_sequence", "")
                identical = int(row.get("identical", 0))
                aligned = int(row.get("alignment_length", 0))
                mismatches = int(row.get("mismatches", 0))
                gaps = int(row.get("gaps", 0))

                # Expect subject_id like stx1c_17_2_A or stx1c_17_2_B
                if subject_id.endswith("_A"):
                    base_id = subject_id[:-2]
                    allele_hits[base_id]["A"] = row
                    allele_stats[base_id]["identical"] += identical
                    allele_stats[base_id]["aligned"] += aligned
                    allele_stats[base_id]["mismatches"] += mismatches
                    allele_stats[base_id]["gaps"] += gaps
                    allele_stats[base_id]["query_seqs"]["A"] = query_sequence
                elif subject_id.endswith("_B"):
                    base_id = subject_id[:-2]
                    allele_hits[base_id]["B"] = row
                    allele_stats[base_id]["identical"] += identical
                    allele_stats[base_id]["aligned"] += aligned
                    allele_stats[base_id]["mismatches"] += mismatches
                    allele_stats[base_id]["gaps"] += gaps
                    allele_stats[base_id]["query_seqs"]["B"] = query_sequence

        best_hits = {"stx1": [], "stx2": []}
        query_sequences = {"stx1": {}, "stx2": {}}
        max_identity = {"stx1": 0.0, "stx2": 0.0}

        for base_id, stats in allele_stats.items():
            # Only consider alleles with both subunits
            a_row = allele_hits[base_id].get("A")
            b_row = allele_hits[base_id].get("B")
            if not (a_row and b_row):
                continue
            # Determine stx type
            stx_type = (
                "stx1"
                if "stx1" in base_id.lower()
                else "stx2"
                if "stx2" in base_id.lower()
                else None
            )

            if not stx_type:
                continue

            # Calculate combined percent identity
            total_identical = stats["identical"]

            # Get subject lengths for A and B
            try:
                subject_length_a = int(a_row.get("subject_length", 0))
                subject_length_b = int(b_row.get("subject_length", 0))
            except (TypeError, ValueError):
                continue

            total_subject_length = subject_length_a + subject_length_b

            if total_subject_length == 0:
                continue

            percent_identity = round(
                100.0 * total_identical / total_subject_length, 2
            )

            # Only consider above cutoff
            if percent_identity < cutoff:
                continue

            # Check for perfect match (no mismatches/gaps in either subunit)
            is_perfect = (
                int(a_row.get("mismatches", 0)) == 0
                and int(a_row.get("gaps", 0)) == 0
                and int(b_row.get("mismatches", 0)) == 0
                and int(b_row.get("gaps", 0)) == 0
            )

            # Store individual subunit sequences for downstream use
            query_sequences[stx_type][base_id + "_A"] = \
                stats["query_seqs"]["A"]
            query_sequences[stx_type][base_id + "_B"] = \
                stats["query_seqs"]["B"]

            # Store combined sequence for legacy compatibility
            query_sequences[stx_type][base_id] = (
                stats["query_seqs"]["A"] + stats["query_seqs"]["B"]
            )

            # If perfect, prefer over others
            if is_perfect and percent_identity == 100.0:
                if max_identity[stx_type] < 100.0:
                    best_hits[stx_type] = []
                max_identity[stx_type] = 100.0
                best_hits[stx_type].append((base_id, 100.0))
            elif max_identity[stx_type] < 100.0:
                # If not perfect, keep best by percent identity
                if percent_identity > max_identity[stx_type]:
                    max_identity[stx_type] = percent_identity
                    best_hits[stx_type] = [(base_id, percent_identity)]
                elif percent_identity == max_identity[stx_type]:
                    best_hits[stx_type].append((base_id, percent_identity))

        return best_hits, query_sequences

    @staticmethod
    def _write_novel_alleles_to_db(
        *,  # Enforce keyword arguments
        db_fasta_file: str,
        report_path: str,
    ):
        """
        Write all unique novel alleles to a FASTA file using SeqIO.

        :param db_fasta_file: Path to the output FASTA file
        :param report_path: Path to the report folder
        """
        # Pull novel alleles from the registry
        registry = CombinedSubunits._novel_allele_registry()

        # Read all existing records from the db_fasta_file (if it exists)
        existing_records = []
        if os.path.isfile(db_fasta_file):
            with open(db_fasta_file, "r", encoding="utf-8") as handle:
                existing_records = list(SeqIO.parse(handle, "fasta"))

        # Prepare new novel records
        novel_records = []
        for entry in registry.values():
            seq_length = len(entry["sequence"])
            record = SeqRecord(
                Seq(entry["sequence"]),
                id=f"{entry['name']}|{seq_length}nt",
                description="",
            )
            novel_records.append(record)

        # Combine all records and overwrite the db_fasta_file
        all_records = existing_records + novel_records
        if all_records:
            # Check if the existing file ends with a newline
            needs_newline = False
            if os.path.isfile(db_fasta_file):
                with open(db_fasta_file, "rb") as f:
                    f.seek(-1, os.SEEK_END)
                    last_char = f.read(1)
                    if last_char not in [b"\n", b"\r"]:
                        needs_newline = True

            with open(db_fasta_file, "w", encoding="utf-8") as handle:
                # Write all records using SeqIO (this will always start each
                # record on a new line)
                SeqIO.write(all_records, handle, "fasta")

                # If the original file did not end with a newline, add one
                if needs_newline:
                    handle.write("\n")
            logging.info(f"Novel alleles written to {db_fasta_file}")

            # Write the records to the report path as well
            novel_allele_file = os.path.join(
                report_path, "novel_alleles.fasta"
            )
            with open(novel_allele_file, "w", encoding="utf-8") as handle:
                SeqIO.write(novel_records, handle, "fasta")

            # Delete all the previous blast database files
            db_path = os.path.dirname(db_fasta_file)
            for database_file in glob(os.path.join(db_path, "*.n*")):
                os.remove(database_file)

    @staticmethod
    def _export_novel_alleles(
        *,  # Enforce keyword arguments
        blastx: bool,
        cutoff: float,
        metadata: dict,
        molecule: str,
        report_path: str,
    ) -> dict:
        """
        Export the FASTA sequence of any best hits that are >cutoff% and <100%
        percent identity to strain-specific files in the report folder.

        :param blastx: Boolean indicating if BLASTX is being used
        :param metadata: Dictionary containing BLAST hit metadata for each
        strain.
        :param molecule: Molecule type (e.g., "nt" or "aa")
        :param report_path: Path to the folder where novel allele FASTA files
        will be saved.
        :return: Updated metadata dictionary.
        """
        logging.info("Exporting novel alleles")

        for strain_name, info in metadata.items():

            # Initialise a dictionary to store novel allele paths
            info[f"novel_{molecule}_allele_paths"] = {}
            novel_allele_paths = {}

            for stx_type, best_hits in info.get(
                f"{molecule}_best_hits", {}
            ).items():
                for allele_id, percent_identity in best_hits:
                    if percent_identity < cutoff or percent_identity == 100:
                        continue

                    # Get the query sequence for the novel allele
                    query_sequence = (
                        info.get(f"{molecule}_query_sequences", {})
                        .get(stx_type, {})
                        .get(allele_id, "")
                    )
                    if not query_sequence:
                        continue

                    allele_subtype = allele_id.split("_")[0]
                    novel_allele_name = f"{strain_name}_{allele_subtype}"

                    if blastx:
                        # Extract A and B subunit sequences from the query
                        subunit_a_seq = (
                            info.get(f"{molecule}_query_sequences", {})
                            .get(stx_type, {})
                            .get(allele_id + "_A", "")
                        )
                        subunit_b_seq = (
                            info.get(f"{molecule}_query_sequences", {})
                            .get(stx_type, {})
                            .get(allele_id + "_B", "")
                        )

                        if not (subunit_a_seq and subunit_b_seq):
                            continue
                        novel_allele_paths = \
                            CombinedSubunits._write_novel_allele_subunit_files(
                                allele_id=allele_id,
                                molecule=molecule,
                                novel_allele_name=novel_allele_name,
                                path=info["path"],
                                subunit_a_seq=subunit_a_seq,
                                subunit_b_seq=subunit_b_seq,
                                report_path=report_path,
                            )
                    else:
                        novel_allele_paths = \
                            CombinedSubunits._write_novel_allele_files(
                                allele_id=allele_id,
                                molecule=molecule,
                                novel_allele_name=novel_allele_name,
                                path=info["path"],
                                query_sequence=query_sequence,
                                report_path=report_path,
                            )

            # Store the novel allele paths in the info dictionary
            info[f"novel_{molecule}_allele_paths"] = novel_allele_paths

        return metadata

    @staticmethod
    def _write_novel_allele_files(
        *,  # Enforce keyword arguments
        allele_id: str,
        molecule: str,
        novel_allele_name: str,
        path: str,
        query_sequence: str,
        report_path: str
    ) -> dict:
        """
        Write the novel allele FASTA file.

        :param allele_id: ID of the allele
        :param molecule: Molecule type (e.g., "nt" or "aa")
        :param novel_allele_name: Name of the novel allele
        :param path: Path to the strain-specific output folder
        :param query_sequence: Query sequence for the novel allele
        :param report_path: Path to the reports folder
        :return: Dictionary containing the paths to the novel alleles
        """
        # Initialise a dictionary to store the path to the novel alleles
        novel_allele_paths = {}

        # Set the path of the output file to be within the strain-
        # specific output folder
        strain_fasta_path = os.path.join(
            path, f"{novel_allele_name}_{molecule}.fasta"
        )

        # Update the novel allele paths dictionary
        novel_allele_paths[allele_id] = strain_fasta_path

        # Set the path of the FASTA file in the reports folder
        report_fasta_path = os.path.join(
            report_path, f"{novel_allele_name}_{molecule}.fasta"
        )

        # Write the novel allele to a FASTA file
        with open(strain_fasta_path, "w", encoding="utf-8") as fasta_file:
            new_record = SeqRecord(
                seq=Seq(query_sequence),
                id=novel_allele_name,
                name="",
                description="",
            )
            SeqIO.write(new_record, fasta_file, "fasta")

        # Copy the FASTA file to the report folder
        shutil.copy(strain_fasta_path, report_fasta_path)

        return novel_allele_paths

    @staticmethod
    def _write_novel_allele_subunit_files(
        *,
        allele_id: str,
        molecule: str,
        novel_allele_name: str,
        path: str,
        subunit_a_seq: str,
        subunit_b_seq: str,
        report_path: str,
        linker: str = "XXXX",
    ) -> dict:
        """
        Write FASTA files for the A subunit, B subunit, and combined
        A+linker+B sequence.

        :param allele_id: ID of the allele
        :param molecule: Molecule type (e.g., "nt" or "aa")
        :param novel_allele_name: Name of the novel allele
        :param path: Path to the strain-specific output folder
        :param subunit_a_seq: Query sequence for the A subunit
        :param subunit_b_seq: Query sequence for the B subunit
        :param report_path: Path to the reports folder
        :param linker: Linker sequence to use between A and B (default: "XXXX")
        :return: Dictionary containing the paths to the novel allele files
        """
        # Initialise a dictionary to store the paths to the novel allele files
        file_paths = {}

        # Write A subunit
        a_filename = f"{novel_allele_name}_A_{molecule}.fasta"
        a_path = os.path.join(path, a_filename)
        a_record = SeqRecord(
            Seq(subunit_a_seq), id=f"{novel_allele_name}_A", description=""
        )
        with open(a_path, "w") as f:
            SeqIO.write(a_record, f, "fasta")
        file_paths["A"] = a_path

        # Write B subunit
        b_filename = f"{novel_allele_name}_B_{molecule}.fasta"
        b_path = os.path.join(path, b_filename)
        b_record = SeqRecord(
            Seq(subunit_b_seq), id=f"{novel_allele_name}_B", description=""
        )
        with open(b_path, "w") as f:
            SeqIO.write(b_record, f, "fasta")
        file_paths["B"] = b_path

        # Write combined A+linker+B
        combined_seq = subunit_a_seq + linker + subunit_b_seq
        combined_filename = f"{novel_allele_name}_combined_{molecule}.fasta"
        combined_path = os.path.join(path, combined_filename)
        combined_record = SeqRecord(
            Seq(combined_seq),
            id=f"{novel_allele_name}_combined",
            description=f"A+{linker}+B",
        )
        with open(combined_path, "w") as f:
            SeqIO.write(combined_record, f, "fasta")
        file_paths["combined"] = combined_path

        # Also copy to report_path
        for key, src in list(file_paths.items()):
            dst = os.path.join(report_path, os.path.basename(src))
            shutil.copy(src, dst)
            file_paths[key + "_report"] = dst

        return file_paths

    @staticmethod
    def _find_aa_hits_novel_alleles(
        *,  # Enforce keyword arguments
        cpus: int,
        metadata: dict,
        num_alignments: int,
        outfmt: str,
        split_aa_dbs: dict
    ) -> dict:
        """
        Find protein hits for novel alleles using BLAST.

        :param cpus: Number of CPUs to use for BLAST
        :param metadata: Metadata dictionary containing information about the
        strains
        :param num_alignments: Number of alignments to return
        :param outfmt: Output format for BLAST
        :param split_aa_dbs: Dictionary of paths to the split AA databases
        :return: Updated metadata dictionary.
        """
        logging.info("Finding protein hits for novel alleles")
        for _, info in metadata.items():

            # Initialise the BLAST output dictionary
            info["novel_aa_blastx_outputs"] = {}

            # Find protein hits for each novel allele
            for allele_id, fasta_path in info.get(
                "novel_nt_allele_paths", {}
            ).items():

                stx_type = "stx1" if "stx1" in allele_id.lower() else "stx2"
                for subunit in ["A", "B"]:
                    db_key = f"{stx_type}{subunit}"
                    aa_db = os.path.splitext(split_aa_dbs[db_key])[0]
                    blast_output = os.path.join(
                        info["path"],
                        f"{allele_id}_novel_{stx_type}{subunit}_blastx.tsv"
                    )
                    CombinedSubunits._run_blast_cmd(
                        allele_file=aa_db,
                        blast_mode='blastx',
                        blast_output=blast_output,
                        cpus=cpus,
                        outfmt=outfmt,
                        query=fasta_path,
                        num_alignments=num_alignments
                    )
                    allele_subunit = f"{allele_id}_{subunit}"
                    if allele_subunit not in info["novel_aa_blastx_outputs"]:
                        info["novel_aa_blastx_outputs"][allele_subunit] = []
                        info["novel_aa_blastx_outputs"][
                            f"{allele_id}_{subunit}"
                        ].append(blast_output)

        return metadata

    @staticmethod
    def _add_headers_novel_aa_reports(
        *,  # Enforce keyword arguments
        blastx: bool,
        extended_fieldnames: list[str],
        fieldnames: list[str],
        metadata: dict,
        cutoff: float
    ):
        """
        Add headers to the novel reports.
        :param blastx: Boolean indicating if this is a BLASTx report
        :param extended_fieldnames: List of extended field names for the
        BLAST report
        :param fieldnames: List of field names for the BLAST report
        :param metadata: Metadata dictionary containing information about the
        strains
        :param cutoff: Cutoff value for the BLAST report
        """
        for _, info in metadata.items():

            # Add headers to the novel allele BLAST outputs
            for _, blast_outputs in info.get(
                "novel_aa_blastx_outputs", {}
            ).items():
                for blast_output in blast_outputs:
                    # Add the headers
                    CombinedSubunits._add_headers(
                        blast_output=blast_output,
                        blastx=blastx,
                        fieldnames=fieldnames,
                        extended_fieldnames=extended_fieldnames,
                        cutoff=cutoff,
                        molecule='nt'
                    )

    @staticmethod
    def _extract_novel_aa_best_hits(
        *,  # Enforce keyword arguments
        metadata: dict,
    ):
        """
        Extract novel best hits and query sequences from the metadata.

        :param metadata: Metadata dictionary containing information about the
        strains
        :return: Updated metadata dictionary.
        """
        logging.info("Extracting novel best hits and query sequences")
        for name, info in metadata.items():
            # Initialise the necessary keys in the metadata
            info["novel_aa_best_hits"] = {}
            info["novel_aa_query_sequences"] = {}
            info["novel_aa_query_ids"] = {}

            # Find protein hits for each novel allele
            for allele_id, blast_outputs in info.get(
                "novel_aa_blastx_outputs", {}
            ).items():
                # Ensure nested dicts are initialized
                if allele_id not in info["novel_aa_best_hits"]:
                    info["novel_aa_best_hits"][allele_id] = {}
                if allele_id not in info["novel_aa_query_sequences"]:
                    info["novel_aa_query_sequences"][allele_id] = {}
                if allele_id not in info["novel_aa_query_ids"]:
                    info["novel_aa_query_ids"][allele_id] = {}

                for blast_output in blast_outputs:
                    # Extract the gene and subunit from the blast output file
                    gene_subunit = CombinedSubunits.extract_gene_subunit(
                        filepath=blast_output
                    )

                    # Extract the best hits from the BLAST output
                    best_hits, query_seq, query_id = \
                        CombinedSubunits._extract_best_hits_gene_subunits(
                            blast_output=blast_output,
                            gene_subunit=gene_subunit
                        )

                    # Store the best hits in the metadata
                    info["novel_aa_best_hits"][allele_id] = best_hits
                    info["novel_aa_query_sequences"][allele_id] = query_seq
                    info["novel_aa_query_ids"][allele_id] = query_id

            for allele_id, best_hits in info.get(
                "novel_aa_best_hits", {}
            ).items():
                logging.debug(
                    "%s novel %s best hits: %s ", name, "aa", best_hits
                )
        return metadata

    @staticmethod
    def extract_gene_subunit(
        *,  # Enforce keyword arguments
        filepath: str,
    ) -> str:
        """
        Extracts the gene subunit (e.g., Stx1A) from a BLASTX output file path.
        """
        filename = os.path.basename(filepath)
        # Look for e.g. novel_stx1A_blastx.tsv or novel_stx2B_blastx.tsv
        match = re.search(
            r"novel_(stx[12][a-zA-Z]*[AB])_blastx\.tsv",
            filename,
            re.IGNORECASE
        )
        if match:
            gene_subunit = match.group(1)
            # Capitalize as Stx1A, Stx2B, etc.
            return gene_subunit[:3].capitalize() + gene_subunit[3:].upper()
        return ""

    @staticmethod
    def _assign_aa_allele_profiles_and_operons(
        *,  # Enforce keyword arguments
        metadata: dict,
        linker: str = "XXXX",
    ) -> dict:
        """
        For each strain, assign the amino acid allele profile (e.g., 22_38)
        and build the operon sequence (A + linker + B) for downstream use.

        :param metadata: Metadata dictionary containing novel_aa_best_hits and
        novel_aa_query_sequences.
        :param linker: Linker sequence to use between A and B subunits.
        :return: Updated metadata dictionary with 'aa_allele_profile' and
        'aa_operon_sequence' keys.
        """
        for strain, info in metadata.items():
            for stx in ["Stx1", "Stx2"]:
                a_allele = ""
                b_allele = ""
                a_seq = ""
                b_seq = ""
                # Search for A and B hits in the nested dicts
                for allele_id, stx_dict in info.get(
                    "novel_aa_best_hits", {}
                ).items():
                    # stx_dict is like {'Stx2A': [('Stx2A_22', 100.0)]}
                    if f"{stx}A" in stx_dict:
                        a_hits = stx_dict[f"{stx}A"]
                        if a_hits:
                            a_allele = a_hits[0][0].split("_")[-1]
                            # Get sequence for this hit
                            a_seq_dict = (
                                info.get("novel_aa_query_sequences", {})
                                .get(allele_id, {})
                                .get(f"{stx}A", {})
                            )
                            if a_seq_dict:
                                a_seq = next(iter(a_seq_dict.values()), "")
                    if f"{stx}B" in stx_dict:
                        b_hits = stx_dict[f"{stx}B"]
                        if b_hits:
                            b_allele = b_hits[0][0].split("_")[-1]
                            b_seq_dict = (
                                info.get("novel_aa_query_sequences", {})
                                .get(allele_id, {})
                                .get(f"{stx}B", {})
                            )
                            if b_seq_dict:
                                b_seq = next(iter(b_seq_dict.values()), "")
                # Only assign if both alleles found
                profile = ""
                operon_seq = ""
                if a_allele and b_allele:
                    profile = f"{a_allele}_{b_allele}"
                    operon_seq = a_seq + linker + b_seq
                info[f"{stx.lower()}_aa_allele_profile"] = profile
                info[f"{stx.lower()}_aa_operon_sequence"] = operon_seq

                # Debug statement
                if info.get(f"{stx.lower()}_aa_allele_profile"):
                    logging.info(
                        "%s %s amino acid allele profile: %s",
                        strain,
                        stx.lower(),
                        info.get(f"{stx.lower()}_aa_allele_profile", ""),
                    )
        return metadata

    @staticmethod
    def _export_aa_novel_operons(
        *,  # Enforce keyword arguments
        metadata: dict,
        report_path: str,
    ):
        """
        Export the operon sequence (A + linker + B) to a FASTA file for
        each strain/stx if either subunit is <100% identity.

        :param metadata: Metadata dictionary containing allele profiles and
        operon sequences.
        :param report_path: Path to the report directory.
        """
        logging.info("Exporting novel operon sequences")
        for strain_name, info in metadata.items():
            for stx in ["stx1", "stx2"]:
                profile = info.get(f"{stx}_aa_allele_profile", "")
                operon_seq = info.get(f"{stx}_aa_operon_sequence", "")
                if not profile or not operon_seq:
                    continue

                # Parse the alleles from the profile string
                try:
                    a_allele, b_allele = profile.split("_")
                except ValueError:
                    continue

                # Find percent identities for each subunit
                a_perc = 100.0
                b_perc = 100.0
                for _, stx_dict in info.get(
                    "novel_aa_best_hits", {}
                ).items():
                    if f"{stx.capitalize()}A" in stx_dict:
                        hits = stx_dict[f"{stx.capitalize()}A"]
                        if hits and hits[0][0].endswith(f"_{a_allele}"):
                            a_perc = hits[0][1]
                    if f"{stx.capitalize()}B" in stx_dict:
                        hits = stx_dict[f"{stx.capitalize()}B"]
                        if hits and hits[0][0].endswith(f"_{b_allele}"):
                            b_perc = hits[0][1]

                # Only export if either subunit is <100%
                if a_perc == 100.0 and b_perc == 100.0:
                    continue

                # Write the operon sequence to a FASTA file
                fasta_name = (
                    f"{strain_name}_{stx.upper()}_operon_{profile}.fasta"
                )
                fasta_path = os.path.join(report_path, fasta_name)
                record_id = f"{strain_name}_{stx.upper()}_operon_{profile}"
                record = SeqRecord(
                    Seq(operon_seq), id=record_id, description=''
                )
                with open(fasta_path, "w") as f:
                    SeqIO.write(record, f, "fasta")
                logging.info("Wrote novel operon FASTA: %s", fasta_path)

    @staticmethod
    def _export_aa_novel_alleles(
        *,  # Enforce keyword arguments
        cutoff: float,
        metadata: dict,
        molecule: str,
        report_path: str,
    ):
        """
        Export novel alleles to FASTA files.

        :param cutoff: Percent identity cutoff for considering novel alleles
        :param metadata: Metadata dictionary containing information about the
        strains
        :param molecule: Molecule type (e.g., "nt" or "aa")
        :param report_path: Path to the report directory
        """
        for strain_name, info in metadata.items():
            # The structure is:
            # info[f"novel_{molecule}_best_hits"] = {
            #     allele_id: {'stx1': [...], 'stx2': [...]}
            # }
            for _, stx_hits in info.get(
                f"novel_{molecule}_best_hits", {}
            ).items():
                for stx_type, best_hits in stx_hits.items():
                    # Check if there are any best hits for this STX type
                    if not best_hits:
                        continue

                    for best_hit in best_hits:

                        # Iterate over the hits for this STX type
                        hit_allele_id, percent_identity = best_hit

                        # Only export alleles with percent identity >cutoff and
                        # <100
                        if (
                            percent_identity < cutoff
                            or percent_identity == 100
                        ):
                            continue

                        # Retrieve the query sequence for this allele
                        query_sequence = (
                            info.get(f"{molecule}_query_sequences", {})
                            .get(stx_type, {})
                            .get(hit_allele_id, "")
                        )
                        # Check if the query sequence is empty
                        if not query_sequence:
                            continue

                        # Use the first part of the allele ID as subtype
                        # (before '_')
                        allele_subtype = hit_allele_id.split("_")[0]
                        novel_allele_name = (
                            f"{strain_name}_{molecule}_{allele_subtype}"
                        )

                        # Write the novel allele files
                        CombinedSubunits._write_novel_allele_files(
                            allele_id=hit_allele_id,
                            molecule=molecule,
                            novel_allele_name=novel_allele_name,
                            path=info["path"],
                            query_sequence=query_sequence,
                            report_path=report_path,
                        )

    @staticmethod
    def _nt_report(
        *,  # Enforce keyword arguments,
        metadata: dict,
        report_path: str,
    ):
        """
        Generate a report from the metadata with DB_Allele and NovelAllele
        columns.
        """
        logging.info("Generating report")

        # Update header
        header = "Strain\tContig\tDB_Allele\tPercentIdentity\tNovelAllele"
        body = str()

        novel = False

        if not any(
            CombinedSubunits._has_novel_aa_best_hits(info) for
            info in metadata.values()
        ):
            novel = True
        else:
            header += "\tAA_Alleles\tA_PercentIdentity\tB_PercentIdentity"

        header += "\tNotes\n"

        for name, info in metadata.items():
            data = False
            for stx in ["stx1", "stx2"]:
                if not info["nt_best_hits"][stx]:
                    continue
                data = True
                # Get the corresponding novel hits for this stx
                novel_hits = info.get("nt_novel_hits", {}).get(stx, [])
                for idx, (allele, perc_ident) in enumerate(
                    info["nt_best_hits"][stx]
                ):
                    note = "Partial Match" if perc_ident < 90 else ""
                    allele_clean = allele.split("|")[0]
                    contig = info.get(
                        "nt_query_ids", {}
                    ).get(stx, {}).get(allele, "")

                    # Get the novel allele name if present, else empty string
                    novel_allele = ""
                    if idx < len(novel_hits):
                        # Only include a novel allele if the word "novel" is
                        # present in the name
                        novel_allele = novel_hits[idx][0] if "novel" in \
                            novel_hits[idx][0] else ""

                    body += (
                        f"{name}\t{contig}\t{allele_clean}\t{perc_ident}\t"
                        f"{novel_allele}"
                    )

                    # Existing AA/notes logic unchanged
                    allele_profile = info.get(f"{stx}_aa_allele_profile", "")
                    a_perc = ""
                    b_perc = ""
                    if allele_profile:
                        try:
                            a_allele, b_allele = allele_profile.split("_")
                        except ValueError:
                            a_allele, b_allele = "", ""
                        for _, stx_dict in info.get(
                            "novel_aa_best_hits", {}
                        ).items():
                            if f"{stx.capitalize()}A" in stx_dict:
                                hits = stx_dict[f"{stx.capitalize()}A"]
                                hit = hits[0][0]
                                if hits and hit.endswith(f"_{a_allele}"):
                                    a_perc = f"{hits[0][1]:.2f}"
                                if a_perc and float(a_perc) < 100:
                                    if note:
                                        note += ";"
                                    note += "Novel Subunit A Match"
                            if f"{stx.capitalize()}B" in stx_dict:
                                hits = stx_dict[f"{stx.capitalize()}B"]
                                hit = hits[0][0]
                                if hits and hit.endswith(f"_{b_allele}"):
                                    b_perc = f"{hits[0][1]:.2f}"
                                if b_perc and float(b_perc) < 100:
                                    if note:
                                        note += ";"
                                    note += "Novel Subunit B Match"
                        body += (
                            f"\t{allele_profile}\t{a_perc}\t{b_perc}\t{note}\n"
                        )
                    else:
                        body += f"\t{note}\n"

            if not data:
                if novel:
                    body += f"{name}\t\t\t\t\t\t\n"
                else:
                    body += f"{name}\t\t\n"

        report_path = os.path.join(
            report_path, "stec_combined_nt_report.tsv"
        )

        with open(report_path, "w", encoding="utf-8") as report_file:
            report_file.write(header)
            report_file.write(body)

        logging.info("Report written to %s", report_path)

    @staticmethod
    def _has_novel_aa_best_hits(info: dict) -> bool:
        """
        Helper function to determine if there are any non-empty novel amino
        acid best hits for a given strain's info dictionary.

        :param info: Metadata dictionary for a single strain.
        :return: True if any stx1 or stx2 list in novel_aa_best_hits is
        non-empty, else False.
        """
        for allele_hits in info.get("novel_aa_best_hits", {}).values():
            for stx in ["Stx1A", "Stx1B", "Stx2A", "Stx2B"]:

                # Check if the allele hits dictionary is empty
                if allele_hits.get(stx):

                    # If there are any hits for this stx type, return True
                    if len(allele_hits[stx]) > 0:
                        return True
        return False

    @staticmethod
    def _aa_report(
        *,  # Enforce keyword arguments,
        metadata: dict,
        report_path: str
    ):
        """
        Write the report for blastx/tblastx analyses
        """
        logging.info("Generating report")

        # Initialise the header and body strings
        header = "Strain\tAllele\tPercentIdentity\tNotes\n"
        body = str()

        # Iterate over each strain in the metadata
        for name, info in metadata.items():
            # Initialise data flag to False
            data = False

            # Iterate over each stx type
            for stx in ["stx1", "stx2"]:
                # Check if there are any best hits for this stx type
                if not info["aa_best_hits"][stx]:
                    continue
                # Set data flag to True
                data = True

                # Iterate over each allele in the results
                for allele, perc_ident in info["aa_best_hits"][stx]:
                    # Check if the allele is a partial match
                    note = "Partial Match" if perc_ident < 90 else ""

                    # Remove anything after a | in the allele e.g. |1236nt
                    allele = allele.split("|")[0]

                    # Add the row to the body
                    body += f"{name}\t{allele}\t{perc_ident:.2f}\t{note}\n"

                # Check if the were outputs for this strain
                if not data:
                    body += f"{name}\t\t\n"

        # Set the name and path of the report
        report_path = os.path.join(
            report_path, "stec_combined_aa_report.tsv"
        )

        # Write the report
        with open(report_path, "w", encoding="utf-8") as report_file:
            report_file.write(header)
            report_file.write(body)

        logging.info("Report written to %s", report_path)

    @staticmethod
    def _combined_report(
        *,  # Enforce keyword arguments,
        metadata: dict,
        report_path: str
    ):
        """
        Write the combined report for both blastn and tblastx analyses
        """
        logging.info("Generating combined report")

        # Initialise the header and body strings
        header = (
            "Strain\tnt_Allele\tnt_PercentIdentity\taa_Allele\t"
            "aa_PercentIdentity\tNotes\n"
        )
        body = str()

        # Iterate over each strain in the metadata
        for name, info in metadata.items():

            # Initialise data flag to False
            data = False

            body += f"{name}\t"
            # Iterate over each stx type
            for stx in ["stx1", "stx2"]:
                # Check if there are any best hits for this stx type
                if (
                    not info["nt_best_hits"][stx]
                    and not info["aa_best_hits"][stx]
                ):
                    continue
                # Set data flag to True
                data = True

                # Create a variable to store if there are partial sequences
                partial = False

                # Iterate over both the nt_best_hits and the aa_best_hits
                for analysis in ['nt_best_hits', 'aa_best_hits']:
                    # Iterate over each allele in the results
                    for allele, perc_ident in info[analysis][stx]:

                        # Check if the allele is a partial match
                        partial = \
                            "Partial Match" if perc_ident < 90 else ""

                        # Remove anything after a | in the allele e.g. |1236nt
                        allele = allele.split("|")[0]

                        # Add the row to the body
                        body += f"{allele}\t{perc_ident}\t"

                # Check if there are any partial sequences
                if partial:
                    body += f"{partial}\n"
                else:
                    body += '\n'

                # Check if the were outputs for this strain
                if not data:
                    body += f"{name}\t\t\t\t\n"

        # Set the name and path of the report
        report_path = os.path.join(
            report_path, "stec_combined_nt_aa_report.tsv"
        )

        # Write the report
        with open(report_path, "w", encoding="utf-8") as report_file:
            report_file.write(header)
            report_file.write(body)

        logging.info("Report written to %s", report_path)


def main():
    """
    Main function to parse arguments and run the combined subunits processing.
    """
    parser = ArgumentParser(
        description="Process STEC alleles with both A and B subunits."
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    parser.add_argument(
        "-a",
        "--allele_path",
        metavar="allele_path",
        default=os.path.join(os.getcwd(), "alleles"),
        help="Specify name and path of folder containing allele files. "
        "If not provided, the nt_alleles folder in the current working "
        "directory will be used by default",
    )
    parser.add_argument(
        '-r', '--report_path',
        metavar='report_path',
        default=os.path.join(os.getcwd(), 'reports'),
        help='Specify name and path of folder into which reports are to be '
        'placed. If not provided, the reports folder in the current working '
        'directory will be used'
    )
    parser.add_argument(
        "-q",
        "--query_path",
        metavar="query_path",
        default=os.path.join(os.getcwd(), "query"),
        help="Specify name and path of folder containing query files in FASTA "
        "format. If not provided, the query folder in the current working "
        "directory will be used",
    )
    parser.add_argument(
        "--blast_mode",
        choices=["blastn", "tblastx", "blastx", "blastn+tblastx"],
        default="blastx",
        help="Choose BLAST mode: blastn (default), tblastx, blastx, or both "
        "blastn and tblastx. If both are selected, the program will run "
        "blastn on all samples and tblastx on any novel alleles",
    )
    parser.add_argument(
        '--preliminary',
        action='store_true',
        help='Run a preliminary screen on genomes using blastn (do not run '
        'additional downstream analyses)'
    )
    parser.add_argument(
        "--split_aa_db_dir",
        metavar="split_aa_db_dir",
        help=(
            "Path to directory containing split AA allele databases "
            "(Stx1A_aa_*.fasta, Stx1B_aa_*.fasta, Stx2A_aa_*.fasta, "
            "Stx2B_aa_*.fasta)"
        ),
    )
    parser.add_argument(
        "-n",
        "--num_alignments",
        metavar="num_alignments",
        type=int,
        default=10,
        help="Specify the number of alignments to return. Default is 10."
    )
    parser.add_argument(
        '-c',
        "--cutoff",
        metavar="cutoff",
        type=float,
        default=99.0,
        help="Specify the percent identity cutoff. Default is 99.0."
    )
    parser.add_argument(
        "-v",
        "--verbosity",
        choices=["debug", "info", "warning", "error", "critical"],
        metavar="verbosity",
        default="info",
        help="Set the logging level. Options are debug, info, warning, error, "
        "and critical. Default is info."
    )

    # Get the arguments into an object
    arguments = parser.parse_args()

    # Set up logging
    setup_logging(arguments=arguments)

    # Start timing
    start_time = time.time()

    # Create a CombinedSubunits instance and process the alleles
    combined_subunits = CombinedSubunits(arguments)
    combined_subunits.process()

    # Print elapsed time
    elapsed = time.time() - start_time
    logging.info("Pipeline completed in %.2f seconds.", elapsed)


if __name__ == "__main__":
    main()
