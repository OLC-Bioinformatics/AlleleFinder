#! /usr/bin/env python3

"""
Unit tests for the FindDuplicateAlleles class.
"""

from collections import defaultdict
from pathlib import Path
import sys

# Third-party imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

# Local imports
from allele_tools.find_duplicate_alleles import (
    FindDuplicateAlleles,
    main
)


def make_fasta(
    tmp_path: Path,
    records: list[SeqRecord]
) -> str:
    """
    Write a list of SeqRecord objects to a FASTA file.

    :param tmp_path: Path object representing the temporary directory
    :type tmp_path: Path
    :param records: List of SeqRecord objects to write
    :type records: list[SeqRecord]
    :return: Path to the created FASTA file as a string
    :rtype: str
    """
    fasta_path = tmp_path / "alleles.fasta"

    # Write each record to the FASTA file with specified encoding
    with open(fasta_path, "w", encoding="utf-8") as f:
        for rec in records:
            f.write(f">{rec.id}\n{rec.seq}\n")

    return str(fasta_path)


def test_extract_records(
    tmp_path: Path
) -> None:
    """
    Test the _extract_records method of FindDuplicateAlleles.

    This function writes a set of SeqRecord objects to a FASTA file,
    then verifies that _extract_records correctly maps sequences to IDs
    and IDs to SeqRecord objects.

    :param tmp_path: Path object representing the temporary directory
    :type tmp_path: Path
    :return: None
    """
    # Prepare test records
    records = [
        SeqRecord(Seq("ATGC"), id="allele1"),
        SeqRecord(Seq("ATGC"), id="allele2"),
        SeqRecord(Seq("GGGG"), id="allele3"),
    ]

    # Write records to a FASTA file
    fasta = make_fasta(tmp_path, records)

    # Extract records using the method under test
    seq_to_ids, id_to_record = FindDuplicateAlleles._extract_records(
        allele_file=fasta
    )

    # Check that duplicate sequences are mapped to correct IDs
    assert set(seq_to_ids["ATGC"]) == {"allele1", "allele2"}
    assert seq_to_ids["GGGG"] == ["allele3"]

    # Check that ID to SeqRecord mapping is correct
    assert id_to_record["allele1"].seq == Seq("ATGC")
    assert id_to_record["allele3"].id == "allele3"


def test_find_duplicates() -> None:
    """
    Test the _find_duplicates method of FindDuplicateAlleles.

    This function creates a mapping from sequences to lists of allele IDs,
    then verifies that _find_duplicates correctly identifies groups of
    duplicate alleles (i.e., those with more than one ID per sequence).

    :return: None
    """
    # Prepare a mapping of sequences to allele IDs, with some duplicates
    seq_to_ids: defaultdict[str, list[str]] = defaultdict(
        list,
        {
            "ATGC": ["a1", "a2"],
            "GGGG": ["b1"],
            "CCCC": ["c1", "c2", "c3"],
        },
    )

    # Find duplicates using the method under test
    duplicates: list[list[str]] = FindDuplicateAlleles._find_duplicates(
        seq_to_ids=seq_to_ids
    )

    # Check that duplicate groups are present
    assert ["a1", "a2"] in duplicates
    assert ["c1", "c2", "c3"] in duplicates

    # Ensure all groups are lists
    assert all(isinstance(group, list) for group in duplicates)

    # Ensure singletons are not included
    assert not any("b1" in group for group in duplicates)


def test_write_duplicates(
    tmp_path: Path
) -> None:
    """
    Test the _write_duplicates method of FindDuplicateAlleles.

    This function writes groups of duplicate allele IDs to a text file,
    then verifies that the file contains the correct header and
    semicolon-separated groups.

    :param tmp_path: Path object representing the temporary directory
    :type tmp_path: Path
    :return: None
    """
    # Prepare output file path and duplicate groups
    output_file = tmp_path / "dups.txt"
    duplicates: list[list[str]] = [["a1", "a2"], ["b1", "b2", "b3"]]

    # Write duplicates using the method under test
    FindDuplicateAlleles._write_duplicates(
        duplicates=duplicates,
        output_file=str(output_file)
    )

    # Read and verify the output file content
    lines = output_file.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "Duplicate_Alleles"
    assert "a1;a2" in lines[1]
    assert "b1;b2;b3" in lines[2]


def test_write_unique_alleles(
    tmp_path: Path
) -> None:
    """
    Test the _write_unique_alleles method of FindDuplicateAlleles.

    This function writes unique allele sequences to a FASTA file,
    ensuring that only the first occurrence of each sequence is present.

    :param tmp_path: Path object representing the temporary directory
    :type tmp_path: Path
    :return: None
    """
    # Prepare test records with duplicate and unique sequences
    records = [
        SeqRecord(Seq("ATGC"), id="allele1"),
        SeqRecord(Seq("ATGC"), id="allele2"),
        SeqRecord(Seq("GGGG"), id="allele3"),
    ]

    # Write records to a FASTA file
    fasta = make_fasta(tmp_path, records)

    # Extract records using the method under test
    seq_to_ids, id_to_record = FindDuplicateAlleles._extract_records(
        allele_file=fasta
    )

    # Prepare output file path for unique alleles
    output_file = tmp_path / "unique.fasta"

    # Write unique alleles using the method under test
    FindDuplicateAlleles._write_unique_alleles(
        seq_to_ids=seq_to_ids,
        id_to_record=id_to_record,
        output_file=str(output_file)
    )

    # Read the output FASTA file and collect IDs
    ids = [
        rec.id
        for rec in list(
            SeqIO.parse(str(output_file), "fasta")
        )
    ]

    # Only the first occurrence of each sequence should be present
    assert "allele1" in ids
    assert "allele3" in ids
    assert "allele2" not in ids  # duplicate, not first


def test_main_prints_and_writes(
    tmp_path: Path,
    capsys
) -> None:
    """
    Test that FindDuplicateAlleles.main() prints expected output and writes
    correct files when run with a FASTA containing duplicates.

    :param tmp_path: Path object representing the temporary directory
    :type tmp_path: Path
    :param capsys: Pytest fixture to capture stdout/stderr
    :type capsys: Any
    :return: None
    """
    # Prepare a FASTA file with duplicate and unique sequences
    records = [
        SeqRecord(Seq("ATGC"), id="a1"),
        SeqRecord(Seq("ATGC"), id="a2"),
        SeqRecord(Seq("GGGG"), id="b1"),
    ]

    # Write records to a FASTA file
    fasta = make_fasta(tmp_path, records)

    # Create a dummy args object with the allele_file attribute
    args = type(
        "Args",
        (),
        {"allele_file": fasta}
    )()

    # Instantiate and run the FindDuplicateAlleles main method
    finder = FindDuplicateAlleles(args)
    finder.main()

    # Check printed output for duplicate alleles
    captured = capsys.readouterr()
    assert "Duplicate alleles found:" in captured.out
    assert "a1, a2" in captured.out

    # Check that output files exist
    dups = tmp_path / "duplicate_alleles.txt"
    uniq = tmp_path / "unique_alleles.fasta"
    assert dups.exists()
    assert uniq.exists()

    # Read unique_alleles.fasta and collect IDs
    ids = [
        rec.id
        for rec in list(
            SeqIO.parse(str(uniq), "fasta")
        )
    ]

    # Only one of the ATGC alleles should be present
    assert "a1" in ids or "a2" in ids
    assert not ("a1" in ids and "a2" in ids)
    assert "b1" in ids


def make_fasta_file(
    tmp_path: Path,
    records: list[SeqRecord]
) -> str:
    """
    Write a list of SeqRecord objects to a FASTA file.

    :param tmp_path: Path object representing the temporary directory
    :type tmp_path: Path
    :param records: List of SeqRecord objects to write
    :type records: list[SeqRecord]
    :return: Path to the created FASTA file as a string
    :rtype: str
    """
    fasta_path = tmp_path / "test.fasta"

    # Write each record to the FASTA file with specified encoding
    with open(fasta_path, "w", encoding="utf-8") as f:
        for rec in records:
            f.write(f">{rec.id}\n{rec.seq}\n")

    return str(fasta_path)


def test_main_creates_outputs_and_prints(
    monkeypatch,
    tmp_path: Path,
    capsys
) -> None:
    """
    Test that main() creates output files and prints expected output
    when run with a FASTA containing duplicate and unique alleles.

    :param monkeypatch: Pytest fixture for patching
    :type monkeypatch: Any
    :param tmp_path: Path object for temporary directory
    :type tmp_path: Path
    :param capsys: Pytest fixture to capture stdout/stderr
    :type capsys: Any
    :return: None
    """
    # Prepare a FASTA file with duplicate and unique sequences
    records = [
        SeqRecord(Seq("ATGC"), id="a1"),
        SeqRecord(Seq("ATGC"), id="a2"),
        SeqRecord(Seq("GGGG"), id="b1"),
    ]

    fasta_file: str = make_fasta_file(tmp_path, records)

    # Patch sys.argv to simulate command-line call
    monkeypatch.setattr(
        sys,
        "argv",
        ["prog", "-a", fasta_file]
    )

    # Run main function
    main()

    # Capture and check printed output for duplicate alleles
    captured = capsys.readouterr()
    assert "Duplicate alleles found:" in captured.out
    assert "a1, a2" in captured.out

    # Check that output files exist
    dups = Path(fasta_file).parent / "duplicate_alleles.txt"
    uniq = Path(fasta_file).parent / "unique_alleles.fasta"
    assert dups.exists()
    assert uniq.exists()

    # Read and verify duplicate_alleles.txt content
    with open(dups, encoding="utf-8") as f:
        lines = f.read().splitlines()
    assert lines[0] == "Duplicate_Alleles"
    assert "a1;a2" in lines[1]

    # Read unique_alleles.fasta and collect IDs
    ids = [
        rec.id
        for rec in SeqIO.parse(str(uniq), "fasta")
    ]

    # Only one of the ATGC alleles should be present
    assert "a1" in ids or "a2" in ids
    assert not ("a1" in ids and "a2" in ids)
    assert "b1" in ids


def test_main_requires_argument(
    monkeypatch
) -> None:
    """
    Test that main() exits with SystemExit when required argument is missing.

    This function patches sys.argv to simulate a command-line call with
    missing required arguments, then verifies that main() raises SystemExit.

    :param monkeypatch: Pytest fixture for patching
    :type monkeypatch: Any
    :return: None
    """
    # Patch sys.argv to simulate missing argument
    monkeypatch.setattr(sys, "argv", ["prog"])

    # Expect SystemExit due to missing required argument
    with pytest.raises(SystemExit):
        main()
