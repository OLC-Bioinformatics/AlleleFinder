#! /usr/bin/env python3

"""
Unit and integration tests for allele_tools/stec_combined_subunits.py
"""

# Standard imports
from csv import DictWriter
import os
import pathlib
import shutil
import subprocess
import tempfile
import types
from unittest.mock import (
    patch,
    MagicMock
)

# Third-party imports
from Bio import SeqIO
import pytest

# pylint: disable=import-error,no-name-in-module
# Local imports
from allele_tools.stec_combined_subunits import (
    main,
    CombinedSubunits,
)

# pylint: disable=protected-access
# Import aliases for CombinedSubunits methods
_aa_report = CombinedSubunits._aa_report
_add_headers = CombinedSubunits._add_headers
_add_headers_novel_aa_reports = CombinedSubunits._add_headers_novel_aa_reports
_analyze_sequence_for_gaps_and_ambiguity = \
    CombinedSubunits._analyze_sequence_for_gaps_and_ambiguity
_assign_aa_allele_profiles_and_operons = \
    CombinedSubunits._assign_aa_allele_profiles_and_operons
_combined_report = CombinedSubunits._combined_report
_export_aa_novel_operons = CombinedSubunits._export_aa_novel_operons
_extract_best_hits = CombinedSubunits._extract_best_hits
_extract_best_hits_gene_subunits = \
    CombinedSubunits._extract_best_hits_gene_subunits
_extract_best_hits_subunits = CombinedSubunits._extract_best_hits_subunits
_export_aa_novel_alleles = CombinedSubunits._export_aa_novel_alleles
_export_novel_alleles = CombinedSubunits._export_novel_alleles
_extract_novel_aa_best_hits = CombinedSubunits._extract_novel_aa_best_hits
_find_aa_hits_novel_alleles = CombinedSubunits._find_aa_hits_novel_alleles
_find_best_blast_hits = CombinedSubunits._find_best_blast_hits
_find_split_aa_dbs = CombinedSubunits._find_split_aa_dbs
_has_novel_aa_best_hits = CombinedSubunits._has_novel_aa_best_hits
_locate_fasta_files = CombinedSubunits._locate_fasta_files
_novel_allele_registry = CombinedSubunits._novel_allele_registry
_novel_registry = CombinedSubunits._novel_registry
_nt_report = CombinedSubunits._nt_report
_parseable_blast_outputs = CombinedSubunits._parseable_blast_outputs
_prep_query_files = CombinedSubunits._prep_query_files
_register_novel_allele = CombinedSubunits._register_novel_allele
_run_blast = CombinedSubunits._run_blast
_run_blast_cmd = CombinedSubunits._run_blast_cmd
_write_novel_allele_files = CombinedSubunits._write_novel_allele_files
_write_novel_allele_subunit_files = \
    CombinedSubunits._write_novel_allele_subunit_files
_write_novel_alleles_to_db = CombinedSubunits._write_novel_alleles_to_db
extract_gene_subunit = CombinedSubunits.extract_gene_subunit
make_blast_db = CombinedSubunits.make_blast_db
split_fasta = CombinedSubunits.split_fasta


# pylint: disable=redefined-outer-name
@pytest.fixture(scope="module")
def variables():
    """
    Fixture to provide a Variables object with temporary directories and
    dummy files for testing.

    :return: Variables object with paths and dummy files
    :rtype: Variables
    """
    class Variables:
        """
        Helper class to store paths and dummy files for tests.
        """
        def __init__(self) -> None:
            """
            Initialize temporary directories and dummy files for testing.
            """
            self.tmpdir = tempfile.mkdtemp()
            self.allele_path = os.path.join(self.tmpdir, "alleles")
            self.report_path = os.path.join(self.tmpdir, "reports")
            self.query_path = os.path.join(self.tmpdir, "query")
            self.split_aa_db_dir = os.path.join(
                self.tmpdir, "split_aa_db"
            )

            # Create required directories
            os.makedirs(self.allele_path, exist_ok=True)
            os.makedirs(self.report_path, exist_ok=True)
            os.makedirs(self.query_path, exist_ok=True)
            os.makedirs(self.split_aa_db_dir, exist_ok=True)

            # Create a dummy allele FASTA file
            self.allele_fasta = os.path.join(
                self.allele_path, "allele1.fasta"
            )
            with open(self.allele_fasta, "w", encoding='utf-8') as f:
                f.write(">stx1A_1\nATGCATGCATGC\n")

            # Create a dummy query FASTA file
            self.query_fasta = os.path.join(
                self.query_path, "query1.fasta"
            )
            with open(self.query_fasta, "w", encoding='utf-8') as f:
                f.write(">query1\nATGCATGCATGC\n")

            # Create dummy split AA db FASTA files
            for stx in ["Stx1A", "Stx1B", "Stx2A", "Stx2B"]:
                with open(
                    os.path.join(
                        self.split_aa_db_dir, f"{stx}_aa_test.fasta"
                    ),
                    "w", encoding='utf-8'
                ) as f:
                    f.write(f">{stx}_1\nMTESTSEQ\n")

            # Mimic argparse.Namespace for arguments
            self.args = type(
                "Args",
                (),
                {
                    "allele_path": self.allele_path,
                    "report_path": self.report_path,
                    "query_path": self.query_path,
                    "blast_mode": "blastn",
                    "preliminary": False,
                    "num_alignments": 5,
                    "cutoff": 99.0,
                    "split_aa_db_dir": self.split_aa_db_dir,
                    "verbosity": "info"
                }
            )()

    v = Variables()
    yield v

    # Clean up temporary directory after tests
    shutil.rmtree(v.tmpdir)


def test_locate_fasta_files(variables) -> None:
    """
    Test that _locate_fasta_files finds the dummy allele FASTA file.

    :param variables: Variables fixture with test paths and files
    :type variables: Variables
    :return: None
    """
    files, errors = _locate_fasta_files(
        directory=variables.allele_path,
        errors=[]
    )

    # Check that the dummy allele FASTA file is found and no errors occurred
    assert variables.allele_fasta in files
    assert errors == []


def test_find_split_aa_dbs(variables) -> None:
    """
    Test that _find_split_aa_dbs locates all required split AA db FASTA files.

    :param variables: Variables fixture with test paths and files
    :type variables: Variables
    :return: None
    """
    dbs: dict = _find_split_aa_dbs(
        split_aa_db_dir=variables.split_aa_db_dir
    )

    # Ensure all expected keys are present
    assert set(dbs.keys()) == {"stx1A", "stx1B", "stx2A", "stx2B"}

    # Ensure all values are FASTA file paths
    for v in dbs.values():
        assert v.endswith(".fasta")


def test_split_fasta(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the split_fasta function to ensure that a sequence containing the
    linker is split into two records with correct IDs.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    input_fasta = tmp_path / "input.fasta"
    output_fasta = tmp_path / "output.fasta"

    # Sequence with linker
    with open(input_fasta, "w", encoding="utf-8") as f:
        f.write(">rec1\nATGXXXXTTT\n")

    split_fasta(
        input_fasta=str(input_fasta),
        output_fasta=str(output_fasta)
    )

    # Parse the output FASTA and check the number and IDs of records
    records = list(
        SeqIO.parse(str(output_fasta), "fasta")
    )
    assert len(records) == 2
    assert records[0].id.endswith("_A")
    assert records[1].id.endswith("_B")


def test_register_novel_allele_registry() -> None:
    """
    Test the _register_novel_allele and _novel_allele_registry methods.

    Ensures that the registry correctly stores unique novel alleles by
    sequence and stx_type, and that the same sequence returns the same
    novel allele name.

    :return: None
    """
    # Clear the registry and reset the counter
    _novel_registry.clear()

    # Register a novel allele for stx1
    name1 = _register_novel_allele(
        sequence="ATGC",
        stx_type="stx1"
    )

    # Register the same sequence again for stx1
    name2 = _register_novel_allele(
        sequence="ATGC",
        stx_type="stx1"
    )

    # Register a different sequence for stx2
    name3 = _register_novel_allele(
        sequence="GGGG",
        stx_type="stx2"
    )

    # The same sequence and stx_type should yield the same name
    assert name1 == name2

    # Names should start with the correct prefix
    assert name1.startswith("novel_stx1_")
    assert name3.startswith("novel_stx2_")

    # The registry should contain two unique entries
    reg = _novel_allele_registry()
    assert len(reg) == 2


def test_prep_query_files(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _prep_query_files function to ensure that it creates the
    correct metadata for query files.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    qf = tmp_path / "q1.fasta"
    with open(qf, "w", encoding="utf-8") as f:
        f.write(">q1\nATGC\n")

    query_files: list[str] = [str(qf)]
    query_path: str = str(tmp_path)

    # Call the method under test
    meta: dict = _prep_query_files(
        query_files=query_files,
        query_path=query_path
    )

    # Check that the metadata contains the expected keys and paths
    assert "q1" in meta
    assert os.path.isdir(meta["q1"]["path"])
    assert (
        os.path.islink(meta["q1"]["file"])
        or os.path.exists(meta["q1"]["file"])
    )


def test_add_headers(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _add_headers function to ensure that it adds the correct header
    and percent_match field to the BLAST output file.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output = tmp_path / "blast.tsv"
    fieldnames: list[str] = [
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
    extended_fieldnames: list[str] = fieldnames.copy()
    extended_fieldnames.insert(14, "percent_match")

    # Write a fake BLAST output file with encoding specified
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write(
            "query1\tsubj1\t10\t0\t0\t1e-10\t50\t12\t12\t12\t1\t12\t1\t12\t"
            "ATGCATGCATGC\tATGCATGCATGC\n"
        )

    # Call the method under test
    _add_headers(
        blast_output=str(blast_output),
        blastx=False,
        cutoff=90.0,
        extended_fieldnames=extended_fieldnames,
        fieldnames=fieldnames,
        molecule="nt",
    )

    # Read the output file and check for correct header and percent_match
    with open(blast_output, "r", encoding="utf-8") as f:
        lines = f.readlines()

    assert lines[0].startswith("query_id")
    assert "percent_match" in lines[0]


def test_extract_best_hits(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _extract_best_hits function to ensure that it correctly parses
    BLAST output and identifies the best hits for stx1 and stx2.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output = tmp_path / "blast.tsv"
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write(
            "query_id\tsubject_id\tidentical\tmismatches\tgaps\tevalue\t"
            "bit_score\tquery_length\tsubject_length\talignment_length\t"
            "query_start\tquery_end\tsubject_start\tsubject_end\t"
            "percent_match\tquery_sequence\n"
        )
        f.write(
            "q1\tstx1_1\t10\t0\t0\t1e-10\t50\t12\t12\t12\t1\t12\t1\t12\t100.0"
            "\tATGC\n"
        )

    # Call the method under test
    best, _, __, ___, ____, _____ = _extract_best_hits(
        blast_output=str(blast_output),
        cutoff=99.0,
        molecule="nt"
    )

    # Check that the best hit for stx1 is correctly identified
    assert "stx1" in best
    assert best["stx1"][0][0].startswith("stx1_1")


def test_extract_best_hits_subunits(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _extract_best_hits_subunits function to ensure that it correctly
    parses BLAST output and identifies the best hits for stx1 subunits.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output = tmp_path / "blast.tsv"
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write(
            "query_id\tsubject_id\tidentical\tmismatches\tgaps\tevalue\t"
            "bit_score\tquery_length\tsubject_length\talignment_length\t"
            "query_start\tquery_end\tsubject_start\tsubject_end\t"
            "query_sequence\tpercent_match\n"
        )
        f.write(
            "q1\tstx1A_1_A\t10\t0\t0\t1e-10\t50\t12\t6\t6\t1\t6\t1\t6\t"
            "ATGC\t100.0\n"
        )
        f.write(
            "q1\tstx1A_1_B\t10\t0\t0\t1e-10\t50\t12\t6\t6\t1\t6\t1\t6\t"
            "ATGC\t100.0\n"
        )

    # Call the method under test
    best, _ = _extract_best_hits_subunits(
        blast_output=str(blast_output),
        cutoff=99.0
    )

    # Check that the best hit for stx1 is correctly identified
    assert "stx1" in best
    assert best["stx1"][0][0].startswith("stx1A_1")


def test_write_novel_allele_files(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _write_novel_allele_files function to ensure that it writes
    the novel allele FASTA file to both the specified path and the report
    directory.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Create a separate report directory
    report_dir = tmp_path / "report"
    report_dir.mkdir()

    # Call the method under test
    out: dict[str, str] = _write_novel_allele_files(
        allele_id="stx1_1",
        molecule="nt",
        novel_allele_name="strain1_stx1",
        path=str(tmp_path),
        query_sequence="ATGC",
        report_path=str(report_dir),
    )

    # Check that the output dictionary contains the expected key
    assert "stx1_1" in out

    # Check that the FASTA file was written to the output path
    assert os.path.isfile(out["stx1_1"])

    # Check that the file was also copied to the report directory
    report_fasta = report_dir / "strain1_stx1_nt.fasta"
    assert os.path.isfile(report_fasta)


def test_write_novel_allele_subunit_files(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _write_novel_allele_subunit_files function to ensure that it
    writes the novel allele subunit FASTA files (A, B, combined) to both
    the specified path and the report directory.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Use a separate directory for report_path to avoid SameFileError
    report_dir = tmp_path / "report"
    report_dir.mkdir()

    # Case 1: Sequences too short (should not write files)
    file_paths, notes = _write_novel_allele_subunit_files(
        molecule="nt",
        novel_allele_name="strain1_stx1",
        path=str(tmp_path),
        subunit_a_seq="ATGC",
        subunit_b_seq="GGGG",
        report_path=str(report_dir),
    )
    # Should not contain file paths
    assert "combined" not in file_paths
    assert file_paths == {}
    # Should contain notes for A and B
    assert "A" in notes and "B" in notes
    assert "too short" in notes["A"]
    assert "too short" in notes["B"]

    # Case 2: Sequences long enough (should write files)
    long_a = "A" * 313
    long_b = "B" * 82
    file_paths, notes = _write_novel_allele_subunit_files(
        molecule="nt",
        novel_allele_name="strain1_stx1",
        path=str(tmp_path),
        subunit_a_seq=long_a,
        subunit_b_seq=long_b,
        report_path=str(report_dir),
    )
    assert "A" in file_paths and "B" in file_paths and "combined" in file_paths
    assert notes == {}


def test_extract_gene_subunit() -> None:
    """
    Test the extract_gene_subunit function to ensure it correctly extracts
    the gene subunit from the given file path.

    :return: None
    """
    # Test with a file path containing "stx1A"
    fn: str = "/tmp/novel_stx1A_blastx.tsv"
    out: str = extract_gene_subunit(filepath=fn)
    assert out == "Stx1A"

    # Test with a file path containing "stx2B"
    fn2: str = "/tmp/novel_stx2B_blastx.tsv"
    assert extract_gene_subunit(filepath=fn2) == "Stx2B"

    # Test with a file path that does not match any gene subunit
    assert extract_gene_subunit(
        filepath="/tmp/other.tsv"
    ) == ""


def test_assign_aa_allele_profiles_and_operons() -> None:
    """
    Test the _assign_aa_allele_profiles_and_operons function to ensure it
    correctly assigns allele profiles and operon sequences from the provided
    metadata.

    :return: None
    """
    meta: dict = {
        "strain1": {
            "novel_aa_best_hits": {
                "allele1": {
                    "Stx1A": [("Stx1A_22", 99.0)],
                    "Stx1B": [("Stx1B_38", 98.0)],
                }
            },
            "novel_aa_query_sequences": {
                "allele1": {
                    "Stx1A": {"Stx1A_22": "AAA"},
                    "Stx1B": {"Stx1B_38": "BBB"},
                }
            },
        }
    }

    # Call the method under test
    out: dict = _assign_aa_allele_profiles_and_operons(
        metadata=meta
    )

    # Check that the allele profile and operon sequence are as expected
    assert out["strain1"]["stx1_aa_allele_profile"] == "22_38"
    assert out["strain1"]["stx1_aa_operon_sequence"] == "AAAXXXXBBB"


def test_has_novel_aa_best_hits() -> None:
    """
    Test the _has_novel_aa_best_hits function to ensure it correctly
    identifies when novel amino acid best hits are present.

    :return: None
    """
    info: dict = {
        "novel_aa_best_hits": {
            "a": {"Stx1A": [("Stx1A_1", 99.0)]}
        }
    }
    # Should return True when best hits are present
    assert _has_novel_aa_best_hits(info)

    info = {"novel_aa_best_hits": {"a": {}}}
    # Should return False when no best hits are present
    assert not _has_novel_aa_best_hits(info)


def test_init_and_process_and_run(variables) -> None:
    """
    Test initialization, process, and run methods of CombinedSubunits

    This test patches all methods that perform file or BLAST operations to
    avoid side effects and focuses on the control flow of initialization,
    process, and run.

    :param variables: Fixture providing test paths and dummy files
    :type variables: Variables
    :return: None
    """
    # Patch methods that do file/BLAST work to avoid side effects
    with patch.object(
        CombinedSubunits,
        "make_blast_db",
        return_value=variables.allele_fasta,
    ), patch.object(
        CombinedSubunits,
        "_run_blast",
        return_value={},
    ), patch.object(
        CombinedSubunits,
        "_parseable_blast_outputs",
        return_value=None,
    ), patch.object(
        CombinedSubunits,
        "_find_best_blast_hits",
        return_value={},
    ), patch.object(
        CombinedSubunits,
        "_write_novel_alleles_to_db",
        return_value=None,
    ), patch.object(
        CombinedSubunits,
        "_nt_report",
        return_value=None,
    ), patch.object(
        CombinedSubunits,
        "_aa_report",
        return_value=None,
    ), patch.object(
        CombinedSubunits,
        "_combined_report",
        return_value=None,
    ), patch.object(
        CombinedSubunits,
        "_export_novel_alleles",
        return_value={},
    ), patch.object(
        CombinedSubunits,
        "_find_aa_hits_novel_alleles",
        return_value={},
    ), patch.object(
        CombinedSubunits,
        "_add_headers_novel_aa_reports",
        return_value=None,
    ), patch.object(
        CombinedSubunits,
        "_extract_novel_aa_best_hits",
        return_value={},
    ), patch.object(
        CombinedSubunits,
        "_assign_aa_allele_profiles_and_operons",
        return_value={},
    ), patch.object(
        CombinedSubunits,
        "_export_aa_novel_operons",
        return_value=None,
    ), patch.object(
        CombinedSubunits,
        "_export_aa_novel_alleles",
        return_value=None,
    ):
        # Initialize CombinedSubunits with test arguments
        cs = CombinedSubunits(variables.args)
        cs.query_files = [variables.query_fasta]
        cs.metadata = {
            "query1": {
                "path": variables.query_path,
                "file": variables.query_fasta,
            }
        }

        # Run the process method (should call patched methods)
        cs.process()

        # Run the run method for a specific mode and molecule
        cs.run(
            blast_mode="blastn",
            molecule="nt",
        )


def test_write_novel_alleles_to_db(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _write_novel_alleles_to_db function to ensure that it writes
    the novel alleles to both the database FASTA file and the report
    directory.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare a dummy registry with a novel allele
    _novel_registry.clear()
    _register_novel_allele(
        sequence="ATGC",
        stx_type="stx1"
    )

    db_fasta: pathlib.Path = tmp_path / "db.fasta"
    report_path: pathlib.Path = tmp_path

    # Call the method under test
    _write_novel_alleles_to_db(
        db_fasta_file=str(db_fasta),
        report_path=str(report_path)
    )

    # Check that the FASTA file was written to the output path
    assert os.path.isfile(db_fasta)

    # Check that the file was also copied to the report directory
    assert os.path.isfile(
        os.path.join(report_path, "novel_alleles.fasta")
    )

    # Delete the novel alleles FASTA file
    if os.path.exists(os.path.join(report_path, "novel_alleles.fasta")):
        os.remove(os.path.join(report_path, "novel_alleles.fasta"))


def test_extract_best_hits_gene_subunits(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _extract_best_hits_gene_subunits function to ensure that it
    correctly parses BLAST output and applies tie-breaker logic for best
    hits.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output = tmp_path / "blast.tsv"

    # Write a BLAST output file with two hits for the same gene subunit
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write(
            "query_id\tsubject_id\tidentical\tmismatches\tgaps\tevalue\t"
            "bit_score\tquery_length\tsubject_length\talignment_length\t"
            "query_start\tquery_end\tsubject_start\tsubject_end\t"
            "percent_match\tquery_sequence\n"
        )
        f.write(
            "q1\tStx1A_1\t10\t0\t0\t1e-10\t50\t12\t12\t12\t1\t12\t1\t12\t"
            "99.0\tATGC\n"
        )
        f.write(
            "q2\tStx1A_2\t10\t0\t0\t1e-10\t50\t12\t13\t12\t1\t12\t1\t12\t"
            "99.0\tGGGG\n"
        )

    # Call the method under test
    best, seqs, ids = _extract_best_hits_gene_subunits(
        blast_output=str(blast_output),
        gene_subunit="Stx1A"
    )

    # Check that the best hits, sequences, and ids are as expected
    assert "Stx1A" in best
    assert isinstance(seqs, dict)
    assert isinstance(ids, dict)


def test_add_headers_empty_and_header_present(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _add_headers function for early returns when the BLAST output
    file is empty or already contains a header.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output = tmp_path / "blast.tsv"

    # Create an empty file to test early return on empty file
    with open(blast_output, "w", encoding="utf-8"):
        pass

    _add_headers(
        blast_output=str(blast_output),
        blastx=False,
        cutoff=90.0,
        extended_fieldnames=["query_id"],
        fieldnames=["query_id"],
        molecule="nt",
    )

    # Create a file with a header already present to test early return
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write("query_id\nrow1\n")

    _add_headers(
        blast_output=str(blast_output),
        blastx=False,
        cutoff=90.0,
        extended_fieldnames=["query_id"],
        fieldnames=["query_id"],
        molecule="nt",
    )

    # Check that the header is present in the file
    with open(blast_output, encoding="utf-8") as f:
        assert f.readline().startswith("query_id")


def test_extract_best_hits_no_hits(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _extract_best_hits function when there are no hits in the BLAST
    output file.

    Ensures that the function returns empty dictionaries for best hits,
    novel hits, sequences, and IDs.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output = tmp_path / "blast.tsv"

    # Write an empty BLAST output file with only the header
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write("query_id\tsubject_id\tpercent_match\tquery_sequence\n")

    best, novel, _, __, ___, ____ = _extract_best_hits(
        blast_output=str(blast_output),
        cutoff=99.0,
        molecule="nt"
    )

    # Check that all returned dictionaries are empty for both stx1 and stx2
    assert best == {"stx1": [], "stx2": []}
    assert novel == {"stx1": [], "stx2": []}


def test_extract_best_hits_subunits_partial(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _extract_best_hits_subunits with only one subunit present.

    Ensures that if only one subunit is present in the BLAST output,
    no best hits are returned for either stx1 or stx2.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output = tmp_path / "blast.tsv"

    # Write a BLAST output file with only one subunit present
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write(
            "query_id\tsubject_id\tidentical\tmismatches\tgaps\t"
            "alignment_length\tquery_sequence\n"
        )
        f.write("q1\tstx1A_1_A\t10\t0\t0\t6\tATGC\n")

    # Call the method under test
    best, _ = _extract_best_hits_subunits(
        blast_output=str(blast_output), cutoff=99.0
    )

    # Since only one subunit is present, best hits should be empty
    assert best == {"stx1": [], "stx2": []}


def test_extract_gene_subunit_case_insensitive() -> None:
    """
    Test extract_gene_subunit for case-insensitive matching.

    Ensures that the function correctly identifies gene subunits
    regardless of case in the file path.

    :return: None
    """
    # Covers extract_gene_subunit with lower/upper case in the file path

    fn: str = "/tmp/novel_stx1a_blastx.tsv"
    assert extract_gene_subunit(filepath=fn) == "Stx1A"


def test_assign_aa_allele_profiles_and_operons_missing_keys() -> None:
    """
    Test _assign_aa_allele_profiles_and_operons with missing keys.

    Ensures that when the required keys are missing from the metadata,
    the function returns empty strings for both the allele profile and
    operon sequence.

    :return: None
    """
    meta: dict = {"strain1": {"novel_aa_best_hits": {}}}

    # Call the method under test
    out: dict = _assign_aa_allele_profiles_and_operons(
        metadata=meta
    )

    # Check that the output contains empty strings for missing keys
    assert out["strain1"]["stx1_aa_allele_profile"] == ""
    assert out["strain1"]["stx1_aa_operon_sequence"] == ""


def test_has_novel_aa_best_hits_empty() -> None:
    """
    Test the _has_novel_aa_best_hits function with empty input.

    Ensures that the function returns False when the input dictionary is
    empty.

    :return: None
    """
    info: dict = {}

    # Should return False when input is empty
    assert not _has_novel_aa_best_hits(info)


def test_locate_fasta_files_no_files(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _locate_fasta_files when no FASTA files are present in the directory.

    Ensures that the function returns an empty list of files and an error
    message indicating that no FASTA files were found.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Covers error path for _locate_fasta_files

    errors: list[str] = []
    files, errors = _locate_fasta_files(
        directory=str(tmp_path),
        errors=errors,
    )

    # Check that no files are found and an appropriate error is returned
    assert files == []
    assert errors and "No FASTA files found" in errors[0]


def test_find_split_aa_dbs_missing(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _find_split_aa_dbs for missing required split AA db FASTA files.

    Ensures that a FileNotFoundError is raised when not all required
    split AA db FASTA files are present in the directory.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :raises FileNotFoundError: If any required split AA db FASTA file is
    missing
    :return: None
    """
    # Create the split directory but do not add all required files
    os.makedirs(tmp_path / "split", exist_ok=True)

    # Should raise FileNotFoundError due to missing files
    with pytest.raises(FileNotFoundError):
        _find_split_aa_dbs(
            split_aa_db_dir=str(tmp_path / "split")
        )


def test_make_blast_db_nt(tmp_path: pathlib.Path) -> None:
    """
    Test the make_blast_db function for the 'nt' molecule branch.

    Ensures that the function returns the expected output path when
    creating a nucleotide BLAST database, and that subprocess.run is
    called (simulated here).

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Create a dummy FASTA file for testing
    fasta = tmp_path / "allele.fasta"
    with open(fasta, "w", encoding="utf-8") as f:
        f.write(">a\nATGC\n")

    # Simulate makeblastdb not present by patching subprocess.run
    with patch("subprocess.run") as run:
        run.return_value = None

        out: str = make_blast_db(
            allele_file=str(fasta),
            allele_path=str(tmp_path),
            molecule="nt"
        )

        # The output path should end with 'allele'
        assert out.endswith("allele")


def test_make_blast_db_aa(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test the make_blast_db function for the 'aa' molecule branch.

    Ensures that the function returns the expected output path when
    creating an amino acid BLAST database, and that subprocess.run is
    called (simulated here).

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Create a dummy FASTA file for testing
    fasta = tmp_path / "allele.fasta"
    with open(fasta, "w", encoding="utf-8") as f:
        f.write(">a\nATGC\n")

    # Simulate makeblastdb not present by patching subprocess.run
    with patch("subprocess.run") as run:
        run.return_value = None

        out: str = make_blast_db(
            allele_file=str(fasta),
            allele_path=str(tmp_path),
            molecule="aa"
        )

        # The output path should end with 'allele'
        assert out.endswith("allele")


def test_make_blast_db_blastx(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test the make_blast_db function with blastx=True.

    Ensures that the function returns the expected output path when
    creating a BLASTX split allele database, and that subprocess.run is
    called (simulated here).

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Create a dummy FASTA file for testing
    fasta = tmp_path / "allele.fasta"
    with open(fasta, "w", encoding="utf-8") as f:
        f.write(">a\nATGCXXXXATGC\n")

    # Simulate makeblastdb not present by patching subprocess.run
    with patch("subprocess.run") as run:
        run.return_value = None

        out: str = make_blast_db(
            allele_file=str(fasta),
            allele_path=str(tmp_path),
            blastx=True,
            molecule="aa",
        )

        # The output path should indicate a blastx split allele database
        assert "blastx_split_allele_database" in out


def test_split_fasta_no_linker(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the split_fasta function when no linker is present in the sequence.

    Ensures that the function writes a single record to the output FASTA
    when the input sequence does not contain the linker.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Create input FASTA file without linker
    fasta = tmp_path / "in.fasta"
    out = tmp_path / "out.fasta"

    with open(fasta, "w", encoding="utf-8") as f:
        f.write(">a\nATGCATGC\n")

    # Call the method under test
    split_fasta(
        input_fasta=str(fasta),
        output_fasta=str(out)
    )

    # Parse the output FASTA and check that only one record is present
    records = list(
        SeqIO.parse(str(out), "fasta")
    )
    assert len(records) == 1


def test_run_blast_cmd_blastn_existing(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _run_blast_cmd for the 'blastn' mode when the output file already
    exists.

    Ensures that the function returns early and does not attempt to run
    the BLAST command if the output file is present.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output = tmp_path / "out.tsv"

    # Create the output file to simulate its existence
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write("dummy")

    # Patch the BLAST command to avoid running the actual command
    with patch("Bio.Blast.Applications.NcbiblastnCommandline") as cmd:
        cmd.return_value = lambda: None

        result: str = _run_blast_cmd(
            allele_file="db",
            blast_mode="blastn",
            blast_output=str(blast_output),
            cpus=1,
            outfmt="6",
            query="q",
        )

        # The result should be a string (the output file path)
        assert isinstance(result, str)


def test_run_blast_cmd_tblastx(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test the _run_blast_cmd function for the 'tblastx' mode.

    Ensures that the function returns the expected output path when
    running tblastx, and that the BLAST command is called (simulated).

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Covers _run_blast_cmd tblastx branch

    blast_output: pathlib.Path = tmp_path / "out.tsv"

    # Patch the BLAST command to avoid running the actual command
    with patch(
        "Bio.Blast.Applications.NcbitblastxCommandline"
    ) as cmd:
        cmd.return_value = lambda: None

        result: str = _run_blast_cmd(
            allele_file="db",
            blast_mode="tblastx",
            blast_output=str(blast_output),
            cpus=1,
            outfmt="6",
            query="q",
        )

        # The result should be a string (the output file path)
        assert isinstance(result, str)


def test_run_blast_cmd_blastx(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test the _run_blast_cmd function for the 'blastx' mode.

    Ensures that the function returns the expected output path when
    running blastx, and that the BLAST command is called (simulated).

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Covers _run_blast_cmd blastx branch

    blast_output: pathlib.Path = tmp_path / "out.tsv"

    # Patch the BLAST command to avoid running the actual command
    with patch(
        "Bio.Blast.Applications.NcbiblastxCommandline"
    ) as cmd:
        cmd.return_value = lambda: None

        result: str = _run_blast_cmd(
            allele_file="db",
            blast_mode="blastx",
            blast_output=str(blast_output),
            cpus=1,
            outfmt="6",
            query="q",
        )

        # The result should be a string (the output file path)
        assert isinstance(result, str)


def test_run_blast_cmd_application_error(
    tmp_path: pathlib.Path, caplog
) -> None:
    """
    Test _run_blast_cmd for ApplicationError handling.

    Ensures that if the BLAST command raises an exception, the function
    returns a string and logs an error message.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :param caplog: Pytest fixture for capturing log messages
    :type caplog: Any
    :return: None
    """
    blast_output: pathlib.Path = tmp_path / "out.tsv"

    # Remove the output file if it exists to simulate a fresh run
    if os.path.exists(blast_output):
        os.remove(blast_output)

    class DummyCmd:
        """
        Dummy command class to simulate a BLAST command that fails.
        """
        def __call__(self) -> None:
            """
            Simulate running the command and raise an exception.
            """
            raise OSError("fail")

        def __str__(self) -> str:
            """
            Return a string representation of the dummy command.
            """
            return "dummy-blast-cmd"

    # Patch the BLAST command to use DummyCmd, which raises an exception
    with patch(
        "Bio.Blast.Applications.NcbiblastnCommandline",
        return_value=DummyCmd()
    ):
        with caplog.at_level("ERROR"):
            result: str = _run_blast_cmd(
                allele_file="db",
                blast_mode="blastn",
                blast_output=str(blast_output),
                cpus=1,
                outfmt="6",
                query="q",
            )

        # The result should be a string (the output file path or error)
        assert isinstance(result, str)

        # Check that an error message was logged
        assert any(
            "BLAST application error" in m for m in caplog.messages
        )


def test_run_blast_cmd_invalid_mode(
    tmp_path: pathlib.Path
) -> None:
    """
    Test _run_blast_cmd for an invalid blast_mode.

    Ensures that the function returns an empty string when an unknown
    blast_mode is provided.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Should return '' for unknown blast_mode
    blast_output = tmp_path / "out.tsv"
    result: str = _run_blast_cmd(
        allele_file="dummy_db",
        blast_mode="invalid_mode",
        blast_output=str(blast_output),
        cpus=1,
        outfmt="6",
        query="dummy_query",
    )
    assert result == ""


def test_run_blast_cmd_empty_mode(
    tmp_path: pathlib.Path
) -> None:
    """
    Test _run_blast_cmd for an empty blast_mode.

    Ensures that the function returns an empty string when an empty
    blast_mode is provided.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Should return '' for empty blast_mode

    blast_output: pathlib.Path = tmp_path / "out.tsv"
    result: str = _run_blast_cmd(
        allele_file="dummy_db",
        blast_mode="",
        blast_output=str(blast_output),
        cpus=1,
        outfmt="6",
        query="dummy_query",
    )

    # The result should be an empty string
    assert result == ""


def test_add_headers_bad_row(tmp_path: pathlib.Path) -> None:
    """
    Test the _add_headers function with bad data rows.

    This test covers cases where the BLAST output contains rows with
    invalid data types (ValueError), missing keys (KeyError), or
    division by zero (ZeroDivisionError). It ensures that the function
    skips such rows and writes only valid data to the output file.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output = tmp_path / "blast.tsv"

    # Write BLAST output with bad and incomplete data rows
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write(
            "query_id\tsubject_id\tidentical\tmismatches\tgaps\t"
            "query_sequence\n"
        )
        f.write("q1\tsubj1\tbad\tbad\tbad\tATGC\n")
        # Row missing subject_length
        f.write("q2\tsubj1\t1\t1\t1\tATGC\n")

    # Call the method under test
    _add_headers(
        blast_output=str(blast_output),
        blastx=False,
        cutoff=90.0,
        extended_fieldnames=[
            "query_id",
            "subject_id",
            "identical",
            "mismatches",
            "gaps",
            "query_sequence",
            "percent_match",
        ],
        fieldnames=[
            "query_id",
            "subject_id",
            "identical",
            "mismatches",
            "gaps",
            "query_sequence",
        ],
        molecule="nt",
    )

    # Read the output file and check that the header is present
    with open(blast_output, encoding="utf-8") as f:
        lines = f.readlines()

    assert lines[0].startswith("query_id")


def test_extract_best_hits_gene_subunits_bad_subject_length(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _extract_best_hits_gene_subunits with a bad subject_length value.

    Ensures that the function handles rows where subject_length is not
    an integer and still returns the expected structure.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output = tmp_path / "blast.tsv"

    # Write a BLAST output file with a bad subject_length value
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write(
            "query_id\tsubject_id\tpercent_match\tquery_sequence"
            "\tsubject_length\n"
        )
        f.write("q1\tStx1A_1\t99.0\tATGC\tbad\n")

    # Call the method under test
    best, _, __ = _extract_best_hits_gene_subunits(
        blast_output=str(blast_output),
        gene_subunit="Stx1A"
    )

    # Check that the best dictionary contains the gene subunit key
    assert "Stx1A" in best


def test_extract_best_hits_gene_subunits_empty(
    tmp_path: pathlib.Path
) -> None:
    """
    Test _extract_best_hits_gene_subunits with an empty BLAST output file.

    Ensures that the function returns an empty dictionary for best hits
    when the BLAST output file contains only the header and no data rows.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output = tmp_path / "blast.tsv"

    # Write only the header to the BLAST output file
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write("query_id\tsubject_id\tpercent_match\tquery_sequence\n")

    best, _, __ = _extract_best_hits_gene_subunits(
        blast_output=str(blast_output),
        gene_subunit="Stx1A"
    )

    # The best dictionary should contain the gene_subunit key with an empty
    # list
    assert best == {"Stx1A": []}


def test_extract_best_hits_subunits_bad_subject_length(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _extract_best_hits_subunits with bad subject_length values.

    Ensures that the function handles rows where subject_length is not
    an integer and still returns the expected structure.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output = tmp_path / "blast.tsv"

    # Write a BLAST output file with bad subject_length values for both
    # subunits
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write(
            "subject_id\tidentical\talignment_length\tmismatches\tgaps\t"
            "query_sequence\tsubject_length\n"
        )
        f.write("stx1A_1_A\t10\t6\t0\t0\tATGC\tbad\n")
        f.write("stx1A_1_B\t10\t6\t0\t0\tATGC\tbad\n")

    # Call the method under test
    best, _ = _extract_best_hits_subunits(
        blast_output=str(blast_output),
        cutoff=99.0,
    )

    # Since subject_length is invalid, best_hits should be empty
    assert best == {"stx1": [], "stx2": []}


def test_write_novel_alleles_to_db_existing_file(
    tmp_path: pathlib.Path
) -> None:
    """
    Test _write_novel_alleles_to_db with an existing file and newline logic.

    Ensures that the function appends novel alleles to an existing FASTA
    file, handles newlines correctly, and writes the output to the report
    directory.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Clear the registry and register a novel allele for stx1
    _novel_registry.clear()
    _register_novel_allele(
        sequence="ATGC",
        stx_type="stx1"
    )

    db_fasta: pathlib.Path = tmp_path / "db.fasta"

    # Create an existing FASTA file with encoding specified
    with open(db_fasta, "w", encoding="utf-8") as f:
        f.write(">old\nATGC\n")

    report_path: pathlib.Path = tmp_path

    # Call the method under test
    _write_novel_alleles_to_db(
        db_fasta_file=str(db_fasta),
        report_path=str(report_path)
    )

    # Check that the FASTA file was written to the output path
    assert os.path.isfile(db_fasta)

    # Check that the file was also copied to the report directory
    assert os.path.isfile(
        os.path.join(report_path, "novel_alleles.fasta")
    )


def test_extract_gene_subunit_no_match() -> None:
    """
    Test extract_gene_subunit when no match is found in the file path.

    Ensures that the function returns an empty string when the file path
    does not contain a recognized gene subunit.

    :return: None
    """
    # Covers extract_gene_subunit when no match is found

    fn: str = "/tmp/novel_otherfile.tsv"
    assert extract_gene_subunit(filepath=fn) == ""


def test_init_errors(tmp_path: pathlib.Path) -> None:
    """
    Test CombinedSubunits initialization error handling.

    Ensures that if required directories are missing, the error_print
    function is called and SystemExit is raised.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    class Args:
        """
        Mimics argparse.Namespace for CombinedSubunits arguments.
        """
        allele_path: str = str(tmp_path / "missing_alleles")
        report_path: str = str(tmp_path / "missing_reports")
        query_path: str = str(tmp_path / "missing_query")
        blast_mode: str = "blastn"
        preliminary: bool = False
        num_alignments: int = 1
        cutoff: float = 99.0
        split_aa_db_dir: str = ''
        verbosity: str = "info"

    # Patch error_print to raise SystemExit and check it is called
    with patch(
        "allele_tools.stec_combined_subunits.error_print",
        side_effect=SystemExit
    ) as ep:
        with pytest.raises(SystemExit):
            CombinedSubunits(Args())
        assert ep.called


def test_make_blast_db_nt_error(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the make_blast_db function for the 'nt' molecule branch when a
    CalledProcessError is raised.

    Ensures that the function raises SystemExit if the subprocess running
    makeblastdb fails.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Create a dummy FASTA file for testing
    fasta = tmp_path / "allele.fasta"
    with open(fasta, "w", encoding="utf-8") as f:
        f.write(">a\nATGC\n")

    # Patch subprocess.run to raise CalledProcessError
    with patch(
        "subprocess.run",
        side_effect=subprocess.CalledProcessError(1, "makeblastdb")
    ):
        # Should raise SystemExit due to the error
        with pytest.raises(SystemExit):
            make_blast_db(
                allele_file=str(fasta),
                allele_path=str(tmp_path),
                molecule="nt"
            )


def test_make_blast_db_aa_error(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the make_blast_db function for the 'aa' molecule branch when a
    CalledProcessError is raised.

    Ensures that the function raises SystemExit if the subprocess running
    makeblastdb fails.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Create a dummy FASTA file for testing
    fasta: pathlib.Path = tmp_path / "allele.fasta"
    with open(fasta, "w", encoding="utf-8") as f:
        f.write(">a\nATGC\n")

    # Patch subprocess.run to raise CalledProcessError
    with patch(
        "subprocess.run",
        side_effect=subprocess.CalledProcessError(1, "makeblastdb")
    ):
        # Should raise SystemExit due to the error
        with pytest.raises(SystemExit):
            make_blast_db(
                allele_file=str(fasta),
                allele_path=str(tmp_path),
                molecule="aa"
            )


def test_split_fasta_with_linker(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the split_fasta function when a linker is present in the sequence.

    Ensures that the function splits the sequence into two records with
    correct IDs when the linker is present.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Covers split_fasta with linker present

    fasta: pathlib.Path = tmp_path / "in.fasta"
    out: pathlib.Path = tmp_path / "out.fasta"

    # Write a sequence with a linker to the input FASTA file
    with open(fasta, "w", encoding="utf-8") as f:
        f.write(">a\nATGCXXXXGGGG\n")

    # Call the split_fasta function under test
    split_fasta(
        input_fasta=str(fasta),
        output_fasta=str(out)
    )

    # Parse the output FASTA and check the number and IDs of records
    records = list(
        SeqIO.parse(str(out), "fasta")
    )
    assert len(records) == 2
    assert records[0].id.endswith("_A")
    assert records[1].id.endswith("_B")


def test_nt_report_and_aa_report(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _nt_report and _aa_report functions with various metadata.

    Ensures that the nucleotide and amino acid reports are generated
    correctly and written to the expected files.

    :param tmp_path: Temporary directory for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare metadata for nucleotide report
    meta: dict = {
        "strain1": {
            "nt_best_hits": {"stx1": [("stx1_1", 99.0)], "stx2": []},
            "nt_novel_hits": {"stx1": [("novel_stx1_1", 100)], "stx2": []},
            "nt_query_ids": {"stx1": {"stx1_1": "contig1"}, "stx2": {}},
            "stx1_aa_allele_profile": "22_38",
            "novel_aa_best_hits": {
                "allele1": {
                    "Stx1A": [("Stx1A_22", 99.0)],
                    "Stx1B": [("Stx1B_38", 98.0)],
                }
            },
        }
    }

    # Call the nucleotide report function
    _nt_report(
        metadata=meta,
        report_path=str(tmp_path)
    )

    # Check that the nucleotide report file exists
    assert os.path.isfile(
        os.path.join(tmp_path, "stec_combined_nt_report.tsv")
    )

    # Prepare metadata for amino acid report
    meta2: dict = {
        "strain1": {
            "aa_best_hits": {"stx1": [("stx1_1", 99.0)], "stx2": []}
        }
    }

    # Call the amino acid report function
    _aa_report(
        metadata=meta2,
        report_path=str(tmp_path)
    )

    # Check that the amino acid report file exists
    assert os.path.isfile(
        os.path.join(tmp_path, "stec_combined_aa_report.tsv")
    )


def test_combined_report(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _combined_report function to ensure that it generates the
    combined nucleotide and amino acid report file when both nt and aa
    best hits are present.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare metadata with both nt and aa best hits for stx1
    meta: dict = {
        "strain1": {
            "nt_best_hits": {"stx1": [("stx1_1", 99.0)], "stx2": []},
            "aa_best_hits": {"stx1": [("stx1_1", 99.0)], "stx2": []},
        }
    }

    # Call the combined report function under test
    _combined_report(
        metadata=meta,
        report_path=str(tmp_path)
    )

    # Check that the combined report file exists in the output directory
    assert os.path.isfile(
        os.path.join(tmp_path, "stec_combined_nt_aa_report.tsv")
    )


def test_export_aa_novel_operons(tmp_path: pathlib.Path) -> None:
    """
    Test the _export_aa_novel_operons function to ensure that it writes
    the operon FASTA file when a partial match is present.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare metadata with a partial match for the operon
    meta: dict = {
        "strain1": {
            "stx1_aa_allele_profile": "22_38",
            "stx1_aa_operon_sequence": "A" * 313 + "XXXX" + "B" * 82,
            "stx1_aa_subunit_a_sequence": "A" * 313,
            "stx1_aa_subunit_b_sequence": "B" * 82,
            "novel_aa_best_hits": {
                "allele1": {
                    "Stx1A": [("Stx1A_22", 99.0)],
                    "Stx1B": [("Stx1B_38", 98.0)],
                }
            },
        }
    }

    # Call the method under test
    _export_aa_novel_operons(metadata=meta, report_path=str(tmp_path))

    # Check that the operon FASTA file was written to the output directory
    files = list(tmp_path.glob("*operon*.fasta"))
    assert files


def test_export_aa_novel_alleles(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _export_aa_novel_alleles function to ensure that it handles
    the case where a novel allele is present in the metadata.

    Ensures that the function does not raise an exception even if no files
    are written.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Covers _export_aa_novel_alleles with a novel allele

    # Create a report directory for output
    report_dir: pathlib.Path = tmp_path / "report"
    report_dir.mkdir()

    # Prepare metadata with a novel allele best hit
    meta: dict = {
        "strain1": {
            "novel_aa_best_hits": {
                "allele1": {"stx1": [("stx1_1", 99.0)]}
            },
            "aa_query_sequences": {
                "stx1": {"stx1_1": "ATGC"}
            },
            "path": str(tmp_path),  # Use tmp_path for strain output
        }
    }

    # Should not raise, even if nothing is written
    _export_aa_novel_alleles(
        cutoff=98.0,
        metadata=meta,
        molecule="aa",
        report_path=str(report_dir)
    )


def test_locate_fasta_files_nonexistent_dir(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _locate_fasta_files with a non-existent directory.

    Ensures that the function returns an empty list of files and an error
    message indicating that no FASTA files were found.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Covers error branch for non-existent directory

    non_existent: pathlib.Path = tmp_path / "does_not_exist"
    errors: list[str] = []
    files, errors = _locate_fasta_files(
        directory=str(non_existent),
        errors=errors,
    )

    # Check that no files are found and an appropriate error is returned
    assert files == []
    assert errors and "No FASTA files found" in errors[0]


def test_find_split_aa_dbs_missing_file(
    tmp_path: pathlib.Path
) -> None:
    """
    Test _find_split_aa_dbs when a required split AA db FASTA file is missing.

    Ensures that a FileNotFoundError is raised if not all required split
    AA db FASTA files are present in the directory.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :raises FileNotFoundError: If any required split AA db FASTA file is
        missing
    :return: None
    """
    # Only create Stx1A, not Stx1B, Stx2A, Stx2B

    db_dir: pathlib.Path = tmp_path / "split_aa_db"
    db_dir.mkdir()

    # Create only the Stx1A FASTA file with encoding specified
    (db_dir / "Stx1A_aa_test.fasta").write_text(
        ">Stx1A_1\nMSEQ\n", encoding="utf-8"
    )

    # Should raise FileNotFoundError due to missing files
    with pytest.raises(FileNotFoundError):
        _find_split_aa_dbs(split_aa_db_dir=str(db_dir))


def test_split_fasta_linker_at_start_end(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the split_fasta function when the linker is at the start or end
    of the sequence.

    Ensures that the function splits the sequence into two records with
    correct IDs when the linker is present at the start or end.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    input_fasta: pathlib.Path = tmp_path / "input.fasta"
    output_fasta: pathlib.Path = tmp_path / "output.fasta"

    # Linker at start
    with open(input_fasta, "w", encoding="utf-8") as f:
        f.write(">rec1\nXXXXATGC\n")

    split_fasta(
        input_fasta=str(input_fasta),
        output_fasta=str(output_fasta)
    )

    # Parse the output FASTA and check the number of records
    records = list(
        SeqIO.parse(str(output_fasta), "fasta")
    )
    assert len(records) == 2

    # Linker at end
    with open(input_fasta, "w", encoding="utf-8") as f:
        f.write(">rec2\nATGCXXXX\n")

    split_fasta(
        input_fasta=str(input_fasta),
        output_fasta=str(output_fasta)
    )

    # Parse the output FASTA and check the number of records
    records = list(
        SeqIO.parse(str(output_fasta), "fasta")
    )
    assert len(records) == 2


def test_write_novel_allele_files_same_path(
    tmp_path: pathlib.Path
) -> None:
    """
    Test _write_novel_allele_files when path and report_path are the same.

    Ensures that a SameFileError is raised if the output path and report
    directory are the same, which would cause a file copy conflict.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :raises shutil.SameFileError: If path and report_path are the same
    :return: None
    """
    # Should raise SameFileError if path and report_path are the same
    with pytest.raises(shutil.SameFileError):
        _write_novel_allele_files(
            allele_id="stx1_1",
            molecule="nt",
            novel_allele_name="strain1_stx1",
            path=str(tmp_path),
            query_sequence="ATGC",
            report_path=str(tmp_path),
        )


def test_export_aa_novel_operons_missing_profile(
    tmp_path: pathlib.Path
) -> None:
    """
    Test _export_aa_novel_operons with missing allele profile and operon
    sequence.

    Ensures that the function does not raise an exception when both the
    allele profile and operon sequence are empty.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    meta: dict = {
        "strain1": {
            "stx1_aa_allele_profile": "",
            "stx1_aa_operon_sequence": "",
            "novel_aa_best_hits": {},
        }
    }

    # Should not raise an exception when called with missing profile/sequence
    _export_aa_novel_operons(
        metadata=meta,
        report_path=str(tmp_path)
    )


def test_export_aa_novel_operons_bad_profile(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _export_aa_novel_operons with a bad allele profile string.

    Ensures that the function does not raise an exception when the
    allele profile cannot be split (e.g., missing expected delimiter).

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    meta: dict = {
        "strain1": {
            "stx1_aa_allele_profile": "bad_profile",
            "stx1_aa_operon_sequence": "AAAXXXXBBB",
            "novel_aa_best_hits": {},
        }
    }

    # Should not raise an exception when called with a bad profile
    _export_aa_novel_operons(
        metadata=meta,
        report_path=str(tmp_path)
    )


def test_extract_best_hits_no_percent_match(
    tmp_path: pathlib.Path
) -> None:
    """
    Test _extract_best_hits when percent_match column is missing.

    Ensures that the function returns a best hit with 0.0 percent identity
    if the percent_match column is not present in the BLAST output.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output: pathlib.Path = tmp_path / "blast.tsv"

    # Write a BLAST output file without percent_match column
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write("query_id\tsubject_id\tquery_sequence\n")
        f.write("q1\tstx1_1\tATGC\n")

    # Call the method under test
    best, _, __, ___, ____, _____ = _extract_best_hits(
        blast_output=str(blast_output),
        cutoff=99.0,
        molecule="nt"
    )

    # Should return 0.0 percent identity for the best hit
    assert best == {"stx1": [("stx1_1", 0.0)], "stx2": []}


def test_add_headers_zero_division(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _add_headers function for ZeroDivisionError handling.

    Ensures that if subject_length is zero, the function skips the row
    and only writes the header to the output file.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    blast_output: pathlib.Path = tmp_path / "blast.tsv"

    # Write a BLAST output file with subject_length = 0
    with open(blast_output, "w", encoding="utf-8") as f:
        f.write(
            "query_id\tsubject_id\tidentical\tmismatches\tgaps\t"
            "subject_length\tquery_sequence\n"
        )
        f.write("q1\tsubj1\t1\t0\t0\t0\tATGC\n")  # subject_length=0

    # Call the method under test
    _add_headers(
        blast_output=str(blast_output),
        blastx=False,
        cutoff=90.0,
        extended_fieldnames=[
            "query_id",
            "subject_id",
            "identical",
            "mismatches",
            "gaps",
            "subject_length",
            "query_sequence",
            "percent_match",
        ],
        fieldnames=[
            "query_id",
            "subject_id",
            "identical",
            "mismatches",
            "gaps",
            "subject_length",
            "query_sequence",
        ],
        molecule="nt",
    )

    # Read the output file and check that only the header is present
    with open(blast_output, encoding="utf-8") as f:
        lines = f.readlines()

    assert lines[0].startswith("query_id")


def test_locate_fasta_files_errors(tmp_path: pathlib.Path) -> None:
    """
    Test _locate_fasta_files error handling for non-existent and empty dirs.

    Ensures that the function returns an empty list of files and an error
    message indicating that no FASTA files were found when the directory
    does not exist or is empty.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """

    # Directory does not exist
    errors: list[str] = []
    files: list[str]
    files, errors = _locate_fasta_files(
        directory=str(tmp_path / "nope"),
        errors=errors
    )
    assert files == []
    assert errors and "No FASTA files found" in errors[0]

    # Directory exists but no FASTA files
    empty_dir: pathlib.Path = tmp_path / "empty"
    empty_dir.mkdir()

    errors = []
    files, errors = _locate_fasta_files(
        directory=str(empty_dir),
        errors=errors
    )
    assert files == []
    assert errors and "No FASTA files found" in errors[0]


def test_blast_mode_assignment() -> None:
    """
    Test that the CombinedSubunits class assigns the correct blast_mode
    attribute based on the input arguments.

    This test covers the assignment for 'blastn', 'tblastx', 'blastx',
    and 'blastn+tblastx' modes.

    :return: None
    """
    # Patch make_blast_db to avoid running real BLAST
    with patch.object(
        CombinedSubunits, "make_blast_db", return_value="dummy_db"
    ):
        args = types.SimpleNamespace(
            allele_path=".",
            report_path=".",
            query_path=".",
            blast_mode="blastn",
            preliminary=False,
            num_alignments=1,
            cutoff=99.0,
            split_aa_db_dir=None,
            verbosity="info",
        )

        # Test blastn mode
        cs = CombinedSubunits(args)
        assert cs.blast_mode == "blastn"

        # Test tblastx mode
        args.blast_mode = "tblastx"
        cs = CombinedSubunits(args)
        assert cs.blast_mode == "tblastx"

        # Test blastx mode
        args.blast_mode = "blastx"
        cs = CombinedSubunits(args)
        assert cs.blast_mode == "blastx"

        # Test blastn+tblastx mode
        args.blast_mode = "blastn+tblastx"
        cs = CombinedSubunits(args)
        assert cs.blast_mode == ["blastn", "tblastx"]


def test_report_folder_creation_oserror(monkeypatch) -> None:
    """
    Test CombinedSubunits initialization when report folder creation fails.

    This test patches os.makedirs to raise OSError and checks that the
    error_print function is called with the appropriate error message.

    :param monkeypatch: Pytest monkeypatch fixture for patching functions
    :type monkeypatch: Any
    :return: None
    """
    # Patch os.path.isdir to always return True so we skip the error
    # branches above
    monkeypatch.setattr("os.path.isdir", lambda path: True)

    # Patch _locate_fasta_files to return dummy data
    monkeypatch.setattr(
        "allele_tools.stec_combined_subunits.CombinedSubunits."
        "_locate_fasta_files",
        lambda self, directory, errors: (["dummy.fasta"], errors),
    )

    # Patch os.makedirs to raise OSError
    monkeypatch.setattr(
        "os.makedirs",
        lambda *a, **k: (_ for _ in ()).throw(OSError("fail")),
    )

    # Patch make_blast_db to avoid running BLAST
    monkeypatch.setattr(
        "allele_tools.stec_combined_subunits.CombinedSubunits.make_blast_db",
        lambda *a, **k: "dummy_db",
    )

    # Patch error_print to capture errors
    captured: dict[str, list[str]] = {}

    def fake_error_print(errors: list[str]) -> None:
        """
        Capture errors passed to error_print.

        :param errors: List of error messages
        :type errors: list[str]
        :return: None
        """
        captured["errors"] = errors

    monkeypatch.setattr(
        "allele_tools.stec_combined_subunits.error_print", fake_error_print
    )

    args = types.SimpleNamespace(
        allele_path=".",
        report_path=".",
        query_path=".",
        blast_mode="blastn",
        preliminary=False,
        num_alignments=1,
        cutoff=99.0,
        split_aa_db_dir=None,
        verbosity="info",
    )

    # Attempt to initialize CombinedSubunits, expecting error_print to be
    # called due to OSError from os.makedirs
    CombinedSubunits(args)

    # Check that the captured errors contain the expected message
    assert any(
        "Error creating report folder" in e for e in captured.get("errors", [])
    )


def test_process_dual_molecule() -> None:
    """
    Test CombinedSubunits.process() with dual molecule mode.

    Ensures that when blast_mode is 'blastn+tblastx', the process method
    calls run twice (once for each molecule) and calls _combined_report
    once.

    :return: None
    """
    # Patch methods to track calls and avoid side effects
    with (
        patch.object(
            CombinedSubunits, "run"
        ) as mock_run,
        patch.object(
            CombinedSubunits, "_combined_report"
        ) as mock_report,
        patch.object(
            CombinedSubunits, "_prep_query_files", return_value={}
        ),
        patch.object(
            CombinedSubunits, "make_blast_db", return_value="dummy_db"
        ),
        patch.object(
            CombinedSubunits, "_locate_fasta_files",
            return_value=(["dummy.fasta"], [])
        ),
        patch("os.path.isdir", return_value=True),
        patch("os.makedirs"),
    ):
        args: types.SimpleNamespace = types.SimpleNamespace(
            allele_path=".",
            report_path=".",
            query_path=".",
            blast_mode="blastn+tblastx",
            preliminary=False,
            num_alignments=1,
            cutoff=99.0,
            split_aa_db_dir=None,
            verbosity="info",
        )

        # Initialize CombinedSubunits with patched dependencies
        cs: CombinedSubunits = CombinedSubunits(args)

        # Call the process method under test
        cs.process()

        # Should call run twice (once for each molecule)
        assert mock_run.call_count == 2

        # Should call _combined_report once
        assert mock_report.call_count == 1


def test_run_blastx_sets_blastx_true():
    """
    Test that CombinedSubunits.run sets blastx argument correctly when
    calling _parseable_blast_outputs for both 'blastx' and 'blastn' modes.
    """
    with (
        patch.object(
            CombinedSubunits, "_parseable_blast_outputs"
        ) as mock_parse,
        patch.object(
            CombinedSubunits, "_find_best_blast_hits", return_value={}
        ),
    ):
        args = types.SimpleNamespace(
            allele_path=".",
            report_path=".",
            query_path=".",
            blast_mode="blastx",
            preliminary=False,
            num_alignments=1,
            cutoff=99.0,
            split_aa_db_dir=None,
            verbosity="info",
        )
        # Patch __init__ to not run real code
        with patch.object(
            CombinedSubunits, "__init__", lambda self, args: None
        ):
            cs = CombinedSubunits(args)
            cs.blast_db_path = "dummy_db"
            cs.cpus = 1
            cs.cutoff = 99.0
            cs.fieldnames = []
            cs.extended_fieldnames = []
            cs.metadata = {
                "dummy": {
                    "aa_best_hits": {"stx1": [], "stx2": []},
                    "path": ".",
                    "file": "dummy_query.fasta",
                }
            }
            cs.molecule = "aa"
            cs.num_alignments = 1
            cs.report_path = "."
            cs.preliminary = False
            cs.outfmt = "6"
            cs.allele_file = "dummy_allele_file"
            cs.split_aa_dbs = {}

            # Run for 'blastx' mode and check blastx=True
            cs.run(blast_mode="blastx", molecule="aa")
            assert mock_parse.call_args[1]["blastx"] is True

            # Run for 'blastn' mode and check blastx=False
            cs.run(blast_mode="blastn", molecule="nt")
            assert mock_parse.call_args[1]["blastx"] is False

            # Delete the dummy allele file
            if os.path.exists("dummy_allele_file"):
                os.remove("dummy_allele_file")


def test_run_logs_novel_hits(caplog) -> None:
    """
    Test that CombinedSubunits.run logs best BLAST hits and novel hits.

    This test patches all downstream methods to avoid side effects and
    checks that the appropriate log messages are generated when novel
    hits are present in the metadata.

    :param caplog: Pytest fixture for capturing log messages
    :type caplog: Any
    :return: None
    """
    # Patch all downstream methods to avoid side effects
    with (
        patch.object(CombinedSubunits, "_run_blast"),
        patch.object(CombinedSubunits, "_parseable_blast_outputs"),
        patch.object(
            CombinedSubunits, "_find_best_blast_hits"
        ) as mock_find_best,
        patch.object(CombinedSubunits, "_write_novel_alleles_to_db"),
        patch.object(CombinedSubunits, "_export_novel_alleles"),
        patch.object(CombinedSubunits, "_find_aa_hits_novel_alleles"),
        patch.object(CombinedSubunits, "_add_headers_novel_aa_reports"),
        patch.object(CombinedSubunits, "_extract_novel_aa_best_hits"),
        patch.object(
            CombinedSubunits, "_assign_aa_allele_profiles_and_operons"
        ),
        patch.object(CombinedSubunits, "_export_aa_novel_operons"),
        patch.object(CombinedSubunits, "_export_aa_novel_alleles"),
        patch.object(CombinedSubunits, "_nt_report"),
        patch.object(CombinedSubunits, "_aa_report"),
    ):
        args = types.SimpleNamespace(
            allele_path=".",
            report_path=".",
            query_path=".",
            blast_mode="blastx",
            preliminary=False,
            num_alignments=1,
            cutoff=99.0,
            split_aa_db_dir=None,
            verbosity="info",
        )

        # Patch __init__ to not run real code
        with patch.object(
            CombinedSubunits, "__init__", lambda self, args: None
        ):
            cs = CombinedSubunits(args)
            cs.blast_db_path = "dummy_db"
            cs.cpus = 1
            cs.cutoff = 99.0
            cs.fieldnames = []
            cs.extended_fieldnames = []
            cs.metadata = {
                "strain1": {
                    "aa_novel_hits": ["novel1", "novel2"],
                    "aa_best_hits": {
                        "stx1": [("hit1", 99.5)],
                        "stx2": [],
                    },
                }
            }
            cs.molecule = "aa"
            cs.num_alignments = 1
            cs.report_path = "."
            cs.preliminary = False
            cs.outfmt = "6"
            cs.allele_file = "dummy_allele_file"
            cs.split_aa_dbs = {}

            # Simulate _find_best_blast_hits returning the same metadata
            mock_find_best.return_value = cs.metadata

            # Capture log messages at DEBUG level
            with caplog.at_level("DEBUG"):
                cs.run(blast_mode="blastx", molecule="aa")

                # Check that best BLAST hits are logged
                assert any(
                    "Best BLAST hits for strain1" in m
                    for m in caplog.messages
                )

                # Check that novel hits are logged
                assert any(
                    "novel1" in m or "novel2" in m
                    for m in caplog.messages
                )

                # Delete the dummy allele file
                if os.path.exists("dummy_allele_file"):
                    os.remove("dummy_allele_file")


def test_run_blast_sets_output_and_cmd() -> None:
    """
    Test that _run_blast sets the output and command for each strain.

    Ensures that _run_blast_cmd is called for each strain in the metadata,
    and that the output path and command are set correctly in the result.

    :return: None
    """
    # Prepare dummy metadata for two strains
    metadata: dict[str, dict[str, str]] = {
        "strainA": {
            "path": "/tmp/strainA",
            "file": "/tmp/strainA/query.fa"
        },
        "strainB": {
            "path": "/tmp/strainB",
            "file": "/tmp/strainB/query.fa"
        },
    }

    # Patch _run_blast_cmd to return a dummy command string
    with patch.object(
        CombinedSubunits,
        "_run_blast_cmd",
        return_value="dummy_cmd"
    ) as mock_cmd:

        result: dict[str, dict[str, str]] = _run_blast(
            allele_file="dummy_db",
            blast_mode="blastn",
            cpus=2,
            metadata=metadata,
            molecule="nt",
            outfmt="dummy_fmt",
        )

        # Check that _run_blast_cmd was called for each strain
        assert mock_cmd.call_count == 2

        for strain, info in result.items():
            # Check that the output path is set correctly
            expected_output: str = os.path.join(
                info["path"],
                f"{strain}_nt_blast.tsv"
            )
            assert info["nt_blast_output"] == expected_output

            # Check that the command is set
            assert info["nt_blast_cmd"] == "dummy_cmd"


def test_parseable_blast_outputs_calls_add_headers() -> None:
    """
    Test that _parseable_blast_outputs calls _add_headers for each strain.

    Ensures that _add_headers is called once for each strain in the
    metadata, and that the correct keyword arguments are passed.

    :return: None
    """
    # Prepare dummy metadata for two strains
    metadata: dict[str, dict[str, str]] = {
        "strainA": {"nt_blast_output": "/tmp/strainA/nt_blast.tsv"},
        "strainB": {"nt_blast_output": "/tmp/strainB/nt_blast.tsv"},
    }

    # Patch _add_headers to track calls
    with patch.object(
        CombinedSubunits, "_add_headers"
    ) as mock_add_headers:
        _parseable_blast_outputs(
            blastx=False,
            cutoff=99.0,
            fieldnames=["query_id", "subject_id"],
            extended_fieldnames=[
                "query_id", "subject_id", "percent_match"
            ],
            metadata=metadata,
            molecule="nt",
        )

        # Should call _add_headers once for each strain
        assert mock_add_headers.call_count == 2

        # Check that each call has the expected keyword arguments
        for call in mock_add_headers.call_args_list:
            kwargs = call.kwargs
            assert "blast_output" in kwargs
            assert "fieldnames" in kwargs
            assert "extended_fieldnames" in kwargs
            assert "cutoff" in kwargs
            assert "molecule" in kwargs


def write_blast_output(content: str) -> str:
    """
    Write the provided BLAST output content to a temporary file.

    This function creates a temporary file with a .tsv suffix, writes the
    given content to it, and returns the file path.

    :param content: The BLAST output content to write to the file.
    :type content: str
    :return: The path to the temporary BLAST output file.
    :rtype: str
    """
    # Create a temporary file with a .tsv suffix
    fd, path = tempfile.mkstemp(suffix=".tsv")

    # Write the content to the file with UTF-8 encoding
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        f.write(content)

    return path


def test_add_headers_subject_length_div3() -> None:
    """
    Test the _add_headers function for subject_length division by 3.

    This test covers the branch where subject_length is divided by 3,
    which is relevant for amino acid BLAST output. It ensures that the
    percent_match column is added and the header is present in the output.

    :return: None
    """
    # Prepare BLAST output content with subject_length for aa molecule
    content = (
        "q1\ts1\t30\t0\t0\t0\t0\t100\t90\t90\t1\t90\t1\t90\tACGT\tACGT\n"
    )
    path: str = write_blast_output(content)

    # Call the method under test
    _add_headers(
        blast_output=path,
        blastx=False,
        cutoff=10.0,
        extended_fieldnames=[
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
            "percent_match",
        ],
        fieldnames=[
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
        ],
        molecule="aa",
    )

    # Read the output file and check for correct header
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()

    # The header should contain 'percent_match'
    assert "percent_match" in lines[0]

    # Clean up the temporary file
    os.remove(path)


def test_add_headers_zerodivisionerror() -> None:
    """
    Test the _add_headers function for ZeroDivisionError handling.

    This test triggers a ZeroDivisionError by providing a row with
    subject_length set to zero. It ensures that the function skips the
    problematic row and only writes the header to the output file.

    :return: None
    """
    # Triggers ZeroDivisionError: subject_length is zero
    content: str = (
        "q1\ts1\t30\t0\t0\t0\tACGT\n"  # subject_length=0
    )
    path: str = write_blast_output(content)

    _add_headers(
        blast_output=path,
        blastx=False,
        cutoff=10.0,
        extended_fieldnames=[
            "query_id",
            "subject_id",
            "identical",
            "mismatches",
            "gaps",
            "subject_length",
            "query_sequence",
            "percent_match",
        ],
        fieldnames=[
            "query_id",
            "subject_id",
            "identical",
            "mismatches",
            "gaps",
            "subject_length",
            "query_sequence",
        ],
        molecule="nt",
    )

    # Only header should be present, row is skipped due to ZeroDivisionError
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()

    assert len(lines) == 1
    assert "percent_match" in lines[0]

    # Clean up the temporary file
    os.remove(path)


def test_add_headers_percent_match_written() -> None:
    """
    Test that the _add_headers function writes the percent_match column.

    This test covers the code lines where percent_identity is calculated
    and written to the row as 'percent_match'. It ensures that the
    percent_match column is present in the output and that the value is
    formatted to two decimal places.

    :return: None
    """
    # Prepare BLAST output content with subject_length for nt molecule
    content: str = (
        "q1\ts1\t30\t0\t0\t30\tACGT\n"
    )
    path: str = write_blast_output(content)

    # Call the method under test
    _add_headers(
        blast_output=path,
        blastx=False,
        cutoff=10.0,
        extended_fieldnames=[
            "query_id",
            "subject_id",
            "identical",
            "mismatches",
            "gaps",
            "subject_length",
            "query_sequence",
            "percent_match",
        ],
        fieldnames=[
            "query_id",
            "subject_id",
            "identical",
            "mismatches",
            "gaps",
            "subject_length",
            "query_sequence",
        ],
        molecule="nt",
    )

    # Read the output file and check for correct header and percent_match
    with open(path, encoding="utf-8") as f:
        lines: list[str] = f.readlines()

    # The header should contain 'percent_match'
    assert "percent_match" in lines[0]

    # The data row should end with the formatted percent identity
    assert lines[1].strip().endswith("100.00")

    # Clean up the temporary file
    os.remove(path)


def test_add_headers_write_updated_report() -> None:
    """
    Test the _add_headers function for writing the updated report.

    This test covers the code line where updated_report.write(...) is called,
    ensuring that the header and data row are written to the output file.

    :return: None
    """
    # Prepare BLAST output content with subject_length for nt molecule
    content: str = (
        "query_id\tsubject_id\tidentical\tmismatches\tgaps\t"
        "subject_length\tquery_sequence\n"
        "q1\ts1\t30\t0\t0\t30\tACGT\n"
    )
    path: str = write_blast_output(content)

    # Call the method under test
    _add_headers(
        blast_output=path,
        blastx=False,
        cutoff=10.0,
        extended_fieldnames=[
            "query_id",
            "subject_id",
            "identical",
            "mismatches",
            "gaps",
            "subject_length",
            "query_sequence",
            "percent_match",
        ],
        fieldnames=[
            "query_id",
            "subject_id",
            "identical",
            "mismatches",
            "gaps",
            "subject_length",
            "query_sequence",
        ],
        molecule="nt",
    )

    # Read the output file and check for header and data row
    with open(path, encoding="utf-8") as f:
        lines: list[str] = f.readlines()

    # Should have header and one data row
    assert len(lines) == 2

    # Clean up the temporary file
    os.remove(path)


def test_find_best_blast_hits_blastx_true() -> None:
    """
    Test the _find_best_blast_hits function when blastx=True.

    Ensures that the function correctly calls _extract_best_hits_subunits
    and updates the metadata with the expected keys and values.

    :return: None
    """
    # Prepare dummy metadata for testing
    metadata: dict = {"strainA": {"aa_blast_output": "dummy_path"}}

    # Patch _extract_best_hits_subunits to return dummy values
    with patch.object(
        CombinedSubunits,
        "_extract_best_hits_subunits",
        return_value=(
            {"stx1": [("allele1", 100.0)], "stx2": []},
            {"stx1": {"allele1": "ATGC"}, "stx2": {}},
        ),
    ):
        result: dict = _find_best_blast_hits(
            blastx=True,
            cutoff=99.0,
            metadata=metadata,
            molecule="aa",
        )

        # Check that the expected keys and values are present in the result
        assert "aa_best_hits" in result["strainA"]
        assert "aa_query_sequences" in result["strainA"]
        assert (
            result["strainA"]["aa_best_hits"]["stx1"][0][0] == "allele1"
        )
        assert (
            result["strainA"]["aa_query_sequences"]["stx1"]["allele1"]
            == "ATGC"
        )


def test_find_best_blast_hits_blastx_false() -> None:
    """
    Test the _find_best_blast_hits function when blastx is False.

    Ensures that the function correctly calls _extract_best_hits and updates
    the metadata with the expected keys and values for nucleotide BLAST.

    :return: None
    """
    # Prepare dummy metadata for testing
    metadata: dict = {"strainB": {"nt_blast_output": "dummy_path"}}

    # Patch _extract_best_hits to return dummy values
    with patch.object(
        CombinedSubunits,
        "_extract_best_hits",
        return_value=(
            {"stx1": [("allele2", 99.5)], "stx2": []},
            {"stx1": [("novel1", 100)], "stx2": []},
            {"stx1": {"allele2": "ATGC"}, "stx2": {}},
            {"stx1": {"allele2": "q1"}, "stx2": {}},
            {},
            {}
        ),
    ):
        result: dict = _find_best_blast_hits(
            blastx=False,
            cutoff=99.0,
            metadata=metadata,
            molecule="nt",
        )

        # Check that the expected keys and values are present in the result
        assert "nt_best_hits" in result["strainB"]
        assert "nt_novel_hits" in result["strainB"]
        assert "nt_query_sequences" in result["strainB"]
        assert "nt_query_ids" in result["strainB"]
        assert result["strainB"]["nt_best_hits"]["stx1"][0][0] == "allele2"
        assert result["strainB"]["nt_novel_hits"]["stx1"][0][0] == "novel1"
        assert (
            result["strainB"]["nt_query_sequences"]["stx1"]["allele2"]
            == "ATGC"
        )
        assert (
            result["strainB"]["nt_query_ids"]["stx1"]["allele2"]
            == "q1"
        )


def test_extract_best_hits_novel_branch(tmp_path: pathlib.Path) -> None:
    """
    Test the _extract_best_hits function for the novel branch.

    This test prepares a BLAST output file with percent_match < 100,
    percent_match >= cutoff, and a query_sequence present. It patches
    _register_novel_allele to return a predictable name, then checks
    that the novel branch is taken and the correct values are returned.
    """
    # Prepare a BLAST output file with percent_match < 100, >= cutoff,
    # and query_sequence present
    blast_content: str = (
        "query_id\tsubject_id\tpercent_match\tquery_sequence\n"
        "q1\tstx1_allele1\t99.5\tATGCATGC\n"
    )
    blast_file: pathlib.Path = tmp_path / "blast.tsv"
    blast_file.write_text(blast_content, encoding="utf-8")

    # Patch _register_novel_allele to return a predictable name
    with patch.object(
        CombinedSubunits,
        "_register_novel_allele",
        staticmethod(lambda sequence, stx_type: "novel_stx1_1"),
    ):
        # Call the method under test
        best_hits, novel_hits, query_sequences, query_ids, _, __ = \
            _extract_best_hits(
                blast_output=str(blast_file),
                cutoff=99.0,
                molecule="nt"
            )

    # Check that the novel branch was taken and values are as expected
    assert best_hits["stx1"] == [("stx1_allele1", 99.5)]
    assert novel_hits["stx1"] == [("novel_stx1_1", 100)]
    assert query_sequences["stx1"]["stx1_allele1"] == "ATGCATGC"
    assert query_ids["stx1"]["stx1_allele1"] == "q1"
    assert query_ids["stx1"]["novel_stx1_1"] == "q1"


def test_extract_best_hits_gene_subunits_tie_on_length(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _extract_best_hits_gene_subunits when two hits have the same
    percent_identity and subject_length.

    Ensures that both hits are included in the best_hits dictionary for
    the gene subunit, and that query_sequences and query_ids are set
    correctly for both.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare a BLAST output file with two rows, same percent_identity and
    # subject_length
    blast_content: str = (
        "query_id\tsubject_id\tpercent_match\tquery_sequence"
        "\tsubject_length\n"
        "q1\tStx1A_1\t99.0\tATGCATGC\t100\n"
        "q2\tStx1A_2\t99.0\tGGCCTTAA\t100\n"
    )
    blast_file: pathlib.Path = tmp_path / "blast.tsv"

    # Write the BLAST output file with UTF-8 encoding
    blast_file.write_text(blast_content, encoding="utf-8")

    # Call the method under test
    best_hits: dict
    query_sequences: dict
    query_ids: dict
    best_hits, query_sequences, query_ids = _extract_best_hits_gene_subunits(
        blast_output=str(blast_file),
        gene_subunit="Stx1A"
    )

    # Both hits should be present in best_hits[gene_subunit]
    assert len(best_hits["Stx1A"]) == 2
    assert ("Stx1A_1", 99.0) in best_hits["Stx1A"]
    assert ("Stx1A_2", 99.0) in best_hits["Stx1A"]

    # Both query sequences should be present
    assert query_sequences["Stx1A"]["Stx1A_1"] == "ATGCATGC"
    assert query_sequences["Stx1A"]["Stx1A_2"] == "GGCCTTAA"

    # Both query_ids should be present
    assert query_ids["Stx1A"]["Stx1A_1"] == "q1"
    assert query_ids["Stx1A"]["Stx1A_2"] == "q2"


def test_extract_best_hits_subunits_skips_if_no_stx_type(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _extract_best_hits_subunits skips rows with unknown stx type.

    Ensures that if the subject_id does not contain 'stx1' or 'stx2',
    the function returns empty best_hits and query_sequences.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare a BLAST output file with a subject_id that does not contain
    # stx1 or stx2
    blast_content: str = (
        "subject_id\tquery_sequence\tidentical\talignment_length\t"
        "mismatches\tgaps\tsubject_length\n"
        "unknown_allele_A\tATGC\t10\t10\t0\t0\t10\n"
        "unknown_allele_B\tGCTA\t10\t10\t0\t0\t10\n"
    )
    blast_file: pathlib.Path = tmp_path / "blast.tsv"

    # Write the BLAST output file with UTF-8 encoding
    blast_file.write_text(blast_content, encoding="utf-8")

    # Call the method under test
    best_hits: dict[str, list]
    query_sequences: dict[str, dict]
    best_hits, query_sequences = _extract_best_hits_subunits(
        blast_output=str(blast_file), cutoff=90.0
    )

    # Since stx_type is None for both rows, best_hits and query_sequences
    # should be empty
    assert best_hits == {"stx1": [], "stx2": []}
    assert query_sequences == {"stx1": {}, "stx2": {}}


def test_extract_best_hits_subunits_skips_zero_subject_length(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _extract_best_hits_subunits skips rows with subject_length == 0.

    Ensures that if both subunits have subject_length 0, the function
    returns empty best_hits and query_sequences.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare a BLAST output file with both subunits having subject_length 0
    blast_content: str = (
        "subject_id\tquery_sequence\tidentical\talignment_length\t"
        "mismatches\tgaps\tsubject_length\n"
        "stx1c_17_2_A\tATGC\t10\t10\t0\t0\t0\n"
        "stx1c_17_2_B\tGCTA\t10\t10\t0\t0\t0\n"
    )
    blast_file: pathlib.Path = tmp_path / "blast.tsv"

    # Write the BLAST output file with UTF-8 encoding
    blast_file.write_text(blast_content, encoding="utf-8")

    # Call the method under test
    best_hits: dict[str, list]
    query_sequences: dict[str, dict]
    best_hits, query_sequences = _extract_best_hits_subunits(
        blast_output=str(blast_file), cutoff=90.0
    )

    # Since total_subject_length == 0, best_hits and query_sequences
    # should be empty
    assert best_hits == {"stx1": [], "stx2": []}
    assert query_sequences == {"stx1": {}, "stx2": {}}


def test_extract_best_hits_subunits_below_cutoff(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _extract_best_hits_subunits when percent_identity < cutoff.

    This test prepares a BLAST output file with both subunits present,
    but the percent_identity is below the cutoff. It ensures that the
    function returns empty best_hits and query_sequences.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare a BLAST output file with both subunits, but percent_identity
    # < cutoff. identical = 40, subject_length_a = 100, subject_length_b = 100,
    # so percent_identity = 40.0.
    blast_content: str = (
        "subject_id\tquery_sequence\tidentical\talignment_length\t"
        "mismatches\tgaps\tsubject_length\n"
        "stx1c_17_2_A\tATGC\t40\t50\t10\t0\t100\n"
        "stx1c_17_2_B\tGCTA\t40\t50\t10\t0\t100\n"
    )
    blast_file: pathlib.Path = tmp_path / "blast.tsv"

    # Write the BLAST output file with UTF-8 encoding
    blast_file.write_text(blast_content, encoding="utf-8")

    # Call the method under test with cutoff above percent_identity (40.0)
    best_hits: dict[str, list]
    query_sequences: dict[str, dict]
    best_hits, query_sequences = _extract_best_hits_subunits(
        blast_output=str(blast_file),
        cutoff=90.0,
    )

    # Since percent_identity < cutoff, best_hits and query_sequences
    # should be empty
    assert best_hits == {"stx1": [], "stx2": []}
    assert query_sequences == {"stx1": {}, "stx2": {}}


def test_extract_best_hits_subunits_perfect_match(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _extract_best_hits_subunits for perfect match scenario.

    Ensures that when both subunits are present with 100% identity,
    no mismatches or gaps, the function returns a perfect match.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Both subunits, 100% identity, no mismatches/gaps (perfect match)
    blast_content: str = (
        "subject_id\tquery_sequence\tidentical\talignment_length\t"
        "mismatches\tgaps\tsubject_length\n"
        "stx1c_17_2_A\tATGC\t50\t50\t0\t0\t50\n"
        "stx1c_17_2_B\tGCTA\t50\t50\t0\t0\t50\n"
    )
    blast_file: pathlib.Path = tmp_path / "blast.tsv"

    # Write the BLAST output file with UTF-8 encoding
    blast_file.write_text(blast_content, encoding="utf-8")

    # Call the method under test
    best_hits: dict[str, list]
    query_sequences: dict[str, dict]
    best_hits, query_sequences = _extract_best_hits_subunits(
        blast_output=str(blast_file), cutoff=99.0
    )

    # Should prefer perfect match at 100% identity
    assert best_hits["stx1"] == [("stx1c_17_2", 100.0)]
    assert query_sequences["stx1"]["stx1c_17_2_A"] == "ATGC"
    assert query_sequences["stx1"]["stx1c_17_2_B"] == "GCTA"
    assert query_sequences["stx1"]["stx1c_17_2"] == "ATGC" + "GCTA"


def test_extract_best_hits_subunits_equal_percent_identity(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _extract_best_hits_subunits with two alleles having the same
    percent identity (not perfect match).

    Ensures that both alleles are present in the best_hits dictionary
    for stx1 with the same percent identity.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare BLAST output content with two alleles, both with same percent
    # identity (90.0), not perfect match
    blast_content: str = (
        "subject_id\tquery_sequence\tidentical\talignment_length\t"
        "mismatches\tgaps\tsubject_length\n"
        "stx1c_17_2_A\tATGC\t45\t50\t1\t0\t50\n"
        "stx1c_17_2_B\tGCTA\t45\t50\t1\t0\t50\n"
        "stx1c_18_2_A\tTTTT\t45\t50\t1\t0\t50\n"
        "stx1c_18_2_B\tCCCC\t45\t50\t1\t0\t50\n"
    )
    blast_file: pathlib.Path = tmp_path / "blast.tsv"

    # Write the BLAST output file with UTF-8 encoding
    blast_file.write_text(blast_content, encoding="utf-8")

    # Call the method under test
    best_hits, _ = _extract_best_hits_subunits(
        blast_output=str(blast_file),
        cutoff=80.0,
    )

    # Both alleles should be present with the same percent identity (90.0)
    assert ("stx1c_17_2", 90.0) in best_hits["stx1"]
    assert ("stx1c_18_2", 90.0) in best_hits["stx1"]
    assert len(best_hits["stx1"]) == 2


def setup_registry_with_novel() -> None:
    """
    Set up the _novel_registry with a single novel allele entry for testing.

    This function clears the global _novel_registry and inserts a test
    entry with a fixed hash and sequence.

    :return: None
    """
    # Clear the registry before adding the test entry
    _novel_registry.clear()

    # Add a test novel allele entry with a fixed hash and sequence
    _novel_registry["hash1"] = {
        "name": "novel_stx1_1",
        "sequence": "ATGC",
    }


def test_write_novel_alleles_to_db_adds_newline(
    tmp_path: pathlib.Path
) -> None:
    """
    Test that _write_novel_alleles_to_db adds a newline if missing.

    This test creates a FASTA file that does not end with a newline,
    runs the function, and checks that the resulting file ends with a
    newline character.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Create a FASTA file that does NOT end with a newline
    fasta_path: pathlib.Path = tmp_path / "db.fasta"
    with open(fasta_path, "wb") as f:
        f.write(b">existing\nATGC")  # No newline at end

    # Set up the registry with a novel allele
    setup_registry_with_novel()
    report_path: pathlib.Path = tmp_path

    # Patch glob and os.remove to avoid deleting real files
    with patch(
        "allele_tools.stec_combined_subunits.glob", return_value=[]
    ):
        _write_novel_alleles_to_db(
            db_fasta_file=str(fasta_path),
            report_path=str(report_path)
        )

    # Check that the new file ends with a newline
    with open(fasta_path, "rb") as f:
        content: bytes = f.read()
        assert content.endswith(b"\n")


def test_write_novel_alleles_to_db_no_newline_needed(
    tmp_path: pathlib.Path
) -> None:
    """
    Test that _write_novel_alleles_to_db does not add a newline if one is
    already present at the end of the file.

    This test creates a FASTA file that ends with a newline, runs the
    function, and checks that the resulting file still ends with a newline
    character.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Create a FASTA file that DOES end with a newline
    fasta_path: pathlib.Path = tmp_path / "db.fasta"
    with open(fasta_path, "wb") as f:
        f.write(b">existing\nATGC\n")  # Has newline at end

    # Set up the registry with a novel allele
    setup_registry_with_novel()
    report_path: pathlib.Path = tmp_path

    # Patch glob to avoid deleting real files
    with patch(
        "allele_tools.stec_combined_subunits.glob", return_value=[]
    ):
        _write_novel_alleles_to_db(
            db_fasta_file=str(fasta_path),
            report_path=str(report_path)
        )

    # Should still end with a newline
    with open(fasta_path, "rb") as f:
        content: bytes = f.read()
        assert content.endswith(b"\n")


def test_write_novel_alleles_to_db_deletes_db_files(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test that _write_novel_alleles_to_db deletes associated BLAST db files.

    This test creates a dummy FASTA file and a dummy .n* file, sets up the
    registry with a novel allele, and checks that os.remove is called to
    delete the .n* file after updating the database.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Create a dummy db FASTA file with encoding specified
    fasta_path: pathlib.Path = tmp_path / "db.fasta"
    with open(fasta_path, "w", encoding="utf-8") as f:
        f.write(">existing\nATGC\n")

    # Create a dummy .n* file to simulate BLAST db index
    nfile: pathlib.Path = tmp_path / "db.nhr"
    nfile.write_text("dummy", encoding="utf-8")

    # Set up the registry with a novel allele
    setup_registry_with_novel()
    report_path: pathlib.Path = tmp_path

    # Patch os.remove to track calls for file deletion
    with patch(
        "allele_tools.stec_combined_subunits.os.remove"
    ) as mock_remove:
        _write_novel_alleles_to_db(
            db_fasta_file=str(fasta_path),
            report_path=str(report_path)
        )
        # Should attempt to remove the .n* file
        mock_remove.assert_any_call(str(nfile))


def test_export_novel_alleles_blastx_false(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _export_novel_alleles for blastx=False with a best hit above cutoff
    and <100% identity.

    Ensures that the file-writing method is called and the metadata is
    updated with the returned paths.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Simulate metadata with a best hit above cutoff and <100% identity
    metadata: dict = {
        "strain1": {
            "nt_best_hits": {"stx1": [("stx1c_17_2", 99.5)], "stx2": []},
            "nt_query_sequences": {
                "stx1": {"stx1c_17_2": "ATGC"},
                "stx2": {},
            },
            "path": str(tmp_path),
        }
    }

    # Patch the file-writing method to avoid actual file IO
    with patch.object(
        CombinedSubunits,
        "_write_novel_allele_files",
        return_value={"stx1c_17_2": "dummy_path"},
    ) as mock_write:
        result: dict = _export_novel_alleles(
            blastx=False,
            cutoff=99.0,
            metadata=metadata,
            molecule="nt",
            report_path=str(tmp_path),
        )

        # Should call the file-writing method
        mock_write.assert_called_once()

        # Should update the metadata with the returned paths
        assert "novel_nt_allele_paths" in result["strain1"]
        assert result["strain1"]["novel_nt_allele_paths"] == {
            "stx1c_17_2": "dummy_path"
        }


def test_export_novel_alleles_blastx_true(tmp_path: pathlib.Path) -> None:
    """
    Test _export_novel_alleles for blastx=True with a best hit above cutoff
    and <100% identity.

    Ensures that the subunit file-writing method is called and the metadata
    is updated with the returned paths.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Simulate metadata with a best hit above cutoff and <100% identity
    metadata: dict = {
        "strain2": {
            "aa_best_hits": {
                "stx2": [("stx2c_5_1", 99.2)],
                "stx1": [],
            },
            "aa_query_sequences": {
                "stx2": {
                    "stx2c_5_1": "ATGC",
                    "stx2c_5_1_A": "AAAA",
                    "stx2c_5_1_B": "CCCC",
                },
                "stx1": {},
            },
            "path": str(tmp_path),
        }
    }

    # Patch the subunit file-writing method to avoid actual file IO
    with patch.object(
        CombinedSubunits,
        "_write_novel_allele_subunit_files",
        return_value=(
            {
                "A": "a_path",
                "B": "b_path",
                "combined": "c_path",
            },
            {}
        ),
    ) as mock_write:
        result: dict = _export_novel_alleles(
            blastx=True,
            cutoff=99.0,
            metadata=metadata,
            molecule="aa",
            report_path=str(tmp_path),
        )

        # Should call the file-writing method once
        mock_write.assert_called_once()

        # Should update the metadata with the returned paths
        assert "novel_aa_allele_paths" in result["strain2"]
        assert result["strain2"]["novel_aa_allele_paths"] == {
            "A": "a_path",
            "B": "b_path",
            "combined": "c_path",
        }


def test_export_novel_alleles_below_cutoff(
    tmp_path: pathlib.Path
) -> None:
    """
    Test _export_novel_alleles for a best hit below the cutoff.

    Ensures that the file-writing method is not called when the best hit
    percent identity is below the cutoff, and that the output key is set
    to an empty dictionary.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Simulate metadata with a best hit below cutoff
    metadata: dict = {
        "strain3": {
            "nt_best_hits": {
                "stx1": [("stx1c_17_2", 95.0)],
                "stx2": [],
            },
            "nt_query_sequences": {
                "stx1": {"stx1c_17_2": "ATGC"},
                "stx2": {},
            },
            "path": str(tmp_path),
        }
    }

    # Patch the file-writing method to avoid actual file IO
    with patch.object(
        CombinedSubunits, "_write_novel_allele_files"
    ) as mock_write:
        result: dict = _export_novel_alleles(
            blastx=False,
            cutoff=99.0,
            metadata=metadata,
            molecule="nt",
            report_path=str(tmp_path),
        )

        # Should NOT call the file-writing method
        mock_write.assert_not_called()

        # Should still set the key, but it will be empty
        assert result["strain3"]["novel_nt_allele_paths"] == {}


def test_export_novel_alleles_exactly_100(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _export_novel_alleles for a best hit at exactly 100% identity.

    Ensures that the file-writing method is NOT called when the best hit
    percent identity is exactly 100, and that the output key is set to
    an empty dictionary.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Simulate metadata with a best hit at 100% identity
    metadata: dict = {
        "strain4": {
            "nt_best_hits": {
                "stx1": [("stx1c_17_2", 100.0)],
                "stx2": [],
            },
            "nt_query_sequences": {
                "stx1": {"stx1c_17_2": "ATGC"},
                "stx2": {},
            },
            "path": str(tmp_path),
        }
    }

    # Patch the file-writing method to avoid actual file IO
    with patch.object(
        CombinedSubunits, "_write_novel_allele_files"
    ) as mock_write:
        result: dict = _export_novel_alleles(
            blastx=False,
            cutoff=99.0,
            metadata=metadata,
            molecule="nt",
            report_path=str(tmp_path),
        )

        # Should NOT call the file-writing method
        mock_write.assert_not_called()

        # Should still set the key, but it will be empty
        assert result["strain4"]["novel_nt_allele_paths"] == {}


def test_export_novel_alleles_skips_empty_query_sequence(
    tmp_path: pathlib.Path
) -> None:
    """
    Test _export_novel_alleles skips writing files if query_sequence is empty.

    This test simulates metadata with a best hit above cutoff and less than
    100% identity, but with no query_sequence present for the allele. It
    ensures that the file-writing method is not called and that the output
    key is set to an empty dictionary.

    :param tmp_path: Temporary directory provided by pytest for file IO.
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Simulate metadata with a best hit above cutoff and <100% identity,
    # but no query_sequence for stx1c_17_2
    metadata: dict = {
        "strain1": {
            "nt_best_hits": {"stx1": [("stx1c_17_2", 99.5)], "stx2": []},
            "nt_query_sequences": {
                "stx1": {},  # No sequence for stx1c_17_2
                "stx2": {},
            },
            "path": str(tmp_path),
        }
    }

    # Patch the file-writing method to avoid actual file IO
    with patch.object(
        CombinedSubunits, "_write_novel_allele_files"
    ) as mock_write:
        result: dict = _export_novel_alleles(
            blastx=False,
            cutoff=99.0,
            metadata=metadata,
            molecule="nt",
            report_path=str(tmp_path),
        )

        # Should NOT call the file-writing method
        mock_write.assert_not_called()

        # Should still set the key, but it will be empty
        assert result["strain1"]["novel_nt_allele_paths"] == {}


def test_export_novel_alleles_skips_missing_subunit_sequences(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _export_novel_alleles skips writing files if a subunit sequence
    is missing.

    This test simulates metadata with a best hit above cutoff and less than
    100% identity, but with a missing subunit B sequence for the allele.
    It ensures that the file-writing method is not called and that the
    output key is set to an empty dictionary.

    :param tmp_path: Temporary directory provided by pytest for file IO.
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Simulate metadata with a best hit above cutoff and <100% identity,
    # but missing subunit B sequence for stx1c_17_2
    metadata: dict = {
        "strain1": {
            "aa_best_hits": {"stx1": [("stx1c_17_2", 99.5)], "stx2": []},
            "aa_query_sequences": {
                "stx1": {
                    "stx1c_17_2": "ATGC",  # combined
                    "stx1c_17_2_A": "AAAA",  # subunit A present
                    # "stx1c_17_2_B" missing   # subunit B missing
                },
                "stx2": {},
            },
            "path": str(tmp_path),
        }
    }

    # Patch the file-writing method to avoid actual file IO
    with patch.object(
        CombinedSubunits, "_write_novel_allele_subunit_files"
    ) as mock_write:
        result: dict = _export_novel_alleles(
            blastx=True,
            cutoff=99.0,
            metadata=metadata,
            molecule="aa",
            report_path=str(tmp_path),
        )

        # Should NOT call the file-writing method
        mock_write.assert_not_called()

        # Should still set the key, but it will be empty
        assert result["strain1"]["novel_aa_allele_paths"] == {}


def test_find_aa_hits_novel_alleles_basic(
    tmp_path: pathlib.Path
) -> None:
    """
    Test _find_aa_hits_novel_alleles with two alleles (stx1 and stx2).

    Ensures that the function calls _run_blast_cmd for each allele and
    subunit, and that the output dictionary contains the expected keys.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Simulate metadata with two alleles, one stx1 and one stx2
    metadata: dict = {
        "strainA": {
            "novel_nt_allele_paths": {
                "stx1c_17_2": "file1.fasta",
                "stx2c_5_1": "file2.fasta",
            },
            "path": str(tmp_path),
        }
    }
    split_aa_dbs: dict = {
        "stx1A": str(tmp_path / "stx1A.fasta"),
        "stx1B": str(tmp_path / "stx1B.fasta"),
        "stx2A": str(tmp_path / "stx2A.fasta"),
        "stx2B": str(tmp_path / "stx2B.fasta"),
    }

    # Patch _run_blast_cmd so it doesn't actually run BLAST
    with patch.object(
        CombinedSubunits,
        "_run_blast_cmd",
        return_value="blast_cmd"
    ) as mock_blast:
        result: dict = _find_aa_hits_novel_alleles(
            cpus=2,
            metadata=metadata,
            num_alignments=5,
            outfmt="6 std",
            split_aa_dbs=split_aa_dbs,
        )

        # Should call _run_blast_cmd 4 times (2 alleles x 2 subunits)
        assert mock_blast.call_count == 4

        # Check that the correct keys are set in the output
        outputs: dict = result["strainA"]["novel_aa_blastx_outputs"]
        assert "stx1c_17_2_A" in outputs
        assert "stx1c_17_2_B" in outputs
        assert "stx2c_5_1_A" in outputs
        assert "stx2c_5_1_B" in outputs

        # Each should have a list with the blast output path
        for key, value in outputs.items():
            assert isinstance(value, list)
            # The output file should end with the expected suffix
            assert (
                value[0].endswith(f"{key}_novel_{key[-5:]}_blastx.tsv")
                or value[0].endswith("_blastx.tsv")
            )


def test_find_aa_hits_novel_alleles_empty(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _find_aa_hits_novel_alleles with empty novel_nt_allele_paths.

    Ensures that if there are no alleles in the metadata, the function
    does not call _run_blast_cmd and still sets the output key.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Simulate metadata with no alleles
    metadata: dict = {
        "strainB": {
            "novel_nt_allele_paths": {},
            "path": str(tmp_path)
        }
    }
    split_aa_dbs: dict = {
        "stx1A": str(tmp_path / "stx1A.fasta"),
        "stx1B": str(tmp_path / "stx1B.fasta"),
        "stx2A": str(tmp_path / "stx2A.fasta"),
        "stx2B": str(tmp_path / "stx2B.fasta"),
    }

    # Patch _run_blast_cmd to ensure it is not called
    with patch.object(
        CombinedSubunits, "_run_blast_cmd"
    ) as mock_blast:
        result: dict = _find_aa_hits_novel_alleles(
            cpus=1,
            metadata=metadata,
            num_alignments=1,
            outfmt="6 std",
            split_aa_dbs=split_aa_dbs,
        )

        # Should not call _run_blast_cmd
        mock_blast.assert_not_called()

        # Should still set the output key
        assert "novel_aa_blastx_outputs" in result["strainB"]
        assert result["strainB"]["novel_aa_blastx_outputs"] == {}


def test_add_headers_novel_aa_reports_single_strain_single_output() -> None:
    """
    Test _add_headers_novel_aa_reports with a single strain and single output.

    Ensures that the _add_headers method is called for each output file
    in the metadata for a single strain.

    :return: None
    """
    metadata: dict = {
        "strain1": {
            "novel_aa_blastx_outputs": {
                "alleleA_A": ["outputA.tsv"],
                "alleleA_B": ["outputB.tsv"],
            }
        }
    }

    # Patch _add_headers to track calls for each output file
    with patch.object(CombinedSubunits, "_add_headers") as mock_add_headers:
        _add_headers_novel_aa_reports(
            blastx=True,
            extended_fieldnames=["f1", "f2"],
            fieldnames=["f1", "f2"],
            metadata=metadata,
            cutoff=99.0,
        )

        # Should be called for each output file
        assert mock_add_headers.call_count == 2

        # Collect the output files for which _add_headers was called
        called_files: list[str] = [
            call.kwargs["blast_output"]
            for call in mock_add_headers.call_args_list
        ]
        assert set(called_files) == {"outputA.tsv", "outputB.tsv"}


def test_add_headers_novel_aa_reports_multiple_strains_multiple_outputs(
) -> None:
    """
    Test _add_headers_novel_aa_reports with multiple strains and multiple
    output files per strain.

    Ensures that the _add_headers method is called for each output file
    in all strains' metadata.

    :return: None
    """
    metadata: dict = {
        "strain1": {
            "novel_aa_blastx_outputs": {
                "alleleA_A": ["outputA1.tsv", "outputA2.tsv"],
            }
        },
        "strain2": {
            "novel_aa_blastx_outputs": {
                "alleleB_A": ["outputB1.tsv"],
                "alleleB_B": ["outputB2.tsv", "outputB3.tsv"],
            }
        },
    }

    # Patch _add_headers to track calls for each output file
    with patch.object(CombinedSubunits, "_add_headers") as mock_add_headers:
        _add_headers_novel_aa_reports(
            blastx=False,
            extended_fieldnames=["f1", "f2"],
            fieldnames=["f1", "f2"],
            metadata=metadata,
            cutoff=98.0,
        )

        # Should be called for each output file in all strains
        assert mock_add_headers.call_count == 5

        # Collect the output files for which _add_headers was called
        called_files: list[str] = [
            call.kwargs["blast_output"]
            for call in mock_add_headers.call_args_list
        ]
        assert set(called_files) == {
            "outputA1.tsv",
            "outputA2.tsv",
            "outputB1.tsv",
            "outputB2.tsv",
            "outputB3.tsv",
        }


def test_add_headers_novel_aa_reports_no_outputs() -> None:
    """
    Test _add_headers_novel_aa_reports when there are no output files.

    Ensures that the _add_headers method is not called if the
    metadata for a strain contains an empty novel_aa_blastx_outputs
    dictionary.

    :return: None
    """
    metadata: dict = {"strain1": {"novel_aa_blastx_outputs": {}}}

    # Patch _add_headers to track calls
    with patch.object(
        CombinedSubunits, "_add_headers"
    ) as mock_add_headers:
        _add_headers_novel_aa_reports(
            blastx=False,
            extended_fieldnames=["f1", "f2"],
            fieldnames=["f1", "f2"],
            metadata=metadata,
            cutoff=97.0,
        )

        # Should not call _add_headers if there are no outputs
        mock_add_headers.assert_not_called()


def test_extract_novel_aa_best_hits_basic() -> None:
    """
    Test the _extract_novel_aa_best_hits function with basic input.

    Ensures that the function initializes and updates the expected keys
    in the result dictionary for each strain and allele.

    :return: None
    """
    metadata: dict = {
        "strain1": {
            "novel_aa_blastx_outputs": {
                "alleleA_A": ["outputA1.tsv", "outputA2.tsv"],
                "alleleA_B": ["outputB1.tsv"],
            }
        },
        "strain2": {
            "novel_aa_blastx_outputs": {
                "alleleB_A": ["outputC1.tsv"]
            }
        },
    }

    # Patch extract_gene_subunit and _extract_best_hits_gene_subunits
    with (
        patch.object(
            CombinedSubunits,
            "extract_gene_subunit",
            side_effect=lambda filepath: (
                "Stx1A" if "A" in filepath else "Stx1B"
            ),
        ),
        patch.object(
            CombinedSubunits,
            "_extract_best_hits_gene_subunits",
            return_value=(
                {"Stx1A": [("Stx1A_22", 100.0)]},
                {"Stx1A": {"Stx1A_22": "SEQ"}},
                {"Stx1A": {"Stx1A_22": "QID"}},
            ),
        ),
    ):
        result: dict = _extract_novel_aa_best_hits(metadata=metadata)

        # Check that the keys are initialized and updated for each strain
        for strain in ["strain1", "strain2"]:
            assert "novel_aa_best_hits" in result[strain]
            assert "novel_aa_query_sequences" in result[strain]
            assert "novel_aa_query_ids" in result[strain]

            # Should have the expected structure for each allele
            for allele_id in metadata[strain]["novel_aa_blastx_outputs"]:
                assert (
                    "Stx1A"
                    in result[strain]["novel_aa_best_hits"][allele_id]
                )
                assert (
                    "Stx1A"
                    in result[strain]["novel_aa_query_sequences"][allele_id]
                )
                assert (
                    "Stx1A"
                    in result[strain]["novel_aa_query_ids"][allele_id]
                )


def test_extract_novel_aa_best_hits_empty_metadata() -> None:
    """
    Test _extract_novel_aa_best_hits with empty metadata.

    Ensures that the function returns an empty dictionary when the input
    metadata is empty.

    :return: None
    """
    metadata: dict = {}

    # Call the method under test
    result: dict = _extract_novel_aa_best_hits(metadata=metadata)

    # The result should be an empty dictionary
    assert result == {}


def test_extract_novel_aa_best_hits_no_outputs() -> None:
    """
    Test _extract_novel_aa_best_hits when no outputs are present.

    Ensures that the function initializes the expected keys in the result
    dictionary for a strain with no novel_aa_blastx_outputs.

    :return: None
    """
    metadata: dict = {"strain1": {"novel_aa_blastx_outputs": {}}}

    # Call the method under test
    result: dict = _extract_novel_aa_best_hits(metadata=metadata)

    # Check that the expected keys are present and empty
    assert "novel_aa_best_hits" in result["strain1"]
    assert "novel_aa_query_sequences" in result["strain1"]
    assert "novel_aa_query_ids" in result["strain1"]
    assert result["strain1"]["novel_aa_best_hits"] == {}
    assert result["strain1"]["novel_aa_query_sequences"] == {}
    assert result["strain1"]["novel_aa_query_ids"] == {}


def test_export_aa_novel_operons_skips_if_no_profile_or_operon_seq(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _export_aa_novel_operons skips writing FASTA files if either the
    allele profile or operon sequence is missing.

    Ensures that the function does not attempt to write any FASTA files
    when either the allele profile is empty or the operon sequence is
    empty for a given strain.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Case 1: profile is empty for stx1, operon_sequence is empty for stx2
    metadata: dict = {
        "strain1": {
            "stx1_aa_allele_profile": "",
            "stx1_aa_operon_sequence": "ATGC",
            "stx2_aa_allele_profile": "22_38",
            "stx2_aa_operon_sequence": "",
        }
    }

    # Patch SeqIO.write to track if any FASTA files are written
    with patch("Bio.SeqIO.write") as mock_write:
        _export_aa_novel_operons(
            metadata=metadata,
            report_path=str(tmp_path)
        )

        # Should not write any FASTA files since required data is missing
        mock_write.assert_not_called()


def test_export_aa_novel_operons_skips_if_profile_split_fails(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _export_aa_novel_operons skips writing FASTA files if the allele
    profile cannot be split (ValueError).

    Ensures that the function does not attempt to write any FASTA files
    when the allele profile string cannot be split as expected.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Case 2: profile cannot be split (ValueError)
    metadata: dict = {
        "strain1": {
            "stx1_aa_allele_profile": "badprofile",  # No underscore
            "stx1_aa_operon_sequence": "ATGC",
            "stx2_aa_allele_profile": "anotherbad",  # No underscore
            "stx2_aa_operon_sequence": "GGCC",
        }
    }

    # Patch SeqIO.write to track if any FASTA files are written
    with patch("Bio.SeqIO.write") as mock_write:
        _export_aa_novel_operons(
            metadata=metadata,
            report_path=str(tmp_path)
        )

        # Should not write any FASTA files since profile split fails
        mock_write.assert_not_called()


def test_export_aa_novel_alleles_skips_if_no_best_hits(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _export_aa_novel_alleles skips if best_hits is empty.

    Ensures that the file-writing method is not called when the best_hits
    list is empty for all alleles and stx types.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare metadata with empty best_hits for all alleles and stx types
    metadata: dict = {
        "strain1": {
            "aa_query_sequences": {"stx1": {}, "stx2": {}},
            "novel_aa_best_hits": {
                "allele1": {"stx1": [], "stx2": []}
            },
            "path": str(tmp_path),
        }
    }

    # Patch the file-writing method to ensure it is not called
    with patch.object(
        CombinedSubunits, "_write_novel_allele_files"
    ) as mock_write:
        _export_aa_novel_alleles(
            cutoff=99.0,
            metadata=metadata,
            molecule="aa",
            report_path=str(tmp_path),
        )
        # Should not call the file-writing method
        mock_write.assert_not_called()


def test_export_aa_novel_alleles_skips_if_below_cutoff_or_100(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _export_aa_novel_alleles skips if percent_identity < cutoff or == 100.

    Ensures that the file-writing method is not called when the best hit
    percent identity is below the cutoff or exactly 100, and that the
    output key is set to an empty dictionary.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare metadata with one hit below cutoff and one at 100% identity
    metadata: dict = {
        "strain1": {
            "aa_query_sequences": {
                "stx1": {"alleleA": "ATGC"},
                "stx2": {},
            },
            "novel_aa_best_hits": {
                "allele1": {
                    "stx1": [("alleleA", 98.0)],
                    "stx2": [("alleleB", 100.0)],
                }
            },
            "path": str(tmp_path),
        }
    }

    # Patch the file-writing method to ensure it is not called
    with patch.object(
        CombinedSubunits, "_write_novel_allele_files"
    ) as mock_write:
        _export_aa_novel_alleles(
            cutoff=99.0,
            metadata=metadata,
            molecule="aa",
            report_path=str(tmp_path),
        )

        # Both hits should be skipped (one < cutoff, one == 100)
        mock_write.assert_not_called()


def test_export_aa_novel_alleles_skips_if_no_query_sequence(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test _export_aa_novel_alleles skips if query_sequence is empty.

    Ensures that the file-writing method is not called when the
    query_sequence is empty for the best hit, and that the output
    key is set to an empty dictionary.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare metadata with an empty query_sequence for the best hit
    metadata: dict = {
        "strain1": {
            "aa_query_sequences": {
                "stx1": {"alleleA": ""},
                "stx2": {},
            },
            "novel_aa_best_hits": {
                "allele1": {
                    "stx1": [("alleleA", 99.5)],
                    "stx2": [],
                }
            },
            "path": str(tmp_path),
        }
    }

    # Patch the file-writing method to ensure it is not called
    with patch.object(
        CombinedSubunits, "_write_novel_allele_files"
    ) as mock_write:
        _export_aa_novel_alleles(
            cutoff=99.0,
            metadata=metadata,
            molecule="aa",
            report_path=str(tmp_path),
        )

        # Should NOT call the file-writing method
        mock_write.assert_not_called()


def test_nt_report_novel_flag(tmp_path: pathlib.Path) -> None:
    """
    Test the _nt_report function when novel_aa_best_hits is not present.

    Ensures that when _has_novel_aa_best_hits returns False (novel=True),
    the report file contains the correct number of tabs for the novel case.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # No novel_aa_best_hits present, so novel=True
    metadata: dict = {
        "strain1": {
            "nt_best_hits": {"stx1": [], "stx2": []},
            "nt_novel_hits": {"stx1": [], "stx2": []},
        }
    }
    report_path: pathlib.Path = tmp_path

    # Patch _has_novel_aa_best_hits to return False (novel=True)
    with patch.object(
        CombinedSubunits, "_has_novel_aa_best_hits", return_value=False
    ):
        _nt_report(
            metadata=metadata,
            report_path=str(report_path)
        )

    # Check that the report file contains the correct number of tabs
    # for the novel=True case
    report_file: str = os.path.join(
        report_path, "stec_combined_nt_report.tsv"
    )
    with open(report_file, encoding="utf-8") as f:
        lines: list[str] = f.readlines()
        assert any(
            line.endswith("\t\t\t\t\t\t\n") for line in lines
        )


def test_nt_report_valueerror_in_allele_profile(
    tmp_path: pathlib.Path
) -> None:
    """
    Test _nt_report with an allele profile that triggers ValueError.

    This test provides an allele profile string that cannot be split
    (missing an underscore), which should trigger a ValueError in the
    report generation logic. The function should not raise an exception
    and should still write the report file.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Provide an allele profile that can't be split (triggers ValueError)
    metadata: dict = {
        "strain1": {
            "nt_best_hits": {"stx1": [("allele1", 99.5)], "stx2": []},
            "nt_novel_hits": {"stx1": [("novel_allele", 100)], "stx2": []},
            "nt_query_ids": {"stx1": {"allele1": "contig1"}, "stx2": {}},
            "stx1_aa_allele_profile": "badprofile",  # No underscore
            "novel_aa_best_hits": {},
        }
    }
    report_path: pathlib.Path = tmp_path

    # Patch _has_novel_aa_best_hits to return True
    with patch.object(
        CombinedSubunits, "_has_novel_aa_best_hits", return_value=True
    ):
        _nt_report(
            metadata=metadata,
            report_path=str(report_path)
        )

    # Should not raise, and should still write a report
    report_file: str = os.path.join(
        report_path, "stec_combined_nt_report.tsv"
    )
    assert os.path.isfile(report_file)


def test_nt_report_note_semicolon(
    tmp_path: pathlib.Path
) -> None:
    """
    Test that the _nt_report function appends a semicolon to the note
    if both conditions are met: a partial match and a novel subunit
    match.

    This test prepares metadata with a partial match (perc_ident < 90)
    and a novel subunit match, then checks that the note in the report
    contains a semicolon separating the two notes.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare metadata with a partial match and a novel subunit match
    metadata: dict = {
        "strain1": {
            "nt_best_hits": {"stx1": [("allele1", 80.0)], "stx2": []},
            "nt_novel_hits": {"stx1": [("novel_allele", 100)], "stx2": []},
            "nt_query_ids": {"stx1": {"allele1": "contig1"}, "stx2": {}},
            "stx1_aa_allele_profile": "22_38",
            "novel_aa_best_hits": {
                "allele1": {
                    "Stx1A": [("Stx1A_22", 80.0)],
                    "Stx1B": [("Stx1B_38", 99.0)],
                }
            },
        }
    }
    report_path: pathlib.Path = tmp_path

    # Patch _has_novel_aa_best_hits to return True to trigger the note logic
    with patch.object(
        CombinedSubunits, "_has_novel_aa_best_hits", return_value=True
    ):
        _nt_report(
            metadata=metadata,
            report_path=str(report_path)
        )

    # Check that the note contains a semicolon
    report_file: str = os.path.join(
        report_path, "stec_combined_nt_report.tsv"
    )
    with open(report_file, encoding="utf-8") as f:
        content: str = f.read()
        assert (
            "Partial Match;Novel Subunit A Match" in content
            or "Partial Match;Novel Subunit B Match" in content
        )


def test_nt_report_else_branch_for_note(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _nt_report function for the else branch that adds only the note.

    This test ensures that when the allele profile is empty and a partial
    match is present, the note is written in the last column of the report.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare metadata with a partial match and empty allele profile
    metadata: dict = {
        "strain1": {
            "nt_best_hits": {"stx1": [("allele1", 80.0)], "stx2": []},
            "nt_novel_hits": {"stx1": [("novel_allele", 100)], "stx2": []},
            "nt_query_ids": {"stx1": {"allele1": "contig1"}, "stx2": {}},
            "stx1_aa_allele_profile": "",
            "novel_aa_best_hits": {},
        }
    }
    report_path: pathlib.Path = tmp_path

    # Patch _has_novel_aa_best_hits to return True to trigger the else branch
    with patch.object(
        CombinedSubunits, "_has_novel_aa_best_hits", return_value=True
    ):
        _nt_report(
            metadata=metadata,
            report_path=str(report_path)
        )

    # Should write the note in the last column of the report
    report_file: str = os.path.join(
        report_path, "stec_combined_nt_report.tsv"
    )
    with open(report_file, encoding="utf-8") as f:
        content: str = f.read()
        assert content.strip().endswith("Partial Match")


def test_nt_report_not_data_novel_and_not_novel(
    tmp_path: pathlib.Path
) -> None:
    """
    Test _nt_report for both novel and not novel branches when there is no
    data for best hits or novel hits.

    Ensures that the report file contains the correct number of tab-separated
    columns for both the novel and not novel cases.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare metadata with no best hits or novel hits for two strains
    metadata: dict = {
        "strain1": {
            "nt_best_hits": {"stx1": [], "stx2": []},
            "nt_novel_hits": {"stx1": [], "stx2": []},
        },
        "strain2": {
            "nt_best_hits": {"stx1": [], "stx2": []},
            "nt_novel_hits": {"stx1": [], "stx2": []},
        },
    }
    report_path: pathlib.Path = tmp_path

    # First, test the novel=True branch
    with patch.object(
        CombinedSubunits, "_has_novel_aa_best_hits", return_value=False
    ):
        _nt_report(
            metadata=metadata,
            report_path=str(report_path)
        )
        report_file: str = os.path.join(
            report_path, "stec_combined_nt_report.tsv"
        )
        # Read the report file and check for correct number of columns
        with open(report_file, "r", encoding="utf-8") as f:
            lines: list[str] = f.readlines()
            assert any(
                line.endswith("\t\t\t\t\t\t\n") for line in lines
            ), (
                f"Expected line ending with six tabs for novel=True, "
                f"got: {lines}"
            )

    # Now, test the novel=False branch
    with patch.object(
        CombinedSubunits, "_has_novel_aa_best_hits", return_value=True
    ):
        _nt_report(
            metadata=metadata,
            report_path=str(report_path)
        )
        report_file: str = os.path.join(
            report_path, "stec_combined_nt_report.tsv"
        )
        # Read the report file and check for correct number of columns
        with open(report_file, "r", encoding="utf-8") as f:
            lines: list[str] = f.readlines()
            assert any(
                line.endswith("\t\t\n") for line in lines
            ), (
                f"Expected line ending with two tabs for novel=False, "
                f"got: {lines}"
            )


def test_aa_report_no_data_branch(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _aa_report function when there are no aa_best_hits for either
    stx1 or stx2.

    Ensures that the report file contains a line with the strain name and
    two tab-separated columns, as expected for the no data branch.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Metadata with no aa_best_hits for either stx1 or stx2
    metadata: dict = {
        "strainX": {
            "aa_best_hits": {"stx1": [], "stx2": []}
        }
    }
    report_path: pathlib.Path = tmp_path

    # Call the _aa_report function under test
    _aa_report(
        metadata=metadata,
        report_path=str(report_path)
    )

    # Read the report file and check for the expected line
    report_file: str = os.path.join(
        report_path, "stec_combined_aa_report.tsv"
    )
    with open(report_file, "r", encoding="utf-8") as f:
        lines: list[str] = f.readlines()

    # The second line (after header) should be: "strainX\t\t\n"
    assert any(line == "strainX\t\t\n" for line in lines), (
        f"Expected 'strainX\\t\\t\\n' in report, got: {lines}"
    )


def test_combined_report_partial_match(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _combined_report function for a partial match scenario.

    Ensures that when the percent identity is less than 90, the report
    file contains "Partial Match" in the Notes column.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Metadata with a partial match (perc_ident < 90)
    metadata: dict = {
        "strainA": {
            "nt_best_hits": {"stx1": [("allele1", 85.0)], "stx2": []},
            "aa_best_hits": {"stx1": [], "stx2": []},
        }
    }

    # Call the combined report function under test
    _combined_report(
        metadata=metadata,
        report_path=str(tmp_path)
    )

    # Path to the generated report file
    report_file: str = os.path.join(
        tmp_path, "stec_combined_nt_aa_report.tsv"
    )

    # Read the report file and check for "Partial Match" in the Notes column
    with open(report_file, "r", encoding="utf-8") as f:
        lines: list[str] = f.readlines()

    # Should contain "Partial Match" in the Notes column
    assert any("Partial Match" in line for line in lines), (
        f"Expected 'Partial Match' in report, got: {lines}"
    )


def test_combined_report_no_data_branch(
    tmp_path: pathlib.Path
) -> None:
    """
    Test the _combined_report function when there are no best hits for
    either stx1 or stx2.

    Ensures that the report file contains a line with the strain name
    and four tab-separated columns, as expected for the no data branch.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Metadata with no best hits for either stx1 or stx2
    metadata: dict = {
        "strainB": {
            "nt_best_hits": {"stx1": [], "stx2": []},
            "aa_best_hits": {"stx1": [], "stx2": []},
        }
    }

    # Call the _combined_report function under test
    _combined_report(
        metadata=metadata,
        report_path=str(tmp_path)
    )

    # Path to the generated report file
    report_file: str = os.path.join(
        tmp_path, "stec_combined_nt_aa_report.tsv"
    )

    # Read the report file and check for the expected line
    with open(report_file, "r", encoding="utf-8") as f:
        lines: list[str] = f.readlines()

    # Should contain a line with the strain name and four tabs
    assert any(line == "strainB\t\t\t\t\n" for line in lines), (
        f"Expected 'strainB\\t\\t\\t\\t\\n' in report, got: {lines}"
    )


def make_args(
    blast_mode: str = "blastx",
    preliminary: bool = False,
    split_aa_db_dir: str = '',
    num_alignments: int = 10,
    cutoff: float = 99.0,
    allele_path: str = "alleles",
    report_path: str = "reports",
    query_path: str = "query",
    verbosity: str = "info",
) -> types.SimpleNamespace:
    """
    Create a mock argparse.Namespace object for testing.

    :param blast_mode: The BLAST mode to use (e.g., "blastx", "blastn").
    :type blast_mode: str
    :param preliminary: Whether to run in preliminary mode.
    :type preliminary: bool
    :param split_aa_db_dir: Directory for split amino acid DBs.
    :type split_aa_db_dir: str
    :param num_alignments: Number of alignments to use.
    :type num_alignments: int
    :param cutoff: Percent identity cutoff.
    :type cutoff: float
    :param allele_path: Path to allele files.
    :type allele_path: str
    :param report_path: Path to report output.
    :type report_path: str
    :param query_path: Path to query files.
    :type query_path: str
    :param verbosity: Logging verbosity level.
    :type verbosity: str
    :return: A SimpleNamespace mimicking argparse.Namespace.
    :rtype: types.SimpleNamespace
    """
    # Mimic argparse.Namespace for test arguments
    return types.SimpleNamespace(
        blast_mode=blast_mode,
        preliminary=preliminary,
        split_aa_db_dir=split_aa_db_dir,
        num_alignments=num_alignments,
        cutoff=cutoff,
        allele_path=allele_path,
        report_path=report_path,
        query_path=query_path,
        verbosity=verbosity,
        version="1.0",
    )


@patch("allele_tools.stec_combined_subunits.CombinedSubunits")
@patch("os.path.isdir")
@patch("argparse.ArgumentParser.parse_args")
def test_main_blastn_requires_split_aa_db_dir(
    mock_parse_args: MagicMock,
    mock_isdir: MagicMock,
    mock_combined: MagicMock
) -> None:
    """
    Test that main() requires split_aa_db_dir for blastn mode.

    Ensures that if blast_mode is 'blastn', preliminary is False, and
    split_aa_db_dir is None, the parser.error is called and SystemExit
    is raised. Also checks that CombinedSubunits is not called.

    :param mock_parse_args: Mock for ArgumentParser.parse_args
    :type mock_parse_args: MagicMock
    :param mock_isdir: Mock for os.path.isdir
    :type mock_isdir: MagicMock
    :param mock_combined: Mock for CombinedSubunits
    :type mock_combined: MagicMock
    :return: None
    """
    args = make_args(
        blast_mode="blastn",
        preliminary=False,
        split_aa_db_dir=''
    )
    mock_parse_args.return_value = args
    mock_isdir.return_value = False

    # Patch parser.error to raise SystemExit, as it does in real argparse
    with patch(
        "argparse.ArgumentParser.error",
        side_effect=SystemExit
    ):
        with pytest.raises(SystemExit):
            main()

    # Now, CombinedSubunits should NOT have been called
    mock_combined.assert_not_called()


@patch("allele_tools.stec_combined_subunits.CombinedSubunits")
@patch("argparse.ArgumentParser.parse_args")
def test_main_preliminary_invalid_with_non_blastn(
    mock_parse_args: MagicMock,
    mock_combined: MagicMock
) -> None:
    """
    Test that main() raises SystemExit if preliminary is True and
    blast_mode is not 'blastn'.

    Ensures that if preliminary is True and blast_mode is not 'blastn',
    the parser.error is called and SystemExit is raised. Also checks
    that CombinedSubunits is not called.

    :param mock_parse_args: Mock for ArgumentParser.parse_args
    :type mock_parse_args: MagicMock
    :param mock_combined: Mock for CombinedSubunits
    :type mock_combined: MagicMock
    :return: None
    """
    # Prepare arguments with blast_mode not 'blastn' and preliminary True
    args = make_args(
        blast_mode="blastx",
        preliminary=True
    )
    mock_parse_args.return_value = args

    # Patch parser.error to raise SystemExit, as it does in real argparse
    with patch(
        "argparse.ArgumentParser.error",
        side_effect=SystemExit
    ):
        with pytest.raises(SystemExit):
            main()

    # CombinedSubunits should NOT have been called
    mock_combined.assert_not_called()


@patch("allele_tools.stec_combined_subunits.CombinedSubunits")
@patch("os.path.isdir")
@patch("argparse.ArgumentParser.parse_args")
def test_main_warns_if_split_aa_db_dir_not_needed(
    mock_parse_args: MagicMock,
    mock_isdir: MagicMock,
    mock_combined: MagicMock
) -> None:
    """
    Test that main() warns if split_aa_db_dir is provided but not needed.

    Ensures that if split_aa_db_dir is provided for a mode that does not
    require it (e.g., "blastx"), a warning is logged and the pipeline
    proceeds as normal.

    :param mock_parse_args: Mock for ArgumentParser.parse_args
    :type mock_parse_args: MagicMock
    :param mock_isdir: Mock for os.path.isdir
    :type mock_isdir: MagicMock
    :param mock_combined: Mock for CombinedSubunits
    :type mock_combined: MagicMock
    :return: None
    """
    # Should warn if split_aa_db_dir is provided but not needed
    args = make_args(
        blast_mode="blastx",
        split_aa_db_dir="some_dir"
    )
    mock_parse_args.return_value = args
    mock_isdir.return_value = True

    # Patch logging.warning to track warning calls
    with patch("logging.warning") as mock_warn:
        main()
        mock_warn.assert_called_once()

    # Ensure CombinedSubunits is called once
    mock_combined.assert_called_once()


@patch("allele_tools.stec_combined_subunits.CombinedSubunits")
@patch("os.path.isdir")
@patch("argparse.ArgumentParser.parse_args")
def test_main_runs_pipeline_and_logs_time(
    mock_parse_args: MagicMock,
    mock_isdir: MagicMock,
    mock_combined: MagicMock
) -> None:
    """
    Test that main() runs the pipeline and logs elapsed time.

    Ensures that when the pipeline is run, the process method is called
    on the CombinedSubunits instance and that a log message indicating
    pipeline completion is generated.

    :param mock_parse_args: Mock for ArgumentParser.parse_args
    :type mock_parse_args: MagicMock
    :param mock_isdir: Mock for os.path.isdir
    :type mock_isdir: MagicMock
    :param mock_combined: Mock for CombinedSubunits
    :type mock_combined: MagicMock
    :return: None
    """
    # Should run the pipeline and log elapsed time

    args = make_args()
    mock_parse_args.return_value = args
    mock_isdir.return_value = True

    # Create a mock instance for CombinedSubunits
    mock_instance: MagicMock = MagicMock()
    mock_combined.return_value = mock_instance

    # Patch logging.info to track log messages
    with patch("logging.info") as mock_info:
        main()

        # Should call process on the CombinedSubunits instance
        mock_instance.process.assert_called_once()

        # Should log pipeline completion message
        assert any(
            "Pipeline completed in" in str(call.args[0])
            for call in mock_info.call_args_list
        )


def test_analyze_sequence_for_gaps_and_ambiguity_nt_and_aa():
    """
    Test _analyze_sequence_for_gaps_and_ambiguity for both nt and aa molecules,
    ensuring 1-based gap and ambiguous base reporting.
    """
    # nt test: gap at 3-5 (1-based), ambiguous at 2 (R) and 5 (N)
    seq_nt = "AR--TGN"
    cleaned_nt, gaps_nt, ambig_nt = (
        _analyze_sequence_for_gaps_and_ambiguity(
            seq=seq_nt, molecule="nt"
        )
    )
    assert cleaned_nt == "ARTGN"
    assert gaps_nt == [(3, 4)]
    assert ambig_nt == [(2, "R"), (5, "N")]

    # aa test: gap at 2-3 (1-based), ambiguous at 1 (B) and 3 (Z)
    seq_aa = "B--AZ"
    cleaned_aa, gaps_aa, ambig_aa = (
        _analyze_sequence_for_gaps_and_ambiguity(
            seq=seq_aa, molecule="aa"
        )
    )
    assert cleaned_aa == "BAZ"
    assert gaps_aa == [(2, 3)]
    assert ambig_aa == [(1, "B"), (3, "Z")]

    # no gaps or ambiguity
    seq_clean = "ACGT"
    cleaned, gaps, ambig = \
        _analyze_sequence_for_gaps_and_ambiguity(
            seq=seq_clean, molecule="nt"
        )
    assert cleaned == "ACGT"
    assert gaps == []
    assert ambig == []


def test_extract_best_hits_skips_ambiguous():
    """
    Test that the extraction of best hits skips ambiguous bases.
    """
    # Prepare a BLAST output with an ambiguous base (N)
    rows = [
        {
            "subject_id": "stx1|allele1",
            "percent_match": "99.0",
            "query_sequence": "ATG--CN",
            "query_id": "q1",
        }
    ]
    with tempfile.NamedTemporaryFile("w+", delete=False) as tmp:
        writer = DictWriter(
            tmp,
            fieldnames=rows[0].keys(),
            dialect="excel-tab"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
        tmp.flush()

        # Call the method
        result = _extract_best_hits(
            blast_output=tmp.name,
            cutoff=95.0,
            molecule="nt"
        )
        (
            best_hits,
            novel_hits,
            query_sequences,
            _,
            gap_events_dict,
            ambiguous_dict,
        ) = result

        # The ambiguous base should cause the hit to be skipped for
        # novel_hits/query_sequences
        assert best_hits["stx1"] == [("stx1", 99.0)]
        assert novel_hits["stx1"] == []
        assert query_sequences["stx1"] == {}
        assert "stx1" in gap_events_dict
        assert "stx1" in ambiguous_dict
        assert ambiguous_dict["stx1"] == [(5, "N")]


def test_nt_report_prefers_first_gap_or_ambiguous():
    """
    Test that _nt_report prefers the first key with gap/ambiguous info.
    """

    # Simulate metadata with both allele_clean and novel_allele keys
    metadata = {
        "strain1": {
            "nt_best_hits": {"stx1": [("alleleA|extra", 99.0)], "stx2": []},
            "nt_novel_hits": {"stx1": [("novel_alleleA", 100)], "stx2": []},
            "nt_query_ids": {"stx1": {"alleleA": "contig1"}, "stx2": {}},
            "nt_gap_events": {"alleleA": [(2, 3)], "novel_alleleA": [(5, 6)]},
            "nt_ambiguous": {
                "alleleA": [(4, "N")], "novel_alleleA": [(7, "R")]
            },
        }
    }

    # Open the report file and check its contents
    with tempfile.TemporaryDirectory() as tmpdir:
        _nt_report(metadata=metadata, report_path=tmpdir)
        report_file = os.path.join(tmpdir, "stec_combined_nt_report.tsv")
        with open(report_file, encoding='utf-8') as f:
            lines = f.readlines()
        # Should use the first key (alleleA) for gap/ambiguous info
        assert "Ambiguous bases at 4:N" in lines[1]
        assert "1 gap(s) at 2-3" in lines[1]
        # Should not use novel_alleleA's gap/ambiguous info


def test_nt_report_notes_for_gaps_and_ambiguous():
    """
    Test that ambiguous and gap notes are correctly added to the report.
    """

    metadata = {
        "strain2": {
            "nt_best_hits": {"stx1": [("alleleB", 88.0)], "stx2": []},
            "nt_novel_hits": {"stx1": [("novel_alleleB", 100)], "stx2": []},
            "nt_query_ids": {"stx1": {"alleleB": "contig2"}, "stx2": {}},
            "nt_gap_events": {"alleleB": [(1, 2), (5, 6)]},
            "nt_ambiguous": {"alleleB": [(3, "N"), (7, "R")]},
        }
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        _nt_report(metadata=metadata, report_path=tmpdir)
        report_file = os.path.join(tmpdir, "stec_combined_nt_report.tsv")
        with open(report_file, encoding='utf-8') as f:
            lines = f.readlines()
        # Should include ambiguous and gap notes, with correct formatting
        assert "Partial Match" in lines[1]
        assert "Ambiguous bases at 3:N, 7:R" in lines[1]
        assert "2 gap(s) at 1-2, 5-6" in lines[1]


def test_write_novel_allele_subunit_files_no_copy_if_same_path(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test that shutil.copy is NOT called if a_path and dst_a
    (or b_path and dst_b) are the same file.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare valid sequences
    a_seq = "A" * 313
    b_seq = "B" * 82

    # Use the same directory for path and report_path
    # This will make a_path == dst_a and b_path == dst_b
    with patch("shutil.copy") as mock_copy:
        file_paths, notes = CombinedSubunits._write_novel_allele_subunit_files(
            molecule="aa",
            novel_allele_name="strainX_stx1",
            path=str(tmp_path),
            subunit_a_seq=a_seq,
            subunit_b_seq=b_seq,
            report_path=str(tmp_path),
        )
        # shutil.copy should NOT be called since src == dst
        mock_copy.assert_not_called()
    assert "A" in file_paths and "B" in file_paths and "combined" in file_paths
    assert notes == {}


def test_write_novel_allele_subunit_files_copy_if_different_path(
    tmp_path: pathlib.Path,
) -> None:
    """
    Test that shutil.copy IS called if a_path and dst_a (or b_path and dst_b)
    are different files.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Prepare valid sequences
    a_seq = "A" * 313
    b_seq = "B" * 82

    # Use different directories for path and report_path
    report_dir = tmp_path / "report"
    report_dir.mkdir()

    with patch("shutil.copy") as mock_copy:
        file_paths, notes = CombinedSubunits._write_novel_allele_subunit_files(
            molecule="aa",
            novel_allele_name="strainX_stx1",
            path=str(tmp_path),
            subunit_a_seq=a_seq,
            subunit_b_seq=b_seq,
            report_path=str(report_dir),
        )
        # Only the combined file should be copied
        assert mock_copy.call_count == 1
    assert "A" in file_paths and "B" in file_paths and "combined" in file_paths
    assert notes == {}


def test_write_novel_allele_subunit_files_real_copy(
    tmp_path: pathlib.Path
) -> None:
    """
    Test that shutil.copy is called if a_path and dst_a (or b_path and dst_b)
    are different files.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Use different directories for path and report_path
    report_dir = tmp_path / "report"
    report_dir.mkdir()

    a_seq = "A" * 313
    b_seq = "B" * 82

    _, notes = CombinedSubunits._write_novel_allele_subunit_files(
        molecule="aa",
        novel_allele_name="strainX_stx1",
        path=str(tmp_path),
        subunit_a_seq=a_seq,
        subunit_b_seq=b_seq,
        report_path=str(report_dir),
    )
    # Check that the combined file was copied to the report directory
    combined_report = report_dir / "strainX_stx1_combined_aa.fasta"
    assert combined_report.exists()
    assert notes == {}


def test_export_aa_novel_operons_skips_if_missing_fields(
    tmp_path: pathlib.Path
) -> None:
    """
    Test skips if profile, operon_seq, subunit_a_seq, or subunit_b_seq is
    missing.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    meta = {
        "strain1": {
            "stx1_aa_allele_profile": "",
            "stx1_aa_operon_sequence": "A" * 313 + "XXXX" + "B" * 82,
            "stx1_aa_subunit_a_sequence": "A" * 313,
            "stx1_aa_subunit_b_sequence": "B" * 82,
        }
    }
    with patch("Bio.SeqIO.write") as mock_write:
        CombinedSubunits._export_aa_novel_operons(
            metadata=meta, report_path=str(tmp_path)
        )
        mock_write.assert_not_called()


def test_export_aa_novel_operons_skip_if_profile_split_fails(
    tmp_path: pathlib.Path
) -> None:
    """
    Test skips if profile can't be split (ValueError).

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    meta = {
        "strain1": {
            "stx1_aa_allele_profile": "badprofile",
            "stx1_aa_operon_sequence": "A" * 313 + "XXXX" + "B" * 82,
            "stx1_aa_subunit_a_sequence": "A" * 313,
            "stx1_aa_subunit_b_sequence": "B" * 82,
        }
    }
    with patch("Bio.SeqIO.write") as mock_write:
        CombinedSubunits._export_aa_novel_operons(
            metadata=meta, report_path=str(tmp_path)
        )
        mock_write.assert_not_called()


def fake_save_and_validate_subunit(
    expected_length: int,
    fasta_path: pathlib.Path,
    novel_allele_name: str,
    subunit_name: str,
    subunit_seq: str
) -> tuple[str, bool, str]:
    """
    Simulate saving and validating a subunit FASTA file.

    :param expected_length: Expected length of the subunit sequence.
    :type expected_length: int
    :param fasta_path: Path to the FASTA file to save.
    :type fasta_path: pathlib.Path
    :param novel_allele_name: Name of the novel allele.
    :type novel_allele_name: str
    :param subunit_name: Name of the subunit (e.g., "A" or "B").
    :type subunit_name: str
    :param subunit_seq: Sequence of the subunit.
    :type subunit_seq: str
    :return: Tuple of (sequence, is_valid, note).
        - sequence: The sequence written to the file.
        - is_valid: True if the sequence is valid, False otherwise.
        - note: Any notes about the validation.
    """
    if subunit_name == "A":
        expected_length = 313
        novel_allele_name = ">fake_allele\n"
        subunit_seq = "A" * expected_length
        with open(fasta_path, "w", encoding='utf-8') as f:
            f.write(novel_allele_name + subunit_seq + "\n")
        return (subunit_seq, True, "")

    if subunit_name == "B":
        expected_length = 82
        novel_allele_name = ">fake_allele\n"
        subunit_seq = "B" * expected_length
        # Actually write the file to simulate the real function
        with open(fasta_path, "w", encoding='utf-8') as f:
            f.write(novel_allele_name + subunit_seq + "\n")
        return (subunit_seq, True, "")
    return ("", False, "")


def test_export_aa_novel_operons_a_perc_100_sets_a_seq_pass_note(
    tmp_path: pathlib.Path
) -> None:
    """
    Test a_seq, a_pass, a_note set to '', None, None if a_perc == 100.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    meta = {
        "strain1": {
            "stx1_aa_allele_profile": "22_38",
            "stx1_aa_operon_sequence": "A" * 313 + "XXXX" + "B" * 82,
            "stx1_aa_subunit_a_sequence": "A" * 313,
            "stx1_aa_subunit_b_sequence": "B" * 82,
            "novel_aa_best_hits": {
                "allele1": {
                    "Stx1A": [("Stx1A_22", 100.0)],
                    "Stx1B": [("Stx1B_38", 99.0)],
                }
            },
        }
    }

    with patch.object(
        CombinedSubunits,
        "_save_and_validate_subunit",
        side_effect=fake_save_and_validate_subunit,
    ):
        CombinedSubunits._export_aa_novel_operons(
            metadata=meta, report_path=str(tmp_path)
        )
    files = list(tmp_path.glob("*"))

    assert any("subunit_B" in f.name for f in files)
    assert not any("subunit_A" in f.name for f in files)
    assert not any("operon" in f.name for f in files)


def test_export_aa_novel_operons_b_perc_100_sets_b_seq_pass_note(
    tmp_path: pathlib.Path
) -> None:
    """
    Test b_seq, b_pass, b_note set to '', None, None if b_perc == 100.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    meta = {
        "strain1": {
            "stx1_aa_allele_profile": "22_38",
            "stx1_aa_operon_sequence": "A" * 313 + "XXXX" + "B" * 82,
            "stx1_aa_subunit_a_sequence": "A" * 313,
            "stx1_aa_subunit_b_sequence": "B" * 82,
            "novel_aa_best_hits": {
                "allele1": {
                    "Stx1A": [("Stx1A_22", 99.0)],
                    "Stx1B": [("Stx1B_38", 100.0)],
                }
            },
        }
    }

    with patch.object(
        CombinedSubunits,
        "_save_and_validate_subunit",
        side_effect=fake_save_and_validate_subunit,
    ):
        CombinedSubunits._export_aa_novel_operons(
            metadata=meta, report_path=str(tmp_path)
        )
    files = list(tmp_path.glob("*"))

    assert any("subunit_A" in f.name for f in files)
    assert not any("subunit_B" in f.name for f in files)
    assert not any("operon" in f.name for f in files)


def test_export_aa_novel_operons_notes_and_no_operon(
    tmp_path: pathlib.Path
) -> None:
    """
    Test notes are set and operon is not written if either subunit fails
    validation.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    meta = {
        "strain1": {
            "stx1_aa_allele_profile": "22_38",
            "stx1_aa_operon_sequence": "A" * 313 + "XXXX" + "B" * 82,
            "stx1_aa_subunit_a_sequence": "A" * 10,  # too short
            "stx1_aa_subunit_b_sequence": "B" * 10,  # too short
            "novel_aa_best_hits": {
                "allele1": {
                    "Stx1A": [("Stx1A_22", 99.0)],
                    "Stx1B": [("Stx1B_38", 99.0)],
                }
            },
        }
    }
    with patch.object(
        CombinedSubunits,
        "_save_and_validate_subunit",
        side_effect=[("A" * 10, False, "A fail"), ("B" * 10, False, "B fail")],
    ):
        CombinedSubunits._export_aa_novel_operons(
            metadata=meta, report_path=str(tmp_path)
        )
    # No operon file should be written
    files = list(tmp_path.glob("*"))
    assert not any("operon" in f.name for f in files)
    # Notes should be set in metadata
    notes = meta["strain1"]["novel_subunit_notes"]
    assert "A" in notes["strain1_STX1_22_38"]
    assert "B" in notes["strain1_STX1_22_38"]


def test_export_aa_novel_operons_skips_if_any_field_missing(
    tmp_path: pathlib.Path
) -> None:
    """
    Test that no files are written if any required field is missing.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    # Each dict omits one required field
    test_cases = [
        {
            "stx1_aa_allele_profile": "",
            "stx1_aa_operon_sequence": "A",
            "stx1_aa_subunit_a_sequence": "A",
            "stx1_aa_subunit_b_sequence": "B"
        },
        {
            "stx1_aa_allele_profile": "22_38",
            "stx1_aa_operon_sequence": "",
            "stx1_aa_subunit_a_sequence": "A",
            "stx1_aa_subunit_b_sequence": "B"
        },
        {
            "stx1_aa_allele_profile": "22_38",
            "stx1_aa_operon_sequence": "A",
            "stx1_aa_subunit_a_sequence": "",
            "stx1_aa_subunit_b_sequence": "B"
        },
        {
            "stx1_aa_allele_profile": "22_38",
            "stx1_aa_operon_sequence": "A",
            "stx1_aa_subunit_a_sequence": "A",
            "stx1_aa_subunit_b_sequence": ""
        },
    ]
    for i, fields in enumerate(test_cases):
        meta = {"strain1": fields}
        CombinedSubunits._export_aa_novel_operons(
            metadata=meta,
            report_path=str(tmp_path)
        )
        files = list(tmp_path.glob("*"))
        assert not files, \
            f"Files were written for missing field case {i}: {files}"


def test_export_aa_novel_operons_skips_if_both_subunits_100(
    tmp_path: pathlib.Path
) -> None:
    """
    Test that no files are written if both subunits are 100% identical.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    meta = {
        "strain1": {
            "stx1_aa_allele_profile": "22_38",
            "stx1_aa_operon_sequence": "A" * 313 + "XXXX" + "B" * 82,
            "stx1_aa_subunit_a_sequence": "A" * 313,
            "stx1_aa_subunit_b_sequence": "B" * 82,
            "novel_aa_best_hits": {
                "allele1": {
                    "Stx1A": [("Stx1A_22", 100.0)],
                    "Stx1B": [("Stx1B_38", 100.0)],
                }
            },
        }
    }
    CombinedSubunits._export_aa_novel_operons(
        metadata=meta,
        report_path=str(tmp_path)
    )
    files = list(tmp_path.glob("*"))
    assert not files, \
        f"Files were written when both subunits are 100%: {files}"


def test_nt_report_includes_novel_subunit_notes_by_novel_allele(
    tmp_path: pathlib.Path
) -> None:
    """
    Test that novel_subunit_notes are included when keyed by novel_allele.

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    meta = {
        "strain1": {
            "nt_best_hits": {"stx1": [("novel1|extra", 99.5)], "stx2": []},
            "nt_query_ids": {"stx1": {"novel1|extra": "contig1"}},
            "nt_novel_hits": {"stx1": [("novel1", 99.5)]},
            "novel_subunit_notes": {"novel1": {"A": "A note", "B": "B note"}},
            "novel_aa_best_hits": {},
            "stx1_aa_allele_profile": "",
            "stx2_aa_allele_profile": "",
        }
    }
    CombinedSubunits._nt_report(metadata=meta, report_path=str(tmp_path))
    report = (tmp_path / "stec_combined_nt_report.tsv").read_text()
    # Should include both subunit notes
    assert "A note" in report
    assert "B note" in report


def test_nt_report_appends_subunit_notes_to_existing_note(
    tmp_path: pathlib.Path
) -> None:
    """
    Test that subunit_notes are appended to an existing note
    (e.g., Partial Match).

    :param tmp_path: Temporary directory provided by pytest for file IO
    :type tmp_path: pathlib.Path
    :return: None
    """
    meta = {
        "strain1": {
            "nt_best_hits": {"stx1": [("novel2|extra", 85.0)], "stx2": []},
            "nt_query_ids": {"stx1": {"novel2|extra": "contig2"}},
            "nt_novel_hits": {"stx1": [("novel2", 85.0)]},
            "novel_subunit_notes": {"novel2": {"A": "A note"}},
            "novel_aa_best_hits": {},
            "stx1_aa_allele_profile": "",
            "stx2_aa_allele_profile": "",
        }
    }
    CombinedSubunits._nt_report(metadata=meta, report_path=str(tmp_path))
    report = (tmp_path / "stec_combined_nt_report.tsv").read_text()
    # Should include both "Partial Match" and the subunit note, separated by ;
    assert "Partial Match;A note" in report

    # Delete the orphaned files from the tests
    if os.path.exists("novel_alleles.fasta"):
        os.remove("novel_alleles.fasta")
    if os.path.exists("stec_combined_aa_report.tsv"):
        os.remove("stec_combined_aa_report.tsv")
    if os.path.exists("stec_combined_nt_report.tsv"):
        os.remove("stec_combined_nt_report.tsv")
