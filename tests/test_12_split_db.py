#!/usr/bin/env python

"""
Unit and integration tests for allele_tools/split_db.py
"""

# Standard imports
import os
import tempfile
from unittest.mock import (
    MagicMock,
    mock_open,
    patch
)

# Third party imports
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Local imports
from allele_tools.split_db import (
    cli,
    FastaProcessor,
)


def test_init():
    """
    Test the initialization of the FastaProcessor.
    """
    processor = FastaProcessor(
        filename='test.fasta',
        split_stx=False
    )
    assert processor.filename == 'test.fasta'


@patch('allele_tools.split_db.SeqIO.write')
@patch('builtins.open', new_callable=mock_open)
def test_write_record_stx1(mock_open, mock_write):
    """
    Test the writing of a record for 'stx1'.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_file_path = os.path.join(temp_dir, 'test.fasta')
        processor = FastaProcessor(
            filename=temp_file_path,
            split_stx=False
        )
        record = MagicMock()

        # Test case for 'stx1'
        processor._write_record(folder='stx1', header='header', record=record)
        mock_write.assert_called_with(record, mock_open(), 'fasta')


@patch('allele_tools.split_db.SeqIO.write')
@patch('builtins.open', new_callable=mock_open)
def test_write_record_stx2(mock_open, mock_write):
    """
    Test the writing of a record for 'stx2'.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_file_path = os.path.join(temp_dir, 'test.fasta')
        processor = FastaProcessor(
            filename=temp_file_path,
            split_stx=False
        )
        record = MagicMock()

        # Test case for 'stx2'
        processor._write_record(folder='stx2', header='header', record=record)
        mock_write.assert_called_with(record, mock_open(), 'fasta')


@patch("allele_tools.split_db.SeqIO.parse")
@patch("builtins.open", new_callable=mock_open)
@patch.object(FastaProcessor, "_write_record", autospec=True)
def test_process_records(mock_write_record, mock_open, mock_parse):
    """
    Test the processing of records in the FastaProcessor.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_file_path = os.path.join(temp_dir, "test.fasta")
        processor = FastaProcessor(filename=temp_file_path, split_stx=False)

        # Create mock records
        record1 = SeqRecord(Seq("ATGC"), id="stx1_record", description="stx1")
        record2 = SeqRecord(Seq("ATGC"), id="stx2_record", description="stx2")

        # Mock the return value of SeqIO.parse
        mock_parse.return_value = [record1, record2]

        # Call the method that contains the code to be tested
        processor.process()

        # Both records should be written to the 'split' folder
        mock_write_record.assert_any_call(
            processor, folder="split", header="stx1", record=record1
        )
        mock_write_record.assert_any_call(
            processor, folder="split", header="stx2", record=record2
        )


@patch('allele_tools.split_db.SeqIO.parse')
@patch('builtins.open', new_callable=mock_open)
def test_process(mock_open, mock_parse):
    """
    Test the process method of FastaProcessor.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_file_path = os.path.join(temp_dir, 'test.fasta')
        processor = FastaProcessor(
            filename=temp_file_path,
            split_stx=True
        )
        processor.process()
        mock_parse.assert_called_with(
            mock_open(), 'fasta'
        )


@patch.object(FastaProcessor, "process")
@patch.object(FastaProcessor, "__init__", return_value=None)
@patch("argparse.ArgumentParser.parse_args")
def test_main(mock_parse_args, mock_init, mock_process):
    """
    Test the main function of the CLI.
    """
    # Mock the command line arguments
    mock_parse_args.return_value = MagicMock(
        filename="test.fasta",
        stx=True,
    )

    # Call the main function
    cli()

    # Check that FastaProcessor was initialized with the correct filename
    mock_init.assert_called_once_with(filename="test.fasta", split_stx=True)


@patch.object(FastaProcessor, "_write_record", autospec=True)
def test_process_stx1_branch(
    mock_write_record: MagicMock
) -> None:
    """
    Test that a record with header starting with 'stx1' triggers the stx1
    branch.

    :param mock_write_record: Mock object for the _write_record method.
    :return: None
    """
    processor = FastaProcessor(filename="dummy.fasta", split_stx=True)
    record = SeqRecord(Seq("ATGC"), id="id", description="desc")
    # Patch the header extraction logic if needed, or set id/description
    # accordingly
    record.id = "stx1_allele"
    record.description = "stx1_allele"
    with (
        patch("allele_tools.split_db.SeqIO.parse", return_value=[record]),
        patch("builtins.open", mock_open()),
    ):
        processor.process()
    mock_write_record.assert_any_call(
        processor, folder="stx1", header="stx1_allele", record=record
    )


@patch.object(FastaProcessor, "_write_record", autospec=True)
def test_process_stx2_branch(
    mock_write_record: MagicMock
) -> None:
    """
    Test that a record with header starting with 'stx2' triggers the stx2
    branch.

    :param mock_write_record: Mock object for the _write_record method.
    """
    processor = FastaProcessor(filename="dummy.fasta", split_stx=True)
    record = SeqRecord(Seq("ATGC"), id="id", description="desc")
    record.id = "stx2_allele"
    record.description = "stx2_allele"
    with (
        patch("allele_tools.split_db.SeqIO.parse", return_value=[record]),
        patch("builtins.open", mock_open()),
    ):
        processor.process()
    mock_write_record.assert_any_call(
        processor, folder="stx2", header="stx2_allele", record=record
    )
