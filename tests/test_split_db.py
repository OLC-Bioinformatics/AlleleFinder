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
    processor = FastaProcessor('test.fasta')
    assert processor.filename == 'test.fasta'


@patch('allele_tools.split_db.SeqIO.write')
@patch('builtins.open', new_callable=mock_open)
def test_write_record_stx1(mock_open, mock_write):
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_file_path = os.path.join(temp_dir, 'test.fasta')
        processor = FastaProcessor(temp_file_path)
        record = MagicMock()

        # Test case for 'stx1'
        processor._write_record('stx1', 'header', record)
        mock_write.assert_called_with(record, mock_open(), 'fasta')


@patch('allele_tools.split_db.SeqIO.write')
@patch('builtins.open', new_callable=mock_open)
def test_write_record_stx2(mock_open, mock_write):
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_file_path = os.path.join(temp_dir, 'test.fasta')
        processor = FastaProcessor(temp_file_path)
        record = MagicMock()

        # Test case for 'stx2'
        processor._write_record('stx2', 'header', record)
        mock_write.assert_called_with(record, mock_open(), 'fasta')


@patch('allele_tools.split_db.SeqIO.parse')
@patch('builtins.open', new_callable=mock_open)
@patch.object(FastaProcessor, '_write_record', autospec=True)
def test_process_records(mock_write_record, mock_open, mock_parse):
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_file_path = os.path.join(temp_dir, 'test.fasta')
        processor = FastaProcessor(temp_file_path)

        # Create mock records
        record1 = SeqRecord(Seq("ATGC"), id="stx1_record", description="stx1")
        record2 = SeqRecord(Seq("ATGC"), id="stx2_record", description="stx2")

        # Mock the return value of SeqIO.parse
        mock_parse.return_value = [record1, record2]

        # Call the method that contains the code to be tested
        processor.process()

        # Check that _write_record was called with the correct arguments
        mock_write_record.assert_any_call(processor, 'stx1', 'stx1', record1)
        mock_write_record.assert_any_call(processor, 'stx2', 'stx2', record2)


@patch('allele_tools.split_db.SeqIO.parse')
@patch('builtins.open', new_callable=mock_open)
def test_process(mock_open, mock_parse):
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_file_path = os.path.join(temp_dir, 'test.fasta')
        processor = FastaProcessor(temp_file_path)
        processor.process()
        mock_parse.assert_called_with(mock_open(), 'fasta')


@patch.object(FastaProcessor, 'process')
@patch.object(FastaProcessor, '__init__', return_value=None)
@patch('argparse.ArgumentParser.parse_args')
def test_main(mock_parse_args, mock_init, mock_process):
    # Mock the command line arguments
    mock_parse_args.return_value = MagicMock(filename='test.fasta')

    # Call the main function
    cli()

    # Check that FastaProcessor was initialized with the correct filename
    mock_init.assert_called_once_with('test.fasta')

    # Check that the process method was called
    mock_process.assert_called_once()
