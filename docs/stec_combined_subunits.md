## Table of Contents
1. [Introduction](#introduction)
2. [Prerequisites](#prerequisites)
3. [Inputs](#inputs)
4. [Running the Script](#running-the-script)
5. [Command Line Arguments](#command-line-arguments)
6. [Outputs](#outputs)
7. [Interpreting the Results](#interpreting-the-results)
8. [Troubleshooting](#troubleshooting)
9. [Additional Resources](#additional-resources)

## Introduction <a name="introduction"></a>

This script processes STEC alleles with both A and B subunits. It performs BLAST analyses to identify the best matching alleles for each subunit, detects novel alleles, and generates reports summarizing the results. This is essential for accurate subunit-level typing and characterization of Shiga toxin-producing _Escherichia coli_ (STEC).

## Prerequisites <a name="prerequisites"></a>

Before running the script, ensure you have the following:

1. Query files in FASTA format containing nucleotide or amino acid sequences.
2. Allele database files for the relevant subunits (A and B) in FASTA format.

## Inputs <a name="inputs"></a>

The script requires the following inputs:

1. Query files in FASTA format.
2. Allele database files for the A and B subunits (e.g., `Stx1A_aa_*.fasta`, `Stx1B_aa_*.fasta`, `Stx2A_aa_*.fasta`, `Stx2B_aa_*.fasta`).
3. Optionally, a directory containing split AA allele databases for subunit-level BLAST analyses.

## Running the Script <a name="running-the-script"></a>

To run the script, use the following command:

```bash
stec.py stec_combined_subunits --allele_path /path/to/allele_folder --report_path /path/to/output_folder --query_path /path/to/query_folder --blast_mode blastn --split_aa_db_dir /path/to/split_aa_db_dir
```

## Command Line Arguments <a name="command-line-arguments"></a>

The script accepts the following command line arguments:

- `-h, --help`: Show this help message and exit.
- `--version`: Show program's version number and exit.
- `-a, --allele_path`: Specify the path to the folder containing allele files. Default: `alleles` in the current working directory.
- `-r, --report_path`: Specify the path to the folder for output reports. Default: `reports` in the current working directory.
- `-q, --query_path`: Specify the path to the folder containing query files in FASTA format. Default: `query` in the current working directory.
- `--blast_mode`: Choose BLAST mode: `blastn`, `tblastx`, `blastx`, or `blastn+tblastx`. Default: `blastx`.
- `--preliminary`: Run a preliminary screen using blastn (do not run additional downstream analyses).
- `--split_aa_db_dir`: Path to directory containing split AA allele databases (required for full analysis).
- `-n, --num_alignments`: Number of alignments to return. Default: 10.
- `-c, --cutoff`: Percent identity cutoff. Default: 99.0.
- `-v, --verbosity`: Set the logging level. Options: debug, info, warning, error, critical. Default: info.

## Outputs <a name="outputs"></a>

The script generates the following output files:

- `novel_alleles.fasta`: FASTA file containing all novel alleles discovered.
- `stec_combined_nt_report.tsv`: TSV report summarizing nucleotide-level matches and novel alleles.
- `stec_combined_aa_report.tsv`: TSV report summarizing amino acid-level matches and novel alleles.
- `stec_combined_nt_aa_report.tsv`: Combined report for both nucleotide and amino acid analyses.
- Individual FASTA files for novel subunits and operons, as appropriate.

## Interpreting the Results <a name="interpreting-the-results"></a>

The main reports (`stec_combined_nt_report.tsv`, `stec_combined_aa_report.tsv`, and `stec_combined_nt_aa_report.tsv`) summarize the best matches for each query, percent identity, and any novel alleles detected. Additional notes may indicate partial matches, ambiguous bases, or gap events.

## Troubleshooting <a name="troubleshooting"></a>

If you encounter issues:

- Ensure all input files are in the correct format and location.
- Confirm that BLAST+ is installed and accessible from your command line.
- Check that the split AA database directory contains the required subunit FASTA files.
- Review error messages for missing files or permission issues.

## Additional Resources <a name="additional-resources"></a>

For more information, see the [`stec_combined_subunits`](https://olc-bioinformatics.github.io/AlleleFinder/stec_combined_subunits) documentation or open an issue on [GitHub](https://github.com/OLC-Bioinformatics/AlleleFinder/issues).