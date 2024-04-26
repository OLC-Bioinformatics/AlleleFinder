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

This script performs BLAST analyses on a nucleotide allele database to discover their sequence types. It updates nucleotide and amino acid profiles and allele databases, which is crucial for understanding the genetic diversity of your samples.

## Prerequisites <a name="prerequisites"></a>

Before running the script, ensure you have the following:

1. Nucleotide query files in FASTA format.
2. All outputs from [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce).

## Inputs <a name="inputs"></a>

The script requires the following inputs:

1. Nucleotide query files in FASTA format.
2. All outputs from [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce).

## Running the Script <a name="running-the-script"></a>

To run the script, use the following command:

```bash
stec.py allele_find --nt_profile /path/to/nt_profile_file --aa_profile /path/to_aa_profile_file --nt_alleles /path/to/nt_allele_folder --aa_alleles /path/to/aa_allele_folder -r /path/to/output_folder -q /path/to/query_folder
```

## Command Line Arguments <a name="command-line-arguments"></a>

The script accepts the following command line arguments:

- `-h, --help`: Show this help message and exit.
- `-version, --version`: Show program's version number and exit.
- `-v verbosity, --verbosity verbosity`: Set the logging level. Options are debug, info, warning, error, and critical. Default is info.
- `--nt_profile nt_profile`: Specify name and path of nucleotide profile file. If not provided, `profile.txt` in the `nt_profile` folder in the current working directory will be used by default.
- `--aa_profile aa_profile`: Specify name and path of amino acid profile file. If not provided, `profile.txt` in the `aa_profile` folder in the current working directory will be used by default.
- `--nt_alleles nt_alleles`: Specify name and path of folder containing nucleotide alleles. If not provided, the `nt_allele` folder in the current working directory will be used by default.
- `--aa_alleles aa_alleles`: Specify name and path of folder containing amino acid alleles. If not provided, the `aa_allele` folder in the current working directory will be used by default.
- `-r report_path, --report_path report_path`: Specify name and path of folder into which reports are to be placed. If not provided, the `reports` folder in the current working directory will be used.
- `-q query_path, --query_path query_path`: Specify name and path of folder containing query files in FASTA format. If not provided, the `query` folder in the current working directory will be used.

## Outputs <a name="outputs"></a>

The script generates the following output files:

- `aa_novel_profiles.txt`: A text file containing all novel amino acid profiles generated from query sequences.
- `aa_gene_name_novel_alleles.fasta`: A FASTA-formatted text file containing all novel amino acid alleles from query sequences.
- `nt_novel_profiles.txt`: A text file containing all novel nucleotide profiles generated from query sequences.
- `nt_gene_name_novel_alleles.fasta`: A FASTA-formatted text file containing all novel nucleotide alleles from query sequences.
- `stec_report.tsv`: A TSV file containing results for each query. Includes sample name, nucleotide allele identifiers, nucleotide sequence type, amino acid allele identifiers, amino acid sequence type, and notes.

## Interpreting the Results <a name="interpreting-the-results"></a>

The `stec_report.tsv` file contains the main results of the analysis. Each row corresponds to a query sequence, and the columns provide information about the identified nucleotide and amino acid alleles and their sequence types.

## Troubleshooting <a name="troubleshooting"></a>

If you encounter any issues while running the script, check the following:

- Ensure all input files are in the correct format and location.
- Make sure you have the necessary permissions to read the input files and write to the output directory.
- If the script fails with an error message, try to understand what the message is saying. Often, the error message provides clues about what went wrong.

## Additional Resources <a name="additional-resources"></a>

For more information about the script and its functionality, refer to the [`allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/allele_find) documentation.