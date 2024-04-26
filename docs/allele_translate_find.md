## Table of Contents
1. [Find alleles with translation](#find-alleles-with-translation)
2. [Inputs](#inputs)
3. [Running the script](#running-the-script)
4. [Usage](#usage)
5. [Outputs](#outputs)

## Find alleles with translation <a name="find-alleles-with-translation"></a>

This script performs analyses on nucleotide sequences against an amino acid database prepared by [allele_translate_reduce](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce) to discover their sequence types. It updates nucleotide and amino acid profiles and allele databases.

## Inputs <a name="inputs"></a>

1. Nucleotide query files in FASTA format
2. All outputs from [allele_translate_reduce](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce)

## Running the script <a name="running-the-script"></a>

You can run the script using the following command:

```bash
stec.py allele_translate_find --nt_profile /path/to/nt_profile_file --aa_profile /path/to_aa_profile_file --nt_alleles /path/to/nt_allele_folder --aa_alleles /path/to/aa_allele_folder -r /path/to/output_folder -q /path/to/query_folder
```

## Usage <a name="usage"></a>

The script can be run with several optional arguments:

| Argument | Description |
| --- | --- |
| `-h, --help` | Show this help message and exit |
| `-version, --version` | Show program's version number and exit |
| `-v verbosity, --verbosity verbosity` | Set the logging level. Options are debug, info, warning, error, and critical. Default is info. |
| `--nt_profile nt_profile` | Specify name and path of nucleotide profile file. If not provided, profile.txt in the nt_profile folder in the current working directory will be used by default |
| `--aa_profile aa_profile` | Specify name and path of amino acid profile file. If not provided, profile.txt in the aa_profile folder in the current working directory will be used by default |
| `--nt_alleles nt_alleles` | Specify name and path of folder containing nucleotide alleles. If not provided, the nt_allele folder in the current working directory will be used by default |
| `--aa_alleles aa_alleles` | Specify name and path of folder containing amino acid alleles. If not provided, the aa_allele folder in the current working directory will be used by default |
| `-r report_path, --report_path report_path` | Specify name and path of folder into which reports are to be placed. If not provided, the reports folder in the current working directory will be used |
| `-q query_path, --query_path query_path` | Specify name and path of folder containing query files in FASTA format. If not provided, the query folder in the current working directory will be used |

## Outputs <a name="outputs"></a>

All outputs from [allele_translate_reduce](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce) will be updated. The output includes:

* `reports`
    * `aa_novel_profiles.txt`: text file containing all novel amino acid profiles generated from query sequences
    * `aa_gene_name_novel_alleles.fasta`: FASTA-formatted text file containing all novel amino acid alleles from query sequences
    * `nt_novel_profiles.txt`: text file containing all novel nucleotide profiles generated from query sequences
    * `nt_gene_name_novel_alleles.fasta`: FASTA-formatted text file containing all novel nucleotide alleles from query sequences
    * `stec_report.tsv`: TSV file containing results for each query. Includes sample name, nucleotide allele identifiers, nucleotide sequence type, amino acid allele identifiers, amino acid sequence type, and notes

Additional information regarding this functionality is available in the [allele_translate_find](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_find) documentation.