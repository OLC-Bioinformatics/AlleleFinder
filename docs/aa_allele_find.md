# Table of Contents
1. [Find Alleles from Amino Query Files](#find-alleles-from-amino-query-files)
2. [Inputs](#inputs)
3. [Running the Script](#running-the-script)
4. [Usage](#usage)
5. [Outputs](#outputs)

## Find Alleles from Amino Query Files <a name="find-alleles-from-amino-query-files"></a>

This script performs BLAST analyses on an amino acid database prepared by [`allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/allele_find) against amino acid query sequences to find matching alleles. Updates allele database.

#### Inputs <a name="inputs"></a>

1. Amino acid query files in FASTA format. One query allele per file. Note that the allele naming scheme must match the outputs from the previous scripts
2. Amino acid allele database prepared by [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce) or [`allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/allele_find)

#### Running the Script <a name="running-the-script"></a>

```bash
stec.py aa_allele_find --aa_alleles /path/to/aa_allele_folder -r /path/to/output_folder -q /path/to/query_folder
```

An example with the amino acid alleles in the `aa_alleles` folder, the amino acid query files in the `query` folder, the desired reports folder `reports`, and a cutoff value of `100` (all in the current working directory):

```bash
stec.py aa_allele_find --aa_alleles aa_alleles -q query -r reports -c 100
```

#### Usage <a name="usage"></a>

```bash
usage: stec.py aa_allele_find [-h] [-version] [-v verbosity]
                              [--aa_alleles aa_alleles] [-r report_path]
                              [-q query_path] [-c cutoff]

Analyse amino acid sequences to determine allele complement. Update profiles and databases. Keep notes

optional arguments:
  -h, --help            show this help message and exit
  -version, --version   show program's version number and exit
  -v verbosity, --verbosity verbosity
                        Set the logging level. Options are debug, info, warning, error, and critical. Default is info.
  --aa_alleles aa_alleles
                        Specify name and path of folder containing amino acid alleles. If not provided, the aa_allele folder in the current working directory will be used by default
  -r report_path, --report_path report_path
                        Specify name and path of folder into which reports are to be placed. If not provided, the reports folder in the current working directory will be used
  -q query_path, --query_path query_path
                        Specify name and path of folder containing query files in FASTA format. If not provided, the query folder in the current working directory will be used
  -c cutoff, --cutoff cutoff
                        Specify the percent identity cutoff for matches. Allowed values are between 90 and 100. Default is 100
```

#### Outputs <a name="outputs"></a>

All amino acid allele files will be updated

* `reports`
    * `allele_report.tsv`: TSV file containing results for each query. Includes sample name, matching allele(s), notes
    * `gene_name_filtered_alleles.fasta`: FASTA-formatted text file containing all filtered query alleles
    * `gene_name_novel_alleles.fasta`: FASTA-formatted text file containing all novel alleles