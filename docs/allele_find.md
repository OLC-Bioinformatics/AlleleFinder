## Find alleles

This script performs BLAST analyses on a nucleotide allele database prepared by [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce) against nucleotide query sequences to discover their sequence types. Updates nucleotide and amino acid profiles and allele databases

#### Inputs

1. nucleotide query files in FASTA format
2. all outputs from [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce)

#### Running the script

```
stec.py allele_find --nt_profile /path/to/nt_profile_file --aa_profile /path/to_aa_profile_file --nt_alleles /path/to/nt_allele_folder --aa_alleles /path/to/aa_allele_folder -r /path/to/output_folder -q /path/to/query_folder
```

An example with the nucleotide profile file `profile.txt` in the `nt_profile` folder, the amino acid profile file `profile.txt` in the `aa_profile` folder, nucleotide alleles in the `nt_alleles` folder, the amino acid alleles in the `nt_alleles` folder, the query files in the `query` folder, and the desired output path for reports in the `reports` folder (all in the current working directory):

```
stec.py allele_find --nt_profile nt_profile/profile.txt --aa_profile aa_profile/profile.txt --nt_alleles nt_alleles --aa_alleles aa_alleles -q query -r reports
```

Please note that all those files and folders are the default values, so you can obtain the same command as follows:

```
stec.py allele_find
```

#### Usage

```
usage: stec.py allele_find [-h] [-version] [-v verbosity]
                           [--nt_profile nt_profile] [--aa_profile aa_profile]
                           [--nt_alleles nt_alleles] [--aa_alleles aa_alleles]
                           [-r report_path] [-q query_path]

Analyse sequences to determine allele complement. Update profiles and databases. Keep notes

optional arguments:
  -h, --help            show this help message and exit
  -version, --version   show program's version number and exit
  -v verbosity, --verbosity verbosity
                        Set the logging level. Options are debug, info, warning, error, and critical. Default is info.
  --nt_profile nt_profile
                        Specify name and path of nucleotide profile file. If not provided, profile.txt in the nt_profile folder in the current working directory will be used by default
  --aa_profile aa_profile
                        Specify name and path of amino acid profile file. If not provided, profile.txt in the aa_profile folder in the current working directory will be used by default
  --nt_alleles nt_alleles
                        Specify name and path of folder containing nucleotide alleles. If not provided, the nt_allele folder in the current working directory will be used by default
  --aa_alleles aa_alleles
                        Specify name and path of folder containing amino acid alleles. If not provided, the aa_allele folder in the current working directory will be used by default
  -r report_path, --report_path report_path
                        Specify name and path of folder into which reports are to be placed. If not provided, the reports folder in the current working directory will be used
  -q query_path, --query_path query_path
                        Specify name and path of folder containing query files in FASTA format. If not provided, the query folder in the current working directory will be used
```

#### Outputs

All outputs from [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce) will be updated

* `reports`
    * `aa_novel_profiles.txt`: text file containing all novel amino acid profiles generated from query sequences
    * `aa_gene_name_novel_alleles.fasta`: FASTA-formatted text file containing all novel amino acid alleles from query sequences
    * `nt_novel_profiles.txt`: text file containing all novel nucleotide profiles generated from query sequences
    * `nt_gene_name_novel_alleles.fasta`: FASTA-formatted text file containing all novel nucleotide alleles from query sequences
    * `stec_report.tsv`: TSV file containing results for each query. Includes sample name, nucleotide allele identifiers, nucleotide sequence type, amino acid allele identifiers, amino acid sequence type, and notes

Additional information regarding this functionality is available in the [`allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/allele_find) documentation

