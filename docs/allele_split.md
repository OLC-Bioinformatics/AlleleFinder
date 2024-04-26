## Table of Contents
1. [Split Allele Database](#split-allele-database)
2. [Inputs](#inputs)
3. [Running the Script](#running-the-script)
4. [Usage](#usage)
5. [Outputs](#outputs)

## Split Allele Database <a name="split-allele-database"></a>

This script splits a single allele database file into multiple files; one sequence per file.

## Inputs <a name="inputs"></a>

1. Allele database file.

## Running the Script <a name="running-the-script"></a>

```bash
stec.py allele_split -q /path/to/query_folder -o /path/to_output_folder
```

An example with the query files in the `query` folder, and the desired outputs in the `split_alleles` folder (all in the current working directory).

## Usage <a name="usage"></a>

```bash
usage: stec.py allele_split [-h] [-version] [-v verbosity] [-q query_path]
                            [-o output_path]

Split combined allele files into individual files

optional arguments:
  -h, --help            show this help message and exit
  -version, --version   show program's version number and exit
  -v verbosity, --verbosity verbosity
                        Set the logging level. Options are debug, info, warning, error, and critical. Default is info.
  -q query_path, --query_path query_path
                        Specify name and path of folder containing query files in FASTA format. If not provided, the query folder in the current working directory will be used
  -o output_path, --output_path output_path
                        Specify name and path of folder into which the split allele files are to be written. If not provided, the split_alleles folder in the current working directory will be used
```

## Outputs <a name="outputs"></a>

The provided output folder will contain all sequences in the provided database files, one per file.

Additional information regarding this functionality is available in the [`allele_split`](https://olc-bioinformatics.github.io/AlleleFinder/allele_split) documentation.