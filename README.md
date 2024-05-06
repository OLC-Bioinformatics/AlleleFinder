# AlleleFinder

[![CircleCI](https://circleci.com/gh/OLC-Bioinformatics/AlleleFinder/tree/main.svg?style=shield)](https://circleci.com/gh/OLC-LOC-Bioinformatics/AlleleFinder/tree/main)
[![codecov](https://codecov.io/gh/OLC-Bioinformatics/AlleleFinder/branch/main/graph/badge.svg?token=Z6SSEJV9GU)](https://codecov.io/gh/OLC-Bioinformatics/AlleleFinder)
[![Anaconda-Server Badge](https://img.shields.io/badge/install%20with-conda-brightgreen)](https://anaconda.org/olcbioinformatics/allelefinder)
[![GitHub Release](https://img.shields.io/github/v/release/OLC-Bioinformatics/AlleleFinder?display_name=release&label=version&color=%20dark%20green
)](https://github.com/OLC-Bioinformatics/AlleleFinder/releases)
[![GitHub issues](https://img.shields.io/github/issues/OLC-Bioinformatics/AlleleFinder)](https://github.com/OLC-LOC-Bioinformatics/AlleleFinder/issues)
[![Documentation Status](https://readthedocs.org/projects/pip/badge/?version=stable)](https://OLC-Bioinformatics.github.io/AlleleFinder/?badge=stable)
[![license](https://img.shields.io/badge/license-MIT-brightgreen)](https://github.com/OLC-Bioinformatics/AlleleFinder/blob/main/LICENSE)

## Overview

AlleleFinder is a Python-based suite of tools designed for the discovery, sequence typing, and profiling of _stx_ alleles in Shiga toxin-producing _Escherichia coli_ (STEC). It provides a comprehensive solution for researchers and professionals working in the field of bacterial genomics.

## Features

AlleleFinder offers seven main functionalities, each represented by a separate script:

1. [`profile_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/profile_reduce): Reduces full wgMLST profile from Enterobase using genes of interest.
2. [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce): Translates allele files from Enterobase in nucleotide format to amino acid, performs content and length checks, and removes duplicates.
3. [`allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/allele_find): Performs BLAST analyses on a nucleotide allele database against nucleotide query sequences to discover their sequence types.
4. [`allele_translate_find`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_find): Performs BLAST analyses on an amino acid database against amino acid query sequences to find matching alleles.
5. [`aa_allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/aa_allele_find): Performs BLAST analyses on an amino acid database against amino acid query sequences to find matching alleles.
6. [`allele_split`](https://olc-bioinformatics.github.io/AlleleFinder/allele_split): Splits a single allele database file into multiple files; one sequence per file.
7. [`allele_concatenate`](https://olc-bioinformatics.github.io/AlleleFinder/allele_concatenate): Concatenates alleles of the _stx_ A and B subunits into a single sequence with a linker.

## Documentation

Detailed documentation for each script and general usage of AlleleFinder is available at the [AlleleFinder GitHub pages site](https://olc-bioinformatics.github.io/AlleleFinder/).

## Quick Start

[`Conda`](https://docs.conda.io/en/latest/) is required to install AlleleFinder. See the [documentation](http://bioconda.github.io/) or [AlleleFinder installation](https://olc-bioinformatics.github.io/AlleleFinder/install/) for instructions on getting conda installed on your system.

Create a new conda environment:

```bash
conda create -n allele_finder -c olcbioinformatics allelefinder
```

Additional installation instructions are available [here](https://olc-bioinformatics.github.io/AlleleFinder/installation).

## Usage

Each script has specific input requirements and usage instructions. Please refer to the respective documentation for each script for detailed information:

- [`profile_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/profile_reduce)
- [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce)
- [`allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/allele_find)
- [`allele_translate_find`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_find)
- [`aa_allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/aa_allele_find)
- [`allele_split`](https://olc-bioinformatics.github.io/AlleleFinder/allele_split)
- [`allele_concatenate`](https://olc-bioinformatics.github.io/AlleleFinder/allele_concatenate)

## Feedback

We welcome your feedback. If you encounter any issues installing or running AlleleFinder, have feature requests, or need assistance, please [open an issue on GitHub](https://github.com/OLC-Bioinformatics/AlleleFinder/issues/new/choose).

## License

AlleleFinder is licensed under the MIT License. See the [LICENSE](https://github.com/OLC-Bioinformatics/AlleleFinder/blob/main/LICENSE) file for more details.
