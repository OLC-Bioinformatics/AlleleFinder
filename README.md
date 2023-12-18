# AlleleFinder

[![CircleCI](https://circleci.com/gh/OLC-Bioinformatics/AlleleFinder/tree/main.svg?style=shield)](https://circleci.com/gh/OLC-LOC-Bioinformatics/AlleleFinder/tree/main)
[![codecov](https://codecov.io/gh/OLC-Bioinformatics/AlleleFinder/branch/main/graph/badge.svg?token=Z6SSEJV9GU)](https://codecov.io/gh/OLC-Bioinformatics/AlleleFinder)
[![Anaconda-Server Badge](https://img.shields.io/badge/install%20with-conda-brightgreen)](https://anaconda.org/olcbioinformatics/allelefinder)
[![GitHub version](https://badge.fury.io/gh/olc-bioinformatics%2Fallelefinder.svg)](https://badge.fury.io/gh/olc-bioinformatics%2Fallelefinder)
[![GitHub issues](https://img.shields.io/github/issues/OLC-Bioinformatics/AlleleFinder)](https://github.com/OLC-LOC-Bioinformatics/AlleleFinder/issues)
[![Documentation Status](https://readthedocs.org/projects/pip/badge/?version=stable)](https://OLC-Bioinformatics.github.io/AlleleFinder/?badge=stable)
[![license](https://img.shields.io/badge/license-MIT-brightgreen)](https://github.com/OLC-Bioinformatics/AlleleFinder/blob/main/LICENSE)



### STEC AlleleFinder

A suite of tools, written in Python, designed for the discovery, sequence typing, and profiling the _stx_ alleles in Shiga toxin-producing _Escherichia coli_  

## Scripts

There is a single STEC script with six separate functionalities:

1. [`profile_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/profile_reduce)
2. [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce)
3. [`allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/allele_find)
4. [`aa_allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/aa_allele_find)
5. [`allele_split`](https://olc-bioinformatics.github.io/AlleleFinder/allele_split)
6. [`allele_concatenate`](https://olc-bioinformatics.github.io/AlleleFinder/allele_concatenate)


Full documentation is available at the [AlleleFinder GitHub pages site](https://olc-bioinformatics.github.io/AlleleFinder/)

## Quick Start

[`Conda`](https://docs.conda.io/en/latest/) is required to install AlleleFinder. See the [documentation](http://bioconda.github.io/) or [AlleleFinder installation](https://olc-bioinformatics.github.io/AlleleFinder/install/) for instructions of getting conda installed on your system


Create a new conda environment:

```
conda create -n allele_finder -c olcbioinformatics allelefinder=0.1.5=py_0
```

Additional documentation is available [here](https://olc-bioinformatics.github.io/AlleleFinder/installation)


## Reduce profiles

This script reduces full wgMLST profile from Enterobase using genes of interest. 

The two _stx_ genes, _stx1_ and _stx2_, have the following identifiers in Enterobase:

* stx1 subunit A: **ECs2974**
* stx1 subunit B: **ECs2973**
* stx2 subunit A: **ECs1205** 
* stx2 subunit B: **ECs1206**


#### Inputs
In order to extract all the unique profiles from a full Enterobase wgMLST profile for both _stx1_ subunits, create a text
file containing the two identifiers (one per row) e.g.:

`genes.txt`

```
ECs2974
ECs2973
```

A full _Escherichia_  wgMLST profile file from Enterobase is also required

#### Running the script

```
stec.py profile_reduce -p /path/to/profile_file -g /path/to/genes_file -o /path/to/output_folder
```

Additional information regarding this functionality is available in the [`profile_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/profile_reduce) documentation.


## Translate and reduce alleles

This script translates allele files from Enterobase in nucleotide format to amino acid, performs content and length checks, and removes duplicates.

In order for a translated allele to pass content and length checks, it must:

1. Start with a `Methionine` residue
2. Pass a minimum length threshold after trimming:
       * The length thresholds are:
         * ECs2973: 82 amino acid residues
         * ECs2974: 313 amino acid residues
         * ECs1205: 313 amino acid residues
         * ECs1206: 84 amino acid residues

3. Not be a duplicate of an allele already in the reduced database

#### Inputs

1. nucleotide allele files from Enterobase in FASTA format. 
2. reduced profile file prepared by [`profile_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/profile_reduce). Note that the allele files must contain sequence for the same genes that were used for the reduction of the profile, e.g.:
    * ECs2973
    * ECs2974

#### Running the script

```
stec.py allele_translate_reduce -a /path/to/allele_folder -p /path/to/profile_file -r /path/to/output/folder/aa_profile -t /path/to/output/folder/aa_alleles
```

Additional information regarding this functionality is available in the [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce) documentation


## Find alleles

This script performs BLAST analyses on a nucleotide allele database prepared by [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce) against nucleotide query sequences to discover their sequence types. Updates nucleotide and amino acid profiles and allele databases

#### Inputs

1. nucleotide query files in FASTA format
2. all outputs from [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce)

#### Running the script

```
stec.py allele_find --nt_profile /path/to/nt_profile_file --aa_profile /path/to_aa_profile_file --nt_alleles /path/to/nt_allele_folder --aa_alleles /path/to/aa_allele_folder -r /path/to/output_folder -q /path/to/query_folder
```

Additional information regarding this functionality is available in the [`allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/allele_find) documentation


## Find alleles from amino query files

This script performs BLAST analyses on an amino acid database prepared by [`allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/allele_find) against amino acid query sequences to find matching alleles. Updates allele database.

#### Inputs

1. amino acid query files in FASTA format. One query allele per file. Note that the allele naming scheme must match the outputs from the previous scripts
2. amino acid allele database prepared by [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce) or [`allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/allele_find)

#### Running the script

```
stec.py aa_allele_find --aa_alleles /path/to/aa_allele_folder -r /path/to/output_folder -q /path/to/query_folder
```

Additional information regarding this functionality is available in the [`aa_allele_find`](https://olc-bioinformatics.github.io/AlleleFinder/aa_allele_find) documentation


## Split allele database

This script splits a single allele database file into multiple files; one sequence per file

#### Inputs

1. allele database file

#### Running the script

```
stec.py allele_split -q /path/to/query_folder -o /path/to_output_folder
```

Additional information regarding this functionality is available in the [`allele_split`](https://olc-bioinformatics.github.io/AlleleFinder/allele_split) documentation


## Concatenate allele database

This script concatenates alleles of the _stx_ A and B subunits into a single sequence with a linker

1. nucleotide and amino acid allele files prepare by [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce)
2. nucleotide and amino acid profile files prepared by [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce). Note that the allele files must contain sequence for the same genes that were used for the reduction of the profile, e.g.:

#### Running the script

```
stec.py allele_concatenate --nt_profile /path/to/nt_profile/profile.txt --aa_profile /path/to/aa_profile/profile.txt --nt_alleles /path/to/nt_alleles --aa_alleles /path/to/aa_alleles -c /path/to/outputs
```

## Feedback

If you encounter any issues installing or running AlleleFinder, have feature requests, or need assistance, please [open an issue on GitHub](https://github.com/OLC-Bioinformatics/AlleleFinder/issues/new/choose)


## License

MIT License

Copyright (c) Government of Canada 2023

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: 

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
