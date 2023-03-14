## Translate and reduce alleles

This script translates allele files from Enterobase in nucleotide format to amino acid, performs content and length checks, and removes duplicates.

In order for a translated allele to pass content and length checks, it must:

1. Start with a `Methionine` residue
2. Pass a minimum length threshold after trimming:
       * The length thresholds are:
         * ECs2973: 90 amino acid residues
         * ECs2974: 316 amino acid residues
         * ECs1205: 320 amino acid residues
         * ECs1206: 88 amino acid residues
3. Not be a duplicate of an allele already in the reduced database

#### Inputs

1. nucleotide allele files from Enterobase in FASTA format. [Download instructions.](downloads.md#download-alleles)
2. reduced profile file prepared by [`profile_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/profile_reduce). Note that the allele files must contain sequence for the same genes that were used for the reduction of the profile, e.g.:
    * ECs2973
    * ECs2974

#### Running the script

```
stec.py allele_translate_reduce -a /path/to/allele_folder -p /path/to/profile_file -r /path/to/output/folder/aa_profile -t /path/to/output/folder/aa_alleles
```

An example with the nucleotide allele folder `nt_alleles`, the profile file `profiles.txt` in the `nt_profile` folder, the desired amino acid profile output folder `aa_profile`, and the desired amino acid allele output folder `aa_alleles` all in the current working directory: 


```
stec.py allele_translate_reduce -a nt_alleles -p nt_profile/profile.txt -r aa_profile -t aa_alleles
```

#### Usage
```
usage: stec.py allele_translate_reduce [-h] [-version] [-v verbosity]
                                       [-a allele_path] [-p profile_file]
                                       [-r report_path] [-t translated_path]

Translate allele files in nucleotide format to amino acid. Remove duplicates. Keep notes

optional arguments:
  -h, --help            show this help message and exit
  -version, --version   show program's version number and exit
  -v verbosity, --verbosity verbosity
                        Set the logging level. Options are debug, info, warning, error, and critical. Default is info.
  -a allele_path, --allele_path allele_path
                        Specify name and path of folder containing allele files. If not provided, the nt_alleles folder in the current working directory will be used by default
  -p profile_file, --profile_file profile_file
                        Optionally specify name and path of profile file. Parse the nucleic acid profile, and create the corresponding reduced amino acid profile
  -r report_path, --report_path report_path
                        Specify the name and path of the folder into which outputs are to be placed. If not provided, the aa_profile folder in the current working directory will be used
  -t translated_path, --translated_path translated_path
                        Specify the name and path of the folder into which alleles are to be placed. If not provided, the aa_alleles folder in the current working directory will be used
```

#### Outputs

These represent the folder structure of the directory containing the outputs

* `aa_alleles`
    * `notes`
        * `gene_name_notes.txt`: text file containing the nucleotide allele, the corresponding amino acid allele, and any notes (whether it is a duplicate, trimmed, or filtered)
    * `gene_name.fasta`: FASTA-formatted text file containing sequences of all amino acid alleles that passed quality/content filters
    * `gene_name_filtered.txt`: FASTA-formatted text file containing sequences of all amino acid alleles that failed quality/content filters
* `aa_profile`
    * `profile`
        * `profile.txt`: text file containing all unique amino acid profiles
        * `reducing_notes.txt`: text file containing all amino acid sequence types, the sequence type following duplicate reduction, and any notes
    * `aa_full_profile.txt`: text file containing all amino acid profiles before duplicate removal
    * `aa_nt_profile_links.tsv`: TSV file containing every amino acid sequence type and the corresponding nucleotide sequence type(s)
    * `gene_names.txt`: text file containing the extracted gene names from the analysis
    * `profile.txt`: copy of `profile.txt` from the `profile` folder. Allows for allele discovery without manually moving the file
* `nt_alleles`
    * `gene_name.fasta`: FASTA-formatted text file containing sequences of all nucleotide alleles that passed quality/content filters
    * `gene_name_filtered.txt`: FASTA-formatted text file containing sequences of all nucleotide alleles that failed quality/content filters
* `nt_profile`
    * `profile.txt`: text file containing all unique nucleotide profiles