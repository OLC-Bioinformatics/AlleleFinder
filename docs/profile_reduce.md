# Table of Contents
1. [Reduce Profiles](#reduce-profiles)
2. [Inputs](#inputs)
3. [Running the Script](#running-the-script)
4. [Usage](#usage)
5. [Outputs](#outputs)

## Reduce Profiles <a name="reduce-profiles"></a>

This script reduces full wgMLST profile from Enterobase using genes of interest. 

The two _stx_ genes, _stx1_ and _stx2_, have the following identifiers in Enterobase:

* stx1 subunit A: **ECs2974**
* stx1 subunit B: **ECs2973**
* stx2 subunit A: **ECs1205** 
* stx2 subunit B: **ECs1206**

#### Inputs <a name="inputs"></a>

In order to extract all the unique profiles from a full Enterobase wgMLST profile for both _stx1_ subunits, create a text
file containing the two identifiers (one per row) e.g.:

`genes.txt`

```bash
ECs2974
ECs2973
```

A full _Escherichia_  wgMLST profile file from Enterobase is also required. [Download instructions.](downloads.md#download-profile)

#### Running the Script <a name="running-the-script"></a>

```bash
stec.py profile_reduce -p /path/to/profile_file -g /path/to/genes_file -o /path/to/output_folder
```

An example with the profile file, `profiles.list` in the `profile` folder, the genes file, `genes.txt`, in the `genes` folder, the reduced profile to be written to the `nt_profile` folder (all in the current working directory):

```bash
stec.py profile_reduce -p profiles/profiles.list -g genes/genes.txt -o nt_profile
```

#### Usage <a name="usage"></a>

```bash
usage: stec.py profile_reduce [-h] [-version] [-v verbosity] [-p profile_file]
                              [-g gene_names] [-o output_folder]

Reduce full wgMLST profile from Enterobase using genes of interest

optional arguments:
  -h, --help            show this help message and exit
  -version, --version   show program's version number and exit
  -v verbosity, --verbosity verbosity
                        Set the logging level. Options are debug, info, warning, error, and critical. Default is info.
  -p profile_file, --profile_file profile_file
                        Specify name and path of profile file. If not provided, the default "profiles.list" in the current working directory will be used
  -g gene_names, --gene_names gene_names
                        Name and path of text file containing gene names to use to filter the profile file (one per line). If not provided, the default "genes.txt" in the current working directory will be used. If the file does not exist, the program will attempt to create a file using the .fasta files in the current working directory
  -o output_folder, --output_folder output_folder
                        Name and path of folder into which the reduced profile and notes are to be placed. If not provided, the default "nt_profile" folder in the current working directory will be used
```

#### Outputs <a name="outputs"></a>

The reduced profile will be written to `profile.txt` in the supplied output folder. It contains all the unique profiles extracted from the full profile file

A notes file will be written to `reducing_notes.txt` in the supplied output folder. It contains notes on every sequence type processed from the full profile. 
If a profile is a duplicate of a previous profile, the `ReducedSequenceType` will be `0`, and the `Notes` column will note that the profile is a `duplicate`