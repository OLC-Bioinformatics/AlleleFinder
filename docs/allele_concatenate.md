# Table of Contents
1. [Concatenate Allele Subunits](#concatenate-allele-subunits)
2. [Inputs](#inputs)
3. [Running the Script](#running-the-script)
4. [Usage](#usage)
5. [Outputs](#outputs)

## Concatenate Allele Subunits <a name="concatenate-allele-subunits"></a>

This script concatenates allele subunits. The outputs from the STEC pipeline (by design) have separate STEC A and B 
subunits. This script concatenates the subunits in the correct
order with the appropriate linker sequence (_stx1_: 9 nt / 3 aa, _stx2_: 12 nt / 4 aa). The linker is `N` for nucleotide 
and `X` for amino acid.

#### Inputs <a name="inputs"></a>

1. Nucleotide and amino acid allele files prepared by [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce)
2. Nucleotide and amino acid profile files prepared by [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce). Note that the allele files must contain sequence for the same genes that were used for the reduction of the profile, e.g.:
    * ECs2973
    * ECs2974

#### Running the Script <a name="running-the-script"></a>

```bash
stec.py allele_concatenate 
--nt_profile /path/to/nt_profile/profile.txt
--aa_profile /path/to/aa_profile/profile.txt
--nt_alleles /path/to/nt_alleles
--aa_alleles /path/to/aa_alleles
-c /path/to/outputs
```

#### Usage <a name="usage"></a>

```bash
usage: stec.py allele_concatenate [-h] [-version] [-v verbosity]
                                  [--nt_profile nt_profile]
                                  [--aa_profile aa_profile]
                                  [--nt_alleles nt_alleles]
                                  [--aa_alleles aa_alleles]
                                  [-c concatenate_path]

Concatenate stx toxin subunit alleles with linkers

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
  -c concatenate_path, --concatenate_path concatenate_path
                        Specify name and path of folder into which concatenated subunit files are to be placed. If not provided, the concatenated_alleles folder in the current working directory will be used
```

#### Outputs <a name="outputs"></a>

The concatenated alleles will be located in the provided `concatenate_path`. Nucleotide alleles will be in the `nt` subdirectory. 
Amino acid alleles will be in the `aa` subdirectory. Files are named based on the ordered alleles: _stx1_: `stx1A_stx1B.fasta`, 
_stx2_: `stx2A_stx2B.fasta`