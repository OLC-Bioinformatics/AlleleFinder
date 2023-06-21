## Concatenate allele subunits

This script concatenates allele subunits. The outputs from the STEC pipeline (by design) have separate STEC A and B 
subunits. This script concatenates the subunits (_stx1_: ECs2974/ECs2973 and _stx2_: ECs1205/ECs1206) in the correct
order with the appropriate linker sequence (_stx1_: 9 nt / 3 aa, _stx2_: 12 nt / 4 aa). The linker is `N` for nucleotide 
and `X` for amino acid.

#### Inputs

1. nucleotide and amino acid allele files prepare by [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce)
2. nucleotide and amino acid profile files prepared by [`allele_translate_reduce`](https://olc-bioinformatics.github.io/AlleleFinder/allele_translate_reduce). Note that the allele files must contain sequence for the same genes that were used for the reduction of the profile, e.g.:
    * ECs2973
    * ECs2974

#### Running the script

```
stec.py allele_concatenate 
--nt_profile /path/to/nt_profile/profile.txt
--aa_profile /path/to/aa_profile/profile.txt
--nt_alleles /path/to/nt_alleles
--aa_alleles /path/to/aa_alleles
-c /path/to/outputs
```


#### Usage

```
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

#### Outputs

The concatenated alleles will be located in the provided `concatenate_path`. Nucleotide alleles will be in the `nt` subdirectory. 
Amino acid alleles will be in the `aa` subdirectory. Files are named based on the ordered alleles: _stx1_: `ECs2974_ECs2973.fasta`, 
_stx2_: `ECs1205_ECs1206.fasta`