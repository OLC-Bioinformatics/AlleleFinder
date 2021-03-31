# AlleleFinder

Allele finding has been simplified in the most recent update!

Simply use `allele_updater.py`, and provide it with the `-p` 
(path) argument. Within the supplied path, you must have two 
folders: `query`, with genome files to profile, and `alleles`, 
with one or more allele files to profile. The allele files must 
be named based on the gene name, and contain no underscores 
e.g. MutS.fasta. The alleles in the gene-specific file must be 
named based on the gene name, and the allele number e.g. 
MutS_0

The alleles created by this can be subtyped as follows:
```
primer_finder.py identity -s /path/to/alleles -pf
/path/to/genemethods/genemethods/assemblypipeline/primers.txt
-m 3 -f fasta -e -cb -a vtyper
```