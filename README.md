# MIG-finder
Metagenomic Integron-associated Gene finder

### Introduction
MIG-finder is a tool developed to parse genomic data searching for integron-associated genes. The method is optimised for parsely assembled metagenomic data but can be used in whole genomes as well. Integrons are genomic mobile elements carry an integrase gene, a promoter and an array of gene cassettes. Because of the fragmented nature of (assembled) metagenomic data, MIG-finder do not look for the whole integron, instead it seaarches for their gene cassettes and require at least two to be found in the same sequence so that the false positive number is reduced.

### Method
MIG-finder is thus a bioinformatics pipeline. It starts by searching for *attC* sites using HattCI, validates these sites checking their secondary structures using Infernal, filters out single hits to reduce number of false positives and finally search for ORFs on the vicinity of the *attC* sites using Prodigal.

### Requirements
MIG-finder is a python script and requires other tools to be installed and available to run properly. You will need:

* [HattCI 1.0b+](http://github.com/maribuon/HattCI)
* [Infernal 1.1.1+](http://eddylab.org/infernal/)
* [Prodigal 2.6.2+](http://github.com/hyattpd/Prodigal)

### Instalation

XXX

### Running
```python
import migfinder as mg
mg.main(fastafile)
```

You can opt to give MIG-finder plenty of arguments:
* fastafile   an input .fasta, .fna, or .fa containing the genomic sequences that may contain integron-associated genes
* both        logical value that indicates if both the ordinary sequences and the complementary sequence should be processed \[default: True\]
* nseq        HattCI reads x sequences at a time and processes them before reading the next x sequences, in order to avoid overextending RAM. This flag gives the option to manually choose number of sequences to read, in the case of large sequences, reduce it \[default: 1000\]
* nthread     HattCI may run a large part of the computations in parallel, i.e. let different threads process a set of sequences, which in turn gives a reduced computation time. The parallelization works best when processing larger chunks of sequences at a time \[default: 6\]
* cm_model=cm_model
, k_cm=20
, k_orf=0
, save_orf=True
, dist_threshold = 4000
, d_CDS_attC = 500


