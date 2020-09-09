# MIG-finder
Metagenomic Integron-associated Gene finder

### Introduction
MIG-finder is a tool developed to parse genomic data searching for integron-associated genes. The method is optimised for partly assembled metagenomic data but can be used in whole genomes as well. Integrons are genomic mobile elements carry an integrase gene, a promoter and an array of gene cassettes. Because of the fragmented nature of (assembled) metagenomic data, MIG-finder do not look for the whole integron, instead it searches for their gene cassettes and require at least two to be found in the same sequence so that the false positive number is reduced.

### Method
MIG-finder is thus a bioinformatics pipeline. It starts by searching for *attC* sites using HattCI, validates these sites checking their secondary structures using Infernal, filters out single hits to reduce number of false positives and finally search for ORFs on the vicinity of the *attC* sites using Prodigal.

### Requirements
MIG-finder is a python script and requires other tools to be installed and available to run properly. You will need:

* [HattCI 1.0b+](http://github.com/maribuon/HattCI)
* [Infernal 1.1.1+](http://eddylab.org/infernal/)
* [Prodigal 2.6.2+](http://github.com/hyattpd/Prodigal)

### Instalation
Download and uncompress the files provided here. In the directory where setup.py and migfinder directory are located run: 
```console
user@machine:~$ pip install -e migfinder
```

Alternatively, clone the repository and then:

1. unzip the folder
2. cd migfinder-master/migfinder
3. run make

### Running
Recommended to create a folder to each fastafile analysed. From this folder, start python and run:
```python
import migfinder as mf
mf.main(fastafile)
```

You can opt to give MIG-finder plenty of arguments:
* **fastafile**   an input .fasta, .fna, or .fa containing the genomic sequences that may contain integron-associated genes
* **both**        logical value that indicates if both the ordinary sequences and the complementary sequence should be processed \[default: True\]
* **nseq**        HattCI reads x sequences at a time and processes them before reading the next x sequences, in order to avoid overextending RAM. This flag gives the option to manually choose number of sequences to read, in the case of large sequences, reduce it \[default: 1000\]
* **nthread**     HattCI may run a large part of the computations in parallel, i.e. let different threads process a set of sequences, which in turn gives a reduced computation time. The parallelization works best when processing larger chunks of sequences at a time \[default: 6\]
* **cm_model**    covariance model used by Infernal to validate the *attC* site secondary structure. A model is provided with the implementation and it is recommended to be used unless you have a good reason use an alternative one.
* **k_cm**        threshold used to filter Infernal results \[default: 20\]
* **k_orf**       threshold used to filter HattCI results \[default: 0\]
* **save_orf**    save ORF results in a separate fasta file \[default: True\]
* **dist_threshold** max distance allowed to consider two adjacent *attC* sites part of the same integron \[default: 4000\]
* **d_CDS_attC**  max distance allowed between ORF and *attC* site in the same gene cassette. Used to validate the first *attC* site in the array \[default: 500\]


### Output
The output is organised in:

* **hmmresults**, dir with HattCI results
* **cmrsults**, dir with Infernal results
* **orfresults**, dir with Prodigal results, only if *attC* sites have been found.

If *attC* sites have been found, the final results will be organised in the **fastaname.results** file. The columns in this file are:
|                |        |
| ------         | ------ |
| id             | sequence id from fastafile |
| element        | if the hit is a CDS (ORF) or *attC* site|
| start          | start position of the element in the original fastafile |
| end            | end position of the element in the original fastafile |
| strand         | strand of the element in the original fastafile |
| score          | Prodigal score for ORFs, Infernal score for *attC* sites | 
| len            | length of the element |
| dist           | distance to previous element. ini0 indicates first in the array|
| array_no       | CDS_# or attC_#, where # is the number of the gene cassette in the array they belong to |
| dist_attC      | distance between *attC* sites, from the end of one to the start of another |
| Vscore         | Viterbi score given by HattCI |
| R', sp', L', loop, L", sp", R" | start position of each one of the *attC* site motifs |
