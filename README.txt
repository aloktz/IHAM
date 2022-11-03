Author: Alon Katz
year: 2022
version: 1.0

General Information
* IHAM infers host-adaptive mutations in human pathogens that emerged from animal hosts
* This script accepts single alignments only, no multiple sequence alignments
* make sure that stop codons and incomplete codons have been removed from alignments before running IHAM
* Your directory should contain alignments of one or more genes (in separate files) from one ancestral animal lineage and alignments of one or more genes from one or more descendant human lineages
* fasta files must follow the naming convention "gene-host-subtype.fasta" (case sensitive) for example env-chimpanzee-cpz. human pathogens must have "human" as the host e.g. env-human-M
* ensure that sequence names in your fasta files begin with the subtype followed by a fullstop e.g. >A.CH.03.HIV_CH_BID_V3538_2003.JQ403028

Dependencies:
* hyphy - https://stevenweaver.github.io/hyphy-site/installation/
* biopython - https://biopython.org/
* ete3 - http://etetoolkit.org/download/
* fasttree - http://www.microbesonline.org/fasttree/#Install
* python libraries



How IHAM works: 
*  alignment files in specified directory are parsed
* terminal gaps and sequences with >= 50% gaps are discarded
* a simple consesnsus sequence is calculated and sequences that aren't similar enough to it (comparison uses a user-specified similarity threshold) are discarded so that fewer than 500 sequences remain
* a tree is constructed using fasttree for use with hyphy
* FUBAR analysis is performed on the remaining sequences 

instructions for use:
1. copy the script into a folder containing one or more files of sequence alignments of a SINGLE gene/genomic region and make sure you are connected to the internet
2. use "python infer_host_adaptive_mutations.py" (TODO: change from ipynb to  py) in the terminal/console followed by args
3. args: 
    -p path: path to the folder of alignments 
    -f format: type the format of your alignment files in lower case
    -e extension: type the file extension of your alignment files in lower case
    -s similarity: type a float between 0 and 1 to set a similarity threshold for if the alignment must be reduced to a smaller set of representative sequences
   you can change how many sequences will be used for FUBAR and BUSTED by changing the int value in line 363 and line 427 respectively
