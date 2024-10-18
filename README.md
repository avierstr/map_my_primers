# map_my_primers

Map_my_primers is a tool for easy visualization of your primers/barcodes/adapters on your sequences.

I have written this tool because I was puzzled why a (large) amount of reads in a barcoded amplicon experiment end up in the "Unknown" file and have barcodes only on one side, or none at all, or sometimes 2 different barcodes, and what are the reads that are longer than expected.  

**Requirements:**

Runs on Linux, Unix, Windows, Mac.

- Python 3
- edlib: Lightweight, super fast C/C++ library for sequence alignment using edit (Levenshtein) distance ([https://pypi.org/project/edlib/#description](https://pypi.org/project/edlib/#description)) (pip: `python3 -m pip install edlib` or conda: `conda install bioconda::python-edlib`) 
(for Win64: pip: `python -m pip install edlib` or conda: `conda install conda-forge::edlib`)
-   biopyton (pip: `pip install biopython` or conda: `conda install biopython` or in linux `sudo apt-get install python3-biopython`)

### Licence
GNU GPL 3.0

### Keywords
amplicon sequencing, MinION, Oxford Nanopore Technologies, barcodes, primers, adapters

### Options:

`-i, --input`: Input **file** in fastq or fasta format.  Also a **folder** can be given as input and will be scanned for .fasta or .fastq files to process.  Make sure the input file(s) is (are) named as .fasta or .fastq because it replaces the extension in parts of the script.

`-o, --outputfolder`: Save the results in the specified outputfolder. Default = "mapped".

`-min, --minlength`: Minimum readlenght to process. 

`-max, --maxlength`: Maximum readlenght to process. Default=No limit

`-pr, --primers`: File in csv format with primers/barcodes used in experiment.

`-er, --error`: Percentage error allowed in editdistance for adapters and barcodes. Default = 0.15 (15%)

### Primer file:

The file with the primers/barcodes/adapters is a csv file (comma or tab separated) with 3 columns: name, forward sequence, reverse sequence (5'-3').  It is best to provide the file not only with the barcode-primers combinations but also with the adapters, barcodes, primers so that the script can find the longest hit or only a partial hit if a complete hit is not found.
```
16S,ACGACGTTGTAGAGAGTTTGATCMTGGCTCAG,GATGGTCGATGACGGTTACCTTGTTACGACTT
LSK114,MMTGTACTTCGTTCAGTTACGTATTGCT,GCAATACGTAACTGAACGAAGTACAGG
Middle,ACTTCGTTCAGTTACGTATTGCT,
BC01_16S,TCGAAGAAAGTTGTCGGTGTCTTTGTGACGACGTTGTAGAGAGTTTGATCMTGGCTCAG,TCGAAGAAAGTTGTCGGTGTCTTTGTGGATGGTCGATGACGGTTACCTTGTTACGACTT
BC02_16S,TCGTCGATTCCGTTTGTAGTCGTCTGTACGACGTTGTAGAGAGTTTGATCMTGGCTCAG,TCGTCGATTCCGTTTGTAGTCGTCTGTGATGGTCGATGACGGTTACCTTGTTACGACTT
BC01,TCGAAGAAAGTTGTCGGTGTCTTTGTG,TCGAAGAAAGTTGTCGGTGTCTTTGTG
BC02,TCGTCGATTCCGTTTGTAGTCGTCTGT,TCGTCGATTCCGTTTGTAGTCGTCTGT
```
### Command example:
*Process "unknown.fasta", only the reads longer than 1800 bp and try to map the primers in the file "primers.csv" :*

`python3 map_my_primers.py -i unknown.fasta -min 1800 -p primers.csv`: 

### Output file:

The output file is a fasta file where the sequences are shown in line lengths of 100 bp and the primers/barcodes/adapters mapped on the reads with the name and an ">" forward or "<" reverse symbol.
```
GAGTGGGCTACACACGTGCTACAATGGTGTCTACAATGGGCTGCAAGGTGCGCAAGCCTAAGCTAATCCCTAAAAGACATCTCAGTTCGGATTGTACTCT
                                                                                                    
                                                                                                    
GCAACTCGAGTACATGAAGTTGGAATCGCTAGTAATCGTGGATCAGCATGCCACGGTGAATACGTTCTCGGGTCTTGTACACACTGCCCGTCACGCCATG
                                                                                                    
                                                                                                    
GGAATTGGTTTCACTCGAAGCTAATGGCCTAACCGCAAGGAAGGAGTTTGTTTGATTGCGATGTGACTGGGGTGAAGTCGTAACAAGGTAACCGTCATCG
                                                                          AAGTCGTAACAAGGTAACCGTCATCG
                                                                          <RV_BC17_16S              
ACCGATCATCAGAGGTACTTTCCTGGAAAG--GAACGAGTCTCTTGGGACCG-TAG-GATGGTCGATGACGGTTACCTTGTTACGACTTCACCCCAGTCA
ACC-ATCATCAGAGGTACTTTCCTGGAGGGTCGAACGAGTCTCTTGGGACCCATAGAGATGGTCGATGACGGTTACCTTGTTACGACTT           
                              FW_BC14_16S>                                                          
CTGATTCCACTTTAAATAACTCCTTCCTTGCGGTTAGGCCATTAGCTTCGAGTGAAACCAATTCCCATGGCGTGACGGGCAGTGTGTACAAGACCCGAGA
                                                                                                    
                                                                                                    
ACGTATTCACCGTGGCATGCTGATCCACGATTACTAGCGATTCCAACTTCATGTACTCGAGTTGCAGAGTACAATCCGAACTGAGATGTCTTCAGGGATT
                                                                                                    
                                                                                                    
AGCTTAGGCTTGCGCACCTTGCAACCCATTGTAGACACCATTGTAGCACGTGTGTAGCCCACTCCATAAAGGCCATGATGACTCGACATCATCCCCACCT
```


> Written with [StackEdit](https://stackedit.io/).
