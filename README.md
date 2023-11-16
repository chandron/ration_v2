This git project serves as a repository for the scripts I have written in the frame of the EU-funded RATION project. For more information, visit https://www.ration-lrp.eu/

# Usage - dsRNA_analyzer.v0.1.py

```
Usage: dsRNA_analyzer.v2.py <fasta file containing the mRNAs> <genome sequence of the target organism (fasta file)> <genome sequences of non-target organisms (fasta file)>
```


## Overview

This script will take in a fasta file containing  mRNA sequences and evaluate their suitability as RNAi targets. More specifically, it will check whether there are short interfering RNAs (siRNAs) contained within this mRNA sequence that:
 1. have specific properties (GC content, asymmetry, no nucleotide runs)
 2. are specific for the particular mRNA in the target organism
 3. don't have off-targets in the target organism
 4. don't have off-targets in non-target organisms (NTOs)

The Python libaries required for this script are usually installed with any Python installation.
 * `os`
 * `getopt`
 * `sys`
 * `re`

In addition, BLAST+ should be installed (v2.2.31+ was used for script development).


This script is based on the DEQOR program, which can be found at
http://144.76.155.9/deqor_new/input.html

DEQOR, however, doesn't search for hits in NTOs and, more importantly,
cannot be downloaded and run locally.


## Output

Currently, there are 3 tab-delimited output files:
 1. siRNAs.all.tsv
 2. siRNAs.good.tsv
 3. dsRNAs_per_gene.tsv

The first two files show details about the siRNAs contained in each input mRNA sequence. The third file contains information about the proposed 500 bp long dsRNA sequence which is predicted for each input mRNA.


## Improvements

A number of improvements are needed even though the script is functional in its current form. Below are some examples:
 * the code would be much more readable if I had functions instead!
 * Give certain parameters (number of threads, deletion of temporary files etc) as command line arguments. Both of them are currently given as constants, which are defined inside the script.
 * Currently, the script runs 3 BLASTN searches for each input sequence. While this made programming easier, it also leads to tens of thousands of BLAST searches in the case of entire gene sets! More specifically it might take 2 days to finish a gene set with 15,000 genes (having only two NTOs). Thus, a great speed improvement would be to run only 3 BLAST searches regardless of the number of input sequences.
 * The length of the dsRNA for each input mRNA shouldn't be fixed to 500 bp (I guess). I guess that it's better for it to be dynamically determined based on the distribution of good siRNAs along the mRNA sequence.
 * A graphical environment would be great! And it wouldn't be terribly difficult to write. But it's so 2000s! It would be so much better if a webapp was written! That would make it easy for anyone to run, but it would be much more challenging to implement...


## Tutorial

To be added!


## Support
You can send me a message here, or email me.
