This git project serves as a repository for the scripts I have written in the frame of the EU-funded RATION project. For more information, visit https://www.ration-lrp.eu/

# *Script 1:* `dsRNA_analyzer.v0.1.py`

```
Usage: dsRNA_analyzer.v0.1.py <transcripts (fasta file)> <target genome (fasta file)> <non-target genome(s) (fasta file)>
```


## Summary

This script reads in a list of transcripts and evaluates their suitability as RNAi targets.


## Dependencies

The Python libaries required for this script are usually installed with any Python installation.
 * `os`
 * `getopt`
 * `sys`
 * `re`

In addition, BLAST+ should be installed (v2.2.31+ was used for script development).

This script is based on the DEQOR program, which can be found at
http://144.76.155.9/deqor_new/input.html


## Novelty

There are reports in the literature about other programs that perform the same task, such as DEQOR (Henschel et al 2004) and OfftargetFinder (Good et al 2016). However, neither one of them seems to be currently working; a website for DEQOR does exist (see below) but it's not really functional (entering a sequence and hitting the button doesn't seem to do something). In addition, DEQOR doesn't check for off targets in non-target organisms (NTOs). As for the latter it's no longer available after personal communication I personally had with one of the authors of the paper.


## Overview

At the technical level, this script will take in a fasta file containing mRNA sequences and calculate a few attributes in order to evaluate their suitability as RNAi targets. More specifically, it will check whether there are short interfering RNAs (siRNAs) contained within this mRNA sequence that:
 1. have specific properties (GC content, asymmetry, no nucleotide runs)
 2. are specific for the particular mRNA in the target organism
 3. don't have off-targets in the target organism
 4. don't have off-targets in non-target organisms (NTOs)

The length of the siRNAs is set to 21-nt and the properties in (1) are calculated as in the DEQOR program.

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

In order to run this tutorial you'll have to:
 1. have access to a Linux machine (workstation/server)
 2. know how to run stuff on the terminal

Even though it's theoretically possible to this in Windows machines, you'll still have to run it in the Windows command line. I haven't tested though on Windows (and won't test it in the forseable future!) so you're on your own there regarding any pecularities of the Windows file system and any dependencies.

The dependencies for Linux are:
 1. the `datasets` tool that can be downloaded from NCBI at https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
 2. the standalone version of BLAST+, which can be downloaded from NCBI at https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
 3. Perl and Python (v3.X). Both of these are pre-installed in all Linux distributions.
 4. the `unzip` tool, which is also pre-installed.

To examine all predicted transcripts in the Colorado Potato Beetle, Leptinotarsa decemlineata do the following:

First you'll have to download the transcripts and genome sequence from NCBI. Even though you can go to the FTP server and download it directly, I'm using the `datasets` tool from NCBI in order to download everything using one command. So create a working directory `cd` into it and issue the following commands

```
mkdir genome

cd genome

datasets download genome accession GCF_000500325.1 --include gff3,rna,cds,protein,genome,seq-report --filename GCF_000500325.1.zip

unzip GCF_000500325.1.zip

mv ncbi_dataset/data/GCF_000500325.1/* .
```

`cd` out of the `genome` directory
```
cd ..
```

Symlink the genome sequence of L. decemlineata
```
ln -s genome/GCF_000500325.1_Ldec_2.0_genomic.fna .
```

Symlink the CDS of the predicted genes
```
ln -s genome/cds_from_genomic.fna .
```

Modify the previous file so that the fasta headers are simpler
```
perl -pe 's/^>lcl\|\S+_(cds_\S+) \S+.+$/>$1/' cds_from_genomic.fna > cds_from_genomic.mod.fna
```

Generate a fasta file containing the genome sequences of any NTOs you want. For this tutorial I'll be using the honey bee (GCF_003254395.2) and the human (GCF_000001405.40) genomes. Using the Refseq accession numbers download them to the `NTO_genomes` directory and merge the genome sequences into one fasta 
```
cd NTO_genome

datasets download genome accession GCF_003254395.2 --include genome --filename GCF_003254395.2.zip

datasets download genome accession GCF_000001405.40 --include genome --filename GCF_000001405.40.zip

unzip GCF_000001405.40.zip
unzip GCF_000001405.40.zip

mv ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna .
mv ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna .
```

Delete the `ncbi_dataset` directory after moving the above two files to the current directory.

The genomes of the two NTOs should be put into one fasta file. The same should be done if more NTOs are used (all NTO genomes into one fasta).

```
cat GCF*fna > NTOs.fna
```

Move to the root directory and symlink the NTO genomes
```
ln -s NTO_genome/NTOs.fna .
```

Finally, run `dsRNA_analyzer.py` like this
```
dsRNA_analyzer.py cds_from_genomic.mod.fna GCF_000500325.1_Ldec_2.0_genomic.fna NTOs.fna > stdout_stderr 2>&1 &
```

The program will take about 2 days to complete! This is mostly due to the large number of BLAST runs (as many as the input transcripts - also see Improvements above).

## Support
You can send me a message here, or email me.
