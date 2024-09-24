![](RATION_logo.png)

### `dsRNA_analyzer_2.py` ###

### This script reads in a list of transcripts and evaluates their suitability as RNAi targets in a target organism. It will also search for off-targets in the target organism as well as potential hits in non-target organisms (NTOs). ####

This tool has been developed under the framework of the EU-funded project [**RATION**](https://www.ration-lrp.eu/).

**NOTE**: Development of `dsRNA_analyzer_2.py` has been based on `dsRNA_analyzer.v0.1.py`, created by Panos Ioannidis. `dsRNA_analyzer.v0.1.py` is available [here](https://gitlab.com/pioannidis/ration).
---

```
Usage: dsRNA_analyzer_2.py [-h] -i INPUT -o ORGANISM_DB -t TO -f TO_GENOME -n NTO -a NTO_GENOMES [-s SIRNA_LENGTH] [-d DSRNA_LENGTH] [-m MISMATCHES] [-p THREADS]

This script reads in a list of transcripts and evaluates their suitability as RNAi targets in a target organism. It will also search for off-targets in the target organism as well as for
other non-target organisms (NTOs).

Options:
  -h, --help  show this help message and exit
  -p THREADS, --threads THREADS
            number of threads to use

required named arguments:
  -i INPUT, --input INPUT
			Path to input transcripts. This has to be a valid fasta/multi-fasta file
  -o ORGANISM_DB, --Organism_db ORGANISM_DB
            Mapping of organisms - one per line ID and scientific name, tab-separated
  -t TO, --TO TO        Target organism; Use scientific name
  -f TO_GENOME, --TO_genome TO_GENOME
            Path to target organism genome and GFF files
  -n NTO, --NTO NTO     List of NTO(s) considered in this run; Scientific names separated by semicolon(;)
  -a NTO_GENOMES, --NTO_genomes NTO_GENOMES
            Path to non-target organism genome(s)
  -s SIRNA_LENGTH, --siRNA_length SIRNA_LENGTH
            length of siRNA
  -d DSRNA_LENGTH, --dsRNA_length DSRNA_LENGTH
            length of dsRNA
  -m MISMATCHES, --mismatches MISMATCHES
            number of mismatches allowed when matching siRNAs to target and NTO genome
```
### Dependencies ###
1. `bowtie`

The Python libaries required for this script are usually installed with any Python installation.
 * `os`
 * `getopt`
 * `sys`
 * `re`

In addition, BLAST+ should be installed (v2.2.31+ was used for script development).