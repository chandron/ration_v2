#!/usr/bin/python3

import os
import sys
import argparse
import re
import pandas as pd
import subprocess
import pysam

parser = argparse.ArgumentParser(description='This script reads in a list of transcripts and evaluates their suitability as RNAi targets in a target organism. It will also search for off-targets in the target organism as well as for other non-target organisms (NTOs).')
# https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument('-i', '--input', help='path to input transcripts', required=True)
requiredNamed.add_argument('-g', '--genome', help='path to target genome', required=True)
requiredNamed.add_argument('-t', '--gff', help='GFF file of target genome', required=True)
requiredNamed.add_argument('-n', '--NTO_genome', help='path to NTO genome(s)', required=True)
requiredNamed.add_argument('-a', '--NTO_ids', help='path to NTO ID(s)', required=True)
requiredNamed.add_argument('-s', '--siRNA_length', type=int, default=20, help='length of siRNA')
requiredNamed.add_argument('-d', '--dsRNA_length', type=int, default=500, help='length of dsRNA')
requiredNamed.add_argument('-m', '--mismatches', type=int, default=2, help='number of mismatches allowed when matching siRNAs to genome')
parser.add_argument('-p', '--threads', type=int, default=8, help='number of threads to use')

try:
	args = parser.parse_args()
except argparse.ArgumentError:
	print("missing options")
	sys.exit(1)

mRNAs = args.input
genome = args.genome
NTOs   = args.NTO_genome
NTO_ids = args.NTO_ids
siRNA_len= args.siRNA_length
GFF = args.gff
ds_len = args.dsRNA_length
mis = args.mismatches
CPUS = args.threads

#########################################
def generate_siRNAs(sequence, si_length):
    siRNA_sequences = []
    for i in range(0, len(sequence) - si_length + 1):
        kmer = sequence[i:i+si_length]
        siRNA_sequences.append(kmer)
    return siRNA_sequences

def run_bowtie1(mis, index_prefix, siRNA_length, sam_out):

    bowtie1_command = [
        "bowtie",
		"-f",
		"-S",
        "-n", str(mis),
        "-l", str(siRNA_length),
        "-a",
		"--no-unal",
		# "--best",
		# "--strata",
        "-x", index_prefix,
        "siRNAs.fa",
        os.path.join(os.getcwd(), sam_out)
    ]

    try:
        subprocess.run(bowtie1_command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {str(e)}")

    return os.path.join(os.getcwd(), sam_out)
#########################################
DELETE_TMP = True  # delete temporary file
#########################################
# Import the genomic.gff file
# can use that for finding specificity on target gene rather than blasting
gff = pd.read_csv(GFF, sep='\t', header=None, comment='#')
gff_cds = gff[gff[2] == "CDS"].copy()
gff_cds[9] = gff_cds[8].replace(to_replace=r'^.+Name=([^;]+);.+$', value=r'\1',regex=True)
gff_cds = gff_cds[[9, 3, 4, 6, 0]].rename(columns={9:'CDS_name', 3:'start', 4:'end', 6:'strand', 0:'chromosome'}).reset_index(drop=True)

#############################

# open the fasta file containing the dsRNA sequences
fh1 = open ( mRNAs )

line = fh1.readline()

fasta = {} # this will hold the sequences of the fasta file

a = True

sys.stderr.write( "Reading the fasta file containing the mRNAs\n" )

while a:
	line = line.strip()
	if line.startswith( ">" ):
		header = line[1:] # remove the leading ">"
		
		sequence = ""
		line = fh1.readline()
		while not line.startswith( ">" ):
			line = line.strip()
			sequence += line
			line = fh1.readline()
			if not line:
				a = False
				break
		fasta[header] = sequence

fh1.close()

sys.stderr.write( "Finished reading the dsRNA fasta file\n" )

# open the file where you'll write the dsRNA sequence (and other details)
fhout_dsRNA = open( "dsRNAs_per_gene.tsv", "w")
fhout_dsRNA.write( "TranscriptID\tbest_dsRNA_start\tbest_dsRNA_stop\tTranscript_length\tsiRNAs_specific_to_target_gene\tsiRNAs_targeting_NTOs\tdsRNA_sequence\n" )

# open the output files
fhgood = open ( "siRNAs.good.tsv", "w" )
fhall  = open ( "siRNAs.all.tsv", "w" )
fhbad  = open ( "siRNAs.bad.tsv", "w" )
fhplainbad  = open ( "siRNAs.plainbad.tsv", "w" )

# print the header of the output
fhgood.write( "siRNA_name\tsequence\tQC_asymmetry\tQC_nucleotide_runs\tQC_GC_content\tQC_specificity\tQC_offtargets\tQC_NTO_offtargets\n" )
fhall.write( "siRNA_name\tsequence\tQC_asymmetry\tQC_nucleotide_runs\tQC_GC_content\tQC_specificity\tQC_offtargets\tQC_NTO_offtargets\n" )
fhbad.write( "siRNA_name\tsequence\tQC_asymmetry\tQC_nucleotide_runs\tQC_GC_content\tQC_specificity\tQC_offtargets\tQC_NTO_offtargets\n" )
fhplainbad.write( "siRNA_name\tsequence\tQC_asymmetry\tQC_nucleotide_runs\tQC_GC_content\tQC_specificity\tQC_offtargets\tQC_NTO_offtargets\n" )

# Now loop through the fasta dictionary and analyze each sequence
#sys.stderr.write( "\nAnalyzing each gene separately\n" )
for gene in fasta:
	sys.stderr.write( "\nAnalyzing gene: " + gene + "...\n" )
	
	tmp_file = gene + ".fa"
	fhout = open ( tmp_file, "w" )
	fhout.write( ">" + gene + "\n" )
	fhout.write( fasta[gene] + "\n" )
	fhout.close()
	
	s_coords = []  # these are the subject (genome) coordinates
	
	q_len = len ( fasta[gene] )      # that's the length of the input mRNA
	
	hit_len = 0   
	
	
	SIRNA_LEN = int(siRNA_len)  # The default was intially set to 21. Think twice before you change it!
	siRNA_sequences = generate_siRNAs(fasta[gene], SIRNA_LEN)
	
	sys.stderr.write( "Examining all " + str(SIRNA_LEN) + "-nt sequences...\n" )
	
	properties = {}
	
	fhout = open ( "siRNAs.fa", "w" )

	for i, siRNA in enumerate(siRNA_sequences):
		
		siRNA_name = gene + "_" + str(i)
		
		properties[siRNA_name] = {}
		properties[siRNA_name]["sequence"] = siRNA
		
		properties[siRNA_name]["qc_asymmetry"] = 0 # the 5' should be A/T and the 3' should be G/C
		if ( siRNA[0] == 'A' or siRNA[0] == 'T' ) and ( siRNA[-1] == 'G' or siRNA[-1] == "C" ):
			properties[siRNA_name]["qc_asymmetry"] = 1
		
		properties[siRNA_name]["qc_nt_runs"] = 0 # no more than 3 identical consecutive nucleotides
		res = re.findall( "AAA|TTT|GGG|CCC", siRNA )
		if len(res) > 0:
			properties[siRNA_name]["qc_nt_runs"] = 1
			
		properties[siRNA_name]["qc_gc_content"] = 0  # should be between 20-50%
		gc_percent = ( siRNA.count("G") + siRNA.count("C") ) / SIRNA_LEN
		if gc_percent > 0.2 and gc_percent < 0.5:
			properties[siRNA_name]["qc_gc_content"] = 1
		
		# also initialize the BLAST-related properties for every siRNA
		properties[siRNA_name]["qc_specificity"] = 0
		properties[siRNA_name]["qc_offtargets"] = 0
		properties[siRNA_name]["qc_nto_offtargets"] = 0
		
		# write the sequence of this siRNA in the output
		# fasta file so that you can run the BLAST searches
		fhout.write( ">" + siRNA_name + "\n" + siRNA + "\n" )

	fhout.close()

	###############################
	## Get chromosome and start-end ranges for the particular transcript
	chrom = set()
	ranges = []
	all_ranges = []
	
	chrom = set(gff_cds[gff_cds['CDS_name'] == gene]['chromosome'].unique())
	# Create a list with all transcript positions
	start_end_pos = gff_cds[gff_cds['CDS_name'] == gene][['start', 'end']]
	start_end_pos['range'] = [list(range(i, j+1)) for i, j in start_end_pos.values]
	ranges = list(start_end_pos['range'].values)
	for sublist in ranges:
		all_ranges.extend(sublist)
	##############################
	# bowtie1 target organism index
	BT_IDX = os.path.splitext(genome)[0]+'_bowtie_idx'

	# Align siRNAs to genome
	sys.stderr.write( "\nAligning to target organism(s) ...\n" )

	sam_file_to = run_bowtie1(mis, BT_IDX, SIRNA_LEN, gene + ".sam")
	
	with pysam.AlignmentFile(sam_file_to, "r") as sam_to:
		for read in sam_to.fetch():
			if not read.is_unmapped:
				siRNA_name = read.query_name
				cds_name	= re.sub(r'_\d+$', '', siRNA_name)
				siRNA_coord = read.reference_start
				# overlap = read.get_overlap(0, SIRNA_LEN)
				mismatches = read.get_tag("NM")
				ref_name = read.reference_name
				aln_length = read.reference_length
				query_sequence = read.query_sequence
				perfect_matches = len(query_sequence) - mismatches
				strand = "reverse" if read.is_reverse else "forward"

				if mismatches <= mis:
					## Examine specificity based on GFF file
					if (ref_name in chrom) and (siRNA_coord in all_ranges):
						# is specific for target
						properties[siRNA_name]["qc_specificity"] = 1
					
					if (ref_name in chrom) and (not siRNA_coord in all_ranges):
						# is off target
						properties[siRNA_name]["qc_offtargets"] = 1
					
					if not ref_name in chrom:
						# is off target: alignment on wrong chromosome
						properties[siRNA_name]["qc_offtargets"] = 1


	#######################################
	# NTO bowtie index
	BT_IDX_NTO = os.path.splitext(NTOs)[0]+'_bowtie_idx'

	# Align with bowtie1 against NTOs
	sys.stderr.write( "\nAligning to NTO(s) ...\n" )

	sam_file_nto = run_bowtie1(mis, BT_IDX_NTO, SIRNA_LEN, gene + "_nto.sam")

	with pysam.AlignmentFile(sam_file_nto, "r") as sam_nto:
		for read in sam_nto.fetch():
			if not read.is_unmapped:
				siRNA_name = read.query_name
				cds_name	= re.sub(r'_\d+$', '', siRNA_name)
				siRNA_coord = read.reference_start
				# overlap = read.get_overlap(0, SIRNA_LEN)
				mismatches = read.get_tag("NM")
				ref_name = read.reference_name
				aln_length = read.reference_length
				query_sequence = read.query_sequence
				perfect_matches = len(query_sequence) - mismatches
				strand = "reverse" if read.is_reverse else "forward"

				if mismatches <= mis:
					properties[siRNA_name]["qc_nto_offtargets"] = 1

	###################
	## Create simplified form of NTO sam file to aid siRNA design
	sam_df = pd.read_csv(sam_file_nto, sep='\t', comment='@', usecols=[0,2,3,12,13], header=0, names=['siRNA_name', 'NTO_hit', 'position', 'mismatched_bases', 'mismatches'])

	sam_df['mismatched_bases'] = sam_df['mismatched_bases'].str.replace('MD:Z:','')
	sam_df['mismatches'] = sam_df['mismatches'].str.replace('NM:i:', '')
	sam_df['cds_name'] = sam_df['siRNA_name'].replace(r'_\d+$', '', regex=True)

	# Add species info
	acc_ids = pd.read_csv(NTO_ids, sep='\t', header=0, names=['NTO_hit', 'species'])

	sam_df_merged = pd.merge(sam_df, acc_ids, how='left', on='NTO_hit').drop('cds_name', axis=1)
	# Rearrange columns
	sam_df_merged = sam_df_merged[['siRNA_name', 'NTO_hit', 'species', 'position', 'mismatched_bases', 'mismatches']]
	
	sam_df_merged.to_csv(gene + "_NTO_hits.tsv", index=False, sep='\t')
	###################

	# print the results for all siRNAs of the current gene #
	good_pos = [] # this array will hold the position of the "good" siRNAs
	bad_pos = [] # this array will hold the position of the "bad" siRNAs, i.e. the ones holding siRNAs that hit NTOs
		      # (to be used in the next section)

	for siRNA_name in properties:
		
		pos = int( siRNA_name.split("_")[-1] )

		out = siRNA_name
		out += "\t" + properties[siRNA_name]["sequence"]
		out += "\t" + str(properties[siRNA_name]["qc_asymmetry"])
		out += "\t" + str(properties[siRNA_name]["qc_nt_runs"])
		out += "\t" + str(properties[siRNA_name]["qc_gc_content"])
		out += "\t" + str(properties[siRNA_name]["qc_specificity"])
		out += "\t" + str(properties[siRNA_name]["qc_offtargets"])
		out += "\t" + str(properties[siRNA_name]["qc_nto_offtargets"])
		out += "\n"
		fhall.write( out )
		
		# if this siRNA is satisfying the TO QC criteria then 
		# print it in the "good" siRNA output.
		if properties[siRNA_name]["qc_asymmetry"] == 1 \
			and properties[siRNA_name]["qc_nt_runs"] == 0 \
			and properties[siRNA_name]["qc_gc_content"] == 1 \
			and properties[siRNA_name]["qc_specificity"] == 1 \
			and properties[siRNA_name]["qc_offtargets"] == 0:
			# and properties[siRNA_name]["qc_nto_offtargets"] == 0:
			fhgood.write( out )
			
			good_pos.append( pos )
		
		# if this siRNA is not satisfying the NTO QC then 
		# print it in the "bad" siRNA output.
		if properties[siRNA_name]["qc_asymmetry"] == 1 \
			and properties[siRNA_name]["qc_nt_runs"] == 0 \
			and properties[siRNA_name]["qc_gc_content"] == 1 \
			and properties[siRNA_name]["qc_specificity"] == 1 \
			and properties[siRNA_name]["qc_nto_offtargets"] == 1:
			# or ((properties[siRNA_name]["qc_offtargets"] == 1) \
			fhbad.write( out )
			
			bad_pos.append( pos )
		
		# if this siRNA is not satisfying the NTO QC criteria no matter the TO QC then 
		# print it in the "plainbad" siRNA output.
		elif properties[siRNA_name]["qc_nto_offtargets"] == 1:
			fhplainbad.write( out )
			
	
	# End of results printing #####################

	
	# Finally, find and print the sequence of the dsRNA
	# based on the position of the "good/bad" siRNAs. Ideally,
	# you'd like your dsRNA to contain as many good siRNAs
	# as possible (to maximize the silencing effect)
	DS_LEN = ds_len
	
	best_good_cnt = 0 # this is the highest number of good siRNAs
	best_bad_cnt = 0  # this is the highest number of bad siRNAs
		              # contained in a given dsRNA of length DS_LEN
	best_pos = 0 # and this is where the dsRNA starts

	for i in range( 0, (q_len - int(DS_LEN)) ):
		curr_good_cnt = 0 # count of "good" siRNAs, i.e not hitting NTOs contained in the current dsRNA
		
		for pos in good_pos:
			if pos > i and pos < i + int(DS_LEN):
				curr_good_cnt += 1
		
		if curr_good_cnt > best_good_cnt:
			best_good_cnt = curr_good_cnt
			best_pos = i

	# Find "bad" off-target siRNAs contained in best dsRNA segment
	best_range = list(range(best_pos, best_pos + int(DS_LEN)))
	for pos in bad_pos:
		if pos >= best_range[0] and pos <= best_range[-1]:
			best_bad_cnt += 1
	
	# get the sequence of the dsRNA
	dsRNA_sequence = fasta[gene][best_pos:(best_pos + int(DS_LEN))]
	
	# write the data
	out = gene
	out += "\t" + str(best_pos)
	out += "\t" + str(best_pos + int(DS_LEN))
	out += "\t" + str(q_len)
	out += "\t" + str(best_good_cnt)
	out += "\t" + str(best_bad_cnt)
	out += "\t" + dsRNA_sequence
	out += "\n"
	fhout_dsRNA.write( out )
	
	if DELETE_TMP:
		os.remove( gene + ".fa" )

	## End of dsRNA printing #############################

fhgood.close()
fhbad.close()
fhplainbad.close()
fhall.close()
fhout_dsRNA.close()
