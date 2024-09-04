#!/usr/bin/python3

import os
import sys
import argparse
import re
# import pandas as pd

parser = argparse.ArgumentParser(description='This script reads in a list of transcripts and evaluates their suitability as RNAi targets in a target organism. It will also search for off-targets in the target organism as well as for other non-target organisms (NTOs).')
# https://stackoverflow.com/questions/24180527/argparse-required-arguments-listed-under-optional-arguments
requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument('-i', '--input', help='path to input transcripts', required=True)
requiredNamed.add_argument('-g', '--genome', help='path to target genome', required=True)
requiredNamed.add_argument('-n', '--NTO_genome', help='path to NTO genome(s)', required=True)
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
siRNA_len= args.siRNA_length
# GFF = sys.argv[6]
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
#########################################
DELETE_TMP = True  # delete temporary file

# first, check if the genome is formatted
if not os.path.isfile(genome + ".nhr"):
	sys.stderr.write ( "Genome is not formatted. Running makeblastdb...\n" )
	return_code = os.system( "makeblastdb -in " + genome + " -dbtype nucl" )
	if return_code != 0:
		sys.stderr.write( "makeblastdb exited with a non-zero code: " + return_code + "\n" )
		exit( return_code )
	else:
		sys.stderr.write( "makeblastdb finished successfully!\n" )

# Also check if the genomes of the NTOs are formatted
if not os.path.isfile(NTOs + ".njs"):
	sys.stderr.write ( "NTOs genomes are not formatted. Running makeblastdb...\n" )
	return_code = os.system( "makeblastdb -in " + NTOs + " -dbtype nucl" )
	if return_code != 0:
		sys.stderr.write( "makeblastdb exited with a non-zero code: " + return_code + "\n" )
		exit( return_code )
	else:
		sys.stderr.write( "makeblastdb finished successfully!\n" )

# #############################
# # Import the genomic.gff file
# # can use that for finding specificity on target gene rather than blasting
# gff = pd.read_csv(sys.argv[4], sep='\t', header=None, comment='#')
# gff_cds = gff[gff[2] == "CDS"]
# gff_cds[9] = gff_cds[8].replace(to_replace=r'^.+Name=([^;]+);.+$', value=r'\1',regex=True)
# gff_cds = gff_cds[[9, 3, 4, 6, 0]]. rename(columns={9:'CDS_name', 3:'start', 4:'end', 6:'strand', 0:'chromosome'})
# #############################

# then, open the fasta file containing the dsRNA sequences
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
		#print (header)
		#print (sequence)
		fasta[header] = sequence

fh1.close()

sys.stderr.write( "Finished reading the dsRNA fasta file\n" )

# open the file where you'll write the dsRNA sequence (and other details)
fhout_dsRNA = open( "dsRNAs_per_gene.tsv", "w")
fhout_dsRNA.write( "TranscriptID\tdsRNA_start\tdsRNA_stop\tTranscript_length\tsiRNAs_not_targeting_NTOs\tsiRNAs_also_targeting_NTOs\tdsRNA_sequence\n" )

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
	sys.stderr.write ( "\nAnalyzing gene: " + gene + "...\n" )
	
	tmp_file = gene + ".fa"
	fhout = open ( tmp_file, "w" )
	fhout.write( ">" + gene + "\n" )
	fhout.write( fasta[gene] + "\n" )
	fhout.close()
	sys.stderr.write( "Find the genomic locus from which the dsRNA originates\n" )
	return_code = os.system( 'blastn -query ' + tmp_file + ' -db ' + genome + ' -out ' + gene + '.blastn.fmt6 -num_threads ' + str(CPUS) + ' -evalue 1e-50 -word_size 7 -dust no -outfmt "6 std qlen slen staxids stitle"' )
	if return_code != 0:
		sys.stderr.write( "blastn exited with a non-zero exit code: " + return_code + "\n" )
		exit( return_code )
	else:
		sys.stderr.write( "blastn finished sucessfully!\n" )
	
	# open the blast output and get the coordinates of the locus of origin
	fhbl = open ( gene + ".blastn.fmt6" )
	
	s_coords = []  # and these are the subject (genome) coordinates
	
	q_len = len ( fasta[gene] )      # that's the length of the input mRNA
	
	hit_len = 0                # this is the cumulative length of the blastn hits
				   # Ideally, it should be the same as q_len
	
	for line in fhbl:
		f = line.split( "\t" )
		if float(f[2]) > 98:   # the dsRNA should match perfectly to the corresponding genome
				       # All other hits are disregarded
			if int(f[8]) < int(f[9]):    # if the hit is on the plus strand
				s_coords.append ( f[1] + "\t" + f[8] + "\t" + f[9] )
			else:
				s_coords.append ( f[1] + "\t" + f[9] + "\t" + f[8] )
			
			hit_len += int(f[7]) - int(f[6]) + 1
	
	fhbl.close()
	
	if DELETE_TMP:
		os.remove( gene + ".fa" )
		os.remove( gene + ".blastn.fmt6" )
	
	if hit_len/q_len > 0.9 and hit_len/q_len <= 1:
		sys.stderr.write( "The majority of the dsRNA sequence was found in the genome: " + str(hit_len) + " / " + str(q_len) + " bp\n" )
	elif hit_len/q_len <= 0.9:
		sys.stderr.write( "A significant part of the dsRNA couldn't be found in the genome: " + str(hit_len) + " / " + str(q_len) + " bp\n" )
	elif hit_len/q_len > 1:
		sys.stderr.write( "The sum of hits is greater than the dsRNA length: " + str(hit_len) + " / " + str(q_len) + " bp\n" )
	else:
		sys.stderr.write( "You should never get here!\n" )

	sys.stderr.write( "Found the genomic locus of origin\n" )

	# Now take all possible sequences of length=siRNA_len from the dsRNA sequence and blast against the genome sequence
	
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

	# BLAST the siRNA sequence against the genome of the target organism
	# in order to verify specificity and find any off-targets

	# # original version
	# return_code = os.system ( 'blastn -query siRNAs.fa -db ' + genome + ' -out siRNAs.blastn.fmt6 -num_threads ' + str(CPUS) + ' -evalue 0.1 -word_size 10 -dust no -outfmt "6 std qlen slen staxids stitle"' )

	# blastn-short
	return_code = os.system ( 'blastn -task blastn-short -query siRNAs.fa -db ' + genome + ' -out siRNAs.blastn.fmt6 -num_threads ' + str(CPUS) + ' -evalue 0.1 -word_size 7 -dust no -outfmt "6 std qlen slen staxids stitle"' )
		
	if return_code > 0:
		sys.stderr.write( "blastn of one siRNA exited with a non-zero exit code: " + return_code + "\n" )
		exit (return_code)
	
	fhbl = open ( "siRNAs.blastn.fmt6" )
	
	
	for line in fhbl:
		f = line.split( "\t" )
		
		cds_name	= re.sub(r'_\d+$', '', f[0])
		siRNA_name      = f[0]
		aln_length  = int(f[3])
		mismatches  = int(f[4])
		gapopen     = int(f[5])
		siRNA_coord = int(f[8])
		chromosome  = f[1]

		if aln_length >= (SIRNA_LEN - 2) and mismatches <= 2 and gapopen == 0:
			for coord in s_coords:
				genome_coords = coord.split ( "\t" )
				
				if chromosome == genome_coords[0] and siRNA_coord >= int(genome_coords[1]) and siRNA_coord <= int(genome_coords[2]):
					properties[siRNA_name]["qc_specificity"] = 1

			# if the current good hit wasn't found within the coordinates of the locus of origin
			# then this good hit is an off-target of the siRNA
			if properties[siRNA_name]["qc_specificity"] == 0:
				properties[siRNA_name]["qc_offtargets"] = 1
	fhbl.close()
	
	if DELETE_TMP:
		os.remove( "siRNAs.blastn.fmt6" )
	### End of BLAST vs the target organism ##############################################





	# BLAST the siRNA sequence against the genome of the non-target organisms (NTOs)
	# in order to find off-targets in those organisms
	
	# # original version
	# return_code = os.system ( 'blastn -query siRNAs.fa -db ' + NTOs + ' -out siRNAs_NTOs.blastn.fmt6 -num_threads ' + str(CPUS) + ' -evalue 0.1 -word_size 10 -dust no -outfmt "6 std qlen slen staxids stitle"' )

	# blastn-short
	return_code = os.system ( 'blastn -task blastn-short -query siRNAs.fa -db ' + NTOs + ' -out siRNAs_NTOs.blastn.fmt6 -num_threads ' + str(CPUS) + ' -evalue 0.1 -word_size 7 -dust no -outfmt "6 std qlen slen staxids stitle"' )
	
	if return_code > 0:
		sys.stderr.write( "blastn of one siRNA exited with a non-zero exit code: " + return_code + "\n" )
		exit (return_code)
	
	fhbl = open ( "siRNAs_NTOs.blastn.fmt6" )
	
	
	for line in fhbl:
		f = line.split( "\t" )
		siRNA_name  = f[0]
		mismatches  = int(f[4])
		gapopen     = int(f[5])
		
		if mismatches <= int(mis) and gapopen == 0:
			properties[siRNA_name]["qc_nto_offtargets"] = 1
	fhbl.close()
	# if DELETE_TMP:
	# 	os.remove( "siRNAs_NTOs.blastn.fmt6" )
	# 	os.remove( "siRNAs.fa" )
	## End of BLAST vs the NTOs ################################



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
		
		# if this siRNA is also satisfying the QC criteria then 
		# print it in the "good" siRNA output as well.
		if properties[siRNA_name]["qc_asymmetry"] == 1 \
			and properties[siRNA_name]["qc_nt_runs"] == 0 \
			and properties[siRNA_name]["qc_gc_content"] == 1 \
			and properties[siRNA_name]["qc_specificity"] == 1 \
			and properties[siRNA_name]["qc_offtargets"] == 0 \
			and properties[siRNA_name]["qc_nto_offtargets"] == 0:
			fhgood.write( out )
			
			good_pos.append( pos )
		
		# if this siRNA is satisfying the TO QC criteria but not the NTO QC then 
		# print it in the "bad" siRNA output.
		elif properties[siRNA_name]["qc_asymmetry"] == 1 \
			and properties[siRNA_name]["qc_nt_runs"] == 0 \
			and properties[siRNA_name]["qc_gc_content"] == 1 \
			and properties[siRNA_name]["qc_specificity"] == 1 \
			and properties[siRNA_name]["qc_offtargets"] == 0 \
			and properties[siRNA_name]["qc_nto_offtargets"] == 1:
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
	# DS_LEN = 500 # the length of the dsRNA
	DS_LEN = ds_len
	
	best_good_cnt = 0 # this is the highest number of good siRNAs
	best_bad_cnt = 0  # this is the highest number of bad siRNAs
		              # contained in a given dsRNA of length DS_LEN
	best_pos = -1 # and this is where the dsRNA starts

	for i in range( 0, (q_len - int(DS_LEN) + 1) ):
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
	
	## End of dsRNA printing #############################

fhgood.close()
fhbad.close()
fhplainbad.close()
fhall.close()
fhout_dsRNA.close()
