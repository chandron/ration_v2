#!/usr/bin/python3

import os
import getopt, sys
import re


if len( sys.argv ) != 4 or os.path.exists( sys.argv[1] ) == False or os.path.exists( sys.argv[2] ) == False or os.path.exists( sys.argv[3] ) == False:
	directories = sys.argv[0].split( "/" )
	print(
"""
Usage: %s <fasta file containing the mRNAs> <genome sequence of the target organism (fasta file)> <genome sequences of non-target organisms (fasta file)>

This script will take in a fasta file containing the mRNAs and
evaluate their suitability as RNAi targets in a target organism.
It will also search for off-targets in the target organism
as well as for other non-target organisms (NTOs).

""" % directories[-1] )
	sys.exit(1)

mRNAs = sys.argv[1]
genome = sys.argv[2]
NTOs   = sys.argv[3]

DELETE_TMP = True  # delete temporary file
CPUS = 4

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
if not os.path.isfile(NTOs + ".nhr"):
	sys.stderr.write ( "NTOs genomes are not formatted. Running makeblastdb...\n" )
	return_code = os.system( "makeblastdb -in " + NTOs + " -dbtype nucl" )
	if return_code != 0:
		sys.stderr.write( "makeblastdb exited with a non-zero code: " + return_code + "\n" )
		exit( return_code )
	else:
		sys.stderr.write( "makeblastdb finished successfully!\n" )

#
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
fhout_dsRNA.write( "GeneID\tdsRNA_start\tdsRNA_stop\tGene_length\tCount_of_good_siRNAs\tdsRNA_sequence\n" )

# open the output files
fhgood = open ( "siRNAs.good.tsv", "w" )
fhall  = open ( "siRNAs.all.tsv", "w" )

# print the header of the output
fhgood.write( "siRNA name\tsequence\tQC_asymmetry\tQC_nucleotide_runs\tQC_GC_content\tQC_specificity\tQC_offtargets\tQC_NTO_offtargets\n" )
fhall.write( "siRNA name\tsequence\tQC_asymmetry\tQC_nucleotide_runs\tQC_GC_content\tQC_specificity\tQC_offtargets\tQC_NTO_offtargets\n" )

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
	return_code = os.system( 'blastn -query ' + tmp_file + ' -db ' + genome + ' -out ' + gene + '.blastn.fmt6 -num_threads ' + str(CPUS) + ' -evalue 1e-50 -word_size 10 -dust no -outfmt "6 std qlen slen staxids stitle"' )
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


	# Now take all possible 21-nt sequences from the dsRNA sequence
	# and blast it against the genome sequence
	SIRNA_LEN = 21  # this is the default. Think twice before you change it!

	sys.stderr.write( "Examining all " + str(SIRNA_LEN) + "-nt sequences...\n" )
	
	properties = {}
	
	fhout = open ( "siRNAs.fa", "w" )

	for i in range( 0, len( fasta[gene] ) - SIRNA_LEN + 1 ):
		siRNA = fasta[gene][ i: (i + SIRNA_LEN) ] # this is the siRNA sequence
		
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
	return_code = os.system ( 'blastn -query siRNAs.fa -db ' + genome + ' -out siRNAs.blastn.fmt6 -num_threads ' + str(CPUS) + ' -evalue 0.1 -word_size 10 -dust no -outfmt "6 std qlen slen staxids stitle"' )
	
	if return_code > 0:
		sys.stderr.write( "blastn of one siRNA exited with a non-zero exit code: " + return_code + "\n" )
		exit (return_code)
	
	fhbl = open ( "siRNAs.blastn.fmt6" )
	
	
	for line in fhbl:
		f = line.split( "\t" )
		
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





	# BLAST the siRNA sequence against the genome of the non-target organisms
	# in order to find off-targets in those organisms
	return_code = os.system ( 'blastn -query siRNAs.fa -db ' + NTOs + ' -out siRNAs_NTOs.blastn.fmt6 -num_threads ' + str(CPUS) + ' -evalue 0.1 -word_size 10 -dust no -outfmt "6 std qlen slen staxids stitle"' )
	
	if return_code > 0:
		sys.stderr.write( "blastn of one siRNA exited with a non-zero exit code: " + return_code + "\n" )
		exit (return_code)
	
	fhbl = open ( "siRNAs_NTOs.blastn.fmt6" )
	
	
	for line in fhbl:
		f = line.split( "\t" )
		siRNA_name  = f[0]
		mismatches  = int(f[4])
		gapopen     = int(f[5])
		
		if mismatches <= 2 and gapopen == 0:
			properties[siRNA_name]["qc_nto_offtargets"] = 1
	fhbl.close()
	if DELETE_TMP:
		os.remove( "siRNAs_NTOs.blastn.fmt6" )
		os.remove( "siRNAs.fa" )
	## End of BLAST vs the NTOs ################################



	# print the results for all siRNAs of the current gene #
	good_pos = [] # this array will hold the position of the "good" siRNAs
		      # (to be used in the next section)

	for siRNA_name in properties:
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
			
			pos = int( siRNA_name.split("_")[1] )
			good_pos.append( pos )
	
	# End of results printing #####################

	
	# Finally, find and print the sequence of the dsRNA
	# based on the position of the "good" siRNAs. Ideally,
	# you'd like your dsRNA to contain as many good siRNAs
	# as possible (to maximize the silencing effect)
	DS_LEN = 500 # the length of the dsRNA
	
	best_cnt = 0  # this is the highest number of good siRNAs
		      # contained in a 500-bp dsRNA
	best_pos = -1 # and this is where the dsRNA starts

	for i in range( 0, (q_len - DS_LEN + 1) ):
		curr_cnt = 0 # count of good siRNAs contained in the current dsRNA
		for pos in good_pos:
			if pos > i and pos < i + DS_LEN:
				curr_cnt += 1
		
		if curr_cnt > best_cnt:
			best_cnt = curr_cnt
			best_pos = i
	
	# get the sequence of the dsRNA
	dsRNA_sequence = fasta[gene][best_pos:(best_pos + DS_LEN)]
	
	# write the data
	out = gene
	out += "\t" + str(best_pos)
	out += "\t" + str(best_pos + DS_LEN)
	out += "\t" + str(q_len)
	out += "\t" + str(best_cnt)
	out += "\t" + dsRNA_sequence
	out += "\n"
	fhout_dsRNA.write( out )
	
	## End of dsRNA printing #############################

fhgood.close()
fhall.close()
fhout_dsRNA.close()
