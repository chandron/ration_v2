import os
import subprocess
import tempfile
import pandas as pd
import pysam
import matplotlib.pyplot as plt
import csv
from brokenaxes import brokenaxes

def generate_siRNAs(sequence, si_length):
    siRNA_sequences = []
    for i in range(0, len(sequence) - si_length + 1):
        kmer = sequence[i:i+si_length]
        siRNA_sequences.append(kmer)
    return siRNA_sequences

def run_bowtie1(siRNA_names, siRNA_sequences, index_prefix, siRNA_length, mis_input="2"):
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fna') as tmp:
        for siRNA_name, siRNA_sequence in zip(siRNA_names, siRNA_sequences):
            tmp.write(f">siRNA_{siRNA_name}\n{siRNA_sequence}\n")
        tmp_fasta_name = tmp.name

    bowtie1_command = [
        "bowtie-1.3.1-macos-x86_64/bowtie-align-s",
        "-n", "2",
        "-l", str(siRNA_length),
        "-a",
        "-x", index_prefix,
        "-f", tmp_fasta_name,
        "-S", os.path.join("off_target_output/combined_siRNAs.sam")
    ]

    try:
        subprocess.run(bowtie1_command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {str(e)}")
    finally:
        os.unlink(tmp_fasta_name)

    return os.path.join("off_target_output/combined_siRNAs.sam")


def read_lethal_genes(species_list):
    all_lethals = {}
    for species in species_list:
        file_path = os.path.join("all_lethals", species + "_all_lethals")
        try:
            with open(file_path, 'r') as file:
                reader = csv.DictReader(file, delimiter='\t')
                all_lethals[species] = [row[species] for row in reader if row[species]]

        except FileNotFoundError:
            print(f"No lethal genes file found for {species}.")
        except KeyError:
            print(f"Species {species} not found in the file headers.")
    return all_lethals

def predict_off_target(input_sequence, species_list, kmer_length, all_lethals, mis_input="2"):
    alignment_data = {}
    siRNA_sequences = generate_siRNAs(input_sequence, kmer_length)

    filtered_siRNA_sequences = []
    siRNA_names = []
    for idx, kmer in enumerate(siRNA_sequences):
        polyAT_count = kmer.count('A') + kmer.count('T')
        if polyAT_count / kmer_length < 0.92:
            filtered_siRNA_sequences.append(kmer)
            siRNA_names.append(f"siRNA_{idx + 1}")

    siRNA_names = [f"siRNA_{idx + 1}" for idx in range(len(filtered_siRNA_sequences))]

    for species in species_list:
        index_prefix = os.path.join("indexes", species, species)
        sam_file_path = run_bowtie1(siRNA_names, filtered_siRNA_sequences, index_prefix, kmer_length, mis_input=mis_input)

        with pysam.AlignmentFile(sam_file_path, "r") as samfile:
            for read in samfile.fetch():
                if not read.is_unmapped:
                    siRNA_name = read.query_name
                    mismatches = read.get_tag("NM")
                    ref_name = read.reference_name
                    query_sequence = read.query_sequence
                    perfect_matches = len(query_sequence) - mismatches
                    strand = "reverse" if read.is_reverse else "forward"

                    # Check if the gene is lethal
                    gene_name = ref_name
                    if gene_name in all_lethals.get(species, []):
                        gene_name += "_LETHAL"

                    alignment_result = {
                        "query_sequence": query_sequence,
                        "fasta_file": os.path.basename(index_prefix).replace("_index", ".fna"),
                        "gene_name": gene_name,
                        "perfect_matches": perfect_matches,
                        "mismatches": mismatches,
                        "strand": strand
                    }

                    if species not in alignment_data:
                        alignment_data[species] = {}
                    if siRNA_name not in alignment_data[species]:
                        alignment_data[species][siRNA_name] = []

                    alignment_data[species][siRNA_name].append(alignment_result)

    return alignment_data

def main(seq_name, seq, species_input="Aphis_mellifera,Bombus_terrestris,Coccinella_septempunctata,Daphnia_magna,Chrysoperla_carnea", siRNA_length=21, mismatch_input=2):
    input_sequence = seq
    species_list = species_input.split(',')
    kmer_length = siRNA_length
    mismatch = str(mismatch_input)
    lethal_plot_data = {}  # Initialize data for lethal genes plot (outside the species loop)


    all_lethals = read_lethal_genes(species_list)
    alignment_data = predict_off_target(input_sequence, species_list, kmer_length, all_lethals, mis_input=mismatch)

    for species, siRNA_results in alignment_data.items():
        if species not in lethal_plot_data:
            lethal_plot_data[species] = {}
        species_df = pd.DataFrame()  # Initialize an empty DataFrame for each species

        # Collect data for each siRNA in the species
        for siRNA_name, data in siRNA_results.items():
            if data:
                df = pd.DataFrame(data)
                df['siRNA_name'] = siRNA_name  # Add a column for siRNA name
                species_df = pd.concat([species_df, df])

                # Collect data for lethal genes plot
                for entry in data:
                    if '_LETHAL' in entry['gene_name']:
                        position = int(siRNA_name.replace('siRNA_', ''))
                        lethal_plot_data[species].setdefault(position, 0)
                        lethal_plot_data[species][position] += 1

        # Write the DataFrame for each species to a tab-delimited text file
        if not species_df.empty:
            species_df.to_csv(f"off_target_output/{seq_name}_{species}_off_target.txt", sep='\t', index=False)

    # Plotting
    plt.figure(figsize=(20, 5))

    # Plot 1: Off-targets for siRNAs
    plt.subplot(1, 2, 1)
    for species, siRNA_results in alignment_data.items():
        plot_data = {}
        for siRNA_name, data in siRNA_results.items():
            for entry in data:
                position = int(siRNA_name.replace('siRNA_', ''))
                plot_data.setdefault(position, 0)
                plot_data[position] += 1
        if plot_data:
            plt.scatter(plot_data.keys(), plot_data.values(), label=species, s=3)

    plt.xlabel('siRNA Position')
    plt.ylabel('Number of Off-targets')
    plt.title(f'Off-targets for {seq_name}')
    #plt.legend(loc='lower left', bbox_to_anchor=(1, 1))  # Legend outside the plot


    plt.tight_layout()  # Adjust layout

    # Plot 2: Lethal genes
    plt.subplot(1, 2, 2)
    for species in lethal_plot_data:
        species_data = lethal_plot_data[species]
        if species_data:
            positions = list(species_data.keys())
            counts = list(species_data.values())
            plt.scatter(positions, counts, label=species, s=3)





    plt.tight_layout()  # Adjust layout



    plt.xlabel('siRNA Position')
    plt.ylabel('Number of Lethal Off-targets')
    plt.title(f'Lethal Off-targets for {seq_name}')
    plt.subplots_adjust(wspace=0.1)

    #
    plt.show()




### Usage ####


main(" Prosβ7", "AACCAACTTGAATCTGAATGACGCTAAAGCTCTTATAAAGAAAGCTGGAACGTCAATTGTCACTACAAACTGTCAATAGTGTCAATTGTCAAATTCTTAATTTTGTATTTTGGAAAATCATCGTGTGTTGTTTGTTCGCCTATTTGTGCCGAATTTGATTATCGTTTTTCTTCAGATACTCACAGAAGTGTTTGAGTTTTGTAATAATGTTTAATCCATTAGAAAGTGCACCACCTTTGTGGCATAACGGTCCAGTTCCAGGTGCCTTTTACAATTTTCCAGGCACTCAAAGAATTCCTGCTGCCAACCCGATAACCCATACTCAATCTCCTGTAACTACATCAACTTCTATAATAGCTATTACTTATGATAAAGGTGTATTAATTGCTGGAGATTTAGTCGCCTCTTATGGATCTTTAGCCAGGTACAGAAATTGTCCCAGGGTAATTGAAGTTAATCCTAATATAATATTAGGAGCTGGTGGAGATTATGCAGATTTTCAGTACGTGAAAAGCTTTATTGAGCAGAAAATAATTGATGAAGACTGTTTAGATGATGGTTTAAAAATGAAACCAAAATCTCTATATTGTTGGCTTACAAGAATCATGTATCAACGTAGAAGTAAACTTGACCCCTTCTGGAATAATCTTGTAGTAGGAGGACTGCAGGATGGTGTACCTTTCTTGGGAACCATAGATAAACTTGGCACAGCTTACACAGATAAAGTTATTTGTACTGGGTATGGAGCTCACATTGCTCTCCCTATTCTTCGGGATGCTCTTGATAAGAAAACCAACTTGAATCTTAATGACGCCAAAGCTCTTATAAAGCGGTGCATGGAAGTGCTTTTTTACAGAGATGCAAGAAGTTTTCCAAAATATCAGTTAGGTATTATTGATAAGGATGATGGCGTGAAAATTGAGGGACCACTTACAGTTGAAGAAAATTGGAATATTGCTTACATGACATAATTTGTCTCTTTAAATTTTGATATTTATAATTTTTGTAAGGTAATTCCCTATAATAAAGGTTTTAATATTTAAAAAAAAAAAAAAAAAAAAAAA")
