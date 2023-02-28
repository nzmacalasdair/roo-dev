import os

from joblib import Parallel, delayed

import numpy as np


from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

#This module contains all functions necessary to read panaroo input, and
# write all outputs


def get_pairwise_differences(str1, str2):
    if len(str1) != len(str2):
        raise ValueError("Sequences are of different lengths!")
    seq1 = np.fromstring(str1.lower(), dtype=np.int8)
    seq2 = np.fromstring(str2.lower(), dtype=np.int8)
    diffs = np.count_nonzero(seq1-seq2)
    length = len(str1)
    result = (np.array([diffs, length]))
    return result

def get_pangenome_pairwise_differences(sequence_files, sequences_to_consider):    
    diffs = []
    names = []
    for sequences in sequence_files:
        seq1 = None
        seq2 = None
        for sequence in sequences[0]:
            if sequences_to_consider[0] in sequence.id:
                seq1=str(sequence.seq)
            elif sequences_to_consider[1] in sequence.id:
                seq2 = str(sequence.seq)
            else:
                continue
        if (type(seq1) == str) and (type(seq2) == str):
            diffs.append(get_pairwise_differences(seq1, seq2))
            names.append(sequences[1].split(".")[0])
    diffs = np.array(diffs)
    return (diffs, names)

def get_pairs(isolate_list):
    pairs_list = [(isolate_list[i], isolate_list[j]) for i in range(len(isolate_list)) for j in range(i+1,len(isolate_list))]
    return pairs_list

def get_all_pairwise_diffs(pairs, filt_genes, alignment_directory, threads):
    pair_diff_len_distributions = {}
    alignment_names = os.listdir(alignment_directory)
    filtered_alignment_names = []
    for file in alignment_names:
        name = file.split(".")[0]
        if name in filt_genes:
            filtered_alignment_names.append(file)
    
    filtered_alignment_paths = [alignment_directory + x for x in filtered_alignment_names]
    sequences = Parallel(n_jobs=threads, prefer="threads")(
        delayed(SeqIO.parse)(x, 'fasta') for x in filtered_alignment_paths)
    sequences = [list(x) for x in sequences]
    
    alignments = [(sequences[x], 
                   filtered_alignment_names[x]) for x in range(len(sequences))]
    ##Legacy single-threaded code
    # alignments = []
    # for alignment in filtered_alignment_names:
    #     alignfile = alignment_directory + alignment
    #     sequences = list(SeqIO.parse(alignfile, 'fasta'))
    #     alignments.append((sequences, alignment))
    #    
    #for pair in pairs:
    #    pairid = "-".join(pair)
    #    pair_diff_len_distributions[pairid] = get_pangenome_pairwise_differences(alignments, pair)
    
    
    diffs_lens = Parallel(n_jobs=threads, prefer="threads")(
        delayed(get_pangenome_pairwise_differences)(pair, alignments)
        for pair in pairs)
    
    pairids = ["-".join(x) for x in pairs]
    pair_diff_len_distributions = {}
    for index in range(len(pairs)):
        pair_diff_len_distributions[pairids[index]] = diffs_lens[index]
    
    return [x.split(".")[0] for x in filtered_alignment_names], pair_diff_len_distributions

def parse_pangenome(output_dir, threads):
    if output_dir[-1] != "/":
        output_dir += "/"
    #Get all the pairwise comparison combinations
    gene_pa_file = output_dir + "gene_presence_absence.csv"
    if not os.path.isfile(gene_pa_file):
        raise ValueError("Panaroo output is missing, is the output directory correct?")
    with open(gene_pa_file) as inhandle:
        firstline = inhandle.readline()
    isolates = firstline.split(",")[3:]
    isolates = [x.strip() for x in isolates]
    pairs = get_pairs(isolates)
    
    #Get the names of the gene alignments
    #use codon alignments if they exist
    gene_alignments_dir = output_dir + "codon_aligned_gene_sequences/"
    if not os.path.isdir(gene_alignments_dir):
        gene_alignments_dir = output_dir + "aligned_gene_sequences/"
    if not os.path.isdir(gene_alignments_dir):
        raise ValueError("aligned_gene_sequences directory is missing!")
        
    gene_alignment_files = os.listdir(gene_alignments_dir)
    #Filter for genes present in >2 isolates    
    gene_alignment_files = [x for x in gene_alignment_files if ".aln.fas" in x]
    
    genes = [x.split(".")[0] for x in gene_alignment_files]
    
    #Filter genes based on entropy scores
    with open(output_dir + "alignment_entropy.csv", 'r') as inhandle:
        lines = inhandle.read().splitlines()
    hc_vals = [x.split(",") for x in lines]
    
    allh = np.array([float(gene[1]) for gene in hc_vals])
    q = np.quantile(allh, [0.25,0.75])
    hc_threshold = max(0.01, q[1] + 1.5*(q[1]-q[0]))
    print(f"Entropy threshold automatically set to {hc_threshold}.")
    
    for gene in hc_vals:
        if float(gene[1]) > hc_threshold:
            name = gene[0].split(".")[0]
            if name in genes:
                genes.remove(name)
    
    #Get all the distributions of pairwise differences
    
    ordered_genes, pairwise_differences = get_all_pairwise_diffs(pairs, genes, 
                                                gene_alignments_dir, threads)

    return(genes, pairwise_differences)

def write_rm_estimate(rm_regression, output_dir):
    outline1 = "Collection r/m estimate: " + str(rm_regression[0])
    outline2 = "Standard Error of Estimate: " + str(rm_regression[1])
    with open(output_dir + "rm_estimate.txt", 'w+') as outhandle:
        outhandle.write(outline1 + '\n')
        outhandle.write(outline2 + '\n')
    return True

def remove_recombinant_seqs(recombinations, out_dir):
    for gene in recombinations:
        alignment = out_dir + "aligned_gene_sequences/" + gene +".aln.fas"
        sequences = list(SeqIO.parse(alignment, 'fasta'))
        sequence_names = [x.id for x in sequences]
        
        for recombinant in recombinations[gene]:
            fasta_ids = [i for i in sequence_names if recombinant in i]
            indexes2remove = []
            for fid in fasta_ids:
                indexes2remove.append(sequence_names.index(fid))
            for index2remove in sorted(indexes2remove, reverse=True):
                del sequences[index2remove]
                del sequence_names[index2remove]
        if len(sequences) > 0:
            outname = out_dir + "recombination_free_aligned_genes/" + gene +".aln.fas"
            SeqIO.write(sequences, outname, 'fasta')
    
    return True

def get_core_gene_nodes(G, threshold, num_isolates):
    # Get the core genes based on percent threshold
    core_nodes = []
    for node in G.nodes():
        if float(G.nodes[node]["size"]) / float(num_isolates) >= threshold:
            core_nodes.append(node)
    return core_nodes

def write_alignment_header(alignment_list, outdir, filename):
    out_entries = []
    # Set the tracking variables for gene positions
    gene_start = 1
    gene_end = 0
    for gene in alignment_list:
        # Get length and name from one sequence in the alignment
        # Set variables that need to be set pre-output
        gene_end += gene[2]
        gene_name = gene[0]
        # Create the 3 line feature entry
        gene_entry1 = (
            "FT   feature         " + str(gene_start) + ".." + str(gene_end) + "\n"
        )
        gene_entry2 = "FT                   /label=" + gene_name + "\n"
        gene_entry3 = "FT                   /locus_tag=" + gene_name + "\n"
        gene_entry = gene_entry1 + gene_entry2 + gene_entry3
        # Add it to the output list
        out_entries.append(gene_entry)
        # Alter the post-output variables
        gene_start += gene[2]
    # Create the header and footer
    header = (
        "ID   Genome standard; DNA; PRO; 1234 BP.\nXX\nFH   Key"
        + "             Location/Qualifiers\nFH\n"
    )
    footer = (
        "XX\nSQ   Sequence 1234 BP; 789 A; 1717 C; 1693 G; 691 T;" + " 0 other;\n//\n"
    )
    # open file and output
    with open(outdir + filename, "w+") as outhandle:
        outhandle.write(header)
        for entry in out_entries:
            outhandle.write(entry)
        outhandle.write(footer)

    return True

def concatenate_core_genome_alignments(core_names, output_dir):

    alignments_dir = output_dir + "/recombination_free_aligned_genes/"
    # Open up each alignment that is associated with a core node
    alignment_filenames = os.listdir(alignments_dir)
    core_filenames = [
        x for x in alignment_filenames if x.split('.')[0] in core_names
    ]
    
    #Read in all these alignments
    gene_alignments = []
    isolates = set()
    for filename in core_filenames:
        gene_name = os.path.splitext(os.path.basename(filename))[0]
        alignment = AlignIO.read(alignments_dir + filename, "fasta")
        gene_dict = {}
        for record in alignment:
            if len(gene_dict)<1:
                gene_length = len(record.seq)

            if record.id[:3] == "_R_":
                record.id = record.id[3:]
            genome_id = record.id.split(";")[0]
            
            if genome_id in gene_dict:
                if str(record.seq).count("-") < str(gene_dict[genome_id][1]).count("-"):
                    gene_dict[genome_id] = (record.id, record.seq)
            else:
                gene_dict[genome_id] = (record.id, record.seq)
            
            isolates.add(genome_id)
        gene_alignments.append((gene_name, gene_dict, gene_length))
    # Combine them
    isolate_aln = []
    for iso in isolates:
        seq = ""
        for gene in gene_alignments:
            if iso in gene[1]:
                seq += gene[1][iso][1]
            else:
                seq += "-" * gene[2]
        isolate_aln.append(SeqRecord(seq, id=iso, description=""))

    # Write out the two output files
    SeqIO.write(isolate_aln, 
                output_dir + "recombination_free_core_gene_alignment.aln", 
                "fasta")
    write_alignment_header(gene_alignments, 
                           output_dir, 
                           "recombination_free_core_alignment_header.embl")

    return core_filenames