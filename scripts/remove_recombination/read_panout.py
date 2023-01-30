import os


import numpy as np


from Bio import SeqIO

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

def get_all_pairwise_diffs(pairs, alignment_directory):
    pair_diff_len_distributions = {}
    alignment_names = os.listdir(alignment_directory)
    alignments = []
    for alignment in alignment_names:
        alignfile = alignment_directory + alignment
        sequences = list(SeqIO.parse(alignfile, 'fasta'))
        alignments.append((sequences, alignment))
        
    for pair in pairs:
        pairid = "-".join(pair)
        pair_diff_len_distributions[pairid] = get_pangenome_pairwise_differences(alignments, pair)
    return pair_diff_len_distributions

def parse_pangenome(output_dir):
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
    genes = [x.split(".")[0] for x in gene_alignment_files]
    
    #Get all the distributions of pairwise differences
    
    pairwise_differences = get_all_pairwise_diffs(pairs, gene_alignments_dir)

    return(genes, pairwise_differences)

def write_rm_estimate(rm_regression, output_dir):
    outline1 = "Collection r/m estimate: " + str(rm_regression.slope)
    outline2 = "R^2 goodness-of-fit value: " + str(rm_regression.rvalue)
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
    