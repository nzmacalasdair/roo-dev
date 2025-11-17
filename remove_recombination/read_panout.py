import os

from joblib import Parallel, delayed

from tqdm import tqdm

import numpy as np

from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

#This module contains all functions necessary to read panaroo input, and
# write all outputs

from remove_recombination.pairwise_comparisons import *

def read_and_close_fasta(filename):
    #This is necessary because SeqIO generaotrs and joblib Parallel do not play
    #nicely, the generators containing sequence require files to remain open,
    #meaning that every gene alignment in the pangenome needs to stay open!
    with open(filename, 'r') as inhandle:
       seq_generator = SeqIO.parse(inhandle, 'fasta')
       return list(seq_generator)

def get_all_pairwise_diffs(pairs, filt_genes, alignment_directory, threads, gpu):
    pair_diff_len_distributions = {}
    alignment_names = os.listdir(alignment_directory)
    print("Reading alignments...")
    filtered_alignment_names = []
    for file in alignment_names:
        name = file.split(".")[0]
        if name in filt_genes:
            filtered_alignment_names.append(file)
    gene_names = np.aray([x.split(".")[0] for x in filtered_alignment_names])
    
    filtered_alignment_paths = [alignment_directory + x for x in filtered_alignment_names]

    print("Calculating all pairwise differences...")
    allval_pairwise = Parallel(n_jobs=threads, prefer="processes")(
         delayed(get_pangenome_pairwise_differences)(alignment)
         for alignment in tqdm(filtered_alignment_paths))
    
    #reformat this data to the expected output format
    #first get dics for fast lookup
    isolate_gene_indices = {}
    alignment_isolate_row_lookup_dics = []
    for gene_idx in range(len(allval_pairwise)):
        #paiwirse_diffs = allval_pairwise[gene_idx][0]
        #comparison_lens = allval_pairwise[gene_idx][1]
        seqids = allval_pairwise[gene_idx][2]
        isolate_ids = [x.split(";")[0] for x in seqids]
        for isolate in isolate_ids:
            isolate_gene_indices[isolate] = isolate_gene_indices.get(isolate, 
                                                        set()) | {gene_idx}
        row_lookup_dic = {value: index for index, value in enumerate(isolate_ids)}
        alignment_isolate_row_lookup_dics.append(row_lookup_dic)
    
    #Go through all pairs and get distances/lengths
    print("Collating genes for each isolate pair...")
    pairids = ["-".join(x) for x in pairs]
    
    ##Multithreaded code slower due to overhead!
    # pair_diff_len_distributions = Parallel(n_jobs=threads, prefer="processes")(
    #     delayed(collate_pair_comparisons)(pairs[index], isolate_gene_indices, 
    #                              alignment_isolate_row_lookup_dics, 
    #                              allval_pairwise, gene_names) 
    #                         for index in tqdm(range(len(pairs))))
    ##Single threaded code
    pair_diff_len_distributions = []
    for iso1, iso2 in tqdm(pairs):
        #iso1 = pairs[index][0]
        #iso2 = pairs[index][1]
        shared_genes = list(isolate_gene_indices[iso1] & isolate_gene_indices[iso2])
        isolate_gene_names = gene_names[shared_genes]
        gene_dists =[]
        gene_lens = []
        
        for gene in shared_genes:            
            generows = alignment_isolate_row_lookup_dics[gene]
            iso1_row = generows[iso1]
            iso2_row = generows[iso2]
            
            dist_matrix = allval_pairwise[gene][0]
            dist = dist_matrix[iso1_row,iso2_row]
            if dist != dist_matrix[iso2_row,iso1_row]:
                raise ValueError("Reverse pairwise dists not equal!")
            
            length_matrix = allval_pairwise[gene][1]
            length1 = length_matrix[iso1_row]
            length2 = length_matrix[iso2_row]
            comparisonlen = min(length1, length2)
            
            gene_dists.append(dist)
            gene_lens.append(comparisonlen)
        pair_diff_len_distributions.append((isolate_gene_names, 
                                           gene_dists, gene_lens))
    
    return zip(pairids, pair_diff_len_distributions)    

def parse_pangenome(output_dir, threads, use_gpu):
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
    print(f"Entropy threshold set to {hc_threshold}.")
    
    for gene in hc_vals:
        if float(gene[1]) > hc_threshold:
            name = gene[0].split(".")[0]
            if name in genes:
                genes.remove(name)
    #Get all the distributions of pairwise differences
    
    pairs, pairwise_differences = get_all_pairwise_diffs(pairs, genes, 
                                                gene_alignments_dir, threads, use_gpu)

    return(pairs, pairwise_differences)

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