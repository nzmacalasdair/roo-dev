import os
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

from joblib import Parallel, delayed

from tqdm import tqdm

import numpy as np

from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

#This module contains all functions necessary to read panaroo input

from remove_recombination.pairwise_comparisons import *

def read_and_close_fasta(filename):
    #This is necessary because SeqIO generaotrs and joblib Parallel do not play
    #nicely, the generators containing sequence require files to remain open,
    #meaning that every gene alignment in the pangenome needs to stay open!
    with open(filename, 'r') as inhandle:
       seq_generator = SeqIO.parse(inhandle, 'fasta')
       return list(seq_generator)

# Global variables to avoid pickling
ISOLATE_GENE_INDICES = None
ALIGN_ISO_ROW_LOOKUPS = None
ALLVAL_PAIRWISE = None
#GENE_NAMES = None


def collate_pair(iso1, iso2):
    
    global ISOLATE_GENE_INDICES, ALIGN_ISO_ROW_LOOKUPS, ALLVAL_PAIRWISE
    
    shared_genes = list(ISOLATE_GENE_INDICES[iso1] & ISOLATE_GENE_INDICES[iso2])
    
    # Precompute the row indices for iso1 and iso2 in a single pass
    iso1_alnrows = np.array([ALIGN_ISO_ROW_LOOKUPS[gene][iso1] for gene in shared_genes])
    iso2_alnrows = np.array([ALIGN_ISO_ROW_LOOKUPS[gene][iso2] for gene in shared_genes])

    # Precompute the matrices for the shared genes
    dist_matrices = [ALLVAL_PAIRWISE[gene][0] for gene in shared_genes]
    length_matrices = [ALLVAL_PAIRWISE[gene][1] for gene in shared_genes]
    
    # Retrieve gene names directly as a NumPy array
    isolate_gene_names = np.array(shared_genes)
    
    # Use advanced indexing and avoid creating extra NumPy arrays
    gene_dists = np.array([dist_matrices[i][r1, r2] for i, r1, r2 in zip(range(len(shared_genes)), iso1_alnrows, iso2_alnrows)])
    gene_lens = np.minimum(
        np.array([length_matrices[i][r1] for i, r1 in zip(range(len(shared_genes)), iso1_alnrows)]),
        np.array([length_matrices[i][r2] for i, r2 in zip(range(len(shared_genes)), iso2_alnrows)])
    )
    
    return (isolate_gene_names, gene_dists, gene_lens)


def _init_worker(isolate_gene_indices,
                 alignment_isolate_row_lookup_dics,
                 allval_pairwise):

    global ISOLATE_GENE_INDICES
    global ALIGN_ISO_ROW_LOOKUPS
    global ALLVAL_PAIRWISE
    #global GENE_NAMES

    ISOLATE_GENE_INDICES = isolate_gene_indices
    ALIGN_ISO_ROW_LOOKUPS = alignment_isolate_row_lookup_dics
    ALLVAL_PAIRWISE = allval_pairwise
    #GENE_NAMES = gene_names


def parallel_collate_pairs(
    pairs,
    isolate_gene_indices,
    alignment_isolate_row_lookup_dics,
    allval_pairwise,
    n_cpu
):
    # Ensure fork mode (only works on POSIX systems)
    mp.set_start_method("fork", force=True)

    # Only small objects will be sent (the iso1, iso2 pairs)
    iso1s = [iso1 for iso1, iso2 in pairs]
    iso2s = [iso2 for iso1, iso2 in pairs]

    with ProcessPoolExecutor(
        max_workers=n_cpu,
        initializer=_init_worker,
        initargs=(
            isolate_gene_indices,
            alignment_isolate_row_lookup_dics,
            allval_pairwise,
        )
    ) as executor:

        results_iter = executor.map(collate_pair, iso1s, iso2s)

        return list(tqdm(results_iter, total=len(pairs)))


def get_all_pairwise_diffs(pairs, filt_genes, alignment_directory, threads):
    pair_diff_len_distributions = {}
    alignment_names = os.listdir(alignment_directory)
    print("Reading alignments...")
    filtered_alignment_names = []
    for file in alignment_names:
        name = file.split(".")[0]
        if name in filt_genes:
            filtered_alignment_names.append(file)
    gene_names = np.array([x.split(".")[0] for x in filtered_alignment_names])
    
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
    
    #Multithread with new fork backend
    pair_diff_len_distributions = parallel_collate_pairs(pairs, 
                                                          isolate_gene_indices, 
                                                          alignment_isolate_row_lookup_dics, 
                                                          allval_pairwise, 
                                                          threads)
        
    #Single threaded code is faster, again
    # pair_diff_len_distributions = []
    # for iso1, iso2 in tqdm(pairs):
     
    #     shared_genes = list(isolate_gene_indices[iso1] & isolate_gene_indices[iso2])
        
    #     iso1_alnrows = np.array([alignment_isolate_row_lookup_dics[gene][iso1] for gene in shared_genes])
    #     iso2_alnrows = np.array([alignment_isolate_row_lookup_dics[gene][iso2] for gene in shared_genes])
        
    #     dist_matrices = [allval_pairwise[gene][0] for gene in shared_genes]
    #     length_matrices = [allval_pairwise[gene][1] for gene in shared_genes]
        
               
    #     isolate_gene_names = gene_names[shared_genes]
    #     gene_dists = np.array([
    #                 dist_matrices[i][r1, r2]
    #                 for i, r1, r2 in zip(range(len(shared_genes)),
    #                                       iso1_alnrows, iso2_alnrows)
    #             ])
    #     gene_lens = np.minimum(np.array([length_matrices[i][r1] 
    #                                       for i, r1 in 
    #                                       zip(range(len(shared_genes)), iso1_alnrows)]), 
    #                         np.array([length_matrices[i][r2] 
    #                                   for i, r2 in 
    #                                   zip(range(len(shared_genes)), iso2_alnrows)]))
        
        # Do I need this debug check?
        # if not np.all(gene_dists == dist_matrices[np.arange(len(shared_genes)), 
        #                                           iso2_alnrows, iso1_alnrows]):
        #     raise ValueError("Reverse pairwise distances are not equal!")
        
        # avoid this loop if I can
        # for gene in shared_genes:            
            
        #     dist_matrix = allval_pairwise[gene][0]
        #     dist = dist_matrix[iso1_row,iso2_row]
        #     if dist != dist_matrix[iso2_row,iso1_row]:
        #         raise ValueError("Reverse pairwise dists not equal!")
            
        #     length_matrix = allval_pairwise[gene][1]
        #     length1 = length_matrix[iso1_row]
        #     length2 = length_matrix[iso2_row]
        #     comparisonlen = min(length1, length2)
            
        #     gene_dists.append(dist)
        #     gene_lens.append(comparisonlen)
        # pair_diff_len_distributions.append((isolate_gene_names, 
        #                                   gene_dists, gene_lens))
    
    
    if len(pairids) != len(pair_diff_len_distributions):
        raise ValueError("Pairwise comparisons not equal to number of pairs!")
    return pairids, pair_diff_len_distributions, gene_names 

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
    print(f"Entropy threshold set to {hc_threshold}.")
    
    for gene in hc_vals:
        if float(gene[1]) > hc_threshold:
            name = gene[0].split(".")[0]
            if name in genes:
                genes.remove(name)
    #Get all the distributions of pairwise differences
    
    pairs, pairwise_differences, gene_names, = get_all_pairwise_diffs(pairs, genes, gene_alignments_dir, threads)

    return pairs, pairwise_differences, gene_names