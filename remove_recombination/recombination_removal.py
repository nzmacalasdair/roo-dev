import os
import math
import decimal

from joblib import Parallel, delayed

from tqdm import tqdm

import numpy as np
from scipy import sparse
from scipy import stats
import networkx as nx
from Bio import SeqIO

from remove_recombination.read_panout import parse_pangenome 
from remove_recombination.read_panout import remove_recombinant_seqs
from remove_recombination.read_panout import write_rm_estimate
from remove_recombination.read_panout import get_core_gene_nodes
from remove_recombination.read_panout import concatenate_core_genome_alignments
from remove_recombination.read_panout import write_alignment_header

from remove_recombination.recomb_model_functions import *


def main():
    import argparse
    #Get arguments, output directory w/ aligned pangenome, and bayesian/frequentist
    parser = argparse.ArgumentParser(
        "Identify and remove recombinant gene sequences from core or pan alignment"
    )
    parser.add_argument("outdir",
                        help='Location of panaroo output directory with aligned genes')
    parser.add_argument('--method',
                        default="frequentist",
                        choices=['frequentist', 'bayesian'],
                        help=("""Which probability framework to use when 
                              detecting recombinations, pairwise. Empirical
                              testing suggests the Bayesian framework is 
                              more sensitive, but has a higher false-positive
                              rate."""))
    parser.add_argument("--plot_rm",
                        action="store_true",
                        help="""Plot fitted regression used to estimate 
                        collection's r/m""")
    parser.add_argument("--core_threshold",
                        dest="core",
                        help="Core-genome sample threshold, recombination free (default=0.95)",
                        type=float,
                        default=0.95)
    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)
    parser.add_argument("--gpu",
                        action="store_true",
                        help="""Use CuPy to speed up pairwise distance estimates. 
                        Highly recomended on datasets of >=10^3 isolates""")
    args = parser.parse_args()
    
    #Make sure formatting is correct for panaroo dir, and create new out dir
    args.outdir = os.path.join(args.outdir, "")
    if not os.path.isdir(args.outdir + "recombination_free_aligned_genes/"):
        os.mkdir(args.outdir + "recombination_free_aligned_genes/")

    #Check to make sure args.method is accurate
    if args.method not in ["bayesian", "frequentist"]:
        raise ValueError("Method must be one of [bayesian, frequentist]")
        
    #Load in relevant info from genes
    pairs, pairwise_differences, gene_names = parse_pangenome(args.outdir, args.n_cpu, args.gpu)
    
    #Order genes from least snps/length to greatest snps/length
    ordered_pairs = order_pairwise_diffs(pairwise_differences)
    
    #Set up some empty dics for results
    gene_recombination_dic = {}
    pairwise_rm_estimates = {}
    cleaned_dists = {}
    total_dists = {}

    #Do analysis, either bayesian or frequentist to identify recomb. gene pairs
    ##single-threaded code, for now    
    # for pair in ordered_pairs:
    #     if args.method == "bayesian":
    #         recombinants, dists = recombination_analysis_bayesian(ordered_pairs[pair])
    #         print(recombinants)
    #         print(dists)
    #     elif args.method == "frequentist":
    #         recombinants, dists = recombination_analysis_frequentist(ordered_pairs[pair])            
    #         print(recombinants)
    #         print(dists)
    #     for gene in recombinants:
    #             gene_recombination_dic[gene] = gene_recombination_dic.get(gene,
    #                                                                   []) + [pair]
    #     total_dists[pair] = dists[0]
    #     cleaned_dists[pair] = dists[1]
    #     pairwise_rm_estimates = dists[2]/dists[1]
    
    #multithreading
    print("Identifying recombinants...")
    
    #
    
    #if args.method =="bayesian":
         # results = Parallel(n_jobs=args.n_cpu, 
         #                                                      prefer="process")(
         #    delayed(recombination_analysis_bayesian)(ordered_pairs[index]) for index in tqdm(range(len(ordered_pairs))))
         # pairwise_recombinant_genes, mean_distances = zip(*results)                                                         
    #elif args.method == "frequentist":
        # results = Parallel(n_jobs=args.n_cpu, 
        #                                                       prefer="process")(
        #     delayed(recombination_analysis_frequentist)(ordered_pairs[index]) for index in tqdm(range(len(ordered_pairs)))) 
        # pairwise_recombinant_genes, mean_distances = zip(*results)
        
    results = do_recombination_analysis(ordered_pairs, args.method, args.n_cpu)    
    pairwise_recombinant_genes, mean_distances = zip(*results) 
    
    #Reformat pairwise results
    for index in range(len(ordered_pairs)):
        pair_recombinants = pairwise_recombinant_genes[index]
        pair_dists = mean_distances[index]
        pair = pairs[index] #pair is first position in the tuple
        for gene in pair_recombinants:
            gene_name = gene_names[gene]
            gene_recombination_dic[gene_name] = gene_recombination_dic.get(gene_name,
                                                                  []) + [pair]

        total_dists[pair] = pair_dists[0]
        cleaned_dists[pair] = pair_dists[1]
        pairwise_rm_estimates = pair_dists[2]/pair_dists[1]
    
    if gene_recombination_dic == {}:
        print("ERROR: Gene recombination dictionary is empty")
        print("Latest pair_recombinants:")
        print(pair_recombinants)
        print("All pairwise recombinants:")
        print(pairwise_recombinant_genes)
        print("results")
        print(results)
        print("Ordered pairs[0]:")
        print(ordered_pairs[0])

    #Reduce recombinant pairs to only isolates where recombination is present
    #Do this by making a network and taking only isolates of degree > 2
    if not os.path.isdir(args.outdir + "pairwise_recombination_networks/"):
        os.mkdir(args.outdir + "pairwise_recombination_networks/")
    actual_recombinants_to_remove = {}
    print("Integrating pairwise results...")
    for gene in tqdm(gene_recombination_dic):
        if len(gene_recombination_dic[gene]) > 1:
            gene_network = nx.Graph()
            for recombinant_pair in gene_recombination_dic[gene]:
                recombination_list = recombinant_pair.split("-")
                gene_network.add_edge(*recombination_list)
            nx.write_gml(gene_network, 
                         args.outdir + "pairwise_recombination_networks/" + gene + ".gml")
            if len(gene_network.nodes) >= 4:
                to_remove = []
                min_degree = min([x[1] for x in gene_network.degree])
                if min_degree == max([x[1] for x in gene_network.degree]):
                    to_remove = list(gene_network.nodes)
                else:
                    for gene_degree in gene_network.degree:
                        if gene_degree[1] > min_degree:
                            to_remove.append(gene_degree[0])
            else:
                to_remove = list(gene_network.nodes)
        else:
            to_remove = gene_recombination_dic[gene][0].split("-")
        actual_recombinants_to_remove[gene] = to_remove
    #Output list of recombinant gene sequence names
    with open(args.outdir + "recombinant_gene_ids.csv", 'w+') as outhandle:
       outhandle.write("Gene,Recombinant_Isolates\n")
       for gene in actual_recombinants_to_remove:
           outline = gene +','+ ";".join(actual_recombinants_to_remove[gene])
           outhandle.write(outline + '\n')
    
    #Remove recombinant sequences and write new alignments to file
    remove_recombinant_seqs(actual_recombinants_to_remove, args.outdir)
    #Write new core genome alignment
    G = nx.read_gml(args.outdir + "final_graph.gml")
    with open(args.outdir + "gene_presence_absence.Rtab", 'r') as inhandle:
        header = inhandle.readline()
    isolate_no = len(header.split()) - 1
    core_nodes = get_core_gene_nodes(G, isolate_no, args.core)
    core_names = [G.nodes[x]["name"] for x in core_nodes]
    concatenate_core_genome_alignments(core_names, args.outdir)
    
    #Estimate the collection r/m by pooling ratios and estimating the slope
    
    rm, stderr, dist_lists = estimate_collection_rm(cleaned_dists, total_dists)
    
    print("Estimated r/m for collection: %s" %rm)
    
    write_rm_estimate((rm, stderr), args.outdir)
    
    if args.plot_rm == True:
        import matplotlib
        import matplotlib.pyplot as plt
        
        plt.scatter(dist_lists[0], dist_lists[1])
        plt.plot(np.arange(max(dist_lists[0])), 
                 rm*np.arange(max(dist_lists[0])), 
                 'r', label='fitted line')
        plt.savefig(args.outdir + "collection_rm_estimate_regression.png")


if __name__ == "__main__":
    main()
    
    
    
    
    
    
    
    
