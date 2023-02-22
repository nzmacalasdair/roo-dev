import os
import math
import decimal

import numpy as np
from scipy import sparse
from scipy import stats

import networkx as nx
from Bio import SeqIO

from scripts.remove_recombination.read_panout import parse_pangenome 
from scripts.remove_recombination.read_panout import remove_recombinant_seqs
from scripts.remove_recombination.read_panout import write_rm_estimate

from scripts.remove_recombination.recomb_model_functions import *

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
    args = parser.parse_args()
    
    #Make sure formatting is correct for panaroo dir, and create new out dir
    args.outdir = os.path.join(args.outdir, "")
    if not os.path.isdir(args.outdir + "recombination_free_aligned_genes/"):
        os.mkdir(args.outdir + "recombination_free_aligned_genes/")

    #Load in relevant info from genes
    gene_names, pairwise_differences = parse_pangenome(args.outdir)
    #Order genes from least snps/length to greatest snps/length
    ordered_pairs = order_pairwise_diffs(pairwise_differences)
    #Set up some empty dics for results
    gene_recombination_dic = {}
    relative_cleaned_dists = {}
    cleaned_dists = {}
    total_dists = {}
    proportions = {}
    #Do analysis, either bayesian or frequentist to identify recomb. gene pairs
    if args.method == "bayesian":
    
        for pair in ordered_pairs:
            model_probabilities, mean_distance = analyse_pair(ordered_pairs[pair][0][:,0], 
                                                        ordered_pairs[pair][0][:,1])
            
            threshold, recombinants = find_threshold(model_probabilities, 
                                                     ordered_pairs[pair][1])
            
            genes = ordered_pairs[pair][1]
            for gene in recombinants:
                    gene_recombination_dic[gene] = gene_recombination_dic.get(gene,
                                                                          []) + [pair]
            
            total_dist = sum(ordered_pairs[pair][0][:,0])
            cleaned_dist = total_dist - sum(ordered_pairs[pair][0][:threshold,0])
            total_dists[pair] = total_dist
            cleaned_dists[pair] = cleaned_dist
            proportion_distance = ordered_pairs[pair][0][:,0]/ordered_pairs[pair][0][:,1]
            proportions[pair] = proportion_distance
            relative_distance = np.mean(proportion_distance)/mean_distance
            relative_cleaned_dists[pair] = relative_distance
    
    elif args.method == "frequentist":
        
        for pair in ordered_pairs:
            pair_proportions = ordered_pairs[pair][0]
            dists = pair_proportions[:,0]
            lens = pair_proportions[:,1]
            
            genes = ordered_pairs[pair][1]
            
            threshold, recombinants = analyse_pair_frequentist(dists, lens, genes)
            
            mean_distance = sum(dists[:threshold]) / sum(lens[:threshold])
            
            genes = ordered_pairs[pair][1]
            
            for gene in recombinants:
                gene_recombination_dic[gene] = gene_recombination_dic.get(gene,
                                                                          []) + [pair]
            total_dist = sum(ordered_pairs[pair][0][:,0])
            cleaned_dist = total_dist - sum(ordered_pairs[pair][0][threshold:,0])
            total_dists[pair] = total_dist
            cleand_dists[pair] = cleaned_dist
            proportion_distance = ordered_pairs[pair][0][:,0]/ordered_pairs[pair][0][:,1]
            proportions[pair] = proportion_distance
            relative_distance = np.mean(proportion_distance)/mean_distance
            relative_cleaned_dists[pair] = relative_distance
    else:
        raise ValueError("Method must be one of [bayesian, frequentist]")
    
    #Reduce recombinant pairs to only isolates where recombination is present
    #Do this by making a network and taking only isolates of degree > 2
    actual_recombinants_to_remove = {}

    for gene in gene_recombination_dic:
        if len(gene_recombination_dic[gene]) > 1:
            gene_network = nx.Graph()
            for recombinant_pair in gene_recombination_dic[gene]:
                recombination_list = recombinant_pair.split("-")
                gene_network.add_edge(*recombination_list)
            if len(gene_network.nodes) > 4:
                to_remove = []
                min_degree = min([x[1] for x in gene_network.degree])
                for gene_degree in gene_network.degree:
                    if gene_degree[1] > min_degree:
                        to_remove.append(gene_degree[0])
            else:
                to_remove = list(gene_network.nodes)
        else:
            to_remove = gene_recombination_dic[gene][0].split("-")
        actual_recombinants_to_remove[gene] = to_remove
    #Remove recombinant sequences and write new alignments to file
    remove_recombinant_seqs(actual_recombinants_to_remove, args.outdir)
    #Write new core genome alignment
    
    #Take all the pairwise r/m estimates to estimate r/m for the collection
    #Use only the middle 80 percent of data (exclude >90th and <10th)
    mean_proportions = {}
    for pair in proportions:
        mean_proportions[pair] = np.mean(proportions[pair])
    
    rmregression, pairwise, xaxis = estimate_collection_rm(mean_proportions)
    
    print("Estimated r/m for collection: %s" %rmregression.slope)
    
    write_rm_estimate(rmregression, args.outdir)
    
    if args.plot_rm == True:
        import matplotlib
        import matplotlib.pyplot as plt
        plt.style.use('ggplot')
        
        plt.scatter(range(len(pairwise)), pairwise)
        plt.plot(range(len(pairwise)), 
                 rmregression.intercept + rmregression.slope*range(len(pairwise)), 
                 'r', label='fitted line')
        plt.savefig(args.outdir + "collection_rm_estimate_regression.png")


if __name__ == "__main__":
    main()
    
    
    
    
    
    
    
    
