import os
from collections import defaultdict

from tqdm import tqdm

import numpy as np
import networkx as nx

from remove_recombination.read_panout import parse_pangenome 

from remove_recombination.write_output import remove_recombinant_seqs
from remove_recombination.write_output import write_rm_estimate
from remove_recombination.write_output import get_core_gene_nodes
from remove_recombination.write_output import concatenate_core_genome_alignments
from remove_recombination.recomb_model_functions import *
from remove_recombination.recombination_networks import *

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
    parser.add_argument("--write_data",
                        action="store_true",
                        help="""Output pairwise distributions and per-gene
                        reconbimation networks used to identify recombinants""")
    args = parser.parse_args()
    
    #Make sure formatting is correct for panaroo dir, and create new out dir
    args.outdir = os.path.join(args.outdir, "")
    output_alignment_dir = os.path.join(args.outdir, "recombination_free_aligned_genes")
    if not os.path.isdir(output_alignment_dir):
        os.mkdir(output_alignment_dir)

    #Check to make sure args.method is accurate
    if args.method not in ["bayesian", "frequentist"]:
        raise ValueError("Method must be one of [bayesian, frequentist]")
    if not 0 < args.core <= 1:
        raise ValueError("Core threshold must be in the range (0, 1].")
        
    #Load in relevant info from genes
    pairs, pairwise_differences, gene_names, alignment_dir = parse_pangenome(
        args.outdir,
        args.n_cpu,
    )
    
    #Order genes from least snps/length to greatest snps/length
    ordered_pairs = order_pairwise_diffs(pairwise_differences)
    if len(ordered_pairs) != len(pairs):
        raise ValueError("Pairwise analysis inputs do not match the pair list.")
    pair_index_to_isolates = [pair.split("-", 1) for pair in pairs]
    
    #Set up some empty dics for results
    gene_recombination_dic = defaultdict(list)
    total_dists = [0] * len(pairs)
    recombinant_gene_pair_dist = defaultdict(dict)

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
    
    #output debug files
    if args.write_data:
        with open(args.outdir + "pairwise_difference_distributions.csv", 
                  'w+') as outhandle:
            outhandle.write("pair,diffs,lens,gene_names\n")
            for pairidx in range(len(pairs)):
                outline = pairs[pairidx] +','
                dists = ';'.join(ordered_pairs[pairidx][0][:,0].astype(str)) +','
                lens = ";".join(ordered_pairs[pairidx][0][:,1].astype(str)) +','
                geneids = ordered_pairs[pairidx][1].astype(str)
                genes = ';'.join(gene_names[int(x)] for x in geneids)
                outline += dists
                outline += lens
                outline += genes
                outhandle.write(outline + '\n')
                
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
    if len(results) != len(pairs):
        raise ValueError("Pairwise analysis results do not match the pair list.")
    pairwise_recombinant_genes, mean_distances = zip(*results) 
    
    #Reformat pairwise results
    no_of_pairs = len(ordered_pairs)
    for pair_idx in range(no_of_pairs):
        pair_recombinants = pairwise_recombinant_genes[pair_idx]
        pair_dists = mean_distances[pair_idx]
        gene_idx_to_dist = dict(
            zip(
                ordered_pairs[pair_idx][1],
                ordered_pairs[pair_idx][0][:, 0],
            )
        )
        for gene in pair_recombinants:
            gene_name = gene_names[gene]
            gene_recombination_dic[gene_name].append(pair_idx)
            recombinant_gene_pair_dist[gene_name][pair_idx] = int(gene_idx_to_dist[gene])

        total_dists[pair_idx] = pair_dists[0]
    
    if not gene_recombination_dic:
        print("No recombinant genes identified.")

    #Reduce recombinant pairs to only isolates where recombination is present
    #Do this by making a network and taking only isolates of degree > 2
    if args.write_data:
        if not os.path.isdir(args.outdir + "pairwise_recombination_networks/"):
            os.mkdir(args.outdir + "pairwise_recombination_networks/")
    actual_recombinants_to_remove = {}
    print("Integrating pairwise results...")
    for gene in tqdm(gene_recombination_dic):
        if len(gene_recombination_dic[gene]) > 1:
            gene_network = build_recombination_network(
                gene_recombination_dic[gene],
                pair_index_to_isolates,
            )
            if args.write_data:
                nx.write_gml(gene_network, 
                         args.outdir + "pairwise_recombination_networks/" + gene + ".gml")
            to_remove = identify_genuine_recombinants(gene_network)
            actual_recombinants_to_remove[gene] = to_remove
        else:
            #Skip genes where there is only one pair of isolates
            continue
    #Output list of recombinant gene sequence names
    with open(args.outdir + "recombinant_gene_ids.csv", 'w+') as outhandle:
       outhandle.write("Gene,Recombinant_Isolates\n")
       for gene in actual_recombinants_to_remove:
           outline = gene +','+ ";".join(actual_recombinants_to_remove[gene])
           outhandle.write(outline + '\n')

    cleaned_dists, recombinant_dists, retained_pairs_by_gene = reconcile_cleaned_distances(
        total_dists,
        recombinant_gene_pair_dist,
        gene_recombination_dic,
        actual_recombinants_to_remove,
        pair_index_to_isolates,
    )

    if args.write_data:
        with open(args.outdir + "retained_recombinant_pairs.csv", "w+") as outhandle:
            outhandle.write("Gene,Recombinant_Isolates,Retained_Pairs\n")
            for gene in actual_recombinants_to_remove:
                retained_pairs = [
                    pairs[pair_idx] for pair_idx in retained_pairs_by_gene.get(gene, [])
                ]
                outline = (
                    gene
                    + ","
                    + ";".join(actual_recombinants_to_remove[gene])
                    + ","
                    + ";".join(retained_pairs)
                )
                outhandle.write(outline + "\n")

        with open(args.outdir + "pairwise_rm_components.csv", "w+") as outhandle:
            outhandle.write("Pair,Total_SNPs,Recombinant_SNPs,Cleaned_SNPs,Pairwise_r_m\n")
            for pair_idx, pair in enumerate(pairs):
                cleaned = cleaned_dists[pair_idx]
                recombinant = recombinant_dists[pair_idx]
                pairwise_rm = recombinant / cleaned if cleaned > 0 else ""
                outline = (
                    f"{pair},{total_dists[pair_idx]},{recombinant},{cleaned},{pairwise_rm}"
                )
                outhandle.write(outline + "\n")
    
    #Remove recombinant sequences and write new alignments to file
    remove_recombinant_seqs(
        actual_recombinants_to_remove,
        args.outdir,
        alignment_dir,
        gene_names,
    )
    #Write new core genome alignment
    G = nx.read_gml(args.outdir + "final_graph.gml")
    with open(args.outdir + "gene_presence_absence.Rtab", 'r') as inhandle:
        header = inhandle.readline()
    isolate_no = len(header.split()) - 1
    core_nodes = get_core_gene_nodes(G, args.core, isolate_no)
    core_names = [G.nodes[x]["name"] for x in core_nodes]
    concatenate_core_genome_alignments(core_names, args.outdir)
    
    #Estimate the collection r/m from retained recombinant and non-recombinant SNPs
    rm, stderr, dist_lists = estimate_collection_rm(cleaned_dists, recombinant_dists)
    
    print("Estimated r/m for collection: %s" %rm)
    
    write_rm_estimate((rm, stderr), args.outdir)
    
    if args.plot_rm == True:
        import matplotlib
        import matplotlib.pyplot as plt
        
        plt.scatter(dist_lists[0], dist_lists[1])
        plt.plot(np.arange(max(dist_lists[0]) + 1),
                 rm*np.arange(max(dist_lists[0]) + 1),
                 'r', label='r/m')
        plt.savefig(args.outdir + "collection_rm_estimate_regression.png")


if __name__ == "__main__":
    main()
    
    
    
    
    
    
    
    
