import os
import multiprocessing as mp
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

from tqdm import tqdm

import numpy as np
import networkx as nx

from remove_recombination.read_panout import collate_pair_chunk
from remove_recombination.read_panout import parse_pangenome_preload

from remove_recombination.write_output import remove_recombinant_seqs
from remove_recombination.write_output import write_rm_estimate
from remove_recombination.write_output import get_core_gene_nodes
from remove_recombination.write_output import concatenate_core_genome_alignments
from remove_recombination.recomb_model_functions import *
from remove_recombination.recombination_networks import *

PAIR_CHUNK_SIZE = 10000
PARALLEL_ANALYSIS_MIN_BLOCKS = 2

WORKER_PAIR_TUPLES = None
WORKER_PAIR_NAMES = None
WORKER_GENE_NAMES = None
WORKER_PRELOAD = None
WORKER_METHOD = None
WORKER_WRITE_DATA = None


def _init_chunk_worker(pair_tuples, pair_names, gene_names, preload, method, write_data):
    global WORKER_PAIR_TUPLES
    global WORKER_PAIR_NAMES
    global WORKER_GENE_NAMES
    global WORKER_PRELOAD
    global WORKER_METHOD
    global WORKER_WRITE_DATA

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"

    WORKER_PAIR_TUPLES = pair_tuples
    WORKER_PAIR_NAMES = pair_names
    WORKER_GENE_NAMES = gene_names
    WORKER_PRELOAD = preload
    WORKER_METHOD = method
    WORKER_WRITE_DATA = write_data


def iter_pair_chunks(pairs, chunk_size):
    for chunk_start in range(0, len(pairs), chunk_size):
        chunk_end = min(chunk_start + chunk_size, len(pairs))
        yield chunk_start, chunk_end, pairs[chunk_start:chunk_end]


def write_pairwise_distribution_rows(outhandle, pairs, gene_names, chunk_start, ordered_chunk):
    for local_idx, (ordered_dist_len, ordered_gene_ids) in enumerate(ordered_chunk):
        pair_idx = chunk_start + local_idx
        outline = pairs[pair_idx] + ","
        dists = ";".join(ordered_dist_len[:, 0].astype(str)) + ","
        lens = ";".join(ordered_dist_len[:, 1].astype(str)) + ","
        genes = ";".join(gene_names[int(x)] for x in ordered_gene_ids)
        outline += dists
        outline += lens
        outline += genes
        outhandle.write(outline + "\n")


def build_pairwise_distribution_rows(pairs, gene_names, chunk_start, ordered_chunk):
    rows = []
    for local_idx, (ordered_dist_len, ordered_gene_ids) in enumerate(ordered_chunk):
        pair_idx = chunk_start + local_idx
        outline = pairs[pair_idx] + ","
        dists = ";".join(ordered_dist_len[:, 0].astype(str)) + ","
        lens = ";".join(ordered_dist_len[:, 1].astype(str)) + ","
        genes = ";".join(gene_names[int(x)] for x in ordered_gene_ids)
        outline += dists
        outline += lens
        outline += genes
        rows.append(outline + "\n")
    return rows


def fold_chunk_results(
    chunk_start,
    ordered_chunk,
    chunk_results,
    gene_names,
    total_dists,
    gene_recombination_dic,
    recombinant_gene_pair_dist,
):
    for local_idx, ((ordered_dist_len, ordered_gene_ids), (pair_recombinants, pair_dists)) in enumerate(
        zip(ordered_chunk, chunk_results)
    ):
        pair_idx = chunk_start + local_idx
        gene_idx_to_dist = dict(zip(ordered_gene_ids, ordered_dist_len[:, 0]))
        for gene in pair_recombinants:
            gene_name = gene_names[gene]
            gene_recombination_dic[gene_name].append(pair_idx)
            recombinant_gene_pair_dist[gene_name][pair_idx] = int(gene_idx_to_dist[gene])

        total_dists[pair_idx] = pair_dists[0]


def process_pair_block(chunk_start, chunk_end):
    chunk_pairs = WORKER_PAIR_TUPLES[chunk_start:chunk_end]
    collated_chunk = collate_pair_chunk(chunk_pairs, WORKER_PRELOAD)
    ordered_chunk = order_pairwise_diffs_chunk(collated_chunk)
    chunk_results = analyse_pair_chunk(ordered_chunk, WORKER_METHOD)

    block_total_dists = []
    recombinant_entries = []
    for local_idx, ((ordered_dist_len, ordered_gene_ids), (pair_recombinants, pair_dists)) in enumerate(
        zip(ordered_chunk, chunk_results)
    ):
        pair_idx = chunk_start + local_idx
        block_total_dists.append(pair_dists[0])
        gene_idx_to_dist = dict(zip(ordered_gene_ids, ordered_dist_len[:, 0]))
        for gene_idx in pair_recombinants:
            recombinant_entries.append(
                (
                    WORKER_GENE_NAMES[gene_idx],
                    pair_idx,
                    int(gene_idx_to_dist[gene_idx]),
                )
            )

    distribution_rows = None
    if WORKER_WRITE_DATA:
        distribution_rows = build_pairwise_distribution_rows(
            WORKER_PAIR_NAMES,
            WORKER_GENE_NAMES,
            chunk_start,
            ordered_chunk,
        )

    return {
        "chunk_start": chunk_start,
        "chunk_end": chunk_end,
        "total_dists": block_total_dists,
        "recombinant_entries": recombinant_entries,
        "distribution_rows": distribution_rows,
    }


def merge_block_result(
    block_result,
    total_dists,
    gene_recombination_dic,
    recombinant_gene_pair_dist,
):
    chunk_start = block_result["chunk_start"]
    for local_idx, total_dist in enumerate(block_result["total_dists"]):
        total_dists[chunk_start + local_idx] = total_dist

    for gene_name, pair_idx, gene_dist in block_result["recombinant_entries"]:
        gene_recombination_dic[gene_name].append(pair_idx)
        recombinant_gene_pair_dist[gene_name][pair_idx] = gene_dist


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
    pair_tuples, gene_names, alignment_dir, preload = parse_pangenome_preload(
        args.outdir,
        args.n_cpu,
    )
    pairs = ["-".join(pair) for pair in pair_tuples]
    pair_index_to_isolates = pair_tuples
    
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
    
    print("Identifying recombinants...")

    block_ranges = [
        (chunk_start, min(chunk_start + PAIR_CHUNK_SIZE, len(pair_tuples)))
        for chunk_start in range(0, len(pair_tuples), PAIR_CHUNK_SIZE)
    ]

    pairwise_out = None
    if args.write_data:
        pairwise_out = open(
            args.outdir + "pairwise_difference_distributions.csv",
            "w+",
        )
        pairwise_out.write("pair,diffs,lens,gene_names\n")

    try:
        with tqdm(total=len(pair_tuples), desc="Identifying recombinants") as pbar:
            if args.n_cpu == 1 or len(block_ranges) < PARALLEL_ANALYSIS_MIN_BLOCKS:
                _init_chunk_worker(
                    pair_tuples,
                    pairs,
                    gene_names,
                    preload,
                    args.method,
                    args.write_data,
                )
                for chunk_start, chunk_end in block_ranges:
                    block_result = process_pair_block(chunk_start, chunk_end)
                    merge_block_result(
                        block_result,
                        total_dists,
                        gene_recombination_dic,
                        recombinant_gene_pair_dist,
                    )
                    if pairwise_out is not None:
                        pairwise_out.writelines(block_result["distribution_rows"])
                    pbar.update(chunk_end - chunk_start)
            else:
                ctx = mp.get_context("fork")
                pending_rows = {}
                next_chunk_idx_to_write = 0
                block_start_order = [chunk_start for chunk_start, _ in block_ranges]

                with ProcessPoolExecutor(
                    max_workers=args.n_cpu,
                    mp_context=ctx,
                    initializer=_init_chunk_worker,
                    initargs=(
                        pair_tuples,
                        pairs,
                        gene_names,
                        preload,
                        args.method,
                        args.write_data,
                    ),
                ) as executor:
                    future_to_range = {
                        executor.submit(process_pair_block, chunk_start, chunk_end): (
                            chunk_start,
                            chunk_end,
                        )
                        for chunk_start, chunk_end in block_ranges
                    }

                    for future in as_completed(future_to_range):
                        chunk_start, chunk_end = future_to_range[future]
                        block_result = future.result()
                        merge_block_result(
                            block_result,
                            total_dists,
                            gene_recombination_dic,
                            recombinant_gene_pair_dist,
                        )
                        if pairwise_out is not None:
                            pending_rows[chunk_start] = block_result["distribution_rows"]
                            while (
                                next_chunk_idx_to_write < len(block_start_order)
                                and block_start_order[next_chunk_idx_to_write] in pending_rows
                            ):
                                next_chunk_start = block_start_order[next_chunk_idx_to_write]
                                pairwise_out.writelines(pending_rows.pop(next_chunk_start))
                                next_chunk_idx_to_write += 1
                        pbar.update(chunk_end - chunk_start)
    finally:
        if pairwise_out is not None:
            pairwise_out.close()
    
    if not gene_recombination_dic:
        print("No recombinant genes identified.")
    else:
        for gene in gene_recombination_dic:
            gene_recombination_dic[gene] = sorted(set(gene_recombination_dic[gene]))

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
    
    
    
    
    
    
    
    
