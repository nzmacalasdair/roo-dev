import os

import networkx as nx
from Bio import SeqIO

from pairsnp import calculate_snp_matrix, calculate_distance_matrix

from concatenate_core_genome import get_core_gene_nodes
from concatenate_core_genome import concatenate_core_genome_alignments

def check_gappedness(genelist, threshold, outdir):
    passed = []
    for gene in genelist:
        filename = outdir + "aligned_gene_sequences/" + gene +".aln.fas"
        sequences = SeqIO.parse(filename, 'fasta')
        seqlen = 0
        gaps = []
        for sequence_record in sequences:
            sequence = sequence_record.seq
            seqlen = len(sequence)
            gaps.append(str(sequence).count("-"))
        gappedness = [x/float(seqlen) for x in gaps]
        average_gappedness = sum(gappedness)/float(len(gappedness))
        if average_gappedness < threshold:
            passed.append(gene)
    return passed
        

def check_diversity(genelist, threshold, outdir):
    passed = []
    for gene in genelist:
        filename = outdir + "aligned_gene_sequences/" + gene +".aln.fas"
        sparse_matrix, consensus, seq_names = calculate_snp_matrix(filename)
        pairwise_dists = calculate_distance_matrix(sparse_matrix, consensus, 
                                                 "dist", False)
        mean_dist = pairwise_dists.mean()
        diversity_stat = mean_dist/len(consensus)
        if diversity_stat < threshold:
            passed.append(gene)
    return passed
        

        
        


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description="Filter core genes for more accurate tree-building")
    
    parser.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of the Panaroo output directory")
    
    parser.add_argument("--core_threshold",
                      dest="core",
                      help="Core-genome sample threshold (default=0.95)",
                      type=float,
                      default=0.95)  
    
    parser.add_argument("-g",
                        "--gap-threshold",
                        dest="gap",
                        help="""Threshold of gaps/isolate at which genes are 
                             discarded from core""",
                        default=0.7)
    
    parser.add_argument("-d",
                        "--distance-threshold",
                        dest="dist",
                        help="""Average pairwise distance threshold 
                        [snps/length] at which genes are discarded from core""",
                        default = 0.05)

    args = parser.parse_args()
    
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    
    # Load isolate names
    seen = set()
    isolate_names = []
    with open(args.output_dir + "gene_data.csv", 'r') as infile:
        next(infile)
        for line in infile:
            iso = line.split(",")[0]
            if iso not in seen:
                isolate_names.append(iso)
                seen.add(iso)
    #load graph
    G = nx.read_gml(args.output_dir + "final_graph.gml")
    
    #Identify core
    core_nodes = get_core_gene_nodes(G, args.core, len(isolate_names))
    core_names = [G.nodes[x]["name"] for x in core_nodes]
    
    #check gapppedness
    passed_gaps = check_gappedness(core_names, args.gap, args.output_dir)
    #check distance
    passed_distance = check_diversity(passed_gaps, args.dist, args.output_dir)
    #output filtered core genome
    concatenate_core_genome_alignments(args.output_dir + 'aligned_gene_sequences/',
                                       passed_distance, args.output_dir, 
                                       "filtered_core_gene")
                                      
    