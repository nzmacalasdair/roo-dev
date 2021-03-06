import os
import tempfile
import shutil
import argparse

import networkx as nx
from tqdm import tqdm
from joblib import Parallel, delayed
from collections import defaultdict, Counter
import math
import numpy as np

from .isvalid import *
from .__init__ import __version__
from .cdhit import run_cdhit
from .clean_network import collapse_families, collapse_paralogs
from .generate_output import *


def make_list(inp):
    if not type(inp) == list:
        inp = [inp]
    return inp


def load_graphs(graph_files, n_cpu=1):
    graphs = []

    for graph_file in graph_files:
        if not os.path.isfile(graph_file):
            print("Missing:", graph_file)
            raise RuntimeError("Missing graph file!")

    graphs = [nx.read_gml(graph_file) for graph_file in tqdm(graph_files)]

    return graphs


def cluster_centroids(graphs,
                      outdir,
                      len_dif_percent=0.95,
                      identity_threshold=0.95,
                      n_cpu=1):

    # create the files we will need
    temp_input_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_input_file.close()
    temp_output_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_output_file.close()

    # create input for cdhit
    with open(temp_input_file.name, 'w') as outfile:
        for i, G in enumerate(graphs):
            for node in G.nodes():
                outfile.write(">" + str(i) + "_" + node + '\n')
                seqs = G.nodes[node]["protein"].split(";")
                seqs = [s for s in seqs if "*" not in s]
                outfile.write(max(seqs, key=len) + "\n")

    # Run cd-hit
    run_cdhit(temp_input_file.name,
              temp_output_file.name,
              id=identity_threshold,
              s=len_dif_percent,
              accurate=True,
              n_cpu=n_cpu)

    # Process output
    clusters = []
    with open(temp_output_file.name + ".clstr", 'rU') as infile:
        c = []
        for line in infile:
            if line[0] == ">":
                clusters.append(c)
                c = []
            else:
                temp = line.split(">")[1].split("...")[0].split("_")
                centroid = (int(temp[0]), temp[1])
                c.append(centroid)
        clusters.append(c)
    clusters = clusters[1:]

    # remove temporary files
    os.remove(temp_input_file.name)
    os.remove(temp_output_file.name)
    os.remove(temp_output_file.name + ".clstr")

    # need to add new centroidsssss!!!

    return clusters


def simple_merge_graphs(graphs, clusters):
    # Here, we only merge nodes that don't conflict

    # first rename each graphs nodes in preperation for merge
    # get mapping
    mapping = defaultdict(dict)
    nnodes = 0
    reverse_mapping = defaultdict(list)

    merge_centroids = {}
    centroid_context = defaultdict(list)
    for c, cluster in enumerate(clusters):
        c = str(c)
        nnodes += 1
        graph_clust_count = Counter()
        for n in cluster:
            graph_clust_count[n[0]] += 1
        non_conflicting_nodes = [
            n for n in cluster if graph_clust_count[n[0]] < 2
        ]
        conflicting_nodes = [n for n in cluster if graph_clust_count[n[0]] > 1]
        for n in non_conflicting_nodes:
            mapping[n[0]][n[1]] = nnodes
            merge_centroids[nnodes] = c
            reverse_mapping[nnodes].append((n[0], nnodes))
        for n in conflicting_nodes:
            nnodes += 1
            mapping[n[0]][n[1]] = nnodes
            merge_centroids[nnodes] = c
            centroid_context[c].append([nnodes, n[0]])
            reverse_mapping[nnodes].append((n[0], nnodes))

    # rename
    for i, G in enumerate(graphs):
        nx.relabel_nodes(G, mapping[i], copy=False)

    # merge graphs
    merged_G = nx.compose_all(graphs)

    # fix up node attributes
    for node in merged_G.nodes():
        size = 0
        members = set()
        lengths = []
        centroid = []
        seqIDs = set()
        protein = []
        dna = []
        annotation = []
        description = []
        paralog = False
        hasEnd = False
        mergedDNA = False

        for prev in reverse_mapping[node]:
            size += graphs[prev[0]].nodes[prev[1]]['size']
            members |= set([
                str(prev[0]) + "_" + str(m)
                for m in make_list(graphs[prev[0]].nodes[prev[1]]['members'])
            ])
            lengths += make_list(graphs[prev[0]].nodes[prev[1]]['lengths'])
            centroid += [
                str(prev[0]) + "_" + str(m) for m in make_list(graphs[
                    prev[0]].nodes[prev[1]]['centroid'].split(";"))
            ]
            seqIDs |= set([
                str(prev[0]) + "_" + d
                for d in make_list(graphs[prev[0]].nodes[prev[1]]['seqIDs'])
            ])
            protein += make_list(
                graphs[prev[0]].nodes[prev[1]]['protein'].split(";"))
            dna += make_list(graphs[prev[0]].nodes[prev[1]]['dna'].split(";"))
            annotation += make_list(
                graphs[prev[0]].nodes[prev[1]]['annotation'])
            description += make_list(
                graphs[prev[0]].nodes[prev[1]]['description'])
            paralog = (paralog or graphs[prev[0]].nodes[prev[1]]['paralog'])
            hasEnd = (paralog or graphs[prev[0]].nodes[prev[1]]['hasEnd'])
            mergedDNA = (paralog
                         or graphs[prev[0]].nodes[prev[1]]['mergedDNA'])

        merged_G.nodes[node]['size'] = size
        merged_G.nodes[node]['members'] = set(members)
        merged_G.nodes[node]['lengths'] = lengths
        merged_G.nodes[node]['prevCentroids'] = str(
            merge_centroids[node])  #";".join(centroid)
        merged_G.nodes[node]['seqIDs'] = set(seqIDs)
        merged_G.nodes[node]['hasEnd'] = hasEnd
        merged_G.nodes[node]['dna'] = [max(dna, key=len)]
        merged_G.nodes[node]['protein'] = [max(protein, key=len)]
        merged_G.nodes[node]['annotation'] = ";".join(annotation)
        merged_G.nodes[node]['description'] = ";".join(description)
        merged_G.nodes[node]['paralog'] = paralog
        merged_G.nodes[node]['mergedDNA'] = mergedDNA
        merged_G.nodes[node]['centroid'] = [str(merge_centroids[node])]

    # fix longcentroid
    if len(merged_G.nodes[node]['centroid']) != len(
            merged_G.nodes[node]['protein']):
        print(merged_G.nodes[node]['protein'])
        print(merged_G.nodes[node]['centroid'])
        raise RuntimeError("protein/centroid count mismatch!")

    for node in merged_G.nodes():
        merged_G.nodes[node]['longCentroidID'] = max([
            (len(s), sid) for s, sid in zip(merged_G.nodes[node]['protein'],
                                            merged_G.nodes[node]['centroid'])
        ])
        merged_G.nodes[node]['maxLenId'] = max([
            (len(s), index)
            for s, index in zip(merged_G.nodes[node]['dna'],
                                range(len(merged_G.nodes[node]['dna'])))
        ])[1]

    # fix up edge attributes
    for edge in merged_G.edges():
        merged_G[edge[0]][edge[1]]['weight'] = 0
        merged_G[edge[0]][edge[1]]['members'] = set()

        for prev1 in reverse_mapping[edge[0]]:
            for prev2 in reverse_mapping[edge[1]]:
                if prev1[0] == prev2[0]:  #same graph
                    if graphs[prev1[0]].has_edge(prev1[1], prev2[1]):
                        merged_G[edge[0]][edge[1]]['weight'] += 1
                        merged_G[edge[0]][edge[1]]['members'] |= set([
                            str(prev1[0]) + "_" + str(m) for m in graphs[
                                prev1[0]][prev1[1]][prev2[1]]['members']
                        ])

    return merged_G, centroid_context


def get_options():
    import argparse

    description = 'Merge independent runs of Panaroo'
    parser = argparse.ArgumentParser(description=description,
                                     prog='panaroo_merge_graphs')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "-d",
        "--directories",
        dest="directories",
        required=True,
        help="Location of seperate Panaroo output directories",
        nargs='+')

    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of a new output directory",
                         type=lambda x: is_valid_folder(parser, x))

    matching = parser.add_argument_group('Matching')
    matching.add_argument("-c",
                          "--threshold",
                          dest="id",
                          help="sequence identity threshold (default=0.95)",
                          default=0.95,
                          type=float)
    matching.add_argument(
        "-f",
        "--family_threshold",
        dest="family_threshold",
        help="protein family sequence identity threshold (default=0.7)",
        default=0.7,
        type=float)
    matching.add_argument("--len_dif_percent",
                          dest="len_dif_percent",
                          help="length difference cutoff (default=0.95)",
                          default=0.95,
                          type=float)

    parser.add_argument(
        "--min_edge_support_sv",
        dest="min_edge_support_sv",
        help=(
            "minimum edge support required to call structural variants" +
            " in the presence/absence sv file (default=max(2, 0.01*n_samples))"
        ),
        type=int)

    # Other options
    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)
    parser.add_argument("--verbose",
                        dest="verbose",
                        help="print additional output",
                        action='store_true',
                        default=False)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()
    return (args)


def main():
    args = get_options()

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    args.directories = [os.path.join(d, "") for d in args.directories]

    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")

    print(
        "Merging graphs is still under active development and may change frequently!"
    )

    # Load graphs
    print("Loading graphs...")
    graphs = load_graphs([d + "final_graph.gml" for d in args.directories],
                         n_cpu=args.n_cpu)

    # cluster centroids
    print("Clustering centroids...")
    clusters = cluster_centroids(graphs=graphs,
                                 outdir=temp_dir,
                                 len_dif_percent=args.len_dif_percent,
                                 identity_threshold=args.id,
                                 n_cpu=args.n_cpu)

    # perform initial merge
    print("Performing inital merge...")
    G, centroid_contexts = simple_merge_graphs(graphs, clusters)

    print("Number of nodes in merged graph: ", G.number_of_nodes())

    # collapse gene families/paralogs at successively lower thresholds
    print("Collapsing paralogs...")
    G = collapse_paralogs(G, centroid_contexts)

    print("Collapsing at DNA...")
    G = collapse_families(G,
                          outdir=temp_dir,
                          dna_error_threshold=0.98,
                          correct_mistranslations=True,
                          n_cpu=args.n_cpu,
                          quiet=(not args.verbose))[0]

    print("Collapsing at families...")
    G = collapse_families(G,
                          outdir=temp_dir,
                          family_threshold=args.family_threshold,
                          correct_mistranslations=False,
                          n_cpu=args.n_cpu,
                          quiet=(not args.verbose))[0]

    print("Number of nodes in merged graph: ", G.number_of_nodes())

    # Generate output
    print("Generating output...")
    # write out roary like gene_presence_absence.csv
    mems_to_isolates = {}
    for i, sub_G in enumerate(graphs):
        for j, iso in enumerate(sub_G.graph['isolateNames']):
            mems_to_isolates[str(i) + "_" + str(j)] = iso

    n_samples = len(mems_to_isolates)
    args.min_edge_support_sv = max(2, math.ceil(0.01 * n_samples))

    # get original annotaiton IDs
    # get original annotaiton IDs, lengts and whether or
    # not an internal stop codon is present
    orig_ids = {}
    ids_len_stop = {}
    for i, d in enumerate(args.directories):
        with open(d + "gene_data.csv", 'r') as infile:
            next(infile)
            for line in infile:
                line = line.split(",")
                orig_ids[str(i) + "_" + line[2]] = line[3]
                ids_len_stop[str(i) + "_" + line[2]] = (len(line[4]),
                                                        "*" in line[4][1:-3])

    G = generate_roary_gene_presence_absence(G,
                                             mems_to_isolates=mems_to_isolates,
                                             orig_ids=orig_ids,
                                             ids_len_stop=ids_len_stop,
                                             output_dir=args.output_dir)

    # write pan genome reference fasta file
    generate_pan_genome_reference(G,
                                  output_dir=args.output_dir,
                                  split_paralogs=False)

    # write out common structural differences in a matrix format
    generate_common_struct_presence_absence(
        G,
        output_dir=args.output_dir,
        mems_to_isolates=mems_to_isolates,
        min_variant_support=args.min_edge_support_sv)

    # add helpful attributes and write out graph in GML format
    for node in G.nodes():
        G.nodes[node]['size'] = len(G.nodes[node]['members'])
        G.nodes[node]['genomeIDs'] = ";".join(G.nodes[node]['members'])
        G.nodes[node]['members'] = list(G.nodes[node]['members'])
        G.nodes[node]['centroid'] = G.nodes[node]['prevCentroids']
        del G.nodes[node]['prevCentroids']
        G.nodes[node]['geneIDs'] = ";".join(G.nodes[node]['seqIDs'])
        G.nodes[node]['seqIDs'] = list(G.nodes[node]['seqIDs'])
        G.nodes[node]['degrees'] = G.degree[node]
        sub_graphs = list(
            set([m.split("_")[0] for m in G.nodes[node]['members']]))
        G.nodes[node]['subGraphs'] = ";".join(conv_list(sub_graphs))

    for edge in G.edges():
        G.edges[edge[0],
                edge[1]]['genomeIDs'] = ";".join(G.edges[edge[0],
                                                         edge[1]]['members'])
        G.edges[edge[0],
                edge[1]]['members'] = list(G.edges[edge[0],
                                                   edge[1]]['members'])
        sub_graphs = list(
            set([
                m.split("_")[0] for m in G.edges[edge[0], edge[1]]['members']
            ]))
        G.edges[edge[0],
                edge[1]]['subGraphs'] = ";".join(conv_list(sub_graphs))

    nx.write_gml(G, args.output_dir + "merged_final_graph.gml")

    # remove temporary directory
    shutil.rmtree(temp_dir)

    return


if __name__ == '__main__':
    main()
