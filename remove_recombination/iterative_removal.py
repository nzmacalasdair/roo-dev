import subprocess

from Bio import Phylo

from scripts.remove_recombination.read_panout import get_all_pairwise_diffs

def get_closest_pair(tree):
    clades = list(tree.find_clades(terminal=True))
    all_pairs = get_pairs(clades)
    min_dist = 9999999
    min_pair = None
    for pair in all_pairs:
        pair_distance = tree.distance(pair[0], pair[1])
        if pair_distance < min_dist:
            min_dist = pair_distance
            min_pair = pair
    return min_pair

def run_iqtree(model, output_dir, threads):
    alignmentfile = output_dir + "core_gene_alignment_filtered.aln"
    command = "iqtree -st DNA -m " + model + " -nt " + threads + " -s "
    command += alignmentfile
    result = subprocess.Popen(command)
    return True
    
    

def get_gpa_subset(gpa, isolates):
    for x in range(len(isolates)):
        if x == 0:
            continue
        elif x ==1:
            genes_in_common = (gpa[isolates[x-1]==1) & (gpa[isolates[x]==1)
        else:
            genes_in_common = genes_in_common | (gpa[isolates[x]==1)

    gene_names_in_common= gpa[genes_in_common].index

    return gene_names_in_common


