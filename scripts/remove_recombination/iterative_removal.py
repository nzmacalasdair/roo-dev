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
    
    

def get_core_subset(isolates, gpa):
    None 