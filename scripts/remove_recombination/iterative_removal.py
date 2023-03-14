from Bio import Phylo

def get_closest_pair(tree):
    clades = list(tree.find_clades(terminal=true))
    all_pairs = get_pairs(clades)
    min_dist = 9999999
    min_pair = None
    for pair in all_pairs:
        pair_distance = tree.distance([pair[0], pair[1])
        if pair_distance < min_dist:
            min_dist = pair_distance
            min_pair = pair
    return min_pair


