import networkx as nx

import sys

def build_recombination_network(recombinant_pairs, pair_index_to_isolates):
    gene_network = nx.Graph()
    for pair_idx in recombinant_pairs:
        gene_network.add_edge(*pair_index_to_isolates[pair_idx])
    return gene_network

def identify_genuine_recombinants(G):
    genuine_isolates= []
    #get disconnected componenets 
    for component in nx.connected_components(G):
        H = G.subgraph(component)
        #Only look at hits with a reasonable amount of support
        if H.number_of_nodes() < 4:
            continue
        #Find the k-clustering threshold
        core = nx.core_number(H)
        max_k = max(core.values())
        
        #Look at max degrees
        degrees = dict(H.degree())
        sorted_degs = sorted(degrees.values(), reverse=True)
        
        #k-core cannot find singletons
        #Is there only one isolate? Check max degree node(s)
        
        #get first non-max node
        second_degree = next((d for d in sorted_degs if d < sorted_degs[0]), None)
        
        if second_degree == None:
            #No non-max node, entire graph is equally connected
            #Gene is suspicious, get rid of all isolates in this component
            genuine_isolates += list(H.nodes)
            continue
        
        if (sorted_degs[0] / sorted_degs[second_degree]) >= 3:
            # one or more obvious hubs
            top_nodes = [n for n, d in degrees.items() if d == sorted_degs[0]]
            genuine_isolates += top_nodes
        else:
            #grab the k-core nodes
            top_nodes = [n for n, k in core.items() if k == max_k]
            #if the k-core is more than 50% of the isolates, disregard
            if len(top_nodes) > (0.5 * H.number_of_nodes()): 
                genuine_isolates += top_nodes
            else:
                continue
    return genuine_isolates
