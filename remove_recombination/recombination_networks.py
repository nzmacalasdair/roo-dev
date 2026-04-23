import networkx as nx

def build_recombination_network(gene_recombination_dic):
    gene_network = nx.Graph()
    for recombinant_pair in gene_recombination_dic:
        recombination_list = recombinant_pair.split("-")
        gene_network.add_edge(*recombination_list)
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