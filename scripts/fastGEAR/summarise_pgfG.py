#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 16:57:00 2019

@author: nm12
"""

import os

class Recomb_event(object):
    # A fastGEAR recombination event class, just for fun
    
    def __init__(self, gene, age, position, donor_recipient, dg_rg, score):
        self.gene = gene
        self.age = age
        self.position = position
        self.donor_recipient = donor_recipient
        self.donorgroup_recipientgroup = dg_rg
        self.score = score
        
    def start(self):
        return self.position[0]
    def end(self):
        return self.position[1]
    def donor_group(self):
        return self.donorgroup_recipientgroup[0]
    def recipient_group(self):
        return self.donorgroup_recipientgroup[1]


def parse_lineage_information(filehandle):
    linenumber = 0
    lineages = {}
    with open(filehandle) as lineage_info:
        for line  in lineage_info:
            linenumber +=1
            if linenumber > 1:
                splitline = line.split()                
                lineages[splitline[1]] = lineages.get(splitline[1], []) + [splitline[-1]]
    return lineages

def parse_cluster_remap(filehandle):
    clustermap = {}
    with open(filehandle, 'r') as temp:
        lines = temp.read().splitlines()
    for line in lines:
        splitline = line.split(",")
        clustermap[splitline[1]] = splitline[0]
    return clustermap

def parse_recent_recombinations(gene_path, clustering):
    gene_recombinations = []
    if os.path.isfile(gene_path + "/output/recombinations_recent.txt") == False:
        errmsg = " fastGEAR recent recombination detection failed"
        print(gene_path.split("/")[-1] + errmsg)
        raise(ValueError)
        return [None]
    with open(gene_path + "/output/recombinations_recent.txt") as handle:
        recombinations = handle.read().splitlines() 
        linenumber = 0
        if len(recombinations) == 2:
            gene_recombinations.append(gene_path.split("/")[-1])
            return gene_recombinations
        for line in recombinations:
            linenumber += 1
            if linenumber > 2:
                splitline = line.split()
                lineages = parse_lineage_information(gene_path + "/output/lineage_information.txt")
                if type(clustering) == dict:
                    cluster_map = parse_cluster_remap(gene_path + "/cluster_mapping.txt")
                position = (splitline[0], splitline[1])
                clean_name = "_".join(splitline[5].split("_")[:-1])
                clean_name= clean_name.split(";")[0]
                donor_recipient = (lineages.get(splitline[2], "Ukwn"), splitline[5])
                if type(clustering) == dict:
                    dg_rg = (cluster_map.get(splitline[2], "Ukwn"), clustering.get(clean_name))
                    recombination = Recomb_event(gene_path.split("/")[-1], 
                                                 "recent", position , donor_recipient, 
                                                 dg_rg, float(splitline[4]))
                else:
                    recombination = Recomb_event(gene_path.split("/")[-1], 
                                                 "recent", position , donor_recipient, 
                                                 ["NA", "NA"], float(splitline[4]))
                gene_recombinations.append(recombination)  
    return gene_recombinations

def parse_ancestral_recombinations(gene_path, clustering):
    gene_recombinations = []
    if os.path.isfile(gene_path + "/output/recombinations_ancestral.txt") == False:
        errmsg = " fastGEAR ancestral recombination detection failed!"
        print(gene_path.split("/")[-1] + errmsg)
        raise(ValueError)
        return [None]
    with open(gene_path + "/output/recombinations_ancestral.txt") as handle:
        recombinations = handle.read().splitlines()
        linenumber = 0
        if len(recombinations) == 2:
            gene_recombinations.append(gene_path.split("/")[-1])
            return gene_recombinations
        for line in recombinations:
            linenumber += 1
            if linenumber > 2:
                splitline = line.split()
                lineages = parse_lineage_information(gene_path + "/output/lineage_information.txt")
                if type(clustering) == dict:
                    cluster_map = parse_cluster_remap(gene_path + "/cluster_mapping.txt")
                position = (splitline[0], splitline[1])
                donor_recipient = (lineages.get(splitline[2], "Ukwn"), lineages.get(splitline[3]))
                if type(clustering) == dict:
                    dg_rg = (cluster_map.get(splitline[2], "Ukwn"), cluster_map.get(splitline[3], "Ukwn"))
                    recombination = Recomb_event(gene_path.split("/")[-1], 
                                                 "ancestral", position , donor_recipient, 
                                                 dg_rg, float(splitline[4]))
                else:
                    recombination = Recomb_event(gene_path.split("/")[-1], 
                                                 "ancestral", position , donor_recipient, 
                                                 ["NA", "NA"], float(splitline[4]))                    
                gene_recombinations.append(recombination)  
    return gene_recombinations

def summarise_fastGEAR(result_directory, clustering_file):
    recombinations = []
    
    gene_list = os.listdir(result_directory)
    
    if type(clustering_file) == str:
        clustering = {}    
        with open(clustering_file, 'r') as temphandle:
            parsed = temphandle.read().splitlines()[1:]
            for line in parsed:
                splitline = line.split(',')
                clustering[splitline[0]] = splitline[1]
    else:
        clustering = None
        
    for gene in gene_list:
        #ignore any files in the results directory
        if os.path.isfile(gene):
            continue
        recombinations += parse_recent_recombinations(result_directory +'/'+ gene, clustering)
        recombinations += parse_ancestral_recombinations(result_directory +'/'+ gene, clustering)
    
    return recombinations

def write_output(summary, outname):
    
    recombinations = filter(None, summary)
    at_least_one_no_recombination = []
    no_recombs = open("genes_without_recombination.txt", 'w+')
    
    if outname[-4:] != ".csv":
        outname = outname + ".csv"

    output = open(outname, 'w+')
    
    header = "Gene,Age,RecombStart,RecombEnd,Score,Donor_lineage,"
    header += "Recipient_strain,Donor_group,Strain_group\n"    
    output.write(header)
        
    for index in range(len(recombinations)):
        if type(recombinations[index]) == Recomb_event:
            outline1 = recombinations[index].gene +','+ recombinations[index].age +','+ recombinations[index].start() +','+ recombinations[index].end()+','
            outline2 = str(recombinations[index].score) +','+ "~".join(recombinations[index].donor_recipient[0])+','   
            outline3 = "~".join(recombinations[index].donor_recipient[1]) +','+ str(recombinations[index].donor_group()) +','+ str(recombinations[index].recipient_group())
            outline = outline1 + outline2 + outline3 + '\n'
            output.write(outline)
    
        else:
            at_least_one_no_recombination.append(recombinations[index])

    output.close()
    neither_recombination = list(set([x for x in at_least_one_no_recombination if at_least_one_no_recombination.count(x) > 1]))
    no_recombs = open(outname.split(".")[0] + "_genes_without_recombination.txt", 'w+')
    for gene in neither_recombination:
        no_recombs.write(gene + '\n')
    no_recombs.close()
    return None

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('result_dir', 
                        help="Directory containing fastGEAR results on a pan-genome")
    parser.add_argument("-o",
                        "--output",
                        dest="output",
                        default="fastGEAR_recombinations",
                        help="Name of output csv, default: fastGEAR_recombinations")
    parser.add_argument('--clustering',
                        default=None,
                        help="""Path to a .csv containing cluster 
                        information if a separate genome-level clustering (such 
                        as PopPUNK) was used to run fastGEAR on every gene""")
    args = parser.parse_args()
    
    
    results = os.path.dirname(os.path.abspath(args.result_dir)) + '/' + args.result_dir.strip('/')
    fG_summary = summarise_fastGEAR(results, args.clustering)
    write_output(fG_summary, args.output)
