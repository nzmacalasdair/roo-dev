import os
import shutil

from Bio import SeqIO
from Bio.Seq import Seq

def parse_recombination_file(filename):
    gene_recombinations = {}
    with open(filename, 'r') as inhandle:
        lines = inhandle.read().splitlines()
    header = "Gene,Age,RecombStart,RecombEnd,Score,Donor_lineage,Recipient_strain,Donor_group,Strain_group"
    for line in lines:
        if line == header:
            continue
        splitline = line.split(",")
        #fastGEAR defaults to 1-indexing, subtract 1 for 0-indexing
        start = int(splitline[2]) - 1
        stop = int(splitline[3]) - 1
        affected_isolates = splitline[6].split("~")
        gene = splitline[0]
        recomb_info = {"start":start, "stop": stop, "isolates" : affected_isolates}
        gene_recombinations[gene] = gene_recombinations.get(gene, []) + [recomb_info]
    
    return gene_recombinations
    
def remove_recombinations(recombinations, alignments, out_dir):
    for gene in recombinations:
        alignment = out_dir + "aligned_gene_sequences/" + gene +".aln.fas"
        sequences = list(SeqIO.parse(alignment, 'fasta'))
        sequence_names = [x.id for x in sequences]
        for recombination in recombinations[gene]:
            length = recombination["stop"] - recombination["start"]
            if type(recombination["isolates"]) != list:
                recombination["isolates"] = [recombination["isolates"]]
            for isolate in recombination["isolates"]:
                isolate_index = sequence_names.index(isolate)
                recombinant_sequence = sequences[isolate_index].seq
                new_sequence = str(recombinant_sequence[:recombination["start"]])
                new_sequence += length * "-"
                new_sequence += str(recombinant_sequence[recombination["stop"]:])
                new_sequence = Seq(new_sequence)
                
                sequences[isolate_index].seq = new_sequence
        outname = out_dir + "recombination_free_aligned_genes/" + gene +".aln.fas"
        SeqIO.write(sequences, outname, 'fasta')
    
    return True
                
        

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description="")
    
    io_opts = parser.add_argument_group('Input/output')

    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of the Panaroo output directory")
    parser.add_argument('--recombinations',
                        default=None,
                        dest="recombs",
                        help="""Path to a .csv containing the recombination data
                        from a fastGEAR run summarised by summarise_pgfG.py""")
    parser.add_argument('--filtered',
                        dest="filt",
                        action="store_true",
                        help="Filtered or unfiltered alignments")
    args = parser.parse_args()

    #prepare files, output folder    
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    os.mkdir(args.output_dir + "recombination_free_aligned_genes/")
    
    #read in data
    if args.filt:
        alignment_list = os.listdir(args.output_dir + "filtered_gene_alignments/")
    else:        
        alignment_list = os.listdir(args.output_dir + "aligned_gene_sequences/")
    gene_recombinations = parse_recombination_file(args.recombs)
    
    #remove recombination regions from alignments
    remove_recombinations(gene_recombinations, alignment_list, args.output_dir)
    
    #Copy over genes with no recombination 
    recomb_free_alignments = os.listdir(args.output_dir + "recombination_free_aligned_genes/")
    no_recomb_alignments = set(alignment_list) - set (recomb_free_alignments)
    
    for gene in no_recomb_alignments:
        shutil.copy(args.output_dir + "aligned_gene_sequences/" + gene, 
                    args.output_dir + "recombination_free_aligned_genes/")
    