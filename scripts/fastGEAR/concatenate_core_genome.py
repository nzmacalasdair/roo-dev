import os

import numpy as np
import networkx as nx

from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def calc_hc(col_counts):
    with np.errstate(divide='ignore', invalid='ignore'):
        col_counts = col_counts/np.sum(col_counts,0)
        hc = -np.nansum(col_counts[0:4,:]*np.log(col_counts[0:4,:]), 0)
    return(np.sum((1-col_counts[4,:]) * hc)/np.sum(1-col_counts[4,:]))

def update_col_counts(col_counts, s):
    s = np.fromstring(s.lower(), dtype=np.int8)
    s[(s!=97) & (s!=99) & (s!=103) & (s!=116)] = 110
    col_counts[0,s==97] += 1
    col_counts[1,s==99] += 1
    col_counts[2,s==103] += 1
    col_counts[3,s==116] += 1
    col_counts[4,s==110] += 1
    return (col_counts)

def get_core_gene_nodes(G, threshold, num_isolates):
    #Get the core genes based on percent threshold
    core_nodes = []
    for node in G.nodes():
        if float(G.nodes[node]["size"]) / float(num_isolates) > threshold:
            core_nodes.append(node)
    return core_nodes


def concatenate_core_genome_alignments(core_names, output_dir, output_prefix, 
                                       hc_threshold):
    #Open up each alignment that is assosciated with a core 
    alignments_dir = output_dir + "fastGEAR_recombination_free_aligned_genes/"
    alignment_filenames = os.listdir(alignments_dir)
    core_filenames = [
        x for x in alignment_filenames if x.split('.')[0] in core_names
    ]
    #Read in all these alginemnts
    gene_alignments = []
    isolates = set()
    for filename in core_filenames:
        gene_name = os.path.splitext(os.path.basename(filename))[0]
        alignment = AlignIO.read(alignments_dir + filename, 'fasta')
        gene_dict = {}
        for record in alignment:
            if len(gene_dict)<1:
                gene_length = len(record.seq)
                col_counts = np.zeros((5,gene_length), dtype=float)
            col_counts = update_col_counts(col_counts, str(record.seq))

            if record.id[:3] == "_R_":
                record.id = record.id[3:]
            genome_id = record.id.split(";")[0]
            
            if genome_id in gene_dict:
                if str(record.seq).count("-") < str(gene_dict[genome_id][1]).count("-"):
                    gene_dict[genome_id] = (record.id, record.seq)
            else:
                gene_dict[genome_id] = (record.id, record.seq)

            gene_length = len(record.seq)
            isolates.add(genome_id)
        gene_alignments.append((gene_name, gene_dict, gene_length,
                                calc_hc(col_counts)))
    #Combine them
    isolate_aln = []
    for iso in isolates:
        seq = ""
        for gene in gene_alignments:
            if iso in gene[1]:
                seq += gene[1][iso][1]
            else:
                seq += "-" * gene[2]
        isolate_aln.append(SeqRecord(seq, id=iso, description=""))

    #Write out the two output files
    SeqIO.write(isolate_aln, output_dir + output_prefix + '_alignment.aln', 'fasta')

    write_alignment_header(gene_alignments, output_dir, output_prefix)
    
    # Calculate threshold for h.
    if hc_threshold is None:
        allh = np.array([gene[3] for gene in gene_alignments])
        q = np.quantile(allh, [0.25,0.75])
        hc_threshold = max(0.01, q[1] + 1.5*(q[1]-q[0]))
        print(f"Entropy threshold automatically set to {hc_threshold}.")

    isolate_aln = []
    keep_count = 0 
    for iso in isolates:
        seq = ""
        for gene in gene_alignments:
            if gene[3]<=hc_threshold:
                keep_count += 1
                if iso in gene[1]:
                    seq += gene[1][iso][1]
                else:
                    seq += "-" * gene[2]
        isolate_aln.append(SeqRecord(seq, id=iso, description=""))

    with open(output_dir + 'fastGEAR_recombination_free_alignment_entropy.csv',
              'w') as outfile:
        for g in gene_alignments:
            outfile.write(str(g[0]) + ',' + str(g[3]) + '\n')

    # Write out the two output files
    SeqIO.write(isolate_aln, output_dir + output_prefix + "_filtered.aln", "fasta")
    write_alignment_header(
        [g for g in gene_alignments if g[3]<=hc_threshold],
        output_dir, output_prefix + "_filtered",
    )

    print(f"""{keep_count/len(isolates)} out of {len(gene_alignments)} genes kept
          in fastGEAR recombination free filtered core genome""")

    return core_filenames

def write_alignment_header(alignment_list, outdir, prefix):
    out_entries = []
    #Set the tracking variables for gene positions
    gene_start = 1
    gene_end = 0
    for gene in alignment_list:
        #Get length and name from one sequence in the alignment
        #Set variables that need to be set pre-output
        gene_end += gene[2]
        gene_name = gene[0]
        #Create the 3 line feature entry
        gene_entry1 = "FT   feature         " + str(gene_start) + ".." + str(
            gene_end) + '\n'
        gene_entry2 = "FT                   /label=" + gene_name + '\n'
        gene_entry3 = "FT                   /locus_tag=" + gene_name + '\n'
        gene_entry = gene_entry1 + gene_entry2 + gene_entry3
        #Add it to the output list
        out_entries.append(gene_entry)
        #Alter the post-output variables
        gene_start += gene[2]
    #Create the header and footer
    header = ("ID   Genome standard; DNA; PRO; 1234 BP.\nXX\nFH   Key" +
              "             Location/Qualifiers\nFH\n")
    footer = ("XX\nSQ   Sequence 1234 BP; 789 A; 1717 C; 1693 G; 691 T;" +
              " 0 other;\n//\n")
    #open file and output
    with open(outdir + prefix+"_alignment_header.embl", "w+") as outhandle:
        outhandle.write(header)
        for entry in out_entries:
            outhandle.write(entry)
        outhandle.write(footer)
    return True



if __name__ == '__main__':
    import argparse
    description = 'Concatenate (recombination-free) core gene alignmets'
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group('Input/output')

    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of the Panaroo output directory",
                         )
    parser.add_argument("--core_threshold",
                      dest="core",
                      help="Core-genome sample threshold (default=0.95)",
                      type=float,
                      default=0.95)
    parser.add_argument("--core_entropy_filter",
                      dest="hc_threshold",
                      help=("""Manually set the Block Mapping and Gathering with  
                            Entropy (BMGE) filter. Can be between 0.0 and 1.0. By   
                            "default this is set using the Tukey outlier method."""),
                      type=float,
                      default=None)
    parser.add_argument("--prefix",
                        dest="prefix",
                        default="fastGEAR_recombination_free_core_gene",
                        help="Provide custom outprefix for core alignment")
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
    concatenate_core_genome_alignments(core_names, args.output_dir,
                                       args.prefix, args.hc_threshold)