import os

from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

#This module contains all functions necessary to write output
def write_rm_estimate(rm_regression, output_dir):
    outline1 = "Collection r/m estimate: " + str(rm_regression[0])
    outline2 = "Standard Error of Estimate: " + str(rm_regression[1])
    with open(output_dir + "rm_estimate.txt", 'w+') as outhandle:
        outhandle.write(outline1 + '\n')
        outhandle.write(outline2 + '\n')
    return True


def remove_recombinant_seqs(recombinations, out_dir, alignment_dir, all_genes):
    for gene in all_genes:
        alignment = os.path.join(alignment_dir, gene + ".aln.fas")
        sequences = list(SeqIO.parse(alignment, 'fasta'))
        sequence_names = [x.id for x in sequences]
        
        for recombinant in recombinations.get(gene, []):
            fasta_ids = [i for i in sequence_names if recombinant in i]
            indexes2remove = []
            for fid in fasta_ids:
                indexes2remove.append(sequence_names.index(fid))
            for index2remove in sorted(indexes2remove, reverse=True):
                del sequences[index2remove]
                del sequence_names[index2remove]
        if len(sequences) > 0:
            outname = os.path.join(
                out_dir,
                "recombination_free_aligned_genes",
                gene + ".aln.fas",
            )
            SeqIO.write(sequences, outname, 'fasta')
        elif len(sequences) == 0:
            print(f"All sequences removed in {gene}")
            
    
    return True

def get_core_gene_nodes(G, threshold, num_isolates):
    # Get the core genes based on percent threshold
    core_nodes = []
    for node in G.nodes():
        if float(G.nodes[node]["size"]) / float(num_isolates) >= threshold:
            core_nodes.append(node)
    return core_nodes

def write_alignment_header(alignment_list, outdir, filename):
    out_entries = []
    # Set the tracking variables for gene positions
    gene_start = 1
    gene_end = 0
    for gene in alignment_list:
        # Get length and name from one sequence in the alignment
        # Set variables that need to be set pre-output
        gene_end += gene[2]
        gene_name = gene[0]
        # Create the 3 line feature entry
        gene_entry1 = (
            "FT   feature         " + str(gene_start) + ".." + str(gene_end) + "\n"
        )
        gene_entry2 = "FT                   /label=" + gene_name + "\n"
        gene_entry3 = "FT                   /locus_tag=" + gene_name + "\n"
        gene_entry = gene_entry1 + gene_entry2 + gene_entry3
        # Add it to the output list
        out_entries.append(gene_entry)
        # Alter the post-output variables
        gene_start += gene[2]
    # Create the header and footer
    header = (
        "ID   Genome standard; DNA; PRO; 1234 BP.\nXX\nFH   Key"
        + "             Location/Qualifiers\nFH\n"
    )
    footer = (
        "XX\nSQ   Sequence 1234 BP; 789 A; 1717 C; 1693 G; 691 T;" + " 0 other;\n//\n"
    )
    # open file and output
    with open(outdir + filename, "w+") as outhandle:
        outhandle.write(header)
        for entry in out_entries:
            outhandle.write(entry)
        outhandle.write(footer)

    return True

def concatenate_core_genome_alignments(core_names, output_dir):

    alignments_dir = output_dir + "/recombination_free_aligned_genes/"
    # Open up each alignment that is associated with a core node
    alignment_filenames = os.listdir(alignments_dir)
    core_filenames = [
        x for x in alignment_filenames if x.split('.')[0] in core_names
    ]
    
    #Read in all these alignments
    gene_alignments = []
    isolates = set()
    for filename in core_filenames:
        gene_name = os.path.splitext(os.path.basename(filename))[0]
        alignment = AlignIO.read(alignments_dir + filename, "fasta")
        gene_dict = {}
        for record in alignment:
            if len(gene_dict)<1:
                gene_length = len(record.seq)

            if record.id[:3] == "_R_":
                record.id = record.id[3:]
            genome_id = record.id.split(";")[0]
            
            if genome_id in gene_dict:
                if str(record.seq).count("-") < str(gene_dict[genome_id][1]).count("-"):
                    gene_dict[genome_id] = (record.id, record.seq)
            else:
                gene_dict[genome_id] = (record.id, record.seq)
            
            isolates.add(genome_id)
        gene_alignments.append((gene_name, gene_dict, gene_length))
    # Combine them
    isolate_aln = []
    for iso in isolates:
        seq = ""
        for gene in gene_alignments:
            if iso in gene[1]:
                seq += gene[1][iso][1]
            else:
                seq += "-" * gene[2]
        isolate_aln.append(SeqRecord(seq, id=iso, description=""))

    # Write out the two output files
    SeqIO.write(isolate_aln, 
                output_dir + "recombination_free_core_gene_alignment.aln", 
                "fasta")
    write_alignment_header(gene_alignments, 
                           output_dir, 
                           "recombination_free_core_alignment_header.embl")

    return core_filenames
