import numpy as np

#import cupy as cp

from pairsnp import calculate_snp_matrix, calculate_distance_matrix

#Legacy
# def check_for_big_indel(byteseq1, byteseq2):
#     gaps1 = np.sum(byteseq1 == 45) #45 == ord("-") 
#     gaps2 = np.sum(byteseq2 == 45)
#     diff = gaps2 - gaps1
#     if (diff/gaps1.size) > 0.15:
#         return diff
#     else: 
#         return None

# def check_for_big_indel_gpu(byteseq1, byteseq2):
#     gaps1 = cp.sum(cp.asarray(byteseq1) == 45) #45 == ord("-") 
#     gaps2 = cp.sum(cp.asarray(byteseq2) == 45)
#     diff = gaps2 - gaps1
#     if (diff/gaps1.size) > 0.15:
#         return diff
#     else: 
#         return None

# def get_pairwise_differences(seq1, seq2, gpu):
#     if seq1.size != seq2.size:
#         raise ValueError("Two aligned sequences are of different lengths!")
    
#     if gpu:
#         seq1 = cp.asarray(seq1)
#         seq2 = cp.asarray(seq2)
#         if check_for_big_indel_gpu == None:
#             diffs = cp.count_nonzero(seq1^seq2)
#             length = seq1.size
#             result = (np.array([diffs.get(), length]))
#             return result
#         else:
#             mask = seq1 != 45
#             consecutive_bases = cp.diff(cp.concatenate((cp.array([0]), 
#                                                         mask.astype(int), 
#                                                         cp.array([0]))))
#             start_positions = cp.where(consecutive_bases == 1)[0]
#             end_positions = cp.where(consecutive_bases == -1)[0]
#             lengths = end_positions - start_positions
            
#             longest_seq_index = cp.argmax(lengths)
#             longest_seq_start = start_positions[longest_seq_index]
#             longest_seq_end = end_positions[longest_seq_index]
            
#             cropped_seq1 = seq1[longest_seq_start:longest_seq_end]
#             cropped_seq2 = seq2[longest_seq_start:longest_seq_end]
            
#             diffs = cp.count_nonzero(cropped_seq1^cropped_seq2)
#             length = cropped_seq1.size
#             result = (np.array([diffs.get(), length]))
#             return result
            
#     else:    
#         if check_for_big_indel == None:
#             diffs = np.count_nonzero(seq1^seq2)
#             length = seq1.size
#             result = (np.array([diffs, length]))
#             return result
#         else:
#             mask = seq1 != 45
#             consecutive_bases = np.diff(np.concatenate(([0], mask.astype(int), [0])))
#             start_positions = np.where(consecutive_bases == 1)[0]
#             end_positions = np.where(consecutive_bases == -1)[0]
#             lengths = end_positions - start_positions
            
#             longest_seq_index = np.argmax(lengths)
#             longest_seq_start = start_positions[longest_seq_index]
#             longest_seq_end = end_positions[longest_seq_index]
            
#             cropped_seq1 = seq1[longest_seq_start:longest_seq_end]
#             cropped_seq2 = seq2[longest_seq_start:longest_seq_end]
            
#             diffs = np.count_nonzero(cropped_seq1^cropped_seq2)
#             length = cropped_seq1.size
#             result = (np.array([diffs, length]))
#             return result
        
        

def get_pairs(isolate_list):
    pairs_list = [(isolate_list[i], isolate_list[j]) for i in range(len(isolate_list)) for j in range(i+1,len(isolate_list))]
    return pairs_list

#Legacy code which did comparisons pairwise
# def get_pangenome_pairwise_differences(gene_alignments, isolates_to_consider, gpu):    
#     diffs = []
#     gene_lens = []

#     #Legacy code which did comparisons pairwise, did not scale!
#     # for gene in gene_alignments:
#     #     seq1 = gene[0].get(isolates_to_consider[0], None)
#     #     seq2 = gene[0].get(isolates_to_consider[1], None)
#     #     if (seq1 is None) or (seq2 is None):
#     #         continue
#     #     diffs.append(get_pairwise_differences(seq1, seq2, gpu))
#     #     names.append(gene[1].split(".")[0])
#     # diffs = np.array(diffs)
    
    
    
#     return (diffs, gene_lens)

def check_row_for_indel(row, csr_matrix):
    rowstart = csr_matrix.indptr[row]
    rowend = csr_matrix.indptr[row+1]
    
    row_data = csr_matrix.data[rowstart:rowend]
    
    gap_count = np.sum(row_data == 110)
    return gap_count 

def instantiate_csr_row(row, csr_matrix, reference):
    rowstart = csr_matrix.indptr[row]
    rowend = csr_matrix.indptr[row+1]
    
    row_data = csr_matrix.data[rowstart:rowend]
    row_indices = csr_matrix.indices[rowstart:rowend]
    
    row_sequence = reference.copy()
    
    for index in range(len(row_data)):
        row_sequence[row_indices[index]] = row_data[index]
    
    return row_sequence
    
def get_gene_lengths(sparse_matrix, consensus):
    total_aln_len = sparse_matrix.shape[1]
    isolate_count = sparse_matrix.shape[0]
    
    nongap_lengths = np.zeros(isolate_count, dtype = int)
    
    biggap_in_ref = False
    ref_nogap_length = -1
    for row_index in range(isolate_count):       
        
        if row_index == 0: #first row is 'ref' for snps so must be done alone
            gappositions = np.where(consensus == 110)[0]
            gapcount = len(gappositions)
            if (gapcount/total_aln_len) > 0.15:
                biggap_in_ref = True
                ref_nogap_length = len(consensus[consensus!= 110])
                nongap_lengths[row_index] = ref_nogap_length
                continue
            else:
                nongap_lengths[row_index] = total_aln_len
                continue
        
        if biggap_in_ref == False:
            gapcount = check_row_for_indel(row_index, sparse_matrix)
            if gapcount/total_aln_len > 0.15: #more than 15% gaps
                length = total_aln_len - gapcount
                nongap_lengths[row_index] = length
            else: #No big indel
                nongap_lengths[row_index] = total_aln_len
        else:
            #need to fill out row to get row sequence (separate function)
            #then count the gaps in the filled out row
            #this is presumably most expensive so better to do other simpler
            #operations first
            row_sequence = instantiate_csr_row(row_index, sparse_matrix, consensus)
            gappositions = np.where(row_sequence == 110)[0]
            gapcount = len(gappositions)
            if gapcount/total_aln_len > 0.15:
                row_sequence_length = len(row_sequence[row_sequence != 110])
                nongap_lengths[row_index] = row_sequence_length
            else:
                nongap_lengths[row_index] = total_aln_len
            
    return nongap_lengths
    


#Spare matricies -- no gpu
def get_pangenome_pairwise_differences(gene_alignment):    
    sparse_matrix, consensus, seq_names = calculate_snp_matrix(gene_alignment)
    gene_diffs = calculate_distance_matrix(sparse_matrix, consensus, "dist", False)
    gene_lens = get_gene_lengths(sparse_matrix, consensus)
    
    return (gene_diffs, gene_lens, seq_names)
