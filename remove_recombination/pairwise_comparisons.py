import numpy as np

import cupy as cp

def check_for_big_indel(byteseq1, byteseq2):
    gaps1 = np.sum(byteseq1 == 45) #45 == ord("-") 
    gaps2 = np.sum(byteseq2 == 45)
    diff = gaps2 - gaps1
    if (diff/gaps1.size) > 0.15:
        return diff
    else: 
        return None

def check_for_big_indel_gpu(byteseq1, byteseq2):
    gaps1 = cp.sum(byteseq1 == 45) #45 == ord("-") 
    gaps2 = cp.sum(byteseq2 == 45)
    diff = gaps2 - gaps1
    if (diff/gaps1.size) > 0.15:
        return diff
    else: 
        return None

def get_pairwise_differences(seq1, seq2, gpu):
    if seq1.size != seq2.size:
        raise ValueError("Two aligned sequences are of different lengths!")
    
    if gpu:
        if check_for_big_indel_gpu == None:
            diffs = cp.count_nonzero(seq1^seq2)
            length = seq1.size
            result = (np.array([diffs, length]))
            return result
        else:
            mask = seq1 != 45
            consecutive_bases = cp.diff(np.concatenate(([0], mask.astype(int), [0])))
            start_positions = cp.where(consecutive_bases == 1)[0]
            end_positions = cp.where(consecutive_bases == -1)[0]
            lengths = end_positions - start_positions
            
            longest_seq_index = cp.argmax(lengths)
            longest_seq_start = start_positions[longest_seq_index]
            longest_seq_end = end_positions[longest_seq_index]
            
            cropped_seq1 = seq1[longest_seq_start:longest_seq_end]
            cropped_seq2 = seq2[longest_seq_start:longest_seq_end]
            
            diffs = cp.count_nonzero(cropped_seq1^cropped_seq2)
            length = cropped_seq1.size
            result = (np.array([diffs, length]))
            return result
            
    else:    
        if check_for_big_indel == None:
            diffs = np.count_nonzero(seq1^seq2)
            length = seq1.size
            result = (np.array([diffs, length]))
            return result
        else:
            mask = seq1 != 45
            consecutive_bases = np.diff(np.concatenate(([0], mask.astype(int), [0])))
            start_positions = np.where(consecutive_bases == 1)[0]
            end_positions = np.where(consecutive_bases == -1)[0]
            lengths = end_positions - start_positions
            
            longest_seq_index = np.argmax(lengths)
            longest_seq_start = start_positions[longest_seq_index]
            longest_seq_end = end_positions[longest_seq_index]
            
            cropped_seq1 = seq1[longest_seq_start:longest_seq_end]
            cropped_seq2 = seq2[longest_seq_start:longest_seq_end]
            
            diffs = np.count_nonzero(cropped_seq1^cropped_seq2)
            length = cropped_seq1.size
            result = (np.array([diffs, length]))
            return result
        
        

def get_pairs(isolate_list):
    pairs_list = [(isolate_list[i], isolate_list[j]) for i in range(len(isolate_list)) for j in range(i+1,len(isolate_list))]
    return pairs_list

def get_pangenome_pairwise_differences(gene_alignments, isolates_to_consider):    
    #Legacy code  -- extremely slow
    # diffs = []
    # names = []
    # for sequences in sequence_files:
    #     seq1 = None
    #     seq2 = None
    #     for sequence in sequences[0]:
    #         if sequences_to_consider[0] in sequence.id:
    #             seq1=str(sequence.seq)
    #         elif sequences_to_consider[1] in sequence.id:
    #             seq2 = str(sequence.seq)
    #         else:
    #             continue
    #     if (type(seq1) == str) and (type(seq2) == str):
    #         diffs.append(get_pairwise_differences(seq1, seq2))
    #         names.append(sequences[1].split(".")[0])
    # diffs = np.array(diffs)  
    diffs = []
    names = []
    for gene in gene_alignments:
        seq1 = gene[0].get(isolates_to_consider[0], None)
        seq2 = gene[0].get(isolates_to_consider[1], None)
        if (seq1 is None) or (seq2 is None):
            continue
        diffs.append(get_pairwise_differences(seq1, seq2))
        names.append(gene[1].split(".")[0])
    diffs = np.array(diffs)
    return (diffs, names)
