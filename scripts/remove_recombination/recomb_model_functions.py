import math

import numpy as np

import scipy.special as sp
from scipy import stats

#This module contains all the functions required to identify recombinant pairs.
#genes, and estimate r/m for the collection

def order_pairwise_diffs(pairwise_matrices):
    all_ordered_diffs = {}
    for pair in pairwise_matrices:
        genes = pairwise_matrices[pair][1] 
        pairwise = pairwise_matrices[pair][0]
        proportion = pairwise[:,0] / pairwise[:,1]
        ordered = pairwise[proportion.argsort()]
        ordered_genes = np.array(genes)[proportion.argsort()]
        all_ordered_diffs[pair] = (ordered, ordered_genes)
    return all_ordered_diffs

def calc_log_likelihood(lengths, diffs, hyp_par_1, hyp_par_2):
    log_likelihoods = []
    if len(lengths) != len(diffs):
        raise ValueError("Length and SNP difference vectors different lengths!")
    for index in range(len(lengths)):
        t1 = math.lgamma(hyp_par_1+hyp_par_2) - math.lgamma(hyp_par_1+hyp_par_2+lengths[index])
        t2 = math.lgamma(hyp_par_2+diffs[index])-math.lgamma(hyp_par_2)
        t3 = math.lgamma(hyp_par_1 + lengths[index] - diffs[index]) - math.lgamma(hyp_par_1)

        log_ml = t1 + t2 + t3
        log_likelihoods.append(log_ml)
    return(log_likelihoods)

def alt_log_likelihood(a0, a1, n0, n1):
    logml = sp.loggamma(a0+a1) - sp.loggamma(a0+a1+n0+n1) + sp.loggamma(a0+n0) + sp.loggamma(a1+n1) - sp.loggamma(a0) - sp.loggamma(a1)
    return logml

def analyse_pair(ordered_diffs, ordered_lengths, a0=9, a1=1):
    #if len(ordered_diffs) != len(ordered_lengths):
    #    raise ValueError("pairwse lengths and differences not the same!")
    #no_seqs = len(ordered_diffs)
    
    #logml_at_each_threshold = calc_log_likelihood(ordered_lengths, ordered_diffs, a0, a1)
    logml_at_each_threshold = alt_log_likelihood(a1,a0, ordered_lengths-ordered_diffs, ordered_diffs)
    
    cumm_lengths = np.cumsum(ordered_lengths)
    cumm_diffs = np.cumsum(ordered_diffs)
    #joint_threshold_logml = calc_log_likelihood(cumm_lengths, cumm_diffs, a0, a1)
    joint_threshold_logml = alt_log_likelihood(a1,a0,cumm_lengths-cumm_diffs,cumm_diffs)
    
    summed_individual_threshold = np.cumsum(logml_at_each_threshold[::-1])
    summed_individual_threshold = np.append(summed_individual_threshold[::-1][1:], 0)
    
    threshold_logmls = joint_threshold_logml + summed_individual_threshold
    normalised_logmls = threshold_logmls - max(threshold_logmls)
    
    threshold_model_probs = np.exp(normalised_logmls)/sum(np.exp(normalised_logmls))
    
    model_distances = (a1 + cumm_diffs) / (cumm_lengths + a0 + a1)
    mean_model_distance = sum(model_distances * threshold_model_probs)
    
    
    return (threshold_model_probs, mean_model_distance)

def find_threshold(model_probs, ordered_genes):
    no_of_genes = len(ordered_genes)
    model_dists = model_probs * np.arange(1, no_of_genes+1)
    
    threshold = np.sum(model_dists)
    rounded_threshold = round(threshold)
    recombinant_genes = ordered_genes[rounded_threshold:]
    return (rounded_threshold, recombinant_genes)


def analyse_pair_frequentist(ordered_diffs, ordered_lengths, ordered_genes):
    
    average_proportion = sum(ordered_diffs) / sum(ordered_lengths)
    threshold = 0
    
    for cutoff in range(len(ordered_diffs), 0, -1):
        pvalue = 1 - stats.binom.cdf(max(ordered_diffs[cutoff-1], 0), 
                                     ordered_lengths[cutoff-1],
                                     average_proportion)
        multtest_alpha = 0.05/len(ordered_diffs)
        if pvalue < multtest_alpha:
            threshold = cutoff
            break
        else:
            continue
    recombinant_genes = ordered_genes[threshold:]
    
    return (threshold, recombinant_genes)
    

def estimate_collection_rm(cleaned_dists, uncleaned_dists):
    raw_proportions = np.array(sorted(list(proportion_dictionary.values())))
    
    tenth_percentile = np.percentile(raw_proportions, 10)
    nintieth_percentile = np.percentile(raw_proportions, 90)
    
    middle_80 = [(tenth_percentile < raw_proportions) & 
                 (raw_proportions < nintieth_percentile)][0]
    
    xaxis = range(len(raw_proportions[middle_80]))
    regression = stats.linregress(xaxis,raw_proportions[middle_80])
    
    return(regression, raw_proportions[middle_80], xaxis)
    
        