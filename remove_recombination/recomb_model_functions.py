import math

import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

from tqdm import tqdm
import numpy as np

import scipy.special as sp
from scipy import stats




#This module contains Fthe functions required to identify recombinant pairs.
#genes, and estimate r/m for the collection

def order_pairwise_diffs(pairwise_matrices):
    all_ordered_diffs = []
    print("Ordering pairwise gene differences...")

    for pair in tqdm(pairwise_matrices):
        genes = pair[0]              # 1D array of gene names
        dist = pair[1]               # 1D array
        length = pair[2]             # 1D array

        # Compute proportion using a reusable array to avoid creating temporaries
        proportion = dist / length

        # argsort once
        order = np.argsort(proportion)

        # Reindex *without allocating multiple subarrays*
        # Use np.empty and fill in-place to reduce peak memory
        ordered = np.empty((len(dist), 2), dtype=dist.dtype)
        ordered[:, 0] = dist[order]
        ordered[:, 1] = length[order]

        # Reindex gene names
        ordered_genes = genes[order]

        all_ordered_diffs.append((ordered, ordered_genes))

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
    
    for cutoff in range(1, len(ordered_diffs)):
        pvalue = 1 - stats.binom.cdf(max(ordered_diffs[cutoff-1], 0), 
                                     ordered_lengths[cutoff-1],
                                     average_proportion)
        multtest_alpha = 0.05/len(ordered_diffs)
        if pvalue < multtest_alpha:
            threshold = cutoff
            break
        else:
            continue
    #Handle case where no genes are recombinant
    if threshold == 0:
        threshold = len(ordered_diffs)
    
    recombinant_genes = ordered_genes[threshold:]
    
    return (threshold, recombinant_genes)
    

def pair_has_retained_recombinant(pair, retained_isolates):
    iso1, iso2 = pair.split("-", 1)
    retained = set(retained_isolates)
    return (iso1 in retained) or (iso2 in retained)


def reconcile_cleaned_distances(total_dists,
                                recombinant_gene_pair_dist,
                                gene_recombination_dic,
                                actual_recombinants_to_remove):
    cleaned_dists = total_dists.copy()
    recombinant_dists = {pair: 0 for pair in total_dists}
    retained_pairs_by_gene = {}

    for gene, retained_isolates in actual_recombinants_to_remove.items():
        if not retained_isolates:
            continue

        retained_pairs = []
        for pair in dict.fromkeys(gene_recombination_dic.get(gene, [])):
            if not pair_has_retained_recombinant(pair, retained_isolates):
                continue

            gene_dist = recombinant_gene_pair_dist.get(gene, {}).get(pair)
            if gene_dist is None:
                raise ValueError(
                    f"Missing pairwise SNP distance for gene '{gene}' in pair '{pair}'."
                )

            cleaned_dists[pair] -= gene_dist
            recombinant_dists[pair] += gene_dist
            retained_pairs.append(pair)

        if retained_pairs:
            retained_pairs_by_gene[gene] = retained_pairs

    for pair, total_dist in total_dists.items():
        if cleaned_dists[pair] < 0:
            raise ValueError(f"Cleaned distance for pair '{pair}' became negative.")
        if cleaned_dists[pair] > total_dist:
            raise ValueError(
                f"Cleaned distance for pair '{pair}' exceeds total distance."
            )

    return cleaned_dists, recombinant_dists, retained_pairs_by_gene


def estimate_collection_rm(cleaned_dists, recombinant_dists):
    flat_cleaned = []
    flat_recombinant = []
    pairwise_rms = []
    
    for x in cleaned_dists:
        flat_cleaned.append(cleaned_dists[x])
        flat_recombinant.append(recombinant_dists[x])
        if cleaned_dists[x] > 0:
            pairwise_rms.append(recombinant_dists[x] / cleaned_dists[x])

    total_cleaned = sum(flat_cleaned)
    total_recombinant = sum(flat_recombinant)
    if total_cleaned <= 0:
        raise ValueError("Cannot estimate r/m with zero non-recombinant SNPs.")

    rm = total_recombinant / total_cleaned
    if len(pairwise_rms) > 1:
        stderr = np.std(pairwise_rms, ddof=1) / np.sqrt(len(pairwise_rms))
    else:
        stderr = 0.0
    
    return (rm, stderr, (flat_cleaned, flat_recombinant))


def do_recombination_analysis(pairs, framework, threads):
    if framework == "frequentist":
        analysis_fn = recombination_analysis_frequentist
    elif framework == "bayesian":
        analysis_fn = recombination_analysis_bayesian
    else:
        raise ValueError("Framework must be one of ['frequentist', 'bayesian']")

    with ProcessPoolExecutor(max_workers=threads) as executor:
        results_iter = executor.map(analysis_fn, pairs)
        return list(tqdm(results_iter, total=len(pairs)))

def recombination_analysis_bayesian(pair):
    #wrapper to help keep recombination_removal neat
    model_probabilities, mean_distance = analyse_pair(pair[0][:,0], 
                                                pair[0][:,1])
            
    threshold, recombinants = find_threshold(model_probabilities, pair[1])
            
    
    expec_pg_muts = mean_distance * sum(pair[0][:,1])
    total_dist = sum(pair[0][:,0])
    
    rec_mutations = total_dist - expec_pg_muts
    cleaned_dist = total_dist - rec_mutations
    
    return (recombinants, [total_dist, cleaned_dist, rec_mutations])

def recombination_analysis_frequentist(pair):
    #wrapper to help keep main script tidy
    pair_proportions = pair[0]
    dists = pair_proportions[:,0]
    lens = pair_proportions[:,1]
    
    genes = pair[1]
    
    threshold, recombinants = analyse_pair_frequentist(dists, lens, genes)
    
    try:
        mean_distance = sum(dists[:threshold]) / sum(lens[:threshold])
    except:
        print(threshold)
        print(dists[:threshold])
        print(dists)
        print(lens[:threshold])
        print(lens)
        import sys
        sys.exit()
    
    expec_pg_muts = mean_distance * sum(pair[0][:,1])
    total_dist = sum(pair[0][:,0])
    rec_mutations = total_dist - expec_pg_muts
    cleaned_dist = total_dist - rec_mutations
    
    return(recombinants, [total_dist, cleaned_dist, rec_mutations])
