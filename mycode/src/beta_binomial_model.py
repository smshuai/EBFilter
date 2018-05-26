#!/usr/bin/env python
import sys
import numpy as np
import scipy as sp
import pandas as pd
from scipy.optimize import fmin_l_bfgs_b
from scipy.special import betaln
from scipy.stats import combine_pvalues


def get_args():
    ''' Process argmuents
    '''
    if len(sys.argv) == 5:
        tumour_dp_path = sys.argv[1]
        normal_panel_path = sys.argv[2]
        ref_map_path = sys.argv[3]
        out_path = sys.argv[4]
    elif len(sys.argv)==4:
        tumour_dp_path = sys.argv[1]
        normal_panel_path = sys.argv[2]
        ref_map_path = sys.argv[3]
        out_path = ''
    else:
        sys.stderr.write('Usage: beta_binomial_model.py tumour_depth normal_panel ref_map [out_path]')
        sys.exit(1)
    return tumour_dp_path, normal_panel_path, ref_map_path, out_path


def beta_binomial_pmf(params, k, n):
    ''' PMF of beta-binomial distribution
    '''
    llh = beta_binomial_loglikelihood(params, k, n)
    return np.exp(llh)


def beta_binom_pvalue(params, k, n):
    ''' Calculate Prob(X>=k|params, n, k)
    '''
    tempPV = 0
    for kk in range(int(k), int(np.floor(n + 1))):
        currentValue = beta_binomial_pmf(params, kk, n)
        tempPV = tempPV + currentValue
    return tempPV


def beta_binomial_loglikelihood(params, Ks, Ns):
    """ Calculating log-likelihood of beta-binomial distribution

    Args:
        params (List[float]): the parameter of beta distribution ([alpha, beta])
        Ks (numpy.array([int])): the counts for success
        Ns (numpy.array([int])): the counts of trials
    Note: This function is used to replace the one used in EBfilter (more efficient and compatiable with python3)
    """
    alpha = params[0]
    beta = params[1]
    p1 = np.log(sp.special.comb(Ns, Ks))
    p2 = betaln(Ks + alpha, Ns - Ks + beta)
    p3 = betaln(alpha, beta)
    # log-likelihood
    llh = np.sum(p1 + p2 - p3)
    return llh


def regularized_cost(params, Ks, Ns, reg=0.5):
    ''' Cost function for fitting beta binomial distribution with regularization
    '''
    alpha = params[0]
    beta = params[1]
    llh = beta_binomial_loglikelihood(params, Ks, Ns)
    cost = reg * np.log(alpha + beta) - llh
    return cost


def fit_beta_binomial(Ks, Ns, cost=regularized_cost):
    """ Obtaining maximum likelihood estimator of beta-binomial distribution

    Args:
        Ks (numpy.array([int])): the counts for success
        Ns (numpy.array([int])): the counts of trials
    """
    result = fmin_l_bfgs_b(cost, [20, 20], args = (Ks, Ns), approx_grad = True, bounds = [(0.1, 10000000), (1, 10000000)])
    return result[0]


def calc_bb_pval(Ks, Ns, k, n):
    fit = fit_beta_binomial(Ks, Ns)
    pval = beta_binom_pvalue(fit, k, n)
    return pval


def test_base(Ks_p, Ns_p, Ks_n, Ns_n, k_p, n_p, k_n, n_n):
    """ Use beta binomial model to call mutation for the base
    Training data based on normal panels: Ks_p, Ns_p, Ks_n, Ns_n
    Test data from tumour: k_p, n_p, k_n, n_n
    """
    # Convert data to int
    arraytoint = lambda x: x.round().astype(np.int)
    Ks_p = arraytoint(Ks_p)
    Ns_p = arraytoint(Ns_p)
    Ks_n = arraytoint(Ks_n)
    Ns_n = arraytoint(Ns_n)
    k_p = np.int(k_p)
    k_n = np.int(k_n)
    n_p = np.int(n_p)
    n_n = np.int(n_n)
    # Filter Ks and Ns by coverage (>=5 reads) and validity (Ks/Ns < 0.5)
    keep_p = np.logical_and(Ns_p>=5, Ks_p / (Ns_p + 0.1) < 0.5)
    keep_n = np.logical_and(Ns_n>=5, Ks_n / (Ns_n + 0.1) < 0.5)
    if k_p > n_p:
        k_p = n_p
    if k_n > n_n:
        k_n = n_n
    pval_p = calc_bb_pval(Ks_p[keep_p], Ns_p[keep_p], k_p, n_p)
    pval_n = calc_bb_pval(Ks_n[keep_n], Ns_n[keep_n], k_n, n_n)
    if pval_p < 1e-60:
        pval_p = 1e-60  # Cap p-val
    if pval_n < 1e-60:
        pval_n = 1e-60
    pval = combine_pvalues([pval_p, pval_n], 'fisher')[1]
    EB_score = 0
    if pval < 1e-60:
        EB_score = 60
    elif pval > 1.0 - 1e-10:
        EB_score = 0
    else:
        EB_score = - round(np.log10(pval), 3)
    return EB_score, np.sum(keep_p), np.sum(keep_n)


if __name__ == '__main__':
    tumour_dp_path, normal_panel_path, ref_map_path, out_path = get_args()
    tumour_depth = pd.read_table(tumour_dp_path)
    normal_panel = pd.read_table(normal_panel_path)
    ref_map = pd.read_table(ref_map_path)  # reference base at each position
    # Check donor ID
    donor_id = tumour_depth['id'][0]
    if donor_id in normal_panel['id'].values:
        # Get paired normal depth
        normal_depth = normal_panel[normal_panel.id == donor_id]
        normal_panel = normal_panel[normal_panel.id != donor_id]
        # Get output path
        if out_path is '':
            out_path = donor_id + "_U1_GT_res.tsv"
    else:
        sys.stderr.write('Paired normal sample cannot be found for this ID: ' + donor_id)
        sys.exit(1)
    normal_panel = normal_panel.set_index(['chrom', 'pos', 'alt']).sort_index()
    # Add ref to tumour_depth and normal_depth
    tumour_depth = pd.merge(ref_map, tumour_depth, on=('chrom', 'pos'))
    tumour_depth = tumour_depth.loc[tumour_depth.ref != tumour_depth.alt]
    tumour_depth['varDP'] = tumour_depth.varP + tumour_depth.varN
    tumour_depth['totDP'] = tumour_depth.dpP + tumour_depth.dpN
    normal_depth = pd.merge(ref_map, normal_depth, on=('chrom', 'pos'))
    normal_depth = normal_depth.loc[normal_depth.ref != normal_depth.alt]
    normal_depth['varDP'] = normal_depth.varP + normal_depth.varN
    normal_depth['totDP'] = normal_depth.dpP + normal_depth.dpN
    # Test for tumour sample
    eb_scores = []
    num_normal_p = []  # number of normal samples used in model for pos strand
    num_normal_n = []  # number of normal samples used in model for neg strand
    for ix, row in tumour_depth.iterrows():
        if row.varDP > 0:
            k_p, n_p, k_n, n_n = row[['varP', 'dpP', 'varN', 'dpN']]
            Ks_p = normal_panel.loc[(row.chrom, row.pos, row.alt), 'varP'].values
            Ns_p = normal_panel.loc[(row.chrom, row.pos, row.alt), 'dpP'].values
            Ks_n = normal_panel.loc[(row.chrom, row.pos, row.alt), 'varN'].values
            Ns_n = normal_panel.loc[(row.chrom, row.pos, row.alt), 'dpP'].values
            eb, num_p, num_n = test_base(Ks_p, Ns_p, Ks_n, Ns_n, k_p, n_p, k_n, n_n)
            eb_scores.append(eb)
            num_normal_p.append(num_p)
            num_normal_n.append(num_n)
        else:
            # No variant read
            eb_scores.append(0)
            num_normal_p.append(np.nan)
            num_normal_n.append(np.nan)            
    tumour_depth['EB'] = eb_scores
    tumour_depth['num_p'] = num_normal_p
    tumour_depth['num_n'] = num_normal_n
    # Repeat the same analysis for paired normal
    eb_scores = []
    num_normal_p = []  # number of normal samples used in model for pos strand
    num_normal_n = []  # number of normal samples used in model for neg strand
    for ix, row in normal_depth.iterrows():
        if row.varDP > 0:
            k_p, n_p, k_n, n_n = row[['varP', 'dpP', 'varN', 'dpN']]
            Ks_p = normal_panel.loc[(row.chrom, row.pos, row.alt), 'varP'].values
            Ns_p = normal_panel.loc[(row.chrom, row.pos, row.alt), 'dpP'].values
            Ks_n = normal_panel.loc[(row.chrom, row.pos, row.alt), 'varN'].values
            Ns_n = normal_panel.loc[(row.chrom, row.pos, row.alt), 'dpP'].values
            eb, num_p, num_n = test_base(Ks_p, Ns_p, Ks_n, Ns_n, k_p, n_p, k_n, n_n)
            eb_scores.append(eb)
            num_normal_p.append(num_p)
            num_normal_n.append(num_n)
        else:
            # No variant read
            eb_scores.append(0)
            num_normal_p.append(np.nan)
            num_normal_n.append(np.nan)
    normal_depth['EB'] = eb_scores
    normal_depth['num_p'] = num_normal_p
    normal_depth['num_n'] = num_normal_n
    res = pd.merge(tumour_depth, normal_depth, on=('chrom', 'pos', 'ref', 'alt', 'id', 'start', 'end', 'gene', 'strand'), suffixes=('_tumour', '_normal'))
    Mut = np.logical_and(res.EB_tumour-res.EB_normal>3, res.EB_tumour>5)
    Mut = np.logical_and(Mut, res.EB_normal<3)
    res['GT'] = 'WT'
    res.loc[Mut, 'GT'] = 'MUT'
    res.loc[res.totDP_tumour < 12, 'GT'] = 'undet'
    res.to_csv(out_path, sep='\t', index=False)