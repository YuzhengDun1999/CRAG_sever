#!/usr/bin/R
# @ annotation: annotation matrix  with the intercept in the last column
# @ h2_anno: the variance of a SNP being effective or not explained by annotation
# @ h2_x: heritability of y/gene expression
# @ pi: the percentage of SNPs being effective/causal
# @ common: the genotype matrix of common SNPs
# @ rare: the genotype matrix of rare SNPs
# @annotation: annotation for rare SNPs
# @h2_anno: the variance of rare SNPs getting selected explained by annnotation

import math
import scipy as sp
def generation_CommonRare(common, rare, annotation, h2_anno, h2_r, h2_c, pi):
    # first deal with NA in the genoytype matrix
    common = common.fillna(common.mean()).values
    rare = rare.fillna(rare.mean()).values
    rare = (rare-rare.mean())/rare.std()
    common = (common-common.mean())/common.std()
    geno = sp.hstack((rare, common))
    
    G = []
    R = rare.shape[1]
    N = common.shape[0]
    C = common.shape[1]
    ###real parameters
    # number of SNPs being causal
    p_causal = math.floor(pi * R)
    sigma_r = h2_r/p_causal
    sigma_b = h2_anno/(annotation.shape[1]-1)
    
    b = sp.random.normal(0, math.sqrt(sigma_b), annotation.shape[1]-1)
    # first determine random noise in the annotation layer
    anno_error = ((1 - h2_anno)/h2_anno)*annotation[:, 0:-1].dot(b).var()
    b = b/math.sqrt(anno_error)
    sigma_b = h2_anno/(annotation.shape[1]-1)/anno_error
    
    # determine the intercept based on how many values larger than 0
    anno_noise = sp.random.normal(0, 1, R)
    anno_val = annotation[:, 0:-1].dot(b) + anno_noise
    anno_intercept = -sorted(anno_val,reverse=True)[p_causal]
    b = sp.append(b, anno_intercept)
    w = annotation.dot(b) + anno_noise
    
    gamma = sp.zeros(shape = R)
    beta_rare = sp.zeros(shape = R)
    gamma[w>0] = 1
    beta_rare[w>0] = sp.random.normal(0, math.sqrt(sigma_r), sum(w>0))
    
    sigma_c = h2_c/C
    beta_common = sp.random.normal(0, math.sqrt(sigma_c), C)
    
    sigma_e = (1-h2_c-h2_r)/(h2_c+h2_r) * (rare.dot(beta_rare)+common.dot(beta_common)).var()
    noise = sp.random.normal(0, math.sqrt(sigma_e), N)
    noise <- noise/noise.std()*math.sqrt(sigma_e)
    
    beta = sp.append(beta_rare, beta_common)
    G = geno.dot(beta) + noise
    
    return G, common, rare, gamma, beta, sigma_r, sigma_c, sigma_e, sigma_b, b, w
