import pandas as pd
import scipy as sp
from pandas_plink import read_plink
from SAME import SAME
import argparse
import datetime
import random
import numpy as np
import math
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from scipy.stats import norm
from scipy import linalg
from scipy.cluster.hierarchy import fcluster

parser = argparse.ArgumentParser()
parser.add_argument('--plink',dest='plink',help='plink file location')
parser.add_argument('--expression',dest='expression',help='gene expression file location')
parser.add_argument('--gamma',dest='gamma',help='gamma initialization file location')
parser.add_argument('--frq',dest='frq',help='frequency file location')
parser.add_argument('--annotation',dest='annotation',help='annotation file location')
parser.add_argument('--outSNP',dest='outputSNP',help='output of SNP file location')
parser.add_argument('--outAnno',dest='outputAnno',help='output of annotation file location')
parser.add_argument('--outPerfer',dest='outputPerfer',help='output of perfermance file location')
parser.add_argument('--outResult',dest='outputResult',help='output of Result file location')

args = parser.parse_args()
input_plink = args.plink
input_frq = args.frq
input_annotation = args.annotation
input_expression = args.expression
input_gamma = args.gamma
output_SNP = args.outputSNP
output_annotation = args.outputAnno
output_perfer = args.outputPerfer
output_result = args.outputResult

(bim, fam, bed) = read_plink(input_plink)
pdTmp = pd.DataFrame(bed.compute())
frq =  pd.read_table(input_frq)
ano =  pd.read_table(input_annotation)
MAF = []
rare_index = []
common_index = []
for i in range(frq.shape[0]):
    maf = float(frq.iloc[i,0].split()[4])
    if (maf>0.005) & (maf<0.05):
        rare_index.append(i)
    elif (maf >= 0.05):
        common_index.append(i)
    MAF.append(maf)
SNP_index = sp.append(rare_index, common_index)
SNP = bim.iloc[SNP_index,:]

rare_geno = pdTmp.iloc[rare_index,:].T
common_geno = pdTmp.iloc[common_index,:].T

annotation_tmp = ano.values
annotation_tmp = sp.delete(annotation_tmp, [0,1,2], axis = 1)
annotation_tmp = annotation_tmp[rare_index]
annotation = sp.ones(shape = (annotation_tmp.shape[0],annotation_tmp.shape[1]+1))
annotation[:,:-1] = annotation_tmp

##prepare the data
common_geno = common_geno.fillna(common_geno.mean()).values
common_geno = (common_geno-common_geno.mean())/common_geno.std()
rare_geno = rare_geno.fillna(rare_geno.mean()).values
rare_geno = (rare_geno-rare_geno.mean())/rare_geno.std()
G = pd.read_table(input_expression)
G = G.values.T[0]
gamma = pd.read_table(input_gamma)
gamma = gamma.values.T[0]
b = sp.zeros(shape = (1, annotation.shape[1]))
b[0,-1] = -5
alpha_e = 0.01
tau_e = 0.01
alpha_r = 0.01
tau_r = 0.01
alpha_c = 0.01
tau_c = 0.01
alpha_b = 0.01
tau_b = 0.01
Ite = 50
rate = 2


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def flatten(lst, inds):
    result = []
    for i in inds:
        result.extend(lst[i])
    return result

def SAME_CV(cutoff, hierarchy, G, rare_geno, common_geno, annotation):
    #p_cutoff = sorted(p_value)[cutoff]
    gamma = sp.zeros(shape = rare_geno.shape[1])
    #gamma[p_value<p_cutoff] = 1
    labels = fcluster(hierarchy, cutoff, criterion='maxclust')
    for i in range(1,max(labels)+1):
        index = sp.where(labels == i)[0]
        if index.shape==0:
            continue
        else:
            gamma[index[sp.random.randint(0,index.shape)]] = 1
    
    b = sp.zeros(shape = (1, annotation.shape[1]))
    b[0,-1] = -5
    alpha_e = 0.01
    tau_e = 0.01 
    alpha_r = 0.01
    tau_r = 0.01 
    alpha_c = 0.01 
    tau_c = 0.01 
    alpha_b = 0.01 
    tau_b = 0.01
    Ite = 50
    rate = 2
    #geno = sp.hstack((rare_geno, common_geno))
    for i in range(5):
        inds=[*range(0,5)]
        inds.remove(i)
        flat = flatten(a,inds)
        G_tmp=G[flat,]
        rare_geno_tmp = rare_geno[flat,]
        common_geno_tmp = common_geno[flat,]
        result_beta, result_b, result_gamma, result_ite, result_residual, result_Sigma_e, result_Sigma_r, result_Sigma_c, result_Sigma_b, result_Gamma= SAME(G_tmp, rare_geno_tmp, common_geno_tmp, annotation, gamma, b, alpha_e, tau_e, alpha_r, tau_r, alpha_c, tau_c, alpha_b, tau_b, Ite, rate)
        rare_geno_test = rare_geno[a[i],]
        common_geno_test = common_geno[a[i],]
        geno_test = sp.hstack((rare_geno_test, common_geno_test))
        out_beta = result_beta.mean(axis=0)
        out_gamma = sp.append(result_Gamma[-1], sp.ones(shape = common_geno_tmp.shape[1]))
        G_pre = geno_test.dot(out_beta*out_gamma)
        if i == 0:
            G_test_true = G[a[i],]
            G_test_pred = G_pre
        else:
            G_test_true = np.concatenate((G_test_true,G[a[i],]))
            G_test_pred = np.concatenate((G_test_pred,G_pre))
    
    cor_pre_cv = sp.stats.pearsonr(G_test_true,G_test_pred)[0]
    SSR = sum((G_test_true-G_test_pred)**2)
    SST = sum((G_test_true-sp.mean(G_test_true))**2)
    R2_pre_cv-SSR/SST
    ###calculate the overall perfermance
    result_beta, result_b, result_gamma, result_ite, result_residual, result_Sigma_e, result_Sigma_r, result_Sigma_c, result_Sigma_b, result_Gamma= SAME(G, rare_geno, common_geno, annotation, gamma, b, alpha_e, tau_e, alpha_r, tau_r, alpha_c, tau_c, alpha_b, tau_b, Ite, rate)
    out_beta = result_beta.mean(axis=0)
    out_gamma = sp.append(result_Gamma[-1], sp.ones(shape = common_geno.shape[1]))
    geno = sp.hstack((rare_geno, common_geno))
    G_pre_all = geno.dot(out_beta*out_gamma)
    cor_pre_all = sp.stats.pearsonr(G, G_pre_all)[0]
    SSR = sum((G-G_pre_all)**2)
    SST = sum((G-sp.mean(G))**2)
    R2_pre_all = 1-SSR/SST
    return cor_pre_cv,R2_pre_cv,cor_pre_all,R2_pre_all,sp.mean(result_beta,axis=0), result_gamma, sp.mean(result_b,axis=0), result_Sigma_b,result_Sigma_c,result_Sigma_e,result_Sigma_r


inds = [*range(1,G.shape[0])]
random.shuffle(inds)
a = list(chunks(inds, int(len(inds)/5)))

corr = sp.corrcoef(rare_geno,rowvar=0)
corr = (corr+corr.T)/2
sp.fill_diagonal(corr,1)
dissimilarity = 1-abs(corr)
hierarchy = linkage(squareform(dissimilarity), method='average')

Cutoff = sp.trunc(sp.exp(sp.linspace(-math.log(rare_geno.shape[1],5),0,18))[0:10]*rare_geno.shape[1]).astype(int)
Per_names = ["cv_cor","cv_R2","overall_cor","overall_R2","init_gamma","result_gamma"]
Perfermance = pd.DataFrame(columns = Per_names)

for j in range(len(Cutoff)):
    starttime1 = datetime.datetime.now()
    cor_pre_cv,R2_pre_cv,cor_pre_all,R2_pre_all,result_beta, result_gamma, result_b,result_Sigma_b,result_Sigma_c,result_Sigma_e,result_Sigma_r = SAME_CV(Cutoff[j], hierarchy, G, rare_geno, common_geno, annotation)
    endtime1 = datetime.datetime.now()
    print (endtime1 - starttime1)
    Perfermance.loc[j] = sp.array([cor_pre_cv,R2_pre_cv,cor_pre_all,R2_pre_all,Cutoff[j],sum(result_gamma)])

Perfermance.to_csv(output_perfer, sep = '\t', index = False)
cv_R2 = Perfermance['cv_R2'].values
flag = sp.zeros(len(Cutoff))
Result = pd.DataFrame(columns = Per_names)
for i  in range(1,len(Cutoff)-1):
    if (cv_R2[i]-cv_R2[i-1]>0)&(cv_R2[i]-cv_R2[i+1]>0):
        flag[i] = 1
if sum(flag)==0:
    Result.loc[0] = Perfermance.iloc[sp.argmax(cv_R2)].values
else:
    tmp = cv_R2*flag
    Result.loc[0] = Perfermance.iloc[sp.argmax(tmp)].values
Result.to_csv(output_result, sep = '\t', index = False)
#geno = sp.hstack((rare_geno, common_geno))
#out_beta = result_beta.mean(axis=0)
#out_b = result_b.mean(axis=0)
#out_gamma = sp.append(result_Gamma[-1], sp.ones(shape = common_geno.shape[1]))
#SNP['beta'] = out_beta
#SNP['gamma'] = out_gamma
#SNP.to_csv(output_SNP, sep = '\t', index = False)

#Ano_names = list(ano)
#del Ano_names[0:3]
#Ano_names.append('intercept')
#Ano_coef = pd.DataFrame(columns = Ano_names)
#Ano_coef.loc[0] = out_b
#Ano_coef.to_csv(output_annotation, sep = '\t', index = False)

#geno = sp.hstack((rare_geno, common_geno))
#G_pre = geno.dot(out_beta*out_gamma)
#h2_pre = sp.var(G_pre)/sp.var(G)
#true_h2 = sp.var(geno.dot(beta))/sp.var(G)
#cor_pre = sp.stats.pearsonr(G, G_pre)[0]
#true_cor = sp.stats.pearsonr(G, geno.dot(beta))[0]
#SSR = sum((G-G_pre)**2)
#SST = sum((G-sp.mean(G))**2)
##SSR1 = sum((G-geno.dot(beta))**2)
##R2_pre = 1-SSR/SST
#true_R = 1-SSR1/SST
#gamma_sum = sum(result_Gamma[-1])
#gamma_true = sum(generate_gamma)
#Per_names = ["gamma","h2","corr","R2"]
#Perfermance = pd.DataFrame(columns = Per_names)
#Perfermance.loc[0] = sp.array([gamma_sum,h2_pre,cor_pre,R2_pre])
##Perfermance.loc[1] = sp.array([gamma_true,true_h2,true_cor,true_R])
