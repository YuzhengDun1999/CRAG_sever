import pandas as pd
import scipy as sp
from pandas_plink import read_plink
from generation import generation_CommonRare
import statsmodels.api as sm
from SAME import Initial_SAME
import argparse
import datetime
import math
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

parser = argparse.ArgumentParser()
parser.add_argument('--plink',dest='plink',help='plink file location')
parser.add_argument('--frq',dest='frq',help='frequency file location')
parser.add_argument('--annotation',dest='annotation',help='annotation file location')
parser.add_argument('--outSNP',dest='outputSNP',help='output of SNP file location')
parser.add_argument('--outAnno',dest='outputAnno',help='output of annotation file location')

args = parser.parse_args()
input_plink = args.plink
input_frq = args.frq
input_annotation = args.annotation
output_SNP = args.outputSNP
output_annotation = args.outputAnno

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

#all_geno=pdTmp.values
#rare_geno = all_geno[rare_index].T
#common_geno = all_geno[common_index].T
rare_geno = pdTmp.iloc[rare_index,:].T
common_geno = pdTmp.iloc[common_index,:].T

annotation_tmp = ano.values
annotation_tmp = sp.delete(annotation_tmp, [0,1,2], axis = 1)
annotation_tmp = annotation_tmp[rare_index]
annotation = sp.ones(shape = (annotation_tmp.shape[0],annotation_tmp.shape[1]+1))
annotation[:,:-1] = annotation_tmp

###test
#rare_geno = rare_geno.T.iloc[0:2000,:].T
#common_geno = common_geno.T.iloc[0:2000,:].T
#annotation = annotation[0:2000,:]

G, common_geno, rare_geno, generate_gamma, beta, sigma_r, sigma_c, sigma_e, sigma_b, b = generation_CommonRare(common = common_geno,rare = rare_geno,annotation = annotation, h2_anno = 0.5,h2_r = 0.2,h2_c = 0.1,pi = 0.03)

#p_value = sp.zeros(shape = rare_geno.shape[1])
#for i in range(rare_geno.shape[1]):
#    mod = sm.OLS(G, rare_geno[:,i])
#    res = mod.fit()
#    p_value[i] = res.pvalues
###cluster
corr = sp.corrcoef(rare_geno,rowvar=0)
corr = (corr+corr.T)/2
sp.fill_diagonal(corr,1)
dissimilarity = 1-abs(corr)
hierarchy = linkage(squareform(dissimilarity), method='average')

### initialization
Log_Likelihood = []
Log_Likelihood1, Log_Likelihood2,Log_Likelihood3,Log_Likelihood4,Log_Likelihood5,Log_Likelihood6,= [],[],[],[],[],[]
Result_Beta = []
Result_B = []
Result_Gamma = []
Cutoff = sp.trunc(sp.exp(sp.linspace(-math.log(rare_geno.shape[1],5),0,18))[0:10]*rare_geno.shape[1]).astype(int)
Result_Sigma_b,Result_Sigma_c,Result_Sigma_e,Result_Sigma_r=[],[],[],[]

for i in range(len(Cutoff)):
    starttime1 = datetime.datetime.now()
    result_beta, result_gamma, result_b, log_likelihood, log_likelihood1, log_likelihood2, log_likelihood3, log_likelihood4, log_likelihood5, log_likelihood6,result_Sigma_b,result_Sigma_c,result_Sigma_e,result_Sigma_r = Initial_SAME(Cutoff[i], hierarchy, G, rare_geno, common_geno, annotation)
    endtime1 = datetime.datetime.now()
    print (endtime1 - starttime1)
    Log_Likelihood.append(log_likelihood)
    Log_Likelihood1.append(log_likelihood1)
    Log_Likelihood2.append(log_likelihood2)
    Log_Likelihood3.append(log_likelihood3)
    Log_Likelihood4.append(log_likelihood4)
    Log_Likelihood5.append(log_likelihood5)
    Log_Likelihood6.append(log_likelihood6)
    Result_Beta.append(result_beta)
    Result_B.append(result_b)
    Result_Gamma.append(result_gamma)
    Result_Sigma_b.append(result_Sigma_b)
    Result_Sigma_c.append(result_Sigma_c)
    Result_Sigma_e.append(result_Sigma_e)
    Result_Sigma_r.append(result_Sigma_r)

out_beta = Result_Beta[sp.argmax(Log_Likelihood)]
out_b = Result_B[sp.argmax(Log_Likelihood)]
out_gamma = sp.append(Result_Gamma[sp.argmax(Log_Likelihood)], sp.ones(shape = common_geno.shape[1]))
SNP['beta'] = out_beta
SNP['gamma'] = out_gamma
SNP.to_csv(output_SNP, sep = '\t', index = False)
Ano_names = list(ano)
del Ano_names[0:3]
Ano_names.append('intercept')
Ano_coef = pd.DataFrame(columns = Ano_names)
Ano_coef.loc[0] = out_b
Ano_coef.to_csv(output_annotation, sep = '\t', index = False)

#evaluationh2_true = sp.var(geno.dot(beta))/sp.var(G)
geno = sp.hstack((rare_geno, common_geno))
H2_Pre_1, H2_Pre_2, G_Pre_1, G_Pre_2, R2_Pre_1, R2_Pre_2, Cor_Pre_1, Cor_Pre_2= [],[],[],[],[],[],[],[]
sum_gamma = sp.zeros(len(Cutoff))
for i in range(len(Cutoff)):
    result_beta = Result_Beta[i]
    result_gamma = sp.append(Result_Gamma[i], sp.ones(shape = common_geno.shape[1]))
    h2_pre_1 = sp.var(geno.dot(result_beta*result_gamma))/sp.var(G)
    h2_pre_2 = sp.var(geno.dot(result_beta))/sp.var(G)
    G_pre_1 = geno.dot(result_beta*result_gamma)
    G_pre_2 = geno.dot(result_beta)
    SSR = sum((G-G_pre_1)**2)
    SST = sum((G-sp.mean(G))**2)
    R2_1 = 1-SSR/SST
    SSR = sum((G-G_pre_2)**2)
    R2_2 = 1-SSR/SST
    cor_pre_1 = sp.stats.pearsonr(G, G_pre_1)
    cor_pre_2 = sp.stats.pearsonr(G, G_pre_2)
    H2_Pre_1.append(h2_pre_1)
    H2_Pre_2.append(h2_pre_2)
    G_Pre_1.append(G_pre_1)
    G_Pre_2.append(G_pre_2)
    R2_Pre_1.append(R2_1)
    R2_Pre_2.append(R2_2)
    Cor_Pre_1.append(cor_pre_1)
    Cor_Pre_2.append(cor_pre_2)
    sum_gamma[i] = sum(Result_Gamma[i])


SSR1 = sum((G-geno.dot(beta))**2)
true_R = 1-SSR1/SST
true_h2 = sp.var(geno.dot(beta))/sp.var(G)
true_cor = sp.stats.pearsonr(G, geno.dot(beta))


