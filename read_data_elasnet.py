import pandas as pd
import scipy as sp
from pandas_plink import read_plink
from generation import generation_CommonRare
from SAME import SAME
import argparse
import datetime
import glmnet_python
from cvglmnet import cvglmnet
from cvglmnetCoef import cvglmnetCoef


parser = argparse.ArgumentParser()
parser.add_argument('--plink',dest='plink',help='plink file location')
parser.add_argument('--frq',dest='frq',help='frequency file location')
parser.add_argument('--annotation',dest='annotation',help='annotation file location')
parser.add_argument('--outSNP',dest='outputSNP',help='output of SNP file location')
parser.add_argument('--outAnno',dest='outputAnno',help='output of annotation file location')
parser.add_argument('--outPerfer',dest='outputPerfer',help='output of perfermance file location')


args = parser.parse_args()
input_plink = args.plink
input_frq = args.frq
input_annotation = args.annotation
output_SNP = args.outputSNP
output_annotation = args.outputAnno
output_perfer = args.outputPerfer

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

##generation
G, common_geno, rare_geno, generate_gamma, beta, sigma_r, sigma_c, sigma_e, sigma_b, b, generate_z = generation_CommonRare(common = common_geno,rare = rare_geno,annotation = annotation, h2_anno = 0.5,h2_r = 0.2,h2_c = 0.1,pi = 0.04)

### initialization by using glmnet
cvfit = cvglmnet(x = rare_geno, y = G, ptype = 'deviance', alpha = 0.5, intr = 0)
cv_coef = cvglmnetCoef(cvfit, s = 'lambda_min')#n*1
gamma = sp.ones(rare_geno.shape[1])
gamma[cv_coef.T[0] == 0]= 0


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

starttime1 = datetime.datetime.now()
result_beta, result_b, result_gamma, result_ite, result_residual, result_Sigma_e, result_Sigma_r, result_Sigma_c, result_Sigma_b, result_Gamma= SAME(G, rare_geno, common_geno, annotation, gamma, b, alpha_e, tau_e, alpha_r, tau_r, alpha_c, tau_c, alpha_b, tau_b, Ite, rate)
endtime1 = datetime.datetime.now()
print (endtime1 - starttime1)

geno = sp.hstack((rare_geno, common_geno))
out_beta = result_beta
out_b = result_b
out_gamma = sp.append(result_gamma, sp.ones(shape = common_geno.shape[1]))
SNP['beta'] = out_beta
SNP['gamma'] = out_gamma
SNP.to_csv(output_SNP, sep = '\t', index = False)

Ano_names = list(ano)
del Ano_names[0:3]
Ano_names.append('intercept')
Ano_coef = pd.DataFrame(columns = Ano_names)
Ano_coef.loc[0] = out_b
Ano_coef.to_csv(output_annotation, sep = '\t', index = False)

geno = sp.hstack((rare_geno, common_geno))
G_pre = geno.dot(result_beta.mean(axis=0)*result_Gamma[-1])
h2_pre = sp.var(G_pre)/sp.var(G)
true_h2 = sp.var(geno.dot(beta))/sp.var(G)
cor_pre = sp.stats.pearsonr(G, G_pre)[0]
true_cor = sp.stats.pearsonr(G, geno.dot(beta))[0]
SSR = sum((G-G_pre)**2)
SST = sum((G-sp.mean(G))**2)
SSR1 = sum((G-geno.dot(beta))**2)
R2_pre = 1-SSR/SST
true_R = 1-SSR1/SST
gamma_sum = sum(result_Gamma[-1])
gamma_true = sum(generate_gamma)
Per_names = ["gamma","h2","corr","R2"]
Perfermance = pd.DataFrame(columns = Per_names)
Perfermance.loc[0] = sp.array([gamma_sum,h2_pre,cor_pre,R2_pre])
Perfermance.loc[1] = sp.array([gamma_true,true_h2,true_cor,true_R])
Perfermance.to_csv(output_perfer, sep = '\t', index = False)



