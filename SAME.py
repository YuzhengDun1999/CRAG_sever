import math
import scipy as sp
from scipy import linalg
from scipy.stats import norm
from scipy.cluster.hierarchy import fcluster


def sample_b(sigma_b, sum_anno, annot, z, number):
    Az = 0
    for i in range(annot.shape[0]):
        Az = Az + annot[i] * z[i]
    
    for j in range(annot.shape[1]-1):
        sum_anno[j,j] = sum_anno[j,j] + 1/sigma_b
    
    Var = linalg.inv(sum_anno)
    E = Var.dot(Az)
    return sp.random.multivariate_normal(E,Var,number)


def sample_beta_chol(G, geno, gamma, sigma_r, sigma_c, sigma_e, number):
    N = len(G)#total number of samples
    R = len(gamma)
    C = len(geno[1])-R
    zero = 10**(-10)
    gamma_temp = sp.append(gamma, sp.ones(C))
    
    ###Sherman–Morrison–Woodbury formula
    temp = sigma_r * sp.ones(int(sum(gamma)))
    temp = sp.append(temp, sp.array([sigma_c]*C))
    Sigma_x = sp.diag(temp) ##sigma_1
    X = geno[:, gamma_temp == 1]/math.sqrt(sigma_e)
    temp1 = Sigma_x.dot(X.T)
    
    Var = Sigma_x-temp1.dot(linalg.inv(X.dot(Sigma_x).dot(X.T)+sp.identity(N))).dot(temp1.T)
    E = Var.dot(geno[:, gamma_temp == 1].T).dot(G)/sigma_e
    CH = linalg.cholesky(Var)
    E_norm = E.T.dot(linalg.inv(CH))#expectation of dependent multivariate normal
    
    new_beta = sp.zeros(shape = (number, R+C))
    for i in range(number):
        new_beta[i, gamma_temp==1] = (E_norm + sp.random.normal(0, 1, int(sum(gamma))+C)).dot(CH)
        new_beta[i, gamma_temp==0] = sp.random.normal(0, math.sqrt(zero), R-int(sum(gamma)))
    return new_beta


def sample_beta_chol_old(G, geno, gamma, sigma_r, sigma_c, sigma_e, number):
    N = len(G)#total number of samples
    R = len(gamma)
    C = len(geno[1])-R
    zero = 10**(-10)
    
    ###Sherman–Morrison–Woodbury formula
    temp = sigma_r**gamma * zero**(1-gamma)
    temp = sp.append(temp, sp.array([sigma_c]*C))
    Sigma_x = sp.diag(temp) ##sigma_1
    X = geno/math.sqrt(sigma_e)
    temp1 = Sigma_x.dot(X.T)
    
    Var = Sigma_x-temp1.dot(linalg.inv(X.dot(Sigma_x).dot(X.T)+sp.identity(N))).dot(temp1.T)
    E = Var.dot(geno.T).dot(G)/sigma_e
    CH = linalg.cholesky(Var)
    E_norm = E.T.dot(linalg.inv(CH))#expectation of dependent multivariate normal
    
    new_beta = sp.zeros(shape = (number, R+C))
    for i in range(number):
        new_beta[i] = (E_norm + sp.random.normal(0, 1, R+C)).dot(CH)
    return new_beta


def sample_gamma(G, rare, b, annotation, beta, gamma, sigma_r, number):
    zero = 10**(-10)
    R = rare.shape[1]
    Probit = norm.pdf(annotation.dot(b.T)).T
    ratio = norm.pdf(beta[:,0:R], 0, math.sqrt(sigma_r)) / (zero + norm.pdf(beta[:,0:R], 0, math.sqrt(zero)))
    odds = (Probit+zero)/(1-Probit+zero)*ratio
    post = odds/(1+odds)
    post[odds == float('inf')] = 1
    post_sum = sum(post)/number
    for i in range(R):
        gamma[i] = sp.random.binomial(1,post_sum[i],1)
    return gamma, odds, ratio, post, Probit

def sample_z(annotation, b, gamma, number):
    R = len(gamma)
    z = sp.zeros(shape = (R))
    sum_b = b.sum(axis = 0)
    Mean = 1/number * annotation.dot(sum_b.T)
    std = math.sqrt(1/number)
    for i in range(R):
        mean = Mean[i]
        if gamma[i] == 1:
            lower, upper = 0, 10**10
        else:
            lower, upper =  -10**10, 0
        z[i] = sp.stats.truncnorm.rvs((lower-mean)/std, (upper-mean)/std, loc = mean, scale = std,size=1)
    return z


def sample_e(G, geno, beta, alpha_e, tau_e, number):
    N = len(G)
    shape = number * alpha_e + N*number/2 + number -1
    scale = number * tau_e
    for m in range(number):
        scale = scale + sum((G - geno.dot(beta[m]))**2)/2
    return 1/sp.random.gamma(shape, 1/scale, 1)


def sample_r(gamma, beta, alpha_r, tau_r, number):
    R = len(gamma)
    shape = number * alpha_r + sum(gamma) * number / 2 + number -1
    scale = number * tau_r
    for i in range(R):
        for m in range(number):
            scale = scale + gamma[i] * beta[m][i]**2 / 2
    return 1/sp.random.gamma(shape, 1/scale, 1)


def sample_c(beta, R, alpha_c, tau_c, number):
    #input of beta should be the effect size of common variants
    C = beta.shape[1] - R
    shape = number * alpha_c + C*number/2 + number -1
    scale = number * tau_c + sum(beta[:,R:(R+C)]**2)/2
    return 1/sp.random.gamma(shape, 1/scale, 1)


def sample_sigma_b(b, alpha_b, tau_b, number):
    K = b.shape[1] - 1
    shape = number * alpha_b + K * number/2 + number -1
    scale = number * tau_b + sum(b[:,0:K]**2)/2
    return 1/sp.random.gamma(shape, 1/scale, 1)


def SAME(G, rare, common, annotation, gamma, b, alpha_e, tau_e, alpha_r, tau_r, 
                alpha_c, tau_c, alpha_b, tau_b, Ite, rate):
    geno = sp.hstack((rare, common))
    R = len(gamma)
    C = common.shape[1]
    #calculate sum A_i^T * A_i
    sum_anno = 0
    for i in range(annotation.shape[0]):
        sum_anno = sum_anno + sp.outer(annotation[i],annotation[i])
    
    ### initialization ###
    z = sample_z(annotation, b, gamma, 1)
    sigma_e = 1
    sigma_r = 0.01
    sigma_c = 0.001
    sigma_b = 1
    beta = sp.zeros(shape = (1, R+C))
    for i in range(R):
        if gamma[i] == 1:
            beta[0,i] = sp.random.normal(0, math.sqrt(sigma_r),1)
    for j in range(C):
        beta[0,R+j] = sp.random.normal(0, math.sqrt(sigma_c),1)

    ###store the result, to break the for loop
    Sigma_e = []
    Sigma_r = []
    Sigma_c = []
    Sigma_b = []
    Gamma = []
    Odds = []
    Ratio = []
    Post = []
    Probit = []
    
    ### run the algorithm ###
    for ite in range(1, Ite+1):
        if(int(math.log(ite, rate)) - math.log(ite, rate) == 0):
            number = ite
        else:
            number = rate**int(math.log(ite, rate))
        b = sample_b(sigma_b, sum_anno, annotation, z, number)
        beta = sample_beta_chol_old(G, geno, gamma, sigma_r, sigma_c, sigma_e, number)
        gamma, odds, ratio, post, probit = sample_gamma(G, rare, b, annotation, beta, gamma, sigma_r, number)
        Odds.append(odds)
        Ratio.append(ratio)
        Post.append(post)
        Probit.append(probit)
        z = sample_z(annotation, b, gamma, number)
        sigma_e = sample_e(G, geno, beta, alpha_e, tau_e, number)
        sigma_r = sample_r(gamma, beta, alpha_r, tau_r, number)
        sigma_c = sample_c(beta, R, alpha_c, tau_c, number)
        #sigma_b = sample_sigma_b(b, alpha_b, tau_b, number)
        sigma_b = 1
        
        Sigma_e.append(sigma_e)
        Sigma_r.append(sigma_r)
        Sigma_c.append(sigma_c)
        Sigma_b.append(sigma_b)
        Gamma.append(gamma)
        residual = sp.zeros(shape = (5))
        
        if ite >= 2:
            #residual[0] = abs(Sigma_e[ite-1] - Sigma_e[ite-2])/Sigma_e[ite-2]
            residual[0] = abs(Sigma_r[ite-1] - Sigma_r[ite-2])/Sigma_r[ite-2]
            residual[1] = abs(Sigma_c[ite-1] - Sigma_c[ite-2])/Sigma_c[ite-2]
            residual[2] = abs(Sigma_b[ite-1] - Sigma_b[ite-2])/Sigma_b[ite-2]
            residual[3] = sum(abs(Gamma[ite-1] - Gamma[ite-2]))/R
            if max(residual)<0.009:
                break

    return beta, b, gamma, ite, residual, Sigma_e, Sigma_r, Sigma_c, Sigma_b, Gamma


def Initial_SAME(cutoff, hierarchy, G, rare_geno, common_geno, annotation):
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
    
    result_beta, result_b, result_gamma, result_ite, result_residual, result_Sigma_e, result_Sigma_r, result_Sigma_c, result_Sigma_b, result_Gamma, result_z = SAME(G, rare_geno, common_geno, annotation, gamma, b, alpha_e, tau_e, alpha_r, tau_r, alpha_c, tau_c, alpha_b, tau_b, Ite, rate)
    
    geno = sp.hstack((rare_geno, common_geno))
    ### G|beta,sigma_e
    log_likelihood1 = sum(sp.log((norm.pdf((G-geno.dot(result_beta.mean(axis=0))) / math.sqrt(result_Sigma_e[-1]))))) - sp.log(result_Sigma_e[-1])/2*G.shape[0]
    zero = 10**(-10)
    ### beta_rare | sigma_r
    log_likelihood2 = sum(-(result_beta.mean(axis=0)[0:rare_geno.shape[1]]/sp.sqrt(result_Gamma[-1]*result_Sigma_r[-1]+zero))**2/2)
    ### beta_common | sigma_c
    log_likelihood3 = sum(sp.log(norm.pdf(result_beta.mean(axis=0)[rare_geno.shape[1]:result_beta.shape[1]]/math.sqrt(result_Sigma_c[-1]))))
    Probit = norm.cdf(annotation.dot(result_b.mean(axis=0)))
    ### gamma bernoulli distribution
    log_likelihood4 = sum(result_Gamma[-1]*sp.log(Probit)+(1-result_Gamma[-1])*sp.log(1-Probit))
    ### likelihood of logistic regression
    log_likelihood5 = sum(sp.log(norm.pdf(result_z-annotation.dot(result_b.mean(axis=0)))))
    ### b | sigma_b
    log_likelihood6 = sum(-(result_b.mean(axis=0)/sp.sqrt(result_Sigma_b[-1]))**2/2)-sp.log(result_Sigma_b[-1])/2*annotation.shape[1]
    
    log_likelihood = log_likelihood1+log_likelihood2+log_likelihood3+log_likelihood4+log_likelihood5 + log_likelihood6
    
    return sp.mean(result_beta,axis=0), result_gamma, sp.mean(result_b,axis=0), log_likelihood, log_likelihood1, log_likelihood2, log_likelihood3, log_likelihood4, log_likelihood5, log_likelihood6,result_Sigma_b,result_Sigma_c,result_Sigma_e,result_Sigma_r
