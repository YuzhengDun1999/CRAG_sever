library("glmnet")
library('genio')
Args <- commandArgs()
plink = Args[6]##plink file input
frq_inp = Args[7]##frq.txt input
anno_inp = Args[8]##annotation file input
gamma_out = Args[9]##initialized gamma file output
expression_out = Args[10]##gene expression file output
perfer_out = Arg[10]##perfermance file output

bim = read_bim(paste(plink,".bim",sep=""))
fam = read_fam(paste("test",".fam",sep=""))
bed = read_bed(paste("test",".bed",sep=""), bim$id, fam$id)
frq = read.table(frq_inp,header = TRUE)
ano = read.table(anno_inp,header = TRUE)
ano = ano[,-c(1:3)]
rare_index = c()
common_index = c()
for(i in 1:length(frq$MAF)){
  if(frq$MAF[i]>0.005 & frq$MAF[i]<0.05){
    rare_index = c(rare_index, i)
  }
  else if(frq$MAF[i]>=0.05){
    common_index = c(common_index, i)
  }
}


rare_geno = t(bed[rare_index,])#rare variants genotype matrix,row=samples,col=rare variant
common_geno = t(bed[common_index,])#common variants genotype matrix,row=samples,col=common variant
annotation = matrix(unlist(ano[rare_index,1:6]), nrow = length(rare_index), ncol = 6)#annotation of rare variants
annotation = cbind(annotation, rep(1,nrow(annotation)))#add intercept
generation_CommonRare <- function(common,rare, annotation,h2_anno,h2_r,h2_c,pi){
  
  G = c()
  #  N = nrow(rare)
  R = ncol(rare)
  N = nrow(common)
  C = ncol(common)
  # first deal with NA in the genoytype matrix
  which(is.na(common),arr.ind=T) -> inds_na
  mean_col <- colMeans(common,na.rm = T)
  common[inds_na] <- mean_col[inds_na[,2]]
  
  which(is.na(rare),arr.ind=T) -> inds_na
  mean_col <- colMeans(rare,na.rm = T)
  rare[inds_na] <- mean_col[inds_na[,2]]
  
  geno = cbind(rare,common)
  # first scale the genotype matrix
  geno <- scale(geno)
  common <- scale(common)
  rare <- scale(rare)
  #  geno = cbind(rare, common)
  
  
  ###real parameters
  # number of SNPs being causal 
  p_causal <- floor(pi*R)
  sigma_r = h2_r/p_causal
  
  sigma_b = h2_anno/(ncol(annotation)-1)
  b = rnorm(ncol(annotation)-1, mean = 0, sd = sqrt(sigma_b))
  #b = c(b,-1)
  # first determine random noise in the annotation layer
  anno_error <- ((1 - h2_anno)/h2_anno)*var(as.matrix(annotation[,-c(ncol(annotation))],ncol=length(b)) %*% as.matrix(b,ncol=1))
  b <- b/sqrt(anno_error[1,1])
  sigma_b <- h2_anno/(ncol(annotation)-1)/anno_error[1,1]
  
  # determine the intercept based on how many values larger than 0
  anno_noise <- rnorm(R,mean=0,sd=1)
  anno_val <- as.matrix(annotation[,-c(ncol(annotation))]) %*% as.matrix(b,ncol=1) + anno_noise
  anno_intercept <- -anno_val[order(anno_val,decreasing = T)[p_causal+1]]
  b <- c(b,anno_intercept)
  w <- as.matrix(annotation) %*% as.matrix(b,ncol=1) + anno_noise
  
  gamma = rep(0,R)
  beta_rare = rep(0,R)
  gamma[which(w>0)] <- 1
  beta_rare[which(w>0)] <- rnorm(sum(w>0),mean=0,sd=sqrt(sigma_r))
  
  sigma_c <- h2_c/C
  beta_common = rnorm(C,mean=0,sd=sqrt(sigma_c))
  
  sigma_e = (1-h2_c-h2_r)/(h2_c+h2_r)* var(rare %*% beta_rare + common %*% beta_common)[1,1]
  noise <- rnorm(N,mean = 0, sd=sqrt(sigma_e))
  noise <- noise/sd(noise)*sqrt(sigma_e)
  
  beta <- c(beta_rare,beta_common)
  G <- geno %*% beta + noise
  
  return(list(G = G, common = common, rare = rare, gamma = gamma, beta = beta, sigma_r = sigma_r, sigma_c = sigma_c, sigma_e = sigma_e, sigma_b = sigma_b, b = b, generate_z = w))
}

# new generation function
generate <- generation_CommonRare(common = common_geno,rare = rare_geno,annotation = annotation, 
                                  h2_anno = 0.5,h2_r = 0.2,h2_c = 0.1,pi = 0.05)

G = generate$G
rare_geno = generate$rare
common_geno = generate$common
geno = cbind(rare_geno,common_geno)
cv.fit1 = cv.glmnet(x=rare_geno,y = G,intercept=0,alpha=0.5,type.measure="deviance")
              #,penalty.factor=c(rep(1,ncol(rare_geno1)),rep(0,ncol(common_geno1))))
coefficients1<-coef(cv.fit1,s=cv.fit1$lambda.min)
temp=rep(1,4610)
rare_coef=coefficients1[1:4610]
temp[rare_coef==0]=0
rownames(G)=NULL
G_generate = geno%*%generate$beta
cor_generate = cor(G,G_generate)
h2_generate = var(G_generate)/var(G)
SSR = sum((G-G_generate)^2)
SST = sum((G-mean(G))^2)
R_square = 1-SSR/SST

write.table(temp,file=gamma_out,sep=" ",row.names = FALSE)
write.table(G,file=expression_out,sep=" ",row.names = FALSE)
write.table(list(gamma=sum(generate$gamma),h2=h2_generate,corr=cor_generate,R2=R_square),file=expression_out,sep=" ")
