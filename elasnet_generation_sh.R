genes <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/build_model/selected_gene_location.txt",header=T,stringsAsFactors=F)
tmp <- strsplit(genes[,1],"/")
chrs <- sapply(tmp,"[[",11)
gene <- sapply(tmp,"[[",12) ##plink prefix and .frq

h2 = c(0.05,0.15,0.25)
rareh = c(0.1,0.3,0.5)
pi = c(0.05,0.1,0.2)
annoh = c(0.1,0.3,0.5)

sink("elasnet_generation.sh")
for (h in 1:length(h2)) {
  for (r in 1:length(rareh)) {
    for (p in 1:length(pi)) {
      for (a in 1:length(annoh)) {
        for(i in 1:length(chrs)){
          # module load RXX;module load XX;
          sent_dir <- paste0("mkdir -p /gpfs/ysm/scratch60/zhao/yd327/CRAG/simul/",chrs[i],"/",gene[i],"/","h_",h2[h],"_rareh_",rareh[r],"_pi_",pi[p],"_hanno_",annoh[a])
          loc = paste0("/gpfs/ysm/scratch60/zhao/yd327/CRAG/simul/",chrs[i],"/",gene[i],"/","h_",h2[h],"_rareh_",rareh[r],"_pi_",pi[p],"_hanno_",annoh[a])
          cat(sent_dir,"\n")
          sent_generate <- paste0("Rscript /gpfs/ysm/scratch60/zhao/yd327/CRAG_server/elasnet_generation.R ",genes[i,2]," ",paste0(genes[i,2],".frq")," ",genes[i,1]," ",paste0(loc,"/","gamma.txt")," ",paste0(loc,"/","expression.txt")," ",paste0(loc,"/","true_evaluation.txt")," ",h2[h]," ",rareh[r]," ",pi[p]," ",annoh[a])
          cat(sent_generate,"\n")
        }
      }
    }
  }
}
sink()