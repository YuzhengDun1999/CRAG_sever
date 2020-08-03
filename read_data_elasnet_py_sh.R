genes <- read.table("/ysm-gpfs/pi/zhao-data/wl382/cWAS/scripts/git_up/build_model/selected_gene_location.txt",header=T,stringsAsFactors=F)
#genes <- read.table("G:/thesis/code/selected_gene_location.txt",header=T,stringsAsFactors=F)

tmp <- strsplit(genes[,1],"/")
chrs <- sapply(tmp,"[[",11)
gene <- sapply(tmp,"[[",12) ##plink prefix and .frq

h2 = c(0.05,0.15,0.25)
rareh = c(0.1,0.3,0.5)
pi = c(0.05,0.1,0.2)
annoh = c(0.1,0.3,0.5)

sink("read_data_elasnet_py.sh")
for (h in 1:length(h2)) {
  for (r in 1:length(rareh)) {
    for (p in 1:length(pi)) {
      for (a in 1:length(annoh)) {
        for(i in 1:length(chrs)){
          # module load RXX;module load XX;
          sent_dir <- paste0("mkdir -p /ysm-gpfs/pi/zhao/yd327/CRAG/simul/",chrs[i],"/",gene[i],"/","h_",h2[h],"_rareh_",rareh[r],"_pi_",pi[p],"_hanno_",annoh[a])
          output_dir = paste0("/ysm-gpfs/pi/zhao/yd327/CRAG/simul/",chrs[i],"/",gene[i],"/","h_",h2[h],"_rareh_",rareh[r],"_pi_",pi[p],"_hanno_",annoh[a])
          
          ## location of input file
          loc = paste0("/gpfs/ysm/scratch60/zhao/yd327/CRAG/simul/",chrs[i],"/",gene[i],"/","h_",h2[h],"_rareh_",rareh[r],"_pi_",pi[p],"_hanno_",annoh[a])
          cat(sent_dir,"\n")
          sent_result <- paste0("python /gpfs/ysm/scratch60/zhao/yd327/CRAG_sever/read_data_elasnet.py --plink ",genes[i,2]," --expression ",paste0(loc,"/","expression.txt")," --gamma ",paste0(loc,"/","gamma.txt")," --frq ",paste0(genes[i,2],".frq")," --annotation ",genes[i,1]," --outSNP ",paste0(output_dir,"/","SNP.txt")," --outAnno ",paste0(output_dir,"/","anno.txt")," --outPerfer ",paste0(output_dir,"/","pre_evaluation.txt"))
          cat(sent_result,"\n")
        }
      }
    }
  }
}
sink()
