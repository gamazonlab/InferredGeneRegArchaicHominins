args<-as.numeric(commandArgs(TRUE))  #chr 1-22

gtex_path<-'/data/coxvgi/zhoud2/data/gtex/geno/v6/dosage/'

#input GTEx v6p sampleid
sample_list<-read.table('/data/coxvgi/zhoud2/projects/cross_tissue/v6p/geno/gtex/sample_all.txt',header = F,stringsAsFactors = F)

#input genotype file in dosage
gtex<-read.table(paste0(gtex_path,'v6p',args,'.txt'),header = F,stringsAsFactors = F)
colnames(gtex)[7:ncol(gtex)]<-sample_list[,1]
colnames(gtex)[1:6]<-c('chr','varID','pos','counted_gtex','ref_gtex','freq')
gtex$varID<-paste0(gtex$chr,':',gtex$pos)

#generate dosage file for each gene
window_size=1000000

#gene pos annotation
anno<-read.table('/data/coxvgi/zhoud2/anno/gencode/gencode.v27',skip = 5,header=F,stringsAsFactors = F)   #GRch37 v6p
anno<-anno[anno$V2=='gene',]
anno[,5]<-sapply(anno[,5], function(x) strsplit(x,"[.]")[[1]][1])
colnames(anno)[5]<-'geneid'
gene_pos=anno[anno$V1==paste0('chr',args),]

#snp and geno
snp<-gtex[,c(1:5)]
geno_data<-gtex[,-c(1,3,4,5,6)]

for (i in 1:nrow(gene_pos)){
  print(paste0(args,' : ',i,' / ',nrow(gene_pos)))
  left<-which(snp$pos>gene_pos[i,3]-window_size)
  right<-which(snp$pos<gene_pos[i,4]+window_size)
  inter<-intersect(left,right)
  
  #write genotype file for each gene
  s_gene<-geno_data[inter,]
  snp_list<-s_gene[,1]; s_gene<-s_gene[,-1]
  s_gene<-as.data.frame(t(s_gene),stringsAsFactors = F)
  s_gene$sampleid<-as.character(sample_list[,1])
  s_gene<-s_gene[,c(ncol(s_gene),1:(ncol(s_gene)-1))]
  colnames(s_gene)<-c('sampleid',snp_list)
  
  geno_gtex=s_gene
  save(geno_gtex,file=paste0('/data/coxvgi/zhoud2/projects/cross_tissue/v6p/geno/gene/',gene_pos[i,5]))
  
  #write snp info for each gene
  snp_info<-snp[inter,c(2,4,5)]
  save(snp_info,file=paste0('/data/coxvgi/zhoud2/projects/cross_tissue/v6p/geno/snp_info/',gene_pos[i,5]))
}
