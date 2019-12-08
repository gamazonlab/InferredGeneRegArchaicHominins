
main_path='/data/coxvgi/zhoud2/projects/neanderthals'
df<-read.table(paste0(main_path,'/data/individual_introgression_by_gene_toFilter.txt'),header = T,stringsAsFactors = F)
df[,1]<-sapply(df[,1],function(x) strsplit(x,"[.]")[[1]][1])
sample_df<-data.frame('sampleid'=unique(df[,1]),count=NA,stringsAsFactors = F)

#output exclude sample list for each gene
for (i in 1:nrow(sample_df)){
  print(i)
  sample_df[i,2]<-length(which(df[,1]==sample_df[i,1]))
  gene_sample<-df[which(df[,1]==sample_df[i,1]),-1]
  write.table(gene_sample,paste0(main_path,'/samples_ex/',sample_df[i,1]),quote = F,sep='\t',row.names = F)
}

# v6p_sample<-read.table((paste0(main_path,'/data/v6.fam'),stringsAsFactors = F)
# v6p_sample$id<-paste0(sapply(v6p_sample[,1], function(x) strsplit(x,"[-]")[[1]][1]),'-',sapply(v6p_sample[,1], function(x) strsplit(x,"[-]")[[1]][2]))

#plot
#png(paste0(main_path,'/plots/rm_sample_distribution.png'))
#hist(sample_df$count)
#dev.off()


#write.table(sample_df,paste0(main_path,'/data/rm_sample_distribution.txt'),quote = F,sep = '\t',row.names = F)

