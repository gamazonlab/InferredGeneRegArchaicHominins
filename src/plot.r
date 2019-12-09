
tissue_list<-dir('c:/Users/zhoud2/Desktop/neanderthals/results/')

for (tissue_i in 1:length(tissue_list)){
  #tissue_i=44
  print(tissue_i)
  
  df<-read.table(paste0('c:/Users/zhoud2/Desktop/neanderthals/results/',tissue_list[tissue_i]),header = T,stringsAsFactors = F)
  
  #only all samples
  #df<-df[df$marker=='all',]
  
  df$prop_sample_rm=df$N_neand/(df$N_neand+df$N_human)
  df$prop_snp_rm=df$N_snps_rm/df$N_snps
  
  #r=0 p=1
  for (na_col in c(9,11,13,15)){
    df[is.na(df[,na_col]),na_col]<-0
    df[is.na(df[,na_col+1]),na_col+1]<-1
  }
  
  #to binary
  df$regular<-ifelse(df$r>0.1 & df$p<0.05,1,0)
  df$rm_sample<-ifelse(df$r_rm_sample>0.1 & df$p_rm_sample<0.05,1,0)
  df$rm_snp<-ifelse(df$r_rm_snp>0.1 & df$p_rm_snp<0.05,1,0)
  df$rm_sample_snp<-ifelse(df$r_rm_sample_snp>0.1 & df$p_rm_sample_snp<0.05,1,0)
  
  #sample distribution
  png(paste0('c:/Users/zhoud2/Desktop/neanderthals/plots/distribution/',tissue_list[tissue_i],'.png'),width = 1200,height = 1200,res=200)
  par(mfcol=c(1,1))
  hist(df$prop_sample_rm,xlab = 'proportion of Samples with Neand segments',ylab = 'Frequency (Genes)',main = paste0(tissue_list[tissue_i],'  N=',max(df[,5],df[,6],na.rm=T)))
  dev.off()
  #hist(df$prop_snp_rm,xlab = 'proportion of Neand SNPs in the 1Mb region',ylab = 'Frequency (Genes)',main = tissue_list[tissue_i])
  
  df$col_reg_sample<-mapply(reg_rm_sample,df$regular,df$rm_sample)
  
  df$col_reg_snp<-mapply(reg_rm_snp,df$regular,df$rm_snp)
  df$col_reg_sample_snp<-mapply(reg_rm_sample_snp,df$regular,df$rm_sample_snp)
  
  df_all<-df[df$marker=='all',]
  df_resample<-df[df$marker=='resampled',]
  
  png(paste0('c:/Users/zhoud2/Desktop/neanderthals/plots/r2/',tissue_list[tissue_i],'.png'),width = 1500,height = 1500,res=150)
  par(mfcol=c(2,2))
  
  plot(0.5,0.5,col='white',xlim=c(0,1),ylim = c(0,1),xlab='',ylab='',xaxt="n",yaxt="n",bty='n',main = paste0(tissue_list[tissue_i],'  N=',max(df[,5],df[,6],na.rm=T)),cex.main=1.5) #
  text(0.0,0.8,'Imputable gene defined by r>0.1 & p<0.05',pos=4)
  text(0.0,0.7,'in test set',pos=4)
  text(0.07,0.6,'Imputable in both status',pos=4)
  text(0.07,0.5,'Not imputable in both status',pos=4)
  text(0.07,0.4,'Only imputable in regular status',pos=4)
  text(0.07,0.3,'Only imputable when remove Neand_Samples',pos=4)
  text(0.07,0.2,'Only imputable when remove Neand_SNPs',pos=4)
  text(0.07,0.1,'Only imputable when remove Neand_Samples',pos=4)
  text(0.07,0.0,'& SNPs',pos=4)
  
  points(0.04,0.6,col='black',pch=20,cex=1.5)
  points(0.04,0.5,col='grey',pch=20,cex=1.5)
  points(0.04,0.4,col='firebrick2',pch=20,cex=1.5)
  points(0.04,0.3,col='dodgerblue3',pch=20,cex=1.5)
  points(0.04,0.2,col='green4',pch=20,cex=1.5)
  points(0.04,0.1,col='goldenrod',pch=20,cex=1.5)
  
  points(0.07,0.6,col='black',pch=20,cex=0.5)
  points(0.07,0.5,col='grey',pch=20,cex=0.5)
  points(0.07,0.4,col='firebrick2',pch=20,cex=0.5)
  points(0.07,0.3,col='dodgerblue3',pch=20,cex=0.5)
  points(0.07,0.2,col='green4',pch=20,cex=0.5)
  points(0.07,0.1,col='goldenrod',pch=20,cex=0.5)

  #---
  #df_all<-df_all[order(df_all$col_reg_sample),]
  plot(-1,-1,xlab = 'r2 (remove_Neand_Samples)',ylab='r2 (regular)',bty='l',xlim = c(0,max(df$r_rm_sample)^2),ylim = c(0,max(df$r)^2))
  points(df_resample$r_rm_sample^2,df_resample$r^2,col=df_resample$col_reg_sample,pch=20,cex=0.1)
  points(df_all$r_rm_sample^2,df_all$r^2,col=df_all$col_reg_sample,pch=20,cex=1,lwd=1)
  segments(-1,-1,1,1,col = 'yellow',lty=2)

  #---
  #df_all<-df_all[order(df_all$col_reg_snp),]
  plot(-1,-1,xlab = 'r2 (remove_Neand_SNPs)',ylab='r2 (regular)',bty='l',xlim = c(0,max(df$r_rm_snp)^2),ylim = c(0,max(df$r)^2))
  points(df_resample$r_rm_snp^2,df_resample$r^2,col=df_resample$col_reg_snp,pch=20,cex=0.1)
  points(df_all$r_rm_snp^2,df_all$r^2,col=df_all$col_reg_snp,pch=20,cex=1,lwd=1)
  segments(-1,-1,1,1,col = 'yellow',lty=2)

  #---
  #df_all<-df_all[order(df_all$col_reg_sample_snp),]
  plot(-1,-1,xlab = 'r2 (remove_Neand_Samples & SNPs)',ylab='r2 (regular)',bty='l',xlim = c(0,max(df$r_rm_sample_snp)^2),ylim = c(0,max(df$r)^2))
  points(df_resample$r_rm_sample_snp^2,df_resample$r^2,col=df_resample$col_reg_sample_snp,pch=20,cex=0.1)
  points(df_all$r_rm_sample_snp^2,df_all$r^2,col=df_all$col_reg_sample_snp,pch=20,lwd=1)
  segments(-1,-1,1,1,col = 'yellow',lty=2)

  dev.off()
  
}






reg_rm_sample<-function(a,b){  #a=regular
  if (a==0 & b==0){
    col='grey'
  }else if(a==1 & b==1){
    col='black'
  }else if(a==1 & b==0){
    col='firebrick2'
  }else if(a==0 & b==1){
    col='dodgerblue3'
  }
  return(col)
}

reg_rm_snp<-function(a,b){  #a=regular
  if (a==0 & b==0){
    col='grey'
  }else if(a==1 & b==1){
    col='black'
  }else if(a==1 & b==0){
    col='firebrick2'
  }else if(a==0 & b==1){
    col='green4'
  }
}

reg_rm_sample_snp<-function(a,b){  #a=regular
  if (a==0 & b==0){
    col='grey'
  }else if(a==1 & b==1){
    col='black'
  }else if(a==1 & b==0){
    col='firebrick2'
  }else if(a==0 & b==1){
    col='goldenrod'
  }
}

