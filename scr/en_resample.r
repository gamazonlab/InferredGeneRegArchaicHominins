
args<-as.numeric(commandArgs(TRUE))
library('glmnet')
main_dir<-'/data/coxvgi/zhoud2/projects/cross_tissue/v6p/'
#'/data/coxvgi/zhoud2/projects/cross_tissue/exp_v6p/'

run_id<-1
run_list<-list()
for (i in c(1,6,18,23,24,28,29,30,36,41,44)){ #1:44
  for (j in 1:50){  #100*200 jobs 5h
    run_list[[run_id]]<-c(i,j)
    run_id=run_id+1
  }
}

folder='neanderthals'
resample_time=100
resample_p=0.5

run_i<-run_list[[args[1]]]
print(run_i)

#tissue
tissue_list<-dir('/data/coxvgi/zhoud2/data/gtex/exp/v6p/exp_residual/')
tissue_list<-as.character(sapply(tissue_list,function(x) sub('..........................$','',x))) 
tissue<-tissue_list[run_i[1]] 

#median_rpkm
median_exp<-read.table('/data/coxvgi/zhoud2/data/gtex/exp/v6p/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct',header = T,stringsAsFactors = F,sep='\t')
median_exp<-median_exp[,c(1:23,25,24,26:ncol(median_exp))]
median_exp<-median_exp[,-c(9,10,21,22,26,27,33,36,39)]
median_exp<-median_exp[,c(1,2,run_i[1]+2)]
median_exp<-median_exp[median_exp[,3]>0.1,]
median_exp[,1]<-sapply(median_exp[,1], function(x) strsplit(x,"[.]")[[1]][1])

#mkdir
com<-paste0('mkdir ',main_dir,'out/',folder,'/ ; mkdir ',main_dir,'out/',folder,'/',tissue)
system(command = com, wait = T)

#load expression
exp_raw<-read.table(paste0('/data/coxvgi/zhoud2/data/gtex/exp/v6p/exp_residual/',tissue,'_Analysis.v6p.exp.residual'),header = T,stringsAsFactors = F)

#gene list
geno_list<-dir('/data/coxvgi/zhoud2/projects/cross_tissue/v6p/geno/gene_1m/')
gene_list<-intersect(exp_raw[,1],geno_list)
gene_list<-gene_list[which(gene_list %in% median_exp[,1])]

#neanderthals snp list
rm_snp_list<-read.table('/data/coxvgi/zhoud2/projects/neanderthals/data/intro_snps_0.1_r2_1kG_SNPs_rsIDs.bed',header = F,stringsAsFactors = F)

#start_id
i_start=(run_i[2]-1)*500+1
i_end=run_i[2]*500
if (i_end>length(gene_list)){
  i_end=length(gene_list)
}

output_all<-output<-as.data.frame(matrix(data=NA,ncol=17,nrow=0))
colnames(output_all)<-colnames(output)<-c('geneid','genename','notes','median_exp','N_neand','N_human','N_snps','N_snps_rm','r_rm_sample','p_rm_sample','r','p','r_rm_snp','p_rm_snp','r_rm_sample_snp','p_rm_sample_snp','marker')

#---------------------------
#en(box_gtex_ct,box_geuvadis,0,gene_list[i],i)

en<-function(cv_df,test_df){  
  out_list<-list()
  y<-as.matrix(cv_df[,2]); x<-as.matrix(cv_df[,c(3:ncol(cv_df))])
  
  set.seed(2019)
  fit<-cv.glmnet(x=x,y=y,alpha=0.5,nfolds = 5,keep = T)
  fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm))
  best.lam <- fit.df[which.min(fit.df[,1]),]
  nrow.best = best.lam[,3] ##position of best lambda in cv.glmnet output
  ret <- as.numeric(fit$glmnet.fit$beta[,nrow.best])
  non_zero<-which(ret!=0)
  if (length(non_zero)!=0){
    ret<-ret[non_zero]
    predicted<-test_df[,c(3:ncol(test_df))]
    predicted<-as.data.frame(predicted[,non_zero])
    
    for (j in 1:length(non_zero)){
      predicted[,j]=predicted[,j]*ret[j]
    }
    
    predicted$pre_exp=NA
    for (k in 1:nrow(predicted)){
      predicted[k,ncol(predicted)]=sum(predicted[k,],na.rm = T)
    }
    rep_test<-cor.test(test_df[,2],predicted$pre_exp)
    
    r=c()
    p=c()
    
    r[1]=rep_test$estimate
    p[1]=rep_test$p.value
    
    out_list$beta<-ret
    out_list$pos<-non_zero
    
    #--resample
    for (sample_i in 2:resample_time){
      set.seed(sample_i)
      resample_id=sample(seq(1:nrow(predicted)),round(nrow(predicted)*resample_p))
      rep_test<-cor.test(test_df[resample_id,2],predicted$pre_exp[resample_id])
      r[sample_i]=rep_test$estimate
      p[sample_i]=rep_test$p.value
    }
    
    out_list$r<-r
    out_list$p<-p
    
    
  }else{
    out_list$r<-NA
    out_list$p<-NA
  }
  return(out_list)
}


#---------------------------
if (i_start<length(gene_list)){
  for (i in i_start:i_end){
    print(i)
    output[1:resample_time,1]<-gene_list[i]
    output[,2]<-median_exp[which(median_exp[,1]==gene_list[i]),2]
    output[,4]<-median_exp[which(median_exp[,1]==gene_list[i]),3]
    output[,17]<-'resampled'
    output[1,17]<-'all'
    
    #input geno
    load(paste0('/data/coxvgi/zhoud2/projects/cross_tissue/v6p/geno/gene/',gene_list[i]))
    #geno_gtex<-read.table(paste0('/data/coxvgi/zhoud2/projects/cross_tissue/v6p/geno/gene/',gene_list[i]),header=T,stringsAsFactors =F)
    if (ncol(geno_gtex)==2){
      output[,3]<-'only_one_snp'
      next
    }
    
    #input expression
    exp_gtex<-exp_raw[which(exp_raw[,1]==gene_list[i]),]
    exp_gtex<-data.frame('sampleid'=colnames(exp_gtex[-1]),'exp'=as.numeric(exp_gtex[1,-1]),stringsAsFactors = F)
    
    #merge geno exp
    df<-merge(exp_gtex,geno_gtex,by=1)
    
    #snp_info
    load(paste0('/data/coxvgi/zhoud2/projects/cross_tissue/v6p/geno/snp_info/',gene_list[i]))
    #snp_info<-read.table(paste0('/data/coxvgi/zhoud2/projects/cross_tissue/v6p/geno/snp_info_1m/',gene_list[i]),header=T,stringsAsFactors =F)
    
    #load sample list with neanderthals haplotypes
    neanderthals_sample<-try(read.table(paste0('/data/coxvgi/zhoud2/projects/neanderthals/samples_ex/',gene_list[i]),header = T,stringsAsFactors = F))
    if('try-error' %in% class(neanderthals_sample)){
      N_neand=0
      N_human=nrow(df)-N_neand
    }else{
      neanderthals_sample[,1]<-paste0('GTEX.',sapply(neanderthals_sample[,1], function(x) strsplit(x,"[-]")[[1]][2]))
      #count N
      N_neand<-length(which(df[,1] %in% neanderthals_sample[,1]))
      N_human<-nrow(df)-N_neand
    }
    
    output[,5]<-N_neand
    output[,6]<-N_human
    
    #rm snp info
    output[,7]<-ncol(df)-2
    hg19_pos<-as.numeric(sapply(colnames(df)[-c(1,2)], function(x) strsplit(x,"[:]")[[1]][2]))
    gene_chr<-as.numeric(sapply(colnames(df)[-c(1,2)], function(x) strsplit(x,"[:]")[[1]][1]))[1]
    rm_snp_list_chr<-rm_snp_list[rm_snp_list[,1]==gene_chr,]
    output[,8]<-length(which(hg19_pos %in% rm_snp_list_chr[,3]))
    
    #consistant sample size for 4 groups
    if (N_neand>50 & N_human>100){
      
      #T1: rm neanderthal samples
      neanderthals_cv_sample<-which(df[,1] %in% neanderthals_sample[,1])
      cv_df=df[-neanderthals_cv_sample,]
      test_df=df[neanderthals_cv_sample,]
      fit<-en(cv_df,test_df)
      output[,9]<-fit$r; output[,10]<-fit$p
      if(!is.na(fit$r)[1]){
        weights<-snp_info[fit$pos,]
        weights$weights=fit$beta
        write.table(weights,paste0('/data/coxvgi/zhoud2/projects/neanderthals/pred_model/T1_rm_samples/',gene_list[i]),quote = F,row.names = F,sep='\t')
      }
      
      #T2: regular
      set.seed(2019)
      cv_sample=sample(seq(1,nrow(df)),nrow(cv_df),replace = F)
      cv_df=df[cv_sample,]
      test_df=df[-cv_sample,]
      fit<-en(cv_df,test_df)
      output[,11]<-fit$r; output[,12]<-fit$p
      if(!is.na(fit$r)[1]){
        weights<-snp_info[fit$pos,]
        weights$weights=fit$beta
        write.table(weights,paste0('/data/coxvgi/zhoud2/projects/neanderthals/pred_model/T2_regular/',gene_list[i]),quote = F,row.names = F,sep='\t')
      }
      
      
      #T3: rm neanderthal snp
      if(output[1,8]>0){
        df_rm_snp<-df[,-c(which(hg19_pos %in% rm_snp_list_chr[,3])+2)]
      }else{
        df_rm_snp<-df
      }
      if (ncol(df_rm_snp)<=3){
        output[,3]<-'no_snp_or_only_one_snp_left_after_rm_neand'
      }else{
        cv_df=df_rm_snp[cv_sample,]
        test_df=df_rm_snp[-cv_sample,]
        fit<-en(cv_df,test_df)
        output[,13]<-fit$r; output[,14]<-fit$p
        if(!is.na(fit$r)[1]){
          weights<-snp_info[fit$pos,]
          weights$weights=fit$beta
          write.table(weights,paste0('/data/coxvgi/zhoud2/projects/neanderthals/pred_model/T3_rm_snps/',gene_list[i]),quote = F,row.names = F,sep='\t')
        }
      }
      
      #T4: rm neanderthal snp and samples
      if (ncol(df_rm_snp)>3){
        cv_df=df_rm_snp[-neanderthals_cv_sample,]
        test_df=df_rm_snp[neanderthals_cv_sample,]
        fit<-en(cv_df,test_df)
        output[,15]<-fit$r; output[,16]<-fit$p
        if(!is.na(fit$r)[1]){
          weights<-snp_info[fit$pos,]
          weights$weights=fit$beta
          write.table(weights,paste0('/data/coxvgi/zhoud2/projects/neanderthals/pred_model/T4_rm_samples_snps/',gene_list[i]),quote = F,row.names = F,sep='\t')
        }
      }
    }
    
    output_all=rbind(output_all,output)
  }
  
  out_path<-paste0(main_dir,'out/',folder ,'/',tissue,'/',i_start,'_',i_end)
  output_all<-output_all[!is.na(output_all[,1]),]
  write.table(output_all,file=out_path,sep='\t',quote=F,row.names=F)
  
}else{
  print('out of range')
}


















