#aggregating batch results
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_post_review_msi_predicted_063016.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_sampleinfo_post_review_msi_predicted_070616.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_sampleinfo_post_review_msi_predicted_070716.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_sampleinfo_post_review_msi_predicted_070516.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch5_sampleinfo_post_review_msi_predicted_070516.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch6_sampleinfo_post_review_msi_predicted_070516.robj')
sampleinfo<-rbind(sampleinfo1n,sampleinfo2n,sampleinfo3n,sampleinfo4n,sampleinfo5n,sampleinfo6n)
save(sampleinfo,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_post_review_070816.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mmr_mutational_data.robj')
mmrmutdat[mmrmutdat==0]<-'ND'
mmrmutdat[mmrmutdat==1]<-'nsSNV'
sampleinfo$EXO1[which(sampleinfo$EXO1 %in% c('','ND','nan','NaN'))]<-mmrmutdat$EXO1[match(sampleinfo$SAMPLE_NAME,rownames(mmrmutdat))][which(sampleinfo$EXO1 %in% c('','ND','nan','NaN'))]
sampleinfo$EXO1[is.na(sampleinfo$EXO1)]<-'ND'
for (i in 10:18){
	sampleinfo[which(sampleinfo[,i] %in% c('','ND','nan','NaN')),i]<-mmrmutdat[,i-8][match(sampleinfo$SAMPLE_NAME,rownames(mmrmutdat))][which(sampleinfo[,i] %in% c('','ND','nan','NaN'))]
sampleinfo[,i][is.na(sampleinfo[,i])]<-'ND'
}

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_results_post_review_msi_predicted_063016.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_results_post_review_msi_predicted_070616.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_results_post_review_msi_predicted_070716.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_results_post_review_msi_predicted_070516.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch5_results_post_review_msi_predicted_070516.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch6_results_post_review_msi_predicted_070516.robj')
sampledat<-rbind(sampledat1n,sampledat2n,sampledat3n,sampledat4n,sampledat5n,sampledat6n)
sampledat$tumor_type<-sampleinfo$TUMOR_TYPE
save(sampledat,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')

a<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_tcga_loci_results_060616.txt',header=T)
b<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_tcga_loci_results_070716.txt',header=T)
c<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_tcga_loci_results_060616.txt',header=T)
d<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_tcga_loci_results_060616.txt',header=T)
e<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch5_tcga_loci_results_060616.txt',header=T)
f<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch6_tcga_loci_results_060616.txt',header=T)
results<-rbind(a,b,c,d,e,f)
results$total_num<-results$num_missing+results$num_present
results<-results[,c(1:7,12,8:11)]
results<-as.data.table(results)
save(results,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_results_post_review_msi_updated_062916.robj')
mosaic_training <- sampledat1n
mosaic_training$defbsite <- sampledat1n$defbsite
save(mosaic_training, file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_training_data_062916.robj')

write.csv(sampledat[,c(1,4,8,9)], '/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/supplementary_table_5_sample_information_raw.csv', row.names=F, quote=F)

#ANALYSIS
#1 - instability landscape across all cancer types
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_post_review_070816.robj')
#516,876 sites, 5930 samples
results_allloci$genomic_class<-factor(results_allloci$genomic_class)
results_allloci2<-results_allloci[which(results_allloci$num_present>(5930/2)),]
results_allloci2$prop<-((results_allloci2$num_unstable/results_allloci2$num_present)*100)
#223,082 sites have sufficient data in at least 2965 samples

#2 - instability landscape of specific cancer types
library(data.table)
results_per_cancer<-results[,list(median_peak_diff=median(median_peak_diff,na.rm=T),num_unstable=sum(num_unstable,na.rm=T),num_missing=sum(num_missing,na.rm=T),num_present=sum(num_present,na.rm=T),total_num=sum(total_num,na.rm=T),genomic_class=genomic_class[1],gene=gene[1],repeat_type=repeat_type[1],repeat_dna_sequence=repeat_dna_sequence[1]),by=list(locus,tumor_type)]
save(results_per_cancer,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_per_cancer_070816.robj')

#3 - core loci associated with MSI
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')
library(data.table)
results_msi<-results[,list(median_peak_diff=median(median_peak_diff,na.rm=T),num_unstable=sum(num_unstable,na.rm=T),num_missing=sum(num_missing,na.rm=T),num_present=sum(num_present,na.rm=T),total_num=sum(total_num,na.rm=T),genomic_class=genomic_class[1],gene=gene[1],repeat_type=repeat_type[1],repeat_dna_sequence=repeat_dna_sequence[1]),by=list(locus,msi_status)]
save(results_msi,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_070816.robj')
#5666 MSS samples and 264 MSI-H samples

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_042515.robj')

library(reshape2)
results_msi2<-results_msi[-which(results_msi$num_present<(results_msi$total_num/2)),]
results_msi3<-split(results_msi2,results_msi2$msi_status)
results_msi4 = Map(function(x, i) setNames(x, ifelse(names(x) %in% "locus",
      names(x), sprintf('%s.%d', names(x), i))), results_msi3, seq_along(results_msi3))
results_msi5<-Reduce(function(...) merge(..., by="locus",all=FALSE),results_msi4)
#204,797 sites have sufficient data for both groups

#fisher.test(matrix(c(85,121-85,631,3879-631),nrow=2,byrow=T))

system.time(fisher_odds<-results_msi5[,fisher.test(matrix(c(num_unstable.1,num_present.1-num_unstable.1,num_unstable.2,num_present.2-num_unstable.2),nrow=2,byrow=T))$estimate,by=locus]) #131
system.time(fisher_pvals<-results_msi5[,fisher.test(matrix(c(num_unstable.1,num_present.1-num_unstable.1,num_unstable.2,num_present.2-num_unstable.2),nrow=2,byrow=T))$p.value,by=locus])
fisher_pvals$V1[fisher_pvals$V1>1]<-1
fisher_pvals[,qvals:=qvalue(fisher_pvals$V1)$qvalues]

output3<-data.frame(locus=results_msi5$locus,msi_peak_unstable=results_msi5$num_unstable.1,msi_locus_calls=results_msi5$num_present.1,mss_peak_unstable=results_msi5$num_unstable.2,mss_locus_calls=results_msi5$num_present.2,pval=fisher_pvals$V1,qval=fisher_pvals$qvals,odds_ratio=fisher_odds$V1,genomic_class=results_msi5$genomic_class.1,gene=results_msi5$gene.1,repeat_type=results_msi5$repeat_type.1,repeat_dna_sequence=results_msi5$repeat_dna_sequence.1)
output3<-output3[order(output3$qval,decreasing=FALSE),]
write.table(output3,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_070816.txt',quote=F,row.names=F)
write.csv(output3,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_070816.csv',row.names=F)

output3<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_070816.txt', header=T)

#output3<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_050315.txt',header=T)
#write.csv(output3,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_050315.csv',row.names=F)

bitmap('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_locus_msi_vs_mss_qqplot_070816.tiff',width=7*1.25,height=7,units="in",type="tifflzw",res=300)
gg_qqplot(output3$pval)
dev.off()