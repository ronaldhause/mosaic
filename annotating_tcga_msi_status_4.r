#batch 4
library(data.table)
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/luad_lusc_blca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_rf_classifier_031115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_rf_classifier_peak_avg_median_scaled_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_tree_classifier_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_tree_classifier_peak_avg_median_scaled_041315.robj')

#load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_sampleinfo_040115.robj')
#sampleinfo4<-sampleinfo4[match(sampledat4$sample_name,sampleinfo4$SAMPLE_NAME),]
#save(sampleinfo4,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_sampleinfo_041515.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_sampleinfo_041715.robj')

data<-data4
#data<-data[TUMOR_TYPE=='',TUMOR_TYPE:='LUAD']
#data[,TUMOR_TYPE:=droplevels(TUMOR_TYPE)]
#data4<-data
#save(data4,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/luad_lusc_blca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric_luad_tumor_type_fixed.robj')
#setkey(data,LOCUS_COORDINATES)
#loc1<-data["10:100008321-100008335"]
#sampleinfo4<-loc1
#sampleinfo4<-sampleinfo4[match(sampledat4$sample_name,sampleinfo4$SAMPLE_NAME),]
#save(sampleinfo4,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_sampleinfo_041715.robj')

setkey(data,SAMPLE_NAME)
system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #69
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #38s
system.time(num_unstable_ks<-data[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=SAMPLE_NAME]) #28s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=SAMPLE_NAME]) #10s

sampledat4<-data.frame(sample_name=peak_avg$SAMPLE_NAME,peak_avg=peak_avg$V1,peak_sd=peak_sd$V1,num_unstable_ks=num_unstable_ks$V1,num_unstable_raw=num_unstable$V1,num_unstable=num_unstable$V1)
treepred <- predict(msimodel,newdata=sampledat4,type='raw')
treepred2 <- predict(rpartfit1,newdata=sampledat4,type='raw')

sampledat4n<-sampledat4
sampledat4n$tumor_type<-sampleinfo4$TUMOR_TYPE
sampledat4n$peak_avg<-ave(sampledat4n$peak_avg,sampledat4n$tumor_type,FUN=function(x)x-median(x,na.rm=T))
treepred3 <- predict(msimodel2,newdata=sampledat4n,type='raw')
treepred4 <- predict(rpartfit2,newdata=sampledat4n,type='raw')

data4[,MSI_STATUS_PRED:=treepred4[match(data4$SAMPLE_NAME,sampledat4$sample_name)]]
save(data4,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/luad_lusc_blca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric_msi_predicted.robj')

sampledat4$msi_status_rf[which(is.na(sampledat4$peak_avg)==F)]<-as.character(treepred)
sampledat4$msi_status_tree[which(is.na(sampledat4$peak_avg)==F)]<-as.character(treepred2)
sampledat4$msi_status_rf_scaled[which(is.na(sampledat4$peak_avg)==F)]<-as.character(treepred3)
sampledat4$msi_status_tree_scaled[which(is.na(sampledat4$peak_avg)==F)]<-as.character(treepred4)
sampledat4$tumor_type<-sampleinfo4$TUMOR_TYPE
sampledat4<-sampledat4[,-6]
#colnames(sampledat4)<-c('sample_name','peak_avg','peak_sd','num_unstable_ks','num_unstable_raw','msi_status_rf','msi_status_tree','msi_status_rf_scaled','msi_status_tree_scaled','tumor_type')
save(sampledat4,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_sampledata_041715.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampledata_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_sampledata_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_sampledata_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_sampledata_041315.robj')
sampledat<-rbind(sampledat1,sampledat2,sampledat4,sampledat4)

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_030115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_sampleinfo_030115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_sampleinfo_030115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_sampleinfo_040115.robj')

library(ggplot2)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/frac_unstable_per_cancer_batches1-4_040915_violin.pdf',width=7*1.75,height=7)
ggplot(sampledat,aes(y=log10(sampledat$num_unstable_raw),x=sampledat$tumor_type)) + geom_violin() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + geom_hline(aes(yintercept=3.5),linetype='dashed',color='red')
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/peak_avg_per_cancer_batches1-4_040915.pdf',width=7*1.5,height=7)
ggplot(sampledat,aes(y=sampledat$peak_avg,x=sampledat$tumor_type)) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + geom_hline(aes(yintercept=0.005),linetype='dashed',color='red') + coord_cartesian(ylim = c(-0.01,0.02))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/peak_avg_per_cancer_batches1-4_040915_unscaled.pdf',width=7*1.5,height=7)
ggplot(sampledat,aes(y=sampledat$peak_avg,x=sampledat$tumor_type)) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + geom_hline(aes(yintercept=0.005),linetype='dashed',color='red')
dev.off()

library(ggplot2)
msifracs<-tapply(sampledat$msi_status_rf,sampledat$tumor_type,function(x)length(which(x=='MSI-H'))/length(x))
bardat<-data.frame(cancer=names(msifracs)[-c(1,8)],prop=msifracs[-c(1,8)]*100)
bardat$cancer<-factor(bardat$cancer,levels=bardat$cancer[order(bardat$prop,decreasing=TRUE)])
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_prop_rf_batch1-4_041315.pdf',width=7*1.5,height=7)
ggplot(data=bardat, aes(x=cancer, y=prop)) + geom_bar(aes(fill=cancer), stat="identity") + guides(fill=FALSE) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('Proportion of MSI individuals (%)') + coord_cartesian(ylim = c(0,50)) + xlab('')
dev.off()

msifracs<-tapply(sampledat$msi_status_rf_scaled,sampledat$tumor_type,function(x)length(which(x=='MSI-H'))/length(x))
bardat<-data.frame(cancer=names(msifracs)[-c(1,8)],prop=msifracs[-c(1,8)]*100)
bardat$cancer<-factor(bardat$cancer,levels=bardat$cancer[order(bardat$prop,decreasing=TRUE)])
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_prop_rf_scaled_batch1-4_041315.pdf',width=7*1.5,height=7)
ggplot(data=bardat, aes(x=cancer, y=prop)) + geom_bar(aes(fill=cancer), stat="identity") + guides(fill=FALSE) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('Proportion of MSI individuals (%)') + coord_cartesian(ylim = c(0,50)) + xlab('')
dev.off()

msifracs<-tapply(sampledat$msi_status_tree,sampledat$tumor_type,function(x)length(which(x=='MSI-H'))/length(x))
bardat<-data.frame(cancer=names(msifracs)[-c(1,8)],prop=msifracs[-c(1,8)]*100)
bardat$cancer<-factor(bardat$cancer,levels=bardat$cancer[order(bardat$prop,decreasing=TRUE)])
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_prop_tree_batch1-4_041315.pdf',width=7*1.5,height=7)
ggplot(data=bardat, aes(x=cancer, y=prop)) + geom_bar(aes(fill=cancer), stat="identity") + guides(fill=FALSE) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('Proportion of MSI individuals (%)') + coord_cartesian(ylim = c(0,50)) + xlab('')
dev.off()

msifracs<-tapply(sampledat$msi_status_tree_scaled,sampledat$tumor_type,function(x)length(which(x=='MSI-H'))/length(x))
bardat<-data.frame(cancer=names(msifracs)[-c(1,8)],prop=msifracs[-c(1,8)]*100)
bardat$cancer<-factor(bardat$cancer,levels=bardat$cancer[order(bardat$prop,decreasing=TRUE)])
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_prop_tree_scaled_batch1-4_041315.pdf',width=7*1.5,height=7)
ggplot(data=bardat, aes(x=cancer, y=prop)) + geom_bar(aes(fill=cancer), stat="identity") + guides(fill=FALSE) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('Proportion of MSI individuals (%)') + coord_cartesian(ylim = c(0,50)) + xlab('')
dev.off()

#3 - core loci associated with MSI

setkey(data,LOCUS_COORDINATES)
loc1<-data["10:100008321-100008335"]

setkey(data,MSI_STATUS)
msidata<-data["MSI-H"]
mssdata<-data["MSS"]

setkey(msidata,LOCUS_COORDINATES)
setkey(mssdata,LOCUS_COORDINATES)
system.time(msi_mean_locus_peaks<-msidata[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=LOCUS_COORDINATES]) #4s
system.time(msi_median_locus_peaks<-msidata[,median(PEAK_DIFFERENCE_VALUE,na.rm=T),by=LOCUS_COORDINATES]) #48s
system.time(mss_mean_locus_peaks<-mssdata[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=LOCUS_COORDINATES]) #12s
system.time(msi_locus_nas<-msidata[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=LOCUS_COORDINATES]) #12s
system.time(mss_locus_nas<-mssdata[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=LOCUS_COORDINATES]) #12s
system.time(msi_sig_locs<-msidata[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=LOCUS_COORDINATES]) #4s
system.time(mss_sig_locs<-mssdata[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=LOCUS_COORDINATES]) #10s
system.time(msi_pos_locs<-msidata[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=LOCUS_COORDINATES]) #4s
system.time(mss_pos_locs<-mssdata[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=LOCUS_COORDINATES]) #10s
#84 MSI-H, 331 MSS

mssnum<-length(which(loc1$MSI_STATUS=='MSS'))
msinum<-length(which(loc1$MSI_STATUS=='MSI-H'))


sites_with_sufficient_data<-intersect(which(msi_locus_nas$V1<=(msinum/2)),which(mss_locus_nas$V1<(mssnum/2)))
locus_peak_diffs<-msi_mean_locus_peaks$V1-mss_mean_locus_peaks$V1
pos_peak_diffs<-which(locus_peak_diffs>0)
sub<-intersect(pos_peak_diffs,sites_with_sufficient_data)
sublocs<-as.character(msi_locus_nas$LOCUS_COORDINATES[sub])

a<-data.table(locus=msi_sig_locs$LOCUS_COORDINATES,msi_sig_locs$V1,mss_sig_locs$V1,msi_pos_locs$V1,mss_pos_locs$V1,msi_locus_nas$V1,mss_locus_nas$V1)[sub,]
setkey(a,locus)
system.time(fisher_odds<-a[,fisher.test(matrix(c(V4,msinum-V4-V6,V5,mssnum-V5-V7),nrow=2,byrow=T))$estimate,by=locus]) #71s
system.time(fisher_pvals<-a[,fisher.test(matrix(c(V4,msinum-V4-V6,V5,mssnum-V5-V7),nrow=2,byrow=T))$p.value,by=locus]) #71s
test<-a[1,]

source('/net/shendure/vol1/home/hauser/Scripts/useful_functions.r')
fisher_pvals[,qvals:=qvalue2(fisher_pvals$V1)$qvalues]

length(intersect(which(fisher_pvals$qvals<1e-20),which(fisher_odds$V1>1))) #21, 27

#positive numbers = msi
length(which(abs(msi_mean_locus_peaks$V1)>0)) #60540/516876, 11.7%; 62455, 12.1%
length(which(abs(mss_mean_locus_peaks$V1)>0)) #67361/516876, 13.0%; 74873, 14.5%
length(which(msi_mean_locus_peaks$V1>0)) #54010/516876, 10.4%; 55911, 10.8%
length(which(mss_mean_locus_peaks$V1>0)) #36103/516876, 7.0%; 40070, 7.8%

setkey(data,LOCUS_COORDINATES)
sig_sites<-a[order(fisher_pvals$qvals),]

a<-data.table(locus=msi_sig_locs$LOCUS_COORDINATES,msi_sig_locs$V1,mss_sig_locs$V1,msi_pos_locs$V1,mss_pos_locs$V1,msi_locus_nas$V1,mss_locus_nas$V1)[sub,]
merged.data.frame = Reduce(function(...) merge(..., by="locus",all=T), list(sig_sites,fisher_pvals,fisher_odds))
msi_mss_locus_output<-merged.data.frame[order(merged.data.frame$qvals),]
msi_mss_locus_output$V6<-msinum-msi_mss_locus_output$V6
msi_mss_locus_output$V7<-mssnum-msi_mss_locus_output$V7
colnames(msi_mss_locus_output)<-c('locus','msi_samples_ks_sig','mss_samples_ks_sig','msi_peak_unstable','mss_peak_unstable','msi_locus_calls','mss_locus_calls','fisher_test_pval','fisher_test_qval','fisher_test_odds')
msi_mss_locus_output$genomic_class<-locusinfo$GENOMIC_CLASS[match(msi_mss_locus_output$locus,locusinfo$LOCUS_COORDINATES)]
msi_mss_locus_output$gene<-locusinfo$GENE[match(msi_mss_locus_output$locus,locusinfo$LOCUS_COORDINATES)]
msi_mss_locus_output$repeat_type<-locusinfo$REPEAT_TYPE[match(msi_mss_locus_output$locus,locusinfo$LOCUS_COORDINATES)]
msi_mss_locus_output$repeat_dna_sequence<-locusinfo$REPEAT_DNA_SEQUENCE[match(msi_mss_locus_output$locus,locusinfo$LOCUS_COORDINATES)]
write.table(msi_mss_locus_output,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_subset_with_msi_calls_loci_results_post_prediction_031115.txt',quote=F,row.names=F)
write.csv(msi_mss_locus_output,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_subset_with_msi_calls_loci_results_post_prediction_031115.csv',quote=F,row.names=F)