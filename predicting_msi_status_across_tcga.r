#batch 2
require(data.table)
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_gbmm_skcm_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_rf_classifier_031115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_rf_classifier_peak_avg_median_scaled_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_tree_classifier_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_tree_classifier_peak_avg_median_scaled_041315.robj')

data<-data2
setkey(data,SAMPLE_NAME)
system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #19s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #30s
system.time(num_unstable_ks<-data[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=SAMPLE_NAME]) #34s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=SAMPLE_NAME]) #10s

#setkey(data,LOCUS_COORDINATES)
#loc1<-data["10:100008321-100008335"]
#sampleinfo2<-loc1
#save(sampleinfo2,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_sampleinfo_030115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_sampleinfo_030115.robj')

sampledat2<-data.frame(sample_name=peak_avg$SAMPLE_NAME,peak_avg=peak_avg$V1,peak_sd=peak_sd$V1,num_unstable_ks=num_unstable_ks$V1,num_unstable_raw=num_unstable$V1,num_unstable=num_unstable$V1)
 treepred <- predict(msimodel,newdata=sampledat2,type='raw')
 treepred2 <- predict(rpartfit1,newdata=sampledat2,type='raw')

sampledat2n<-sampledat2
sampledat2n$tumor_type<-sampleinfo2$TUMOR_TYPE
sampledat2n$peak_avg<-ave(sampledat2n$peak_avg,sampledat2n$tumor_type,FUN=function(x)x-median(x))
 treepred3 <- predict(msimodel2,newdata=sampledat2n,type='raw')
 treepred4 <- predict(rpartfit2,newdata=sampledat2n,type='raw')

data2$MSI_STATUS_PRED<-treepred4[match(data2$SAMPLE_NAME,sampledat2$sample_name)]
save(data2,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_gbmm_skcm_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric_msi_predicted.robj')

sampledat2$msi_status_rf<-treepred
sampledat2$msi_status_tree<-treepred2
sampledat2$msi_status_rf_scaled<-treepred3
sampledat2$msi_status_tree_scaled<-treepred4
sampledat2$tumor_type<-sampleinfo2$TUMOR_TYPE
sampledat2<-sampledat2[,-6]
save(sampledat2,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_sampledata_041315.robj')

#batch 3
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/brca_prad_thca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_sampleinfo_041715.robj')

data<-data3
library(data.table)
setkey(data,SAMPLE_NAME)
system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #55s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #22s
system.time(num_unstable_ks<-data[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=SAMPLE_NAME]) #28s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=SAMPLE_NAME]) #16s

sampledat3<-data.frame(sample_name=peak_avg$SAMPLE_NAME,peak_avg=peak_avg$V1,peak_sd=peak_sd$V1,num_unstable_ks=num_unstable_ks$V1,num_unstable_raw=num_unstable$V1,num_unstable=num_unstable$V1)
treepred <- predict(msimodel,newdata=sampledat3,type='raw')
treepred2 <- predict(rpartfit1,newdata=sampledat3,type='raw')

sampledat3n<-sampledat3
sampledat3n$tumor_type<-sampleinfo3$TUMOR_TYPE
sampledat3n$peak_avg<-ave(sampledat3n$peak_avg,sampledat3n$tumor_type,FUN=function(x)x-median(x,na.rm=T))
treepred3 <- predict(msimodel2,newdata=sampledat3n,type='raw')
treepred4 <- predict(rpartfit2,newdata=sampledat3n,type='raw')

data3[,MSI_STATUS_PRED:=treepred4[match(data3$SAMPLE_NAME,sampledat3$sample_name)]]
save(data3,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/brca_prad_thca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric_msi_predicted.robj')

sampledat3$msi_status_rf[which(is.na(sampledat3$peak_avg)==F)]<-as.character(treepred)
sampledat3$msi_status_tree[which(is.na(sampledat3$peak_avg)==F)]<-as.character(treepred2)
sampledat3$msi_status_rf_scaled[which(is.na(sampledat3$peak_avg)==F)]<-as.character(treepred3)
sampledat3$msi_status_tree_scaled[which(is.na(sampledat3$peak_avg)==F)]<-as.character(treepred4)
sampledat3$tumor_type<-sampleinfo3$TUMOR_TYPE
sampledat3<-sampledat3[,-6]
save(sampledat3,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_sampledata_041715.robj')

#batch 4
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/luad_lusc_blca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_sampleinfo_041715.robj')

data<-data4

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
save(sampledat4,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_sampledata_041715.robj')

#batch 5
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/stad_ovmu_kirc_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch5_sampleinfo_041515.robj')

data<-data5

setkey(data,SAMPLE_NAME)
system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #53s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #21s
system.time(num_unstable_ks<-data[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=SAMPLE_NAME]) #28s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=SAMPLE_NAME]) #10s

sampledat5<-data.frame(sample_name=peak_avg$SAMPLE_NAME,peak_avg=peak_avg$V1,peak_sd=peak_sd$V1,num_unstable_ks=num_unstable_ks$V1,num_unstable_raw=num_unstable$V1,num_unstable=num_unstable$V1)
treepred <- predict(msimodel,newdata=sampledat5,type='raw')
treepred2 <- predict(rpartfit1,newdata=sampledat5,type='raw')

sampledat5n<-sampledat5
sampledat5n$tumor_type<-sampleinfo5$TUMOR_TYPE
sampledat5n$peak_avg<-ave(sampledat5n$peak_avg,sampledat5n$tumor_type,FUN=function(x)x-median(x,na.rm=T))
treepred3 <- predict(msimodel2,newdata=sampledat5n,type='raw')
treepred4 <- predict(rpartfit2,newdata=sampledat5n,type='raw')

data5[,MSI_STATUS_PRED:=treepred4[match(data5$SAMPLE_NAME,sampledat4$sample_name)]]
badsamples<-sampledat5$sample_name[which(is.na(sampledat5$peak_sd)==T)]
save(data5,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/stad_ovmu_kirc_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric_msi_predicted.robj')

sampledat5$msi_status_rf[which(is.na(sampledat5$peak_avg)==F)]<-as.character(treepred)
sampledat5$msi_status_tree[which(is.na(sampledat5$peak_avg)==F)]<-as.character(treepred2)
sampledat5$msi_status_rf_scaled[which(is.na(sampledat5$peak_avg)==F)]<-as.character(treepred3)
sampledat5$msi_status_tree_scaled[which(is.na(sampledat5$peak_avg)==F)]<-as.character(treepred4)
sampledat5$tumor_type<-sampleinfo5$TUMOR_TYPE
sampledat5<-sampledat5[,-6]
save(sampledat5,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch5_sampledata_042215.robj')

#batch 6
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/kirp_hnsc_lggm_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch6_sampleinfo_042215.robj')

data<-data6
setkey(data,SAMPLE_NAME)
system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #53s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #21s
system.time(num_unstable_ks<-data[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=SAMPLE_NAME]) #28s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=SAMPLE_NAME]) #10s

sampledat6<-data.frame(sample_name=peak_avg$SAMPLE_NAME,peak_avg=peak_avg$V1,peak_sd=peak_sd$V1,num_unstable_ks=num_unstable_ks$V1,num_unstable_raw=num_unstable$V1,num_unstable=num_unstable$V1)
treepred <- predict(msimodel,newdata=sampledat6,type='raw')
treepred2 <- predict(rpartfit1,newdata=sampledat6,type='raw')

sampledat6n<-sampledat6
sampledat6n$tumor_type<-sampleinfo6$TUMOR_TYPE
sampledat6n$peak_avg<-ave(sampledat6n$peak_avg,sampledat6n$tumor_type,FUN=function(x)x-median(x,na.rm=T))
treepred3 <- predict(msimodel2,newdata=sampledat6n,type='raw')
treepred4 <- predict(rpartfit2,newdata=sampledat6n,type='raw')

data6$MSI_STATUS_PRED<-treepred4[match(data6$SAMPLE_NAME,sampledat6$sample_name)]
save(data6,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/kirp_hnsc_lggm_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric_msi_predicted.robj')

sampledat6$msi_status_rf[which(is.na(sampledat6$peak_avg)==F)]<-as.character(treepred)
sampledat6$msi_status_tree[which(is.na(sampledat6$peak_avg)==F)]<-as.character(treepred2)
sampledat6$msi_status_rf_scaled[which(is.na(sampledat6$peak_avg)==F)]<-as.character(treepred3)
sampledat6$msi_status_tree_scaled[which(is.na(sampledat6$peak_avg)==F)]<-as.character(treepred4)
sampledat6$tumor_type<-sampleinfo6$TUMOR_TYPE
sampledat6<-sampledat6[,-6]
save(sampledat6,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch6_sampledata_042215.robj')

#aggregating results
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampledata_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_sampledata_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_sampledata_041715.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_sampledata_041715.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch5_sampledata_042215.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch6_sampledata_042215.robj')
sampledat<-rbind(sampledat1,sampledat2,sampledat3,sampledat4,sampledat5,sampledat6)
save(sampledat,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042215.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_030115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_sampleinfo_030115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_sampleinfo_030115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_sampleinfo_040115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch5_sampleinfo_041515.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch6_sampleinfo_040115.robj')
sampleinfo<-rbind(sampleinfo1,sampleinfo2,sampleinfo3,sampleinfo4,sampleinfo5,sampleinfo6)
save(sampleinfo,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042215.robj')