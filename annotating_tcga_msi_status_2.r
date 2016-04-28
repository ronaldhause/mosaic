#batch 2
library(data.table)
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_gbmm_skcm_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_rf_classifier_031115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_rf_classifier_peak_avg_median_scaled_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_tree_classifier_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_tree_classifier_peak_avg_median_scaled_041315.robj')

data<-data2
#library(data.table)
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
#colnames(sampledat2)<-c('sample_name','peak_avg','peak_sd','num_unstable_ks','num_unstable_raw','msi_status_rf','msi_status_tree','msi_status_rf_scaled','msi_status_tree_scaled','tumor_type')
save(sampledat2,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_sampledata_041315.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampledata_030115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_sampledata_030115.robj')
sampledat<-rbind(sampledat1,sampledat2)