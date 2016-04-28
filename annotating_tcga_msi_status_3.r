#batch 3
library(data.table)
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/brca_prad_thca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_rf_classifier_031115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_rf_classifier_peak_avg_median_scaled_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_tree_classifier_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_tree_classifier_peak_avg_median_scaled_041315.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_sampleinfo_041715.robj')
#load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_sampleinfo_030115.robj')

#data<-data[TUMOR_TYPE=='',TUMOR_TYPE:='THCA']
#data[,TUMOR_TYPE:=droplevels(TUMOR_TYPE)]
#data[TUMOR_TYPE!='',]
#data3<-data
#save(data3,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/brca_prad_thca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric_thca_tumor_type_fixed.robj')
#setkey(data,LOCUS_COORDINATES)
#loc1<-data["10:100008321-100008335"]
#sampleinfo3<-loc1
#sampleinfo3[,TUMOR_TYPE:=droplevels(TUMOR_TYPE)]
#sampleinfo3<-sampleinfo3[match(sampledat3$sample_name,sampleinfo3$SAMPLE_NAME),]
#save(sampleinfo3,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_sampleinfo_041715.robj')

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
#colnames(sampledat3)<-c('sample_name','peak_avg','peak_sd','num_unstable_ks','num_unstable_raw','msi_status_rf','msi_status_tree','msi_status_rf_scaled','msi_status_tree_scaled','tumor_type')
save(sampledat3,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_sampledata_041715.robj')

library(ggplot2)
alllocs<-rbind(sampleinfo1,sampleinfo2,sampleinfo3)
msifracs<-tapply(alllocs$MSI_STATUS,alllocs$TUMOR_TYPE,function(x)length(which(x=='MSI-H'))/length(x))
bardat<-data.frame(cancer=names(msifracs)[-c(1,10)],prop=msifracs[-c(1,10)]*100)
bardat$cancer<-factor(bardat$cancer,levels=bardat$cancer[order(bardat$prop,decreasing=TRUE)])
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_prop_batch1-3.pdf',width=7,height=7)
ggplot(data=bardat, aes(x=cancer, y=prop)) + geom_bar(aes(fill=cancer), stat="identity") + guides(fill=FALSE) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('Proportion of MSI individuals (%)') + coord_cartesian(ylim = c(0,50)) + xlab('')
dev.off()