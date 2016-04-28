#batch 1, 215 million rows

#batch 2, 233 million rows
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/LIHC_1-8_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/GBMm_1-4_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/SKCM_1-4_merged_full_analysis.robj')
library(data.table)
data2<-rbindlist(list(lihc,gbmm,skcm))
system.time(save(data2,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_gbmm_skcm_merged_analysis.robj')) #14 min
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_gbmm_skcm_merged_analysis.robj')
colnames(data2)<-c('GENOMIC_CLASS','GENE','REPEAT_TYPE','REPEAT_DNA_SEQUENCE','LOCUS_COORDINATES','SAMPLE_NAME','TUMOR_TYPE','MSI_STATUS','EXO1','MLH1','MLH3','MSH2','MSH3','MSH6','PMS1','PMS2','POLD1','POLE','KS_VALUE','PEAK_DIFFERENCE_VALUE')
data2$PEAK_DIFFERENCE_VALUE<-as.numeric(as.matrix(data2$PEAK_DIFFERENCE_VALUE))
data2$KS_VALUE<-as.numeric(as.matrix(data2$KS_VALUE))
system.time(save(data2,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_gbmm_skcm_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')) #14 min

#batch 3, 544 million rows
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/BRCA_1-12_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/PRAD_1-16_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/THCA_1-17_merged_full_analysis.robj')
library(data.table)
data3<-rbindlist(list(brca,prad,thca))
system.time(save(data3,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/brca_prad_thca_merged_analysis.robj')) #14 min
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/brca_prad_thca_merged_analysis.robj')
data3$PEAK_DIFFERENCE_VALUE<-as.numeric(as.matrix(data3$PEAK_DIFFERENCE_VALUE))
data3$KS_VALUE<-as.numeric(as.matrix(data3$KS_VALUE))
system.time(save(data3,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/brca_prad_thca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')) #14 min

#batch 4
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/LUAD_1-18_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/LUSC_1-14_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/BLCA_1-10_merged_full_analysis.robj')
library(data.table)
data4<-rbindlist(list(luad,lusc,blca))
system.time(save(data4,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/luad_lusc_blca_merged_analysis.robj')) #14 min
data4$PEAK_DIFFERENCE_VALUE<-as.numeric(as.matrix(data4$PEAK_DIFFERENCE_VALUE))
data4$KS_VALUE<-as.numeric(as.matrix(data4$KS_VALUE))
system.time(save(data4,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/luad_lusc_blca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')) #14 min

#batch 5
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/STAD_1-14_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/OVMU_1-3_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/KIRC_1-5_merged_full_analysis.robj')
library(data.table)
data5<-rbindlist(list(stad,ovmu,kirc))
system.time(save(data5,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/stad_ovmu_kirc_merged_analysis.robj')) #14 min
data5$PEAK_DIFFERENCE_VALUE<-as.numeric(as.matrix(data5$PEAK_DIFFERENCE_VALUE))
data5$KS_VALUE<-as.numeric(as.matrix(data5$KS_VALUE))
system.time(save(data5,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/stad_ovmu_kirc_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')) #14 min

#batch 6
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/KIRP_1-7_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/HNSC_1-21_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/LGGM_1-22_merged_full_analysis.robj')
library(data.table)
data6<-rbindlist(list(kirp,hnsc,lggm))
system.time(save(data6,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/kirp_hnsc_lggm_merged_analysis.robj')) #14 min
data6$PEAK_DIFFERENCE_VALUE<-as.numeric(as.matrix(data6$PEAK_DIFFERENCE_VALUE))
data6$KS_VALUE<-as.numeric(as.matrix(data6$KS_VALUE))
system.time(save(data6,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/kirp_hnsc_lggm_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')) #14 min