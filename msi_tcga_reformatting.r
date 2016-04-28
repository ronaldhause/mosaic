column_names<-read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/TCGA_Analysis_header.txt',header=F,sep='\t')
system.time(ucec<-read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/UCEC_1-3_merged_full_analysis.txt.gz',header=F)) #11 min
system.time(read<-read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/READ_1-3_merged_full_analysis.txt.gz',header=F)) #6 min
system.time(lihc<-read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/LIHC_1-8_merged_full_analysis.txt.gz',header=F)) #15 min
system.time(coad<-read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/COAD_1-9_merged_full_analysis.txt.gz',header=F)) #16 min
system.time(gbmm<-read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/GBMm_1-4_merged_full_analysis.txt.gz',header=F)) #28 min
system.time(skcm<-read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/SKCM_1-4_merged_full_analysis.txt.gz',header=F)) #10 min
system.time(ucec2<-read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/UCEC_4-7_merged_full_analysis.txt.gz',header=F)) #12 min
system.time(brca<-read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/BRCA_1-12_merged_full_analysis.txt.gz',header=F)) #130 min
system.time(prad<-read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/PRAD_1-16_merged_full_analysis.txt.gz',header=F)) #71 min
system.time(thca<-read.delim('/net/shendure/vol7/stevesal/MSI_TCGA_project/THCA_1-17_merged_full_analysis.txt.gz',header=F)) #81 min, line 180906601 did not have 20 elements
system.time(luad<-read.delim('/net/shendure/vol7/stevesal/MSI_TCGA_project/LUAD_1-18_merged_full_analysis.txt.gz',header=F)) #73 min, line 57890113 did not have 20 elements
system.time(lusc<-read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/LUSC_1-14_merged_full_analysis.txt.gz',header=F)) #57 min
system.time(blca<-read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/BLCA_1-10_merged_full_analysis.txt.gz',header=F)) #57 min
system.time(stad<-read.delim('/net/shendure/vol7/stevesal/MSI_TCGA_project/STAD_1-14_merged_full_analysis.txt.gz',header=F)) #68 min
system.time(ovmu<-read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/OV_1-3_merged_full_analysis.txt.gz',header=F)) #10 min
system.time(kirc<-read.delim('/net/shendure/vol7/stevesal/MSI_TCGA_project/KIRC_1-5_merged_full_analysis.txt.gz',header=F)) #19 min
system.time(kirp<-read.delim('/net/shendure/vol7/stevesal/MSI_TCGA_project/KIRP_1-7_merged_full_analysis.txt.gz',header=F)) #45 min
system.time(hnsc<-read.delim('/net/shendure/vol7/stevesal/MSI_TCGA_project/HNSC_1-21_merged_full_analysis.txt.gz',header=F)) #83 min
system.time(lggm<-read.delim('/net/shendure/vol7/stevesal/MSI_TCGA_project/LGG_1-22_merged_full_analysis.txt.gz',header=F)) #118 min



#gold standards available for coad, read, and upec; comparing vs. MSI+/-

colnames(ucec)<-as.character(as.matrix(column_names))
system.time(save(ucec,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/UCEC_1-3_merged_full_analysis.robj')) #2 min
colnames(read)<-as.character(as.matrix(column_names))
system.time(save(read,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/READ_1-3_merged_full_analysis.robj')) #1 min
colnames(lihc)<-as.character(as.matrix(column_names))
system.time(save(lihc,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/LIHC_1-8_merged_full_analysis.robj')) #3 min
colnames(coad)<-as.character(as.matrix(column_names))
system.time(save(coad,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/COAD_1-9_merged_full_analysis.robj')) #3 min
colnames(gbmm)<-as.character(as.matrix(column_names))
system.time(save(gbmm,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/GBMm_1-4_merged_full_analysis.robj')) #3 min
colnames(skcm)<-as.character(as.matrix(column_names))
system.time(save(skcm,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/SKCM_1-4_merged_full_analysis.robj')) #2 min
colnames(ucec2)<-as.character(as.matrix(column_names))
system.time(save(ucec2,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/UCEC_4-7_merged_full_analysis.robj')) #2 min
colnames(brca)<-as.character(as.matrix(column_names))
system.time(save(brca,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/BRCA_1-12_merged_full_analysis.robj')) #8 min
colnames(prad)<-as.character(as.matrix(column_names))
system.time(save(prad,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/PRAD_1-16_merged_full_analysis.robj')) #11 min
colnames(thca)<-as.character(as.matrix(column_names))
system.time(save(thca,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/THCA_1-17_merged_full_analysis.robj')) #12 min
colnames(luad)<-as.character(as.matrix(column_names))
system.time(save(luad,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/LUAD_1-18_merged_full_analysis.robj')) #12 min
colnames(lusc)<-as.character(as.matrix(column_names))
system.time(save(lusc,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/LUSC_1-14_merged_full_analysis.robj')) #9 min
colnames(blca)<-as.character(as.matrix(column_names))
system.time(save(blca,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/BLCA_1-10_merged_full_analysis.robj')) #7 min
colnames(stad)<-as.character(as.matrix(column_names))
system.time(save(stad,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/STAD_1-14_merged_full_analysis.robj')) #9 min
colnames(ovmu)<-as.character(as.matrix(column_names))
system.time(save(ovmu,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/OVMU_1-3_merged_full_analysis.robj')) #9 min
colnames(kirc)<-as.character(as.matrix(column_names))
system.time(save(kirc,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/KIRC_1-5_merged_full_analysis.robj')) #3 min
colnames(kirp)<-as.character(as.matrix(column_names))
system.time(save(kirp,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/KIRP_1-7_merged_full_analysis.robj')) #4 min
colnames(hnsc)<-as.character(as.matrix(column_names))
system.time(save(hnsc,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/HNSC_1-21_merged_full_analysis.robj')) #14 min
colnames(lggm)<-as.character(as.matrix(column_names))
system.time(save(lggm,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/LGGM_1-22_merged_full_analysis.robj')) #15 min

library(data.table)
ucecorig<-ucec
ucec<-rbindlist(list(ucec,ucec2))
system.time(save(ucec,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/UCEC_1-7_merged_full_analysis.robj')) #2 min
