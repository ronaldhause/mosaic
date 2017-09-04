require(qvalue)
require(caret)
require(doMC)
require(reshape2)
require(plyr)
require(scales)
require(data.table)

#batch 2, LIHC, GBM, and SKCM
#read in old data
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_gbmm_skcm_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')
fwrite(data2, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_gbmm_skcm_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.txt', sep='\t', quote=F)

#read in new data
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_extended_exomes.robj')
fwrite(lihcnew,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_extended_exomes.txt', sep='\t', quote=F)
gbmnew<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/gbm_extended_exomes.txt', sep='\t')
skcmnew<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/skcm_extended_exomes.txt', sep='\t')

#newdat<-rbindlist(list(lihcnew, gbmnew, skcmnew)) #230 million rows
#fwrite(newdat,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/new_lihc_gbm_skcm_merged_full_analysis_msi_values_converted_to_numeric_062316.txt', sep='\t', quote=F)
#newdat<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/new_lihc_gbm_skcm_merged_full_analysis_msi_values_converted_to_numeric_062316.txt', sep='\t')

#shell script for combining old and new data
cd /net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/
awk 'FNR==1 && NR!=1{next;}{print}' lihc_gbmm_skcm_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.txt lihc_extended_exomes.txt gbm_extended_exomes.txt skcm_extended_exomes.txt > lihc_gbm_skcm_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_063016.txt

#summary statistics
data<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_gbm_skcm_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_063016.txt', sep='\t', showProgress=T) #899 samples
data[,PEAK_DIFFERENCE_VALUE:=as.numeric(as.matrix(PEAK_DIFFERENCE_VALUE))]
data[,KS_VALUE:=as.numeric(as.matrix(KS_VALUE))]
#eliminate duplicates
#data<-unique(data,by=c("SAMPLE_NAME","LOCUS_COORDINATES"))

loc1<-subset(data, LOCUS_COORDINATES == "10:100008321-100008335")
sampleinfo1n<-loc1
sampleinfo1n$MLH1_silencing<-NA

system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #14s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #15s
system.time(num_unstable_ks<-data[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=SAMPLE_NAME]) #15s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=SAMPLE_NAME]) #10s
system.time(num_na<-data[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=SAMPLE_NAME]) #10s

sampledat1n<-data.frame(sample_name=peak_avg$SAMPLE_NAME, peak_avg=peak_avg$V1, peak_sd=peak_sd$V1, num_unstable=num_unstable$V1, num_called = 516876-num_na$V1)
sampledat1n$prop_unstable<-sampledat1n$num_unstable/sampledat1n$num_called
sampledat1n$defbsite<-rep('stable',dim(sampledat1n)[1])
sampledat1n$msi<-rep('MSS',dim(sampledat1n)[1])

emptysamplelocs<-which(sampledat1n$num_called==0)
emptysamples<-sampledat1n$sample_name[emptysamplelocs]
sampledat1n<-sampledat1n[-emptysamplelocs,]
sampleinfo1n<-sampleinfo1n[-emptysamplelocs,]

emptysamples<-as.character(emptysamples)
data<-data[!(SAMPLE_NAME %in% emptysamples)]
fwrite(data, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_gbm_skcm_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_duplicates_empty_samples_removed_063016.txt', sep='\t')

data<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_gbm_skcm_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_duplicates_empty_samples_removed_063016.txt', sep='\t', showProgress=T) 

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_classifier_063016.robj')
#do at the end for all samples
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mmr_mutational_data.robj')

defbsite<-data[LOCUS_COORDINATES == "8:7679723-7679741"]

unstable_discretize<-function(x){
	output<-rep('NA',length(x))
	output[x<=0]<-'stable'
	output[x>0]<-'unstable'
	return(output)
}

defbsite2<-unstable_discretize(defbsite$PEAK_DIFFERENCE_VALUE)
defbsite3<-defbsite2[match(sampledat1n$sample_name,defbsite$SAMPLE_NAME)]
sampledat1n$defbsite<-defbsite3

ucecsmall2<-sampledat1n
ucecsmall2$X8.7679723.7679741<-ucecsmall2$defbsite
sampledat1n$msi<-predict(mosaic, ucecsmall2, type="raw")

sampledat2n<-sampledat1n
sampleinfo2n<-sampleinfo1n
sampleinfo2n$TUMOR_TYPE<-revalue(sampleinfo2n$TUMOR_TYPE,c("GBMm"="GBM"))
sampleinfo2n$MSI_STATUS<-sampledat2n$msi
save(sampledat2n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_results_post_review_msi_predicted_070516.robj')
save(sampleinfo2n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_sampleinfo_post_review_msi_predicted_070516.robj')

data[, MSI_STATUS := sampledat2n$msi[match(data$SAMPLE_NAME,sampledat2n$sample_name)]]
data[, TUMOR_TYPE := sampleinfo2n$TUMOR_TYPE[match(data$SAMPLE_NAME,sampleinfo2n$SAMPLE_NAME)]]

fwrite(data, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lihc_gbm_skcm_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_duplicates_empty_samples_removed_msi_called_063016.txt', sep='\t')

#locus output
system.time(locus_output<-data[,j=list(median_peak_diff = as.double(median(PEAK_DIFFERENCE_VALUE,na.rm=T)), num_unstable = length(which(PEAK_DIFFERENCE_VALUE>0)), num_missing =  length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)), num_present = length(which(is.na(PEAK_DIFFERENCE_VALUE)==F))),by=list(TUMOR_TYPE,MSI_STATUS,LOCUS_COORDINATES)]) #500s

load('/net/shendure/vol12/projects/msi_tcga/tcga_locusinfo.robj')
colnames(locus_output)[1:3]<-c('tumor_type', 'msi_status', 'locus')
locus_output[, genomic_class := locusinfo$GENOMIC_CLASS[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]]
locus_output$gene<-locusinfo$GENE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_type<-locusinfo$REPEAT_TYPE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_dna_sequence<-locusinfo$REPEAT_DNA_SEQUENCE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
write.table(locus_output,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_tcga_loci_results_060616.txt',quote=F,row.names=F)

#batch 3, BRCA, PRAD, THCA
#read in old data
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/brca_prad_thca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')
fwrite(data3, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/brca_prad_thca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.txt', sep='\t', quote=F, showProgress=T)

#read in new data
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/prad_extended_exomes.robj')
fwrite(pradnew,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/prad_extended_exomes.txt', sep='\t', quote=F, verbose=T)
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/thca_extended_exomes.robj')
fwrite(thcanew,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/thca_extended_exomes.txt', sep='\t', quote=F, verbose=T)

#shell script for combining old and new data
cd /net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/
awk 'FNR==1 && NR!=1{next;}{print}' brca_prad_thca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.txt prad_extended_exomes.txt thca_extended_exomes.txt > brca_prad_thca_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_063016.txt

rm(list=ls())
gc()

#summary statistics
data<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/brca_prad_thca_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_063016.txt', sep='\t', showProgress=T) #1226 samples, 38 min
data[,PEAK_DIFFERENCE_VALUE:=as.numeric(as.matrix(PEAK_DIFFERENCE_VALUE))]
data[,KS_VALUE:=as.numeric(as.matrix(KS_VALUE))]
#eliminate 4 duplicates, 1222 samples

data<-unique(data,by=c("SAMPLE_NAME","LOCUS_COORDINATES"))

loc1<-subset(data, LOCUS_COORDINATES == "10:100008321-100008335")
sampleinfo1n<-loc1
sampleinfo1n$MLH1_silencing<-NA

system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #14s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #15s
system.time(num_unstable_ks<-data[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=SAMPLE_NAME]) #15s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=SAMPLE_NAME]) #10s
system.time(num_na<-data[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=SAMPLE_NAME]) #10s

sampledat1n<-data.frame(sample_name=peak_avg$SAMPLE_NAME, peak_avg=peak_avg$V1, peak_sd=peak_sd$V1, num_unstable=num_unstable$V1, num_called = 516876-num_na$V1)
sampledat1n$prop_unstable<-sampledat1n$num_unstable/sampledat1n$num_called
sampledat1n$defbsite<-rep('stable',dim(sampledat1n)[1])
sampledat1n$msi<-rep('MSS',dim(sampledat1n)[1])

emptysamplelocs<-which(sampledat1n$num_called==0) #8 samples
emptysamples<-sampledat1n$sample_name[emptysamplelocs]
sampledat1n<-sampledat1n[-emptysamplelocs,]
sampleinfo1n<-sampleinfo1n[-emptysamplelocs,]

emptysamples<-as.character(emptysamples)
data<-data[!(SAMPLE_NAME %in% emptysamples)] #1214 samples

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_classifier_063016.robj')

defbsite<-data[LOCUS_COORDINATES == "8:7679723-7679741"]

unstable_discretize<-function(x){
	output<-rep('NA',length(x))
	output[x<=0]<-'stable'
	output[x>0]<-'unstable'
	return(output)
}

defbsite2<-unstable_discretize(defbsite$PEAK_DIFFERENCE_VALUE)
defbsite3<-defbsite2[match(sampledat1n$sample_name,defbsite$SAMPLE_NAME)]
sampledat1n$defbsite<-defbsite3

ucecsmall2<-sampledat1n
ucecsmall2$X8.7679723.7679741<-ucecsmall2$defbsite
sampledat1n$msi<-predict(mosaic, ucecsmall2, type="raw")

sampledat3n<-sampledat1n
sampleinfo3n<-sampleinfo1n
sampleinfo3n$MSI_STATUS<-sampledat3n$msi
save(sampledat3n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_results_post_review_msi_predicted_070516.robj')
save(sampleinfo3n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_sampleinfo_post_review_msi_predicted_070516.robj')

data[, MSI_STATUS := sampledat3n$msi[match(data$SAMPLE_NAME,sampledat3n$sample_name)]]

fwrite(data, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/brca_prad_thca_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_duplicates_empty_samples_removed_msi_called_063016.txt', sep='\t')

#locus output
system.time(median_locus_peaks<-data[,median(PEAK_DIFFERENCE_VALUE,na.rm=T),by=list(TUMOR_TYPE,MSI_STATUS,LOCUS_COORDINATES)]) #130s
system.time(pos_locs<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=list(TUMOR_TYPE,MSI_STATUS,LOCUS_COORDINATES)]) #202s
system.time(locus_nas<-data[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=list(TUMOR_TYPE,MSI_STATUS,LOCUS_COORDINATES)]) #10s

system.time(locus_output<-data[,j=list(median_peak_diff = median(PEAK_DIFFERENCE_VALUE,na.rm=T), num_unstable = length(which(PEAK_DIFFERENCE_VALUE>0)), num_missing =  length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)), num_present = length(which(is.na(PEAK_DIFFERENCE_VALUE)==F))),by=list(TUMOR_TYPE,MSI_STATUS,LOCUS_COORDINATES)]) #500s

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_locusinfo.robj')
colnames(locus_output)[1:3]<-c('tumor_type', 'msi_status', 'locus')
locus_output[, genomic_class := locusinfo$GENOMIC_CLASS[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]]
locus_output$gene<-locusinfo$GENE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_type<-locusinfo$REPEAT_TYPE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_dna_sequence<-locusinfo$REPEAT_DNA_SEQUENCE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
write.table(locus_output,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_tcga_loci_results_060616.txt',quote=F,row.names=F)

#batch 4, LUAD, LUSC, BLCA
#read in old data
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/luad_lusc_blca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')
fwrite(data4, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/luad_lusc_blca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.txt', sep='\t', quote=F, verbose=T)

#read in new data
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/luad_extended_exomes.robj')
fwrite(luadnew,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/luad_extended_exomes.txt', sep='\t', quote=F, verbose=T)
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/blca_extended_exomes.robj')
fwrite(blcanew,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/blca_extended_exomes.txt', sep='\t', quote=F, verbose=T)

#shell script for combining old and new data
cd /net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/
awk 'FNR==1 && NR!=1{next;}{print}' luad_lusc_blca_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.txt luad_extended_exomes.txt lusc_extended_exomes.txt blca_extended_exomes.txt > luad_lusc_blca_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_063016.txt

rm(list=ls())
gc()

#summary statistics
data<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/luad_lusc_blca_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_063016.txt', sep='\t', showProgress=T) #1223 samples, 35 min
data[,PEAK_DIFFERENCE_VALUE:=as.numeric(as.matrix(PEAK_DIFFERENCE_VALUE))]
data[,KS_VALUE:=as.numeric(as.matrix(KS_VALUE))]

loc1<-subset(data, LOCUS_COORDINATES == "10:100008321-100008335")
dim(loc1)
dim(unique(loc1,by="SAMPLE_NAME"))
uniquesamples<-as.character(unique(loc1$SAMPLE_NAME))
#eliminate 12 duplicates, 1211 samples

system.time(data<-unique(data,by=c("SAMPLE_NAME","LOCUS_COORDINATES"))) #1214 samples, 220s

loc1<-subset(data, LOCUS_COORDINATES == "10:100008321-100008335")
sampleinfo1n<-loc1
sampleinfo1n$MLH1_silencing<-NA

system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #14s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #15s
system.time(num_unstable_ks<-data[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=SAMPLE_NAME]) #15s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=SAMPLE_NAME]) #10s
system.time(num_na<-data[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=SAMPLE_NAME]) #10s

sampledat1n<-data.frame(sample_name=peak_avg$SAMPLE_NAME, peak_avg=peak_avg$V1, peak_sd=peak_sd$V1, num_unstable=num_unstable$V1, num_called = 516876-num_na$V1)
sampledat1n$prop_unstable<-sampledat1n$num_unstable/sampledat1n$num_called
sampledat1n$defbsite<-rep('stable',dim(sampledat1n)[1])
sampledat1n$msi<-rep('MSS',dim(sampledat1n)[1])

emptysamplelocs<-which(sampledat1n$num_called==0) #35 samples
emptysamples<-sampledat1n$sample_name[emptysamplelocs]
sampledat1n<-sampledat1n[-emptysamplelocs,]
#sampleinfo1n<-sampleinfo1n[-emptysamplelocs,]

emptysamples<-as.character(emptysamples)
data<-data[!(SAMPLE_NAME %in% emptysamples)] #1188 samples

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_classifier_063016.robj')

defbsite<-data[LOCUS_COORDINATES == "8:7679723-7679741"]

unstable_discretize<-function(x){
	output<-rep('NA',length(x))
	output[x<=0]<-'stable'
	output[x>0]<-'unstable'
	return(output)
}

defbsite2<-unstable_discretize(defbsite$PEAK_DIFFERENCE_VALUE)
defbsite3<-defbsite2[match(sampledat1n$sample_name,defbsite$SAMPLE_NAME)]
sampledat1n$defbsite<-defbsite3

ucecsmall2<-sampledat1n
ucecsmall2$X8.7679723.7679741<-ucecsmall2$defbsite
sampledat1n$msi<-predict(mosaic, ucecsmall2, type="raw")

sampledat4n<-sampledat1n
sampleinfo4n<-sampleinfo1n
#sampleinfo2n$TUMOR_TYPE<-revalue(sampleinfo2n$TUMOR_TYPE,c("GBMm"="GBM"))
sampleinfo4n$MSI_STATUS<-sampledat4n$msi
save(sampledat4n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_results_post_review_msi_predicted_070516.robj')
save(sampleinfo4n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_sampleinfo_post_review_msi_predicted_070516.robj')

data[, MSI_STATUS := sampledat4n$msi[match(data$SAMPLE_NAME,sampledat4n$sample_name)]]
#data[, TUMOR_TYPE := sampleinfo2n$TUMOR_TYPE[match(data$SAMPLE_NAME,sampleinfo2n$sample_name)]]

fwrite(data, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/luad_lusc_blca_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_duplicates_empty_samples_removed_msi_called_063016.txt', sep='\t')

#locus output
system.time(locus_output<-data[,j=list(median_peak_diff = median(PEAK_DIFFERENCE_VALUE,na.rm=T), num_unstable = length(which(PEAK_DIFFERENCE_VALUE>0)), num_missing =  length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)), num_present = length(which(is.na(PEAK_DIFFERENCE_VALUE)==F))),by=list(TUMOR_TYPE,MSI_STATUS,LOCUS_COORDINATES)]) #500s

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_locusinfo.robj')
colnames(locus_output)[1:3]<-c('tumor_type', 'msi_status', 'locus')
locus_output[, genomic_class := locusinfo$GENOMIC_CLASS[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]]
locus_output$gene<-locusinfo$GENE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_type<-locusinfo$REPEAT_TYPE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_dna_sequence<-locusinfo$REPEAT_DNA_SEQUENCE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
write.table(locus_output,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_tcga_loci_results_060616.txt',quote=F,row.names=F)

#batch 5, OV, KIRC, KIRP
#read in old data
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/stad_ovmu_kirc_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')
data5<-subset(data5, TUMOR_TYPE != "STAD")
fwrite(data5, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/ov_kirc_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.txt', sep='\t', quote=F, verbose=T)

load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/kirp_hnsc_lggm_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')
kirp<-subset(data6, TUMOR_TYPE == "KIRP")
fwrite(kirp, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/kirp_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.txt', sep='\t', quote=F, verbose=T)

kirp<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/kirp_extended_exomes.txt', sep='\t', showProgress=T) #691 samples, 35 min

#read in new data
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/ov_extended_exomes.robj')
fwrite(ovnew,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/ov_extended_exomes.txt', sep='\t', quote=F, verbose=T)

#shell script for combining old and new data
cd /net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/
awk 'FNR==1 && NR!=1{next;}{print}' ov_kirc_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.txt kirp_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.txt ov_extended_exomes.txt kirc_extended_exomes.txt kirp_extended_exomes.txt > ov_kirc_kirp_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_063016.txt

rm(list=ls())
gc()

#summary statistics
data<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/ov_kirc_kirp_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_063016.txt', sep='\t', showProgress=T) #691 samples, 35 min
data[,PEAK_DIFFERENCE_VALUE:=as.numeric(as.matrix(PEAK_DIFFERENCE_VALUE))]
data[,KS_VALUE:=as.numeric(as.matrix(KS_VALUE))]

loc1<-subset(data, LOCUS_COORDINATES == "10:100008321-100008335")
dim(loc1)
dim(unique(loc1,by="SAMPLE_NAME"))
#eliminate 59 duplicates, 632 samples

system.time(data<-unique(data,by=c("SAMPLE_NAME","LOCUS_COORDINATES"))) #632 samples, 110s

loc1<-subset(data, LOCUS_COORDINATES == "10:100008321-100008335")
sampleinfo1n<-loc1
sampleinfo1n$MLH1_silencing<-NA

system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #14s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #15s
system.time(num_unstable_ks<-data[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=SAMPLE_NAME]) #15s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=SAMPLE_NAME]) #10s
system.time(num_na<-data[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=SAMPLE_NAME]) #10s

sampledat1n<-data.frame(sample_name=peak_avg$SAMPLE_NAME, peak_avg=peak_avg$V1, peak_sd=peak_sd$V1, num_unstable=num_unstable$V1, num_called = 516876-num_na$V1)
sampledat1n$prop_unstable<-sampledat1n$num_unstable/sampledat1n$num_called
sampledat1n$defbsite<-rep('stable',dim(sampledat1n)[1])
sampledat1n$msi<-rep('MSS',dim(sampledat1n)[1])

emptysamplelocs<-which(sampledat1n$num_called==0) #83 samples
emptysamples<-sampledat1n$sample_name[emptysamplelocs]
sampledat1n<-sampledat1n[-emptysamplelocs,]
sampleinfo1n<-sampleinfo1n[-emptysamplelocs,]

emptysamples<-as.character(emptysamples)
data<-data[!(SAMPLE_NAME %in% emptysamples)] #549 samples

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_classifier_063016.robj')

defbsite<-data[LOCUS_COORDINATES == "8:7679723-7679741"]

unstable_discretize<-function(x){
	output<-rep('NA',length(x))
	output[x<=0]<-'stable'
	output[x>0]<-'unstable'
	return(output)
}

defbsite2<-unstable_discretize(defbsite$PEAK_DIFFERENCE_VALUE)
defbsite3<-defbsite2[match(sampledat1n$sample_name,defbsite$SAMPLE_NAME)]
sampledat1n$defbsite<-defbsite3

ucecsmall2<-sampledat1n
ucecsmall2$X8.7679723.7679741<-ucecsmall2$defbsite
sampledat1n$msi<-predict(mosaic, ucecsmall2, type="raw")

#assessing PPV of mosaic on training data
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_results_post_review_msi_predicted_063016.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_post_review_msi_updated_062916.robj')
a<-sampledat1n[sampleinfo1n$MSI_STATUS=='MSI-H' | sampleinfo1n$MSI_STATUS=='MSS' | sampleinfo1n$MSI_STATUS=='MSI-L',]
b<-sampleinfo1n[sampleinfo1n$MSI_STATUS=='MSI-H' | sampleinfo1n$MSI_STATUS=='MSS' | sampleinfo1n$MSI_STATUS=='MSI-L',]
b$MSI_STATUS[b$MSI_STATUS=='MSI-L']<-'MSS'
confusionMatrix(a$msi,b$MSI_STATUS)

sampledat5n<-sampledat1n
sampleinfo5n<-sampleinfo1n
sampleinfo5n$TUMOR_TYPE<-revalue(sampleinfo5n$TUMOR_TYPE,c("OVmu"="OV"))
sampleinfo5n$MSI_STATUS<-sampledat5n$msi
save(sampledat5n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch5_results_post_review_msi_predicted_070516.robj')
save(sampleinfo5n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch5_sampleinfo_post_review_msi_predicted_070516.robj')

data[, MSI_STATUS := sampledat5n$msi[match(data$SAMPLE_NAME,sampledat5n$sample_name)]]
data[, TUMOR_TYPE := sampleinfo5n$TUMOR_TYPE[match(data$SAMPLE_NAME,sampleinfo5n$SAMPLE_NAME)]]

data<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/ov_kirc_kirp_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_duplicates_empty_samples_removed_msi_called_063016.txt', sep='\t', showProgress=T) #691 samples, 35 min

fwrite(data, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/ov_kirc_kirp_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_duplicates_empty_samples_removed_msi_called_063016.txt', sep='\t')

#locus output
system.time(locus_output<-data[,j=list(median_peak_diff = as.double(median(PEAK_DIFFERENCE_VALUE,na.rm=T)), num_unstable = length(which(PEAK_DIFFERENCE_VALUE>0)), num_missing =  length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)), num_present = length(which(is.na(PEAK_DIFFERENCE_VALUE)==F))),by=list(TUMOR_TYPE,MSI_STATUS,LOCUS_COORDINATES)]) #500s

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_locusinfo.robj')
colnames(locus_output)[1:3]<-c('tumor_type', 'msi_status', 'locus')
locus_output[, genomic_class := locusinfo$GENOMIC_CLASS[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]]
locus_output$gene<-locusinfo$GENE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_type<-locusinfo$REPEAT_TYPE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_dna_sequence<-locusinfo$REPEAT_DNA_SEQUENCE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
write.table(locus_output,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch5_tcga_loci_results_060616.txt',quote=F,row.names=F)

#batch 6, HNSC and LGG
#read in old data
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/kirp_hnsc_lggm_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.robj')
data6<-subset(data6, TUMOR_TYPE != "KIRP")
fwrite(data6, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/hnsc_lgg_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.txt', sep='\t', quote=F, verbose=T)

#read in new data
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lgg_extended_exomes.robj')
fwrite(lggnew,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/lgg_extended_exomes.txt', sep='\t', quote=F, verbose=T)

load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/hnsc_extended_exomes.robj')
fwrite(hnscnew,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/hnsc_extended_exomes.txt', sep='\t', quote=F, verbose=T)

#shell script for combining old and new data
cd /net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/
awk 'FNR==1 && NR!=1{next;}{print}' hnsc_lgg_merged_analysis_merged_full_analysis_msi_values_converted_to_numeric.txt lgg_extended_exomes.txt hnsc_extended_exomes.txt > hnsc_lgg_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_063016.txt

rm(list=ls())
gc()

#summary statistics
data<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/hnsc_lgg_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_063016.txt', sep='\t', showProgress=T) #1019 samples, 30 min
data[,PEAK_DIFFERENCE_VALUE:=as.numeric(as.matrix(PEAK_DIFFERENCE_VALUE))]
data[,KS_VALUE:=as.numeric(as.matrix(KS_VALUE))]

loc1<-subset(data, LOCUS_COORDINATES == "10:100008321-100008335")
dim(loc1)
dim(unique(loc1,by="SAMPLE_NAME"))
#eliminate 0 duplicates, 1019 samples

sampleinfo1n<-loc1
sampleinfo1n$MLH1_silencing<-NA

system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #14s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #15s
system.time(num_unstable_ks<-data[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=SAMPLE_NAME]) #15s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=SAMPLE_NAME]) #10s
system.time(num_na<-data[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=SAMPLE_NAME]) #10s

sampledat1n<-data.frame(sample_name=peak_avg$SAMPLE_NAME, peak_avg=peak_avg$V1, peak_sd=peak_sd$V1, num_unstable=num_unstable$V1, num_called = 516876-num_na$V1)
sampledat1n$prop_unstable<-sampledat1n$num_unstable/sampledat1n$num_called
sampledat1n$defbsite<-rep('stable',dim(sampledat1n)[1])
sampledat1n$msi<-rep('MSS',dim(sampledat1n)[1])

emptysamplelocs<-which(sampledat1n$num_called==0) #0 samples
emptysamples<-sampledat1n$sample_name[emptysamplelocs]
sampledat1n<-sampledat1n[-emptysamplelocs,]
sampleinfo1n<-sampleinfo1n[-emptysamplelocs,]

emptysamples<-as.character(emptysamples)
data<-data[!(SAMPLE_NAME %in% emptysamples)] #549 samples

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_classifier_063016.robj')

defbsite<-data[LOCUS_COORDINATES == "8:7679723-7679741"]

unstable_discretize<-function(x){
	output<-rep('NA',length(x))
	output[x<=0]<-'stable'
	output[x>0]<-'unstable'
	return(output)
}

defbsite2<-unstable_discretize(defbsite$PEAK_DIFFERENCE_VALUE)
defbsite3<-defbsite2[match(sampledat1n$sample_name,defbsite$SAMPLE_NAME)]
sampledat1n$defbsite<-defbsite3

ucecsmall2<-sampledat1n
ucecsmall2$X8.7679723.7679741<-ucecsmall2$defbsite
sampledat1n$msi<-predict(mosaic, ucecsmall2, type="raw")

sampledat6n<-sampledat1n
sampleinfo6n<-sampleinfo1n
sampleinfo6n$TUMOR_TYPE<-revalue(sampleinfo6n$TUMOR_TYPE,c("LGGm"="LGG"))
sampleinfo6n$MSI_STATUS<-sampledat6n$msi
save(sampledat6n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch6_results_post_review_msi_predicted_070516.robj')
save(sampleinfo6n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch6_sampleinfo_post_review_msi_predicted_070516.robj')

data[, MSI_STATUS := sampledat6n$msi[match(data$SAMPLE_NAME,sampledat6n$sample_name)]]
data[, TUMOR_TYPE := sampleinfo6n$TUMOR_TYPE[match(data$SAMPLE_NAME,sampleinfo6n$SAMPLE_NAME)]]

fwrite(data, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/hnsc_lgg_merged_full_analysis_msi_values_converted_to_numeric_post_review_exomes_added_duplicates_empty_samples_removed_msi_called_063016.txt', sep='\t')

#locus output
system.time(locus_output<-data[,j=list(median_peak_diff = median(PEAK_DIFFERENCE_VALUE,na.rm=T), num_unstable = length(which(PEAK_DIFFERENCE_VALUE>0)), num_missing =  length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)), num_present = length(which(is.na(PEAK_DIFFERENCE_VALUE)==F))),by=list(TUMOR_TYPE,MSI_STATUS,LOCUS_COORDINATES)]) #500s

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_locusinfo.robj')
colnames(locus_output)[1:3]<-c('tumor_type', 'msi_status', 'locus')
locus_output[, genomic_class := locusinfo$GENOMIC_CLASS[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]]
locus_output$gene<-locusinfo$GENE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_type<-locusinfo$REPEAT_TYPE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_dna_sequence<-locusinfo$REPEAT_DNA_SEQUENCE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
write.table(locus_output,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch6_tcga_loci_results_060616.txt',quote=F,row.names=F)
