load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
#load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
#load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_locusinfo.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/finalsampledata_062216.robj')

column_names<-as.character(as.matrix(read.table('/net/shendure/vol7/stevesal/MSI_TCGA_project/TCGA_Analysis_header.txt',header=F,sep='\t')))

#for combining all individuals
#cd OV_tumor_normal_comparisons
#awk 'NR==1 {print; next} FNR==1 {next} {print FILENAME " \t " $0}' *.txt > all.txt &

source('~/helper_functions.R')
require(qvalue)
require(caret)
require(doMC)
require(reshape2)
require(plyr)
require(scales)

#load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated.robj')

#setkey(data,SAMPLE_NAME)
#system.time(data['TCGA-B5-A1MU',TUMOR_TYPE := 'UCEC'])
#system.time(data['TCGA-FI-A2EW',TUMOR_TYPE := 'UCEC'])
#system.time(data['TCGA-B5-A1MU',MSI_STATUS := 'MSS'])
#system.time(data['TCGA-FI-A2EW',MSI_STATUS := 'ND'])
#data<-droplevels.data.table(data)

install.packages("data.table", type = "source",
    repos = "http://Rdatatable.github.io/data.table")

#save(data, file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated_062316.robj')

#load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated_062316.robj')
#fwrite(data,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated_062316.txt', sep='\t', quote=F)

#load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/ucec_extended_exomes.robj')
#fwrite(ucecnew,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/ucec_extended_exomes.txt', sep='\t', quote=F)
#load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/coad_extended_exomes.robj')
#fwrite(coadnew,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/coad_extended_exomes.txt', sep='\t', quote=F)
#system.time(load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/read_extended_exomes.robj')) #36s, ~100MB
#system.time(fwrite(readnew,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/read_extended_exomes.txt', sep='\t', quote=F)) #15s
#load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/stad_extended_exomes.robj')
#fwrite(stadnew,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/stad_extended_exomes.txt', sep='\t', quote=F)

load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/STAD_1-14_merged_full_analysis.robj')
system.time(fwrite(stad,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/STAD_1-14_merged_full_analysis.txt', sep='\t', quote=F)) 

#system.time(save(readnew, file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/read_extended_exomes_test.robj')) #76s
#system.time(fwrite(readnew,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/read_extended_exomes.txt', sep='\t', quote=F)) #15s
#system.time(readnew2<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/read_extended_exomes.txt', sep='\t')) #21s, ~800 MB

data<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated_062316.txt', sep='\t')

ucecnew<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/ucec_extended_exomes.txt', sep='\t')
coadnew<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/coad_extended_exomes.txt', sep='\t')
readnew<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/read_extended_exomes.txt', sep='\t')
stadnew<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/stad_extended_exomes.txt', sep='\t')

#newdat<-rbindlist(list(ucecnew, coadnew, readnew, stadnew)) #350 million rows
#system.time(fwrite(newdat,'/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/new_TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated_062316.txt', sep='\t', quote=F))
newdat<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/new_TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated_062316.txt', sep='\t')

install.packages("data.table", type = "source", repos = "http://Rdatatable.github.io/data.table")

cd /net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/
awk 'FNR==1 && NR!=1{next;}{print}' TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated_062316.txt STAD_1-14_merged_full_analysis.txt ucec_extended_exomes.txt coad_extended_exomes.txt read_extended_exomes.txt stad_extended_exomes.txt > TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated_post_review_exomes_added_062316.txt &

data2<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_msi_values_converted_to_numeric_duplicates_eliminated_post_review_exomes_added_062416.txt', sep='\t', showProgress=T) #43 min, 38 min after removing duplicates

#data2<-rbindlist(list(data, stad, ucecnew, coadnew, readnew, stadnew))
data2[,PEAK_DIFFERENCE_VALUE:=as.numeric(as.matrix(PEAK_DIFFERENCE_VALUE))]
data2[,KS_VALUE:=as.numeric(as.matrix(KS_VALUE))]
data2<-unique(data2,by=c("SAMPLE_NAME","LOCUS_COORDINATES"))
fwrite(data2, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_msi_values_converted_to_numeric_duplicates_eliminated_post_review_exomes_added_062416.txt', sep='\t')
#went down from 1404 to 1269 samples

loc1<-subset(data2, LOCUS_COORDINATES == "10:100008321-100008335")
sampleinfo1n<-loc1
save(sampleinfo1n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_post_review_062416.robj')

data<-data2
system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #14s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #15s
system.time(num_unstable_ks<-data[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=SAMPLE_NAME]) #15s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=SAMPLE_NAME]) #10s
system.time(num_na<-data[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=SAMPLE_NAME]) #10s

sampledat1n<-data.frame(sample_name=peak_avg$SAMPLE_NAME, peak_avg=peak_avg$V1, peak_sd=peak_sd$V1, num_unstable_ks=num_unstable_ks$V1,num_unstable=num_unstable$V1, num_called = 516876-num_na$V1)
sampledat1n$prop_unstable<-sampledat1n$num_unstable/sampledat1n$num_called

emptysamplelocs<-which(sampledat1n$num_called==0)
emptysamples<-sampledat1n$sample_name[emptysamplelocs]
sampledat1n<-sampledat1n[-emptysamplelocs,]
save(sampledat1n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_results_post_review_msi_updated_062416.robj')

save(emptysamples,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/emptysamples_062416.robj')
emptysamples<-as.character(emptysamples)
data3<-data2[!(SAMPLE_NAME %in% emptysamples)]
fwrite(data3, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_msi_values_converted_to_numeric_duplicates_eliminated_post_review_exomes_added_empty_samples_removed_062416.txt', sep='\t')

data2<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_msi_values_converted_to_numeric_duplicates_eliminated_post_review_exomes_added_empty_samples_removed_062416.txt', sep='\t', showProgress=T) #32 min

coadread_msi_status<-read.delim('/net/shendure/vol7/stevesal/MSI_TCGA_project/coadread_tcga_pub_clinical_data.tsv', header=T, stringsAsFactors=F)
stad_msi_status<-read.delim('/net/shendure/vol7/stevesal/MSI_TCGA_project/stad_tcga_pub_clinical_data.tsv', header=T, stringsAsFactors=F)
ucec_msi_status<-read.delim('/net/shendure/vol7/stevesal/MSI_TCGA_project/ucec_tcga_pub_clinical_data.tsv', header=T, stringsAsFactors=F)
colnames(ucec_msi_status)[14]<-'MSI.Status'
msiclindat<-Reduce(function(...) merge(..., all=T), list(coadread_msi_status,stad_msi_status,ucec_msi_status))

sampleinfo1n<-sampleinfo1n[-emptysamplelocs,]
newmsi<-msiclindat$MSI.Status[match(sampleinfo1n$SAMPLE_NAME,msiclindat$Patient.ID)]
sampleinfo1n$MSI_STATUS[which(is.na(newmsi)==F)]<-newmsi[which(is.na(newmsi)==F)]
sampleinfo1n$MSI_STATUS[sampleinfo1n$MSI_STATUS=='']<-'ND'
sampleinfo1n$MSI_STATUS[sampleinfo1n$MSI_STATUS=='Indeterminant']<-'ND'
sampleinfo1n$MSI_STATUS[is.na(sampleinfo1n$MSI_STATUS)]<-'ND'
sampleinfo1n$MLH1_silencing<-msiclindat$MLH1.Silencing[match(sampleinfo1n$SAMPLE_NAME,msiclindat$Patient.ID)]
save(sampleinfo1n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_post_review_msi_updated_062416.robj')

sampledat1n$msi<-sampleinfo1n$MSI_STATUS
sampledat1n$tumor_type<-sampleinfo1n$TUMOR_TYPE
training<-sampledat1n[which(sampledat1n$msi=='MSI-H' | sampledat1n$msi=='MSS'),]
msil<-sampledat1n[which(sampledat1n$msi=='MSI-L'),]
msil$msi<-rep('MSS',dim(msil)[1])
mlh10<-sampledat1n[which(sampleinfo1n$MLH1_silencing==0),]
mlh11<-sampledat1n[which(sampleinfo1n$MLH1_silencing==1),]
mss<-sampledat1n[which(sampledat1n$msi=='MSS'),]
msih<-sampledat1n[which(sampledat1n$msi=='MSI-H'),]

require(ggplot2)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_msi_subset_prop_unstable_per_cancer.pdf',width=7*1.75,height=7)
ggplot(aes(y=prop_unstable, x=tumor_type, fill=msi), data=training) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_msi_subset_num_unstable_per_cancer.pdf',width=7*1.75,height=7)
ggplot(aes(y=num_unstable, x=tumor_type, fill=msi), data=training) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_msi_subset_peak_average_per_cancer.pdf',width=7*1.75,height=7)
ggplot(aes(y=peak_avg, x=tumor_type, fill=msi), data=training) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + geom_hline(aes(yintercept=0.0053), linetype='dashed')
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_msi_subset_peak_sd_per_cancer.pdf',width=7*1.75,height=7)
ggplot(aes(y=peak_sd, x=tumor_type, fill=msi), data=training) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + geom_hline(aes(yintercept=0.0053), linetype='dashed')
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_msi_subset_peak_average_per_cancer_final.pdf',width=7*1.75,height=7)
ggplot(aes(y=peak_avg, x=tumor_type, fill=msi), data=training) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + geom_hline(aes(yintercept=0.0053), linetype='dashed') + geom_jitter(data=msil,aes(y=msil$peak_avg,x=msil$tumor_type,color='a')) + geom_jitter(data=mlh10,aes(y=mlh10$peak_avg,x=mlh10$tumor_type,color='b')) + geom_jitter(data=mlh11,aes(y=mlh11$peak_avg,x=mlh11$tumor_type,color='c'))
dev.off()

sampledat1ntemp<-sampledat1n

#looking at old model
sampledat1n$old_msi_pred<-rep('MSS',dim(sampledat1n)[1])
sampledat1n$old_msi_pred[which(sampledat1n$peak_avg>=0.0053)]<-'MSI-H'
sampledat1n2<-sampledat1n
sampledat1n2$peak_avg<-ave(sampledat1n2$peak_avg,sampledat1n2$tumor_type,FUN=function(x)x-median(x,na.rm=T))
sampledat1n$old_msi_pred2<-predict(rpartfit2,sampledat1n2,type="raw")

sampledat1n$defbsite<-siglocs2[,93]

#weighted
sampledat1n$new_msi_pred<-rep('MSS',dim(sampledat1n)[1])
sampledat1n$new_msi_pred[which(sampledat1n$peak_avg>=0.0055)]<-'MSI-H'
sampledat1n$new_msi_pred[Reduce(intersect, list(which(sampledat1n$peak_avg<0.0055), which(sampledat1n$peak_avg>0.0029), which(sampledat1n$defbsite=='unstable')))]<-'MSI-H'

#unweighted
sampledat1n$new_msi_predw<-rep('MSS',dim(sampledat1n)[1])
sampledat1n$new_msi_predw[which(sampledat1n$peak_avg>=0.0057)]<-'MSI-H'
sampledat1n$new_msi_predw[Reduce(intersect, list(which(sampledat1n$peak_avg<0.0057), which(sampledat1n$peak_avg>0.0045), which(sampledat1n$defbsite=='unstable')))]<-'MSI-H'

sampledat1n$corrected_msi<-sampledat1n$msi
sampledat1n$corrected_msi[sampledat1n$sample_name=='TCGA-B5-A11E']<-'MSI-H'
sampledat1n$corrected_msi[sampledat1n$sample_name=='TCGA-A5-A1OF']<-'MSS'
sampledat1n$corrected_msi[sampledat1n$sample_name=='TCGA-AP-AOLS']<-'MSS'
sampledat1n$corrected_msi[sampledat1n$sample_name=='TCGA-AP-A1DK']<-'MSS'
sampledat1n$corrected_msi[sampledat1n$sample_name=='TCGA-AX-A1C9']<-'MSS'
sampledat1n$corrected_msi[sampledat1n$sample_name=='TCGA-AX-A2HC']<-'MSS'
sampledat1n$corrected_msi[sampledat1n$sample_name=='TCGA-EC-A24G']<-'MSS'

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_training_data_030315.robj') #old ucecsmall2

sampledat1nold<-sampledat1n[match(ucecsmall2$peak_avg,sampledat1n$peak_avg),]
sampledat1nnew<-sampledat1n[-match(ucecsmall2$peak_avg,sampledat1n$peak_avg),]

msisubnew<-sampledat1nnew[sampledat1nnew$msi=='MSI-H' | sampledat1nnew$msi=='MSS',] #independent cohort of 99 UCEC and 118 STAD samples
msisubold<-sampledat1nold[sampledat1nold$msi=='MSI-H' | sampledat1nold$msi=='MSS',]
msisub<-sampledat1n[sampledat1n$msi=='MSI-H' | sampledat1n$msi=='MSS',]

confusionMatrix(msisubold$old_msi_pred2,msisubold$msi) #97.4% accuracy, 95% sensitivity, 98.2% specificity, scaled peak avg
confusionMatrix(msisubold$old_msi_pred,msisubold$msi) #97.1% accuracy, 95% sensitivity, 97.8% specificity, straight threshold
confusionMatrix(msisubold$new_msi_pred,msisubold$msi) #97.1% accuracy, 95% sensitivity, 97.8% specificity, straight threshold


confusionMatrix(msisubnew$old_msi_pred2,msisubnew$msi) #93.1% accuracy, 81.0% sensitivity, 98.9% specificity, scaled peak avg
confusionMatrix(msisubnew$old_msi_pred,msisubnew$msi) #93.5% accuracy, 83.3% sensitivity, 98.3% specificity, straight threshold

msisubnew<-sampledat1nnew[sampledat1nnew$msi=='MSI-H' | sampledat1nnew$msi=='MSS' | sampledat1nnew$msi=='MSI-L',] #independent cohort of 65 COAD, 20 READ, 135 STAD and 113 UCEC samples
msisubnew$msi[msisubnew$msi=='MSI-L']<-'MSS'

confusionMatrix(msisubnew$old_msi_pred2,msisubnew$msi) #94.3% accuracy, 81.0% sensitivity, 98.8% specificity, scaled peak avg
confusionMatrix(msisubnew$old_msi_pred,msisubnew$msi) #94.6% accuracy, 83.3% sensitivity, 98.4% specificity, straight threshold

confusionMatrix(msisubnew$old_msi_pred,msisubnew$corrected_msi) #94.6% accuracy, 83.3% sensitivity, 98.4% specificity, straight threshold; after correcting MSI-PCR calls, 95.8% accuracy, 88.8% sensitivity, 98.9% specificity

confusionMatrix(msisubnew$old_msi_pred[which(msisubnew$tumor=='STAD')],msisubnew$msi[which(msisubnew$tumor=='STAD')]) #99.3% accuracy, 96.9% sensitivity, 100% specificity, straight threshold
confusionMatrix(msisubnew$new_msi_pred[which(msisubnew$tumor=='STAD')],msisubnew$msi[which(msisubnew$tumor=='STAD')]) #100% accuracy, 100% specificity, 100% sensitivity


confusionMatrix(msisubnew$old_msi_pred[which(msisubnew$tumor=='COAD')],msisubnew$msi[which(msisubnew$tumor=='COAD')]) #98.5% accuracy, 100% sensitivity, 98.3% specificity, straight threshold
confusionMatrix(msisubnew$new_msi_pred[which(msisubnew$tumor=='COAD')],msisubnew$msi[which(msisubnew$tumor=='COAD')]) #100% accuracy, 100% specificity, 100% sensitivity


confusionMatrix(msisubnew$old_msi_pred[which(msisubnew$tumor=='UCEC')],msisubnew$msi[which(msisubnew$tumor=='UCEC')]) #85.8% accuracy, 70.5% sensitivity, 95.7% specificity, straight threshold
confusionMatrix(msisubnew$new_msi_pred[which(msisubnew$tumor=='UCEC')],msisubnew$msi[which(msisubnew$tumor=='UCEC')]) #85.8% accuracy, 70.5% sensitivity, 95.7% specificity, straight threshold; #90.9% accuracy, 84.1% sensitivity, 96.4% specificity


confusionMatrix(msisubnew$old_msi_pred[which(msisubnew$tumor=='UCEC')],msisubnew$corrected_msi[which(msisubnew$tumor=='UCEC')]) #85.8% accuracy, 70.5% sensitivity, 95.7% specificity, straight threshold; after correcting MSI-PCR calls, 90.9% accuracy, 80.0% sensitivity, 98.3% specificity

msiall<-sampledat1n[sampledat1n$msi=='MSI-H' | sampledat1n$msi=='MSS' | sampledat1n$msi=='MSI-L',] #independent cohort of 65 COAD, 20 READ, 135 STAD and 113 UCEC samples
msiall$msi[msiall$msi=='MSI-L']<-'MSS'
msiall$corrected_msi[msiall$corrected_msi=='MSI-L']<-'MSS'
confusionMatrix(msiall$old_msi_pred,msiall$msi) #95.9% accuracy, 89.0% sensitivity, 98.1% specificity, straight threshold; after correcting MSI-PCR calls, 96.8% accuracy, 91.9% sensitivity, 98.3% specificity

confusionMatrix(msiall$new_msi_pred,msiall$msi) #97.3% accuracy, 93.6% sensitivity, 98.5% specificity, straight thresholds (two)
confusionMatrix(msiall$new_msi_predw,msiall$msi) #97.3% accuracy, 93.6% sensitivity, 98.5% specificity, straight thresholds (two)

misclass<-msiall[which(msiall$new_msi_pred!=msiall$msi),]
misclassoutput<-cbind(misclass,sampleinfo1n[match(misclass$sample_name,sampleinfo1n$SAMPLE_NAME),c(9:18,21), with=FALSE])
write.csv(misclassoutput,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_misclassifications_new_model_all_samples.csv',row.names=F)

misclass_samples<-misclass$sample_name
misclass_clindat<-msiclindat[match(misclass_samples,msiclindat$Patient.ID),]
misclassoutput2<-cbind(misclass,sampleinfo1n[match(misclass$sample_name,sampleinfo1n$SAMPLE_NAME),c(9:18,21), with=FALSE],misclass_clindat)
write.csv(misclassoutput2,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_misclassifications_old_model_new_samples_with_clindat.csv',row.names=F)

allmisclass<-msisub[which(msisub$old_msi_pred!=msisub$msi),]

#training a new mosaic, old model said peak_avg >= 0.0053 was MSI-H
require(caret)
trainctrl <- trainControl(summaryFunction=twoClassSummary,classProb=TRUE)
rfFuncs$summary<-defaultSummary
ctrl <- rfeControl(functions=rfFuncs,method = "loocv")
ctrl <- rfeControl(functions=rfFuncs,method = "cv", number = 10)
set.seed(101)

training$msi<-factor(training$msi)

tdat<-training[,-c(1,6)]

results <- rfe(tdat[,-c(6:7)], as.factor(tdat$msi), sizes=seq(1:5), rfeControl=ctrl,trControl=trainctrl, metric="Accuracy")

require(party)
require(partykit)
results <- ctree(tdat$msi ~ peak_avg + prop_unstable, data=tdat, control = ctree_control(minsplit=6, maxdepth=3))

controls = ctree_control(minsplit=6)

weight_test<-training$msi
levels(weight_test)<-rev(table(weight_test))
weight_test<-as.numeric(as.matrix(weight_test))

fitControl <- trainControl(## 10-fold CV
  method = "loocv",
  returnData = TRUE,
  returnResamp = "final",
  repeats =  1)

gbmGrid <-  expand.grid(mtry=2:3)

set.seed(101)
system.time(msifitrfwn <- train(msi ~ ., data = tdat,
                   method = "rf",
                   trControl = fitControl,
                   verbose = FALSE,
                   prox=TRUE,
                   weights=weight_test,
                   n.trees=1000,
                   allowParallel=TRUE))

fitControl <- trainControl(## 10-fold CV
method = "cv",
returnData = TRUE,
returnResamp = "final",
repeats =  1)

system.time(rpartfit1 <- train(msi ~ peak_avg, data = tdat,
                   method = "rpart",
                   trControl = trainControl(method="loocv"),
                   control=rpart.control(minsplit=6),
                   tuneGrid=expand.grid(.cp=c(0,0.001,0.01,0.1,0.45,0.95)), 
                   weights=weight_test
                   ))

gbmGrid <-  expand.grid(mtry=2:3)

require(doMC)
registerDoMC(cores = 4)
  system.time(msirffit <- train(msi ~ peak_avg + prop_unstable, data = tdat,
                   method = "rf",
                   trControl = fitControl,
                   verbose = FALSE,
                   prox=TRUE,
                   #tuneGrid=gbmGrid,
                   n.trees=100,
                   allowParallel=TRUE))

library(rattle)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_classification_tree_weighted_fancy_post_review_10xcv_062816.pdf',width=7*1.5,height=7*1.5)
fancyRpartPlot(rpartfit1$finalModel)
dev.off()

library(rattle)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_classification_tree_weighted_fancy_post_review_10xcv_just_peak_avg_prop_unstable_062716.pdf',width=7*1.5,height=7*1.5)
fancyRpartPlot(rpartfit1_peak_avg$finalModel)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_classification_tree_weighted_post_review_062716.pdf',width=7*1.5,height=7*1.5)
plot(rpartfit1$finalModel, main="Classification Tree for MSI")
text(rpartfit1$finalModel, use.n=TRUE, all=TRUE, cex=.8)
dev.off()

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_post_review_msi_updated_062916.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_results_post_review_msi_updated_062916.robj')

data3[,MSI_STATUS:=sampleinfo1n$MSI_STATUS[match(data3$SAMPLE_NAME,sampleinfo1n$SAMPLE_NAME)]]

#looking for specific sites that can better differentiate MSI-H and MSS tumors
msidata<-data3[MSI_STATUS == "MSI-H"]
mssdata<-data3[MSI_STATUS == "MSS"]
msi_mean_locus_peaks<-msidata[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=LOCUS_COORDINATES]
mss_mean_locus_peaks<-mssdata[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=LOCUS_COORDINATES]
msi_locus_nas<-msidata[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=LOCUS_COORDINATES]
mss_locus_nas<-mssdata[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=LOCUS_COORDINATES]
msi_pos_locs<-msidata[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=LOCUS_COORDINATES]
mss_pos_locs<-mssdata[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=LOCUS_COORDINATES]

#171 MSI-H samples and 446 MSS so filtering for only sites with sufficient data in at least half of the samples in each category
sites_with_sufficient_data<-intersect(which(msi_locus_nas$V1<=86),which(mss_locus_nas$V1<223))
locus_peak_diffs<-msi_mean_locus_peaks$V1-mss_mean_locus_peaks$V1
pos_peak_diffs<-which(locus_peak_diffs>0)
sub<-intersect(pos_peak_diffs,sites_with_sufficient_data)
sublocs<-as.character(msi_locus_nas$LOCUS_COORDINATES[sub])

#Fisher exact test for disproportionate instability
a<-data.table(locus=msi_locus_nas$LOCUS_COORDINATES,msi_pos_locs$V1,mss_pos_locs$V1,msi_locus_nas$V1,mss_locus_nas$V1)[sub,]
system.time(fisher_odds<-a[,fisher.test(matrix(c(V2,171-V2-V4,V3,446-V3-V5),nrow=2,byrow=T))$estimate,by=locus]) #282s
system.time(fisher_pvals<-a[,fisher.test(matrix(c(V2,171-V2-V4,V3,446-V3-V5),nrow=2,byrow=T))$p.value,by=locus]) #282s

require(qvalue)
fisher_pvals$V1[fisher_pvals$V1>1]<-1
#adjust for multiple hypotheses testing using Storey's q-value
fisher_pvals[,qvals:=qvalue(fisher_pvals$V1)$qvalues]

sig_sites<-a[order(fisher_pvals$qvals),]

sigloc1<-data3[LOCUS_COORDINATES == as.character(sig_sites$locus[1])]
sigloc10<-data3[LOCUS_COORDINATES %in% as.character(sig_sites$locus[1:100])]

mssnum<-446
msinum<-171

merged.data.frame = Reduce(function(...) merge(..., by="locus",all=T), list(sig_sites,fisher_pvals,fisher_odds))
msi_mss_locus_output<-merged.data.frame[order(merged.data.frame$qvals),]
msi_mss_locus_output$V4<-msinum-msi_mss_locus_output$V4
msi_mss_locus_output$V5<-mssnum-msi_mss_locus_output$V5
colnames(msi_mss_locus_output)<-c('locus','msi_peak_unstable','mss_peak_unstable','msi_locus_calls','mss_locus_calls','fisher_test_pval','fisher_test_qval','fisher_test_odds')
msi_mss_locus_output$genomic_class<-locusinfo$GENOMIC_CLASS[match(msi_mss_locus_output$locus,locusinfo$LOCUS_COORDINATES)]
msi_mss_locus_output$gene<-locusinfo$GENE[match(msi_mss_locus_output$locus,locusinfo$LOCUS_COORDINATES)]
msi_mss_locus_output$repeat_type<-locusinfo$REPEAT_TYPE[match(msi_mss_locus_output$locus,locusinfo$LOCUS_COORDINATES)]
msi_mss_locus_output$repeat_dna_sequence<-locusinfo$REPEAT_DNA_SEQUENCE[match(msi_mss_locus_output$locus,locusinfo$LOCUS_COORDINATES)]

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_results_post_review_msi_updated_062416.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_post_review_msi_updated_062416.robj')

require(reshape2)
sigloc10_test<-dcast(sigloc10, SAMPLE_NAME  ~ LOCUS_COORDINATES)
siglocs <- sigloc10_test[,-1]
unstable_discretize<-function(x){
	output<-rep('NA',length(x))
	output[x<=0]<-'stable'
	output[x>0]<-'unstable'
	return(output)
}

siglocs2<-apply(siglocs,2,unstable_discretize)
siglocs2<-siglocs2[match(sampledat1n$sample_name,sigloc10_test$SAMPLE_NAME),]

ucecsmall<-data.frame(sampledat1n, siglocs2)
ucecsmall$tumor_type<-sampleinfo1n$TUMOR_TYPE[match(ucecsmall$sample_name,sampleinfo1n$SAMPLE_NAME)]
ucecsmall$msi<-sampleinfo1n$MSI_STATUS[match(ucecsmall$sample_name,sampleinfo1n$SAMPLE_NAME)]

subsets<-c(seq(1,10),seq(20,100,10))
trainctrl <- trainControl(summaryFunction=twoClassSummary,classProb= TRUE)
rfFuncs$summary<-defaultSummary
ctrl <- rfeControl(functions=rfFuncs,method = "cv", number = 10)
set.seed(101)
ucecsmall[is.na(ucecsmall)]<-0
ucecsmall2<-ucecsmall[ucecsmall$msi=='MSI-H' | ucecsmall$msi=='MSS',]
longresults2 <- rfe(ucecsmall2[,c(2:5,7,13:107)], factor(ucecsmall2$msi), sizes=subsets,rfeControl=ctrl,trControl=trainctrl,metric="Accuracy")

longresults3 <- rfe(ucecsmall2[,c(2,13:107)], factor(ucecsmall2$corrected_msi), sizes=subsets,rfeControl=ctrl,trControl=trainctrl,metric="Kappa")


info<-cbind(misclass,ucecsmall[match(misclass$sample_name,ucecsmall$sample_name),c(1:10,104)])

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_classification_tree_weighted_fancy_post_review_with_sig_locs_062716.pdf',width=7*1.5,height=7*1.5)
plot(longresults2, bty='n', ylim=c(0,1))
dev.off()

ucecsmall2$msi<-as.factor(ucecsmall2$msi)
ucecsmall2$corrected_msi<-sampledat1n$corrected_msi[match(ucecsmall2$sample_name,sampledat1n$sample_name)]
weight_test<-as.factor(ucecsmall2$msi)
levels(weight_test)<-rev(table(weight_test))
weight_test<-as.numeric(as.matrix(weight_test))

system.time(rpartfit2 <- train(msi ~ peak_avg + as.factor(ucecsmall2$X8.7679723.7679741), data = ucecsmall2,
                   method = "rpart",
                   trControl = fitControl,
                   control=rpart.control(minbucket=6, maxdepth=3),
                   tuneGrid=expand.grid(.cp=c(0,0.001,0.01,0.1,0.45,0.95)), 
                   weights=weight_test
                   ))

system.time(rpartfit2_unweighted <- train(msi ~ peak_avg + as.factor(ucecsmall2$X8.7679723.7679741), data = ucecsmall2,
                   method = "rpart",
                   trControl = fitControl,
                   control=rpart.control(minbucket=6, maxdepth=3),
                   tuneGrid=expand.grid(.cp=c(0,0.001,0.01,0.1,0.45,0.95)) 
                   #weights=weight_test
                   ))

#LOOCV
fitControl <- trainControl(## 10-fold CV
  method = "loocv",
  returnData = TRUE,
  returnResamp = "final",
  repeats =  1,
  classProbs = TRUE,
  summaryFunction = twoClassSummary)

ucecsmall2$msi[ucecsmall2$msi=="MSI-H"]<-'MSI.H'
ucecsmall2$msi[ucecsmall2$msi=="MSI.H"]<-'MSI-H'

system.time(rpartfit2_loocv <- train(msi ~ peak_avg + as.factor(ucecsmall2$X8.7679723.7679741), data = ucecsmall2,
                   method = "rpart",
                   trControl = fitControl,
                   control=rpart.control(minbucket=6, maxdepth=3),
                   tuneGrid=expand.grid(.cp=c(0.001,0.01,0.1,0.45)), 
                   weights=weight_test
                   ))

mosaic<-rpartfit2_loocv
save(mosaic,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_classifier_063016.robj')

ucecsmall2$corrected_msi<-as.factor(ucecsmall2$corrected_msi)
irisct <- ctree(corrected_msi ~ ., data = ucecsmall2, controls = ctree_control(max_depth=5, minsplit=6))

library(rattle)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_classification_tree_weighted_fancy_post_review_10xcv_with_best_loc_063016.pdf',width=7*1.5,height=7*1.5)
fancyRpartPlot(rpartfit2$finalModel)
dev.off()

library(rattle)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_classification_tree_weighted_fancy_post_review_loocv_with_best_loc_063016.pdf',width=7*1.5,height=7*1.5)
fancyRpartPlot(rpartfit2_loocv$finalModel)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_classification_tree_weighted_fancy_post_review_10xcv_with_best_loc_unweighted_063016.pdf',width=7*1.5,height=7*1.5)
fancyRpartPlot(rpartfit2_unweighted$finalModel)
dev.off()

library(rattle)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_ctree_classification_tree_062916.pdf',width=7*1.5,height=7*1.5)
plot(as.simpleparty(results))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msiseq_sind_vs_prop_unstable.pdf',width=7*1.5,height=7*1.5)
plot()
dev.off()

#seeing whether mutational data better discretizes UCEC
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_clinical_data.robj')
clindat<-Reduce(function(...) merge(..., all=T), clindatlist)
clindat$tumor_type<-revalue(clindat$tumor_type,c('OVmu'='OV'))
clindat$num_unstable<-sampledat$num_unstable_raw[match(clindat$sample_name,sampledat$sample_name)]
clindat$msi_status<-sampledat$msi_status_tree_scaled[match(clindat$sample_name,sampledat$sample_name)]
clindat$yearstobirth<-as.numeric(as.matrix(clindat$yearstobirth))
clindat$peak_avg<-sampledat$peak_avg[match(clindat$sample_name,sampledat$sample_name)]
clindat<-clindat[-which(clindat$num_unstable>10000),]

require(TCGA2STAT)
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/finalsampledata_062216.robj')
cancers<-as.character(unique(sampledat$tumor_type))
mutdatlist<-list()
for (i in 1:length(cancers)){
	subset.mut<-getTCGA(disease=cancers[i],data.type='Mutation',type='somatic')
	mutdatlist[[i]]<-subset.mut
}

mmrgenes<-c('EXO1','MLH1','MLH3','MSH2','MSH3','MSH6','PMS1','PMS2', 'POLD1', 'POLE')

mmrmutdatlist<-lapply(mutdatlist,function(x){
	output <- x$dat[match(mmrgenes,rownames(x$dat)),]
	rownames(output) <- mmrgenes
	colnames(output) <- colnames(x$dat)
	return(as.data.frame(t(output)))})

mmrmutdat<-do.call(rbind, mmrmutdatlist)
save(mmrmutdat,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mmr_mutational_data.robj')

sampleinfo1n[583:1105,9:18]<-mmrmutdat[match(sampleinfo1n$SAMPLE_NAME,rownames(mmrmutdat))[583:1105],]

save(sampleinfo1n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_post_review_msi_updated_062916.robj')
save(sampledat1n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_results_post_review_msi_updated_062916.robj')

sampledat1n_orig
sampledat1n<-sampledat1n[,c(1:3,5:7,13,14),]
colnames(sampledat1n)[8]<-'msi'

save(sampledat1n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_results_post_review_msi_predicted_063016.robj')
sampleinfo1n$MSI_STATUS<-sampledat1n$msi
save(sampleinfo1n,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_post_review_msi_predicted_063016.robj')

data2<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_msi_values_converted_to_numeric_duplicates_eliminated_post_review_exomes_added_empty_samples_removed_062416.txt', sep='\t', showProgress=T) #32 min

data[, MSI_STATUS := sampledat1n$msi[match(data$SAMPLE_NAME,sampledat1n$sample_name)]]
data[, TUMOR_TYPE := sampleinfo1n$TUMOR_TYPE[match(data$SAMPLE_NAME,sampleinfo1n$SAMPLE_NAME)]]

fwrite(data, '/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_msi_values_converted_to_numeric_duplicates_eliminated_post_review_exomes_added_empty_samples_removed_063016.txt', sep='\t')

#locus output
system.time(locus_output<-data[,j=list(median_peak_diff = as.double(median(PEAK_DIFFERENCE_VALUE,na.rm=T)), num_unstable = length(which(PEAK_DIFFERENCE_VALUE>0)), num_missing =  length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)), num_present = length(which(is.na(PEAK_DIFFERENCE_VALUE)==F))),by=list(TUMOR_TYPE,MSI_STATUS,LOCUS_COORDINATES)]) #500s

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_locusinfo.robj')
colnames(locus_output)[1:3]<-c('tumor_type', 'msi_status', 'locus')
locus_output[, genomic_class := locusinfo$GENOMIC_CLASS[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]]
locus_output$gene<-locusinfo$GENE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_type<-locusinfo$REPEAT_TYPE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_dna_sequence<-locusinfo$REPEAT_DNA_SEQUENCE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
write.table(locus_output,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_tcga_loci_results_060616.txt',quote=F,row.names=F)