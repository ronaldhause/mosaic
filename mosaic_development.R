source('~/helper_functions.R')
require(qvalue)
require(caret)
require(doMC)
require(reshape2)

load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/UCEC_1-7_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/READ_1-3_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/COAD_1-9_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/UCEC_4-7_merged_full_analysis.robj')

#cleaning up duplicate sample entries in ucec
ucec2<-unique(ucec,by=c("SAMPLE_NAME","LOCUS_COORDINATES"))

#save exome results from cancers with MSI statuses for MOSAIC training
save(data,file='/net/shendure/vol7/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated.robj')

#fix errors in MSI statuses from TCGA server
setkey(data,SAMPLE_NAME)
system.time(data['TCGA-B5-A1MU',TUMOR_TYPE := 'UCEC'])
system.time(data['TCGA-FI-A2EW',TUMOR_TYPE := 'UCEC'])
system.time(data['TCGA-B5-A1MU',MSI_STATUS := 'MSS'])
system.time(data['TCGA-FI-A2EW',MSI_STATUS := 'ND'])
data<-droplevels.data.table(data)

#calculating average peak_diff etc. for each sample
system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #19s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #8s
system.time(num_unstable<-data[,length(which(KS_VALUE<0.05)),by=SAMPLE_NAME]) #7s

#load info for these cancers
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_030115.robj')

#separating data into MSI-H and MSS
setkey(data,MSI_STATUS)
msidata<-data["MSI-H"]
mssdata<-data["MSS"]


#preliminary examine of sites differentially unstable in MSI-H vs. MSS
setkey(msidata,LOCUS_COORDINATES)
setkey(mssdata,LOCUS_COORDINATES)
msi_mean_locus_peaks<-msidata[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=LOCUS_COORDINATES]
msi_median_locus_peaks<-msidata[,median(PEAK_DIFFERENCE_VALUE,na.rm=T),by=LOCUS_COORDINATES]
mss_mean_locus_peaks<-mssdata[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=LOCUS_COORDINATES]
msi_locus_nas<-msidata[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=LOCUS_COORDINATES]
mss_locus_nas<-mssdata[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=LOCUS_COORDINATES]
msi_sig_locs<-msidata[,length(which(KS_VALUE<0.05)),by=LOCUS_COORDINATES]
mss_sig_locs<-mssdata[,length(which(KS_VALUE<0.05)),by=LOCUS_COORDINATES]
msi_pos_locs<-msidata[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=LOCUS_COORDINATES]
mss_pos_locs<-mssdata[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=LOCUS_COORDINATES]

#90 MSI-H samples and 287 MSS so filtering for only sites with sufficient data in at least half of the samples in each category
sites_with_sufficient_data<-intersect(which(msi_locus_nas$V1<=40),which(mss_locus_nas$V1<135))
locus_peak_diffs<-msi_mean_locus_peaks$V1-mss_mean_locus_peaks$V1
pos_peak_diffs<-which(locus_peak_diffs>0)
sub<-intersect(pos_peak_diffs,sites_with_sufficient_data)
sublocs<-as.character(msi_locus_nas$LOCUS_COORDINATES[sub])

#Fisher exact test for disproportionate instability
a<-data.table(locus=msi_sig_locs$LOCUS_COORDINATES,msi_sig_locs$V1,mss_sig_locs$V1,msi_pos_locs$V1,mss_pos_locs$V1,msi_locus_nas$V1,mss_locus_nas$V1)[sub,]
setkey(a,locus)
system.time(fisher_odds<-a[,fisher.test(matrix(c(V4,80-V4-V6,V5,270-V5-V7),nrow=2,byrow=T))$estimate,by=locus]) #71s
system.time(fisher_pvals<-a[,fisher.test(matrix(c(V4,80-V4-V6,V5,270-V5-V7),nrow=2,byrow=T))$p.value,by=locus]) #71s
test<-a[1,]

fisher_pvals[,qvals:=qvalue2(fisher_pvals$V1)$qvalues]

setkey(data,LOCUS_COORDINATES)
head(fisher_pvals[order(qvals)])
sig_sites<-a[order(fisher_pvals$qvals),]
sigloc1<-data[as.character(sig_sites$locus[1])]
sigloc10<-data[as.character(sig_sites$locus[1:10000])]

sigloc10_test<-tapply(sigloc10$PEAK_DIFFERENCE_VALUE,as.character(sigloc10$LOCUS_COORDINATES),function(x)x[match(peak_avg$SAMPLE_NAME,sigloc1$SAMPLE_NAME)])
siglocs <- matrix(unlist(sigloc10_test), ncol = 10000)
ucecsmall<-data.frame(sample_name=peak_avg$SAMPLE_NAME,peak_avg=peak_avg$V1,peak_sd=peak_sd$V1,num_unstable=num_unstable$V1,msi=loc1$MSI_STATUS[match(peak_avg$SAMPLE_NAME,loc1$SAMPLE_NAME)])
ucecsmall2<-ucecsmall[which(ucecsmall$msi=='MSI-H' | ucecsmall$msi=='MSS'),]
ucecsmall2$msi<-factor(ucecsmall2$msi)
ucecsmall2<-ucecsmall2[,-1]
ucecsmall<-data.frame(sample_name=peak_avg$SAMPLE_NAME,peak_avg=peak_avg$V1,peak_sd=peak_sd$V1,num_unstable=num_unstable$V1,msi=loc1$MSI_STATUS[match(peak_avg$SAMPLE_NAME,loc1$SAMPLE_NAME)],siglocs[,1:1000])
ucecsmall22<-ucecsmall[which(ucecsmall$msi=='MSI-H' | ucecsmall$msi=='MSS'),]
ucecsmall22$msi<-factor(ucecsmall22$msi)
ucecsmall22<-ucecsmall22[,-1]
save(ucecsmall,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_training_data_sig_loc_incl_030315.robj')
save(ucecsmall2,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_training_data_030315.robj')

ucecsmallnew<-ucecsmall
ucecsmall2new<-ucecsmall2

sample_sub<-ucecsmall[which(ucecsmall$msi=='MSI-H' | ucecsmall$msi=='MSS'),]
ucecsmall2<-ucecsmall[which(ucecsmall$msi=='MSI-H' | ucecsmall$msi=='MSS'),]
ucecsmall2$msi<-factor(ucecsmall2$msi)
ucecsmall2$tumor_type<-loc1$TUMOR_TYPE[match(ucecsmall2$sample_name,loc1$SAMPLE_NAME)]
ucecsmall2$num_unstable_raw<-ucecsmall$num_unstable_raw[match(ucecsmall2$sample_name,ucecsmall$sample_name)]

sampledat1<-data.frame(sample_name=ucecsmall2$sample_name,peak_avg=ucecsmall2$peak_avg,peak_sd=ucecsmall2$peak_sd,num_unstable_ks=ucecsmall2$num_unstable,num_unstable_raw=ucecsmall2$num_unstable_raw,msi_status=ucecsmall2$msi,tumor_type=ucecsmall2$tumor_type)
save(sampledat1,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampledata_030115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampledata_030115.robj')

ucecsmall22[is.na(ucecsmall22)]<-0

 subsets<-c(1:25,50,100)

# define the control using a random forest selection function
#optimizing using accuracy
trainctrl <- trainControl(summaryFunction=twoClassSummary,classProb= TRUE)
rfFuncs$summary<-defaultSummary
ctrl <- rfeControl(functions=rfFuncs,method = "cv",number = 10)
ctrl <- rfeControl(functions=rfFuncs,method = "loocv")
set.seed(101)
 #dat<-preProcess(ucecsmall2[,-4], method='knnImpute')
 #variable selection
 system.time(results <- rfe(ucecsmall2[,-4], ucecsmall2$msi, sizes=subsets,rfeControl=ctrl,trControl=trainctrl,metric="Accuracy")) #60s, 7 variables achieve 97.7% accuracy; 3 variables achieves 97.4% accuracy, 92.6% kappa (03/05/15); 
  system.time(longresults <- rfe(ucecsmall22[,-5], ucecsmall22$msi, sizes=subsets,rfeControl=ctrl,trControl=trainctrl,metric="Kappa")) #60s, 7 variables achieve 97.7% accuracy; 3 variables achieves 97.4% accuracy, 92.6% kappa (03/05/15);

 r<-rfe(ucecsmall2[,-4], ucecsmall2$msi, sizes=subsets)

rfFuncs$summary<-twoClassSummary
trainctrl <- trainControl(classProb= TRUE,summaryFunction = twoClassSummary)
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
system.time(results2 <- rfe(ucecsmall2[,-4], ucecsmall2$msi, sizes=subsets, rfeControl=control,metric="Accuracy",trControl=trainctrl)) #60s, 7 variables achieves 98.9% ROC, 92.5% sensitivity, 98.5% specificity; update: 99.26% ROC, 91.25% sens, 98.52% spec with 3 variables, 03/01/15

library(plyr)
library(scales)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_accuracy_vs_features_031115.pdf',width=7,height=7)
#trellis.par.set(caretTheme())
plotdatvar<-data.frame(y=longresults$results$Accuracy,ymin=longresults$results$Accuracy-2*longresults$results$AccuracySD,ymax=longresults$results$Accuracy+2*longresults$results$AccuracySD,variables=longresults$results$Variables,accuracy=longresults$results$Accuracy,kappa=longresults$results$Kappa)
plotdat<-data.frame(variables=longresults$results$Variables,accuracy=longresults$results$Accuracy,kappa=longresults$results$Kappa)
ggplot() + geom_line(aes(x=variables,y=accuracy),data=plotdat) + geom_smooth(aes(x=variables,y=accuracy,ymin=ymin,ymax=ymax),data=plotdatvar,stat="identity") + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylim(0.50,1.05) + scale_x_log10(breaks=c(1,2,10,100,1000),labels=comma(c(1,2,10,100,1000))) + geom_vline(xintercept=2,linetype="dashed")
#plot(results$results$Accuracy,type=c('g','o'),ylim=c(0,1))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/peak_avg_between_mss_msi.pdf',width=7,height=7)
a<-ggplot(ucecsmall2,aes(y=ucecsmall2$peak_avg,x=ucecsmall2$msi)) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20)
#plot(results$results$Accuracy,type=c('g','o'),ylim=c(0,1))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/peak_sd_between_mss_msi.pdf',width=7,height=7)
b<-ggplot(ucecsmall2,aes(y=ucecsmall2$peak_sd,x=ucecsmall2$msi)) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20)
#plot(results$results$Accuracy,type=c('g','o'),ylim=c(0,1))
dev.off()

num_unstable<-ucecsmall2$num_unstable

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/frac_unstable_between_mss_msi.pdf',width=7,height=7)
c<-ggplot(ucecsmall2,aes(y=num_unstable+0.5,x=ucecsmall2$msi)) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + scale_y_log10()
#plot(results$results$Accuracy,type=c('g','o'),ylim=c(0,1))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/best_site_between_mss_msi.pdf',width=7,height=7)
d<-ggplot(ucecsmall2,aes(y=ucecsmall2$X812,x=ucecsmall2$msi)) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20)
#plot(results$results$Accuracy,type=c('g','o'),ylim=c(0,1))
dev.off()

#in tumor relative to paired normal
a<-ggplot(ucecsmall2,aes(y=ucecsmall2$peak_avg,x=ucecsmall2$msi)) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('global average allele number difference') + xlab('') + theme(axis.title=element_text(size=16))
b<-ggplot(ucecsmall2,aes(y=ucecsmall2$peak_sd,x=ucecsmall2$msi)) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('global sd in allele number difference') + xlab('') + theme(axis.title=element_text(size=16))
c<-ggplot(ucecsmall2,aes(y=num_unstable+0.5,x=ucecsmall2$msi)) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + scale_y_log10() + ylab('global number of significantly unstable loci') + xlab('') + theme(axis.title=element_text(size=16))
d<-ggplot(ucecsmall2,aes(y=ucecsmall2$X812,x=ucecsmall2$msi)) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('allele number difference at 1:115226986-115227007') + xlab('') + theme(axis.title=element_text(size=16))

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_comparison_summary.pdf',width=7*1.5,height=7*1.5)
multiplot(a,c,b,d, cols=2)
dev.off()

ucecsmall3<-ucecsmall
ucecsmall3<-ucecsmall[which(ucecsmall$msi=='MSI-L'),]
ucecsmall3$msi<-factor(ucecsmall3$msi)
ucecsmall3$lab<-rep('MSS',dim(ucecsmall3)[1])
ucecsmall3<-ucecsmall3[,-1]
num_unstable3<-ucecsmall3$num_unstable
a<-ggplot(ucecsmall2,aes(y=ucecsmall2$peak_avg,x=ucecsmall2$msi)) + geom_boxplot() + geom_jitter(data=ucecsmall3,aes(y=ucecsmall3$peak_avg,x=ucecsmall3$lab,color='a')) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('global average allele number difference') + xlab('') + theme(axis.title=element_text(size=16)) + theme(legend.position="none")
b<-ggplot(ucecsmall2,aes(y=ucecsmall2$peak_sd,x=ucecsmall2$msi)) + geom_boxplot() + geom_jitter(data=ucecsmall3,aes(y=ucecsmall3$peak_sd,x=ucecsmall3$lab,color='a')) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('global sd in allele number difference') + xlab('') + theme(axis.title=element_text(size=16)) + theme(legend.position="none")
c<-ggplot(ucecsmall2,aes(y=num_unstable+0.5,x=ucecsmall2$msi)) + geom_boxplot() + geom_jitter(data=ucecsmall3,aes(y=num_unstable3+0.5,x=ucthemecsmall3$lab,color='a')) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + scale_y_log10() + ylab('global number of significantly unstable loci') + xlab('') + theme(axis.title=element_text(size=16)) + theme(legend.position="none")
d<-ggplot(ucecsmall2,aes(y=ucecsmall2$X812,x=ucecsmall2$msi)) + geom_boxplot() + geom_jitter(data=ucecsmall3,aes(y=ucecsmall3$X812,x=ucecsmall3$lab,color='a')) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('allele number difference at 1:115226986-115227007') + xlab('') + theme(axis.title=element_text(size=16)) + theme(legend.position="none")

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_comparison_summary_with_msil.pdf',width=7*1.5,height=7*1.5)
multiplot(a,c,b,d, cols=2)
dev.off()

# ucecsmall2[is.na(ucecsmall2)]<-0
# set.seed(10)
#  rf.profileROC.Radial <- rfe(ucecsmall2[,-4], ucecsmall2$msi, 
#                              sizes=subsets,
#                              rfeControl=ctrl,
#                              method="rf",
#                              ## I also added this line to
#                              ## avoid a warning:
#                              metric = "ROC",
#                              trControl = trainctrl)



fitControl <- trainControl(## 10-fold CV
  method = "loocv",
  returnData = TRUE,
  returnResamp = "final",
  #number = 10,
  repeats =  1)

#fitControl <- trainControl(method = "loocv")

gbmGrid <-  expand.grid(mtry=2:3)

set.seed(101)

registerDoMC(cores = 4)
  system.time(msirffitt <- train(msi ~ ., data = ucecsmallt,
                   method = "rf",
                   trControl = fitControl,
                   verbose = FALSE,
                   prox=TRUE,
                   #tuneGrid=gbmGrid,
                   n.trees=1000,
                   allowParallel=TRUE)) #3s, mtry=2, 97.1% accuracy, 91.7% kappa; 96.9% accuracy, 91.1% kappa

  ucecsmallt<-ucecsmall2
  ucecsmallt$num_unstable_ks<-ucecsmall22$num_unstable_ks
  ucecsmallt$loc1<-ucecsmall22$X698
  ucecsmallt$loc2<-ucecsmall22$X333
  ucecsmallt$loc3<-ucecsmall22$X64

#weighted
weight_test<-ucecsmall2$msi
levels(weight_test)<-rev(table(weight_test))
weight_test<-as.numeric(as.matrix(weight_test))

set.seed(101)
  system.time(msirffitw1 <- train(msi ~ peak_avg + peak_sd + num_unstable, data = ucecsmall3,
                   method = "rf",
                   trControl = fitControl,
                   verbose = FALSE,
                   prox=TRUE,
                   #tuneGrid=gbmGrid,
                   weights=weight_test,
                   n.trees=1000,
                   allowParallel=TRUE))

  fitControl <- trainControl(## 10-fold CV
  method = "cv",
  returnData = TRUE,
  returnResamp = "final",
  number = 10,
  repeats =  1)

    system.time(rpartfit1 <- train(msi ~ ., data = ucecsmall2,
                   method = "rpart",
                   trControl = trainControl(method="loocv"),
                   control=rpart.control(minsplit=2),
                   tuneGrid=expand.grid(.cp=c(0,0.001,0.01,0.1,0.45,0.95)), 
                   weights=weight_test
                   ))

msimodel<-rffit2
save(msimodel,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_rf_classifier.robj')
msimodel2<-msirffitw1
save(msimodel2,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_rf_classifier_peak_avg_median_scaled_041315.robj')
save(rpartfit1,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_tree_classifier_041315.robj')
save(rpartfit2,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_tree_classifier_peak_avg_median_scaled_041315.robj')
save(rpartfitt,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_tree_classifier_peak_avg_median_scaled_best_loc_included_041315.robj')


library(rattle)
      pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_classification_tree_test_2_041315.pdf',width=7*1.5,height=7*1.5)
  fancyRpartPlot(rpartfitt2$finalModel)
  dev.off()

#misclassified samples, reference Resample number
misclass<-rpartfit1$resample$Resample[which(rpartfit1$resample$Accuracy==0)]
a<-ucecsmallt[c(8,11,141,20,83,125,24,291),]
b<-sampledat1[c(8,11,141,20,83,125,24,291),]
misclass2<-rpartfitt$resample$Resample[which(rpartfitt$resample$Accuracy==0)]

        pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_classification_tree_1_041215.pdf',width=7*1.5,height=7*1.5)
  plot(rpartfit2$finalModel, main="Classification Tree for MSI")
  text(rpartfit2$finalModel, use.n=TRUE, all=TRUE, cex=.8)
  dev.off()


  fitControl <- trainControl(## 10-fold CV
  method = "loocv",
  returnData = TRUE,
  returnResamp = "final",
  #number = 10,
  repeats =  1)

  fitControl <- trainControl(## 10-fold CV
  method = "cv",
  returnData = TRUE,
  returnResamp = "final",
  #number = 10,
  repeats =  1)

set.seed(101)
    system.time(msirffitwscaled <- train(msi ~ peak_avg + peak_sd + num_unstable, data = ucecsmall3,
                   method = "rf",
                   trControl = fitControl,
                   verbose = FALSE,
                   prox=TRUE,
                   #tuneGrid=gbmGrid,
                   weights=weight_test,
                   n.trees=1000,
                   allowParallel=TRUE))

     system.time(results <- rfe(sampledat1[,2:5], ucecsmall2$msi, sizes=1:4,rfeControl=ctrl,trControl=trainctrl,metric="Accuracy"))

  ucecsmall3<-ucecsmall2
  ucecsmall3$peak_avg<-ave(sampledat1$peak_avg,sampledat1$tumor_type,FUN=function(x)x-median(x))
  ucecsmall3$peak_sd<-ave(sampledat1$peak_sd,sampledat1$tumor_type,FUN=function(x)x-median(x))
  ucecsmall3$num_unstable<-ave(sampledat1$num_unstable_raw,sampledat1$tumor_type,FUN=function(x)x-median(x))

  #when doing LOOCV to identify misclassified samples:
msirffitw$resample
sample_sub[c(8,11,20,24,29,83),1:5]

  library(rpart)
  fit <- rpart(msi ~ peak_avg + peak_sd + num_unstable, data = ucecsmall2,method='class',control=rpart.control(minsplit=10,cp=0.001))

  fit <- rpart(msi ~ ., data = ucecsmall2, method='class',weights=weight_test,control=rpart.control(xval=350))

  pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_classification_tree_041215.pdf',width=7*1.5,height=7*1.5)
  plot(fit, uniform=TRUE, main="Classification Tree for MSI")
  text(fit, use.n=TRUE, all=TRUE, cex=.8)
  dev.off()

library("rpart.plot")
  pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_classification_rpart_tree_041215.pdf',width=7*1.5,height=7*1.5)
  rpart.plot(fit,extra=1)
  dev.off()

  #with 0 loci - 97.4% accuracy, 92.7% kappa with mtry = 2
  #with 1 locus - 97.8% accuracy, 94.1% kappa with mtry = 2
  #with 10 loci - 98.1% accuracy, 93.8% kappa with mtry = 7
  #with 100 loci - 90.8% accuracy with 500 trees and m.try=4
  #with 1000 loci

  registerDoMC(cores = 4)
  system.time(gbmfitloconly <- train(msi ~ ., data = ucecsmall2,
                   method = "gbm",
                   trControl = fitControl,
                   verbose = FALSE)) #3s, n.trees=50, interaction.depth=1, 97.4% accuracy, 92.5% kappa

    registerDoMC(cores = 4)
  system.time(svmfit1 <- train(msi ~ ., data = ucecsmall2,
                   method = "svmRadial",
                   trControl = fitControl,
                   verbose = FALSE)) #5s, C=0.5, 96.6% accuracy, 90.3% kappa

  varImp(rffit1,scale=F)

  treepred <- predict(msirffitw,newdata=ucecsmall2,type='raw')
  treeconfmat<-table(pred = treepred, true = ucecsmall2$msi)

  pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_gbm_var_imp_batch1.pdf',width=7,height=7)
  plot(varImp(gbmFit1,scale=FALSE),top = 40)
  dev.off()
  pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_gbm_accuracy_batch1.pdf',width=7,height=7)
  ggplot(gbmFit1,metric = "Accuracy") + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20)
  dev.off()

  treepred <- predict(msirffit,newdata=ucecsmall2,type='raw')
  treeconfmat<-table(pred = treepred, true = ucecsmall2$msi)

  #loci that have different sample proportions of MSI between tumor types
  #MSI in genomic features per cancer type
  #sites that best distinguish cancer types based on msi frequencies

  peak_diff<-as.numeric(as.matrix(ucec$PEAK_DIFFERENCE_VALUE))
  ucec_dat<-data.frame(MSI_STATUS=ucec$MSI_STATUS,LOCUS_COORDINATES=ucec$LOCUS_COORDINATES,PEAK_DIFF=peak_diff,SAMPLE_NAME=ucec$SAMPLE_NAME)
  ucec_dat<-ucec_dat[-which(ucec_dat$MSI_STATUS=='NaN'),]
  ucec_dat<-ucec_dat[-which(ucec_dat$MSI_STATUS=='MSS'),]
  ucec_dat$MSI_STATUS<-factor(ucec_dat$MSI_STATUS)

  require(rpart)
  require(caret)

  mod<-rpart(MSI_STATUS ~ .,data=ucec_dat)

  #looking to identify MSI-H vs. MSS without class IDs
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_sep_2d_031515.pdf',width=7,height=7)
ggplot(ucecsmall2,aes(y=num_unstable_raw,x=peak_avg)) + geom_point(aes(color=factor(msi))) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + theme(legend.position=c(0.2,0.8)) + theme(legend.title=element_blank())
#plot(results$results$Accuracy,type=c('g','o'),ylim=c(0,1))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_sep_2d_cancer_split_031515.pdf',width=7*3,height=7)
ggplot(ucecsmall2,aes(y=num_unstable_raw,x=peak_avg)) + geom_point(aes(color=factor(msi))) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + theme(legend.position=c(0.2,0.8)) + theme(legend.title=element_blank()) + facet_grid(. ~ tumor_type)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/frac_unstable_per_cancer_between_mss_msi_031515.pdf',width=7*1.5,height=7)
ggplot(ucecsmall2,aes(y=num_unstable_raw,x=ucecsmall2$msi)) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + facet_grid(. ~ tumor_type)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/peak_avg_per_cancer_between_mss_msi_031515.pdf',width=7*1.5,height=7)
ggplot(ucecsmall2,aes(y=peak_avg,x=ucecsmall2$msi)) + geom_boxplot() + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + facet_grid(. ~ tumor_type)
dev.off()

#testing for enriched MMR mutations correlated at specific MSI sites
#first batch of 118 for UCEC
#i = length(unique(ucec$LOCUS_COORDINATES))

pvals<-c()

for (i in 1:10){
  for (j in 3:length(unique(ucec$EXO1))){
    locuscoord<-which(ucec$LOCUS_COORDINATES==unique(ucec$LOCUS_COORDINATES)[i])
    mutcoord<-which(ucec$EXO1==unique(ucec$EXO1)[j])
    fisher.test(matrix(c(length(intersect(locuscoord,mutcoord)),118-length(intersect(locuscoord,mutcoord)),length(mutcoord),length(unique(ucec$LOCUS_COORDINATES))-length(mutcoord)),nrow=2,byrow=T))
  }

}

library(ggplot2)
data<-data.frame(EXO1_status=ucec$EXO1,peak_difference=as.numeric(as.matrix(ucec$PEAK_DIFFERENCE_VALUE)))
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/eco1_peak_difference_value.pdf',width=7*1.5, height=7)
p <- ggplot(data,aes(EXO1_status,peak_difference))
p + geom_boxplot() + theme_bw(base_size=14)
dev.off()

data<-data3
library(data.table)
setkey(data,SAMPLE_NAME)
system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #17s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #9s
system.time(num_unstable_ks<-data[,length(intersect(which(KS_VALUE<0.05),which(PEAK_DIFFERENCE_VALUE>0))),by=SAMPLE_NAME]) #11s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=SAMPLE_NAME]) #7s

sampledat1<-data.frame(sample_name=peak_avg$SAMPLE_NAME,peak_avg=peak_avg$V1,peak_sd=peak_sd$V1,num_unstable_ks=num_unstable_ks$V1,num_unstable_raw=num_unstable$V1,num_unstable=num_unstable$V1)
treepred <- predict(msimodel,newdata=sampledat1,type='raw')
treepred2 <- predict(rpartfit1,newdata=sampledat1,type='raw')

sampledat1n<-sampledat1
sampledat1n$tumor_type<-sampleinfo1$TUMOR_TYPE
sampledat1n$peak_avg<-ave(sampledat1n$peak_avg,sampledat1n$tumor_type,FUN=function(x)x-median(x))
treepred3 <- predict(msimodel2,newdata=sampledat1n,type='raw')
treepred4 <- predict(rpartfit2,newdata=sampledat1n,type='raw')

sampledat1$msi_status_rf<-treepred
sampledat1$msi_status_tree<-treepred2
sampledat1$msi_status_rf_scaled<-treepred3
sampledat1$msi_status_tree_scaled<-treepred4
sampledat1$tumor_type<-sampleinfo1$TUMOR_TYPE
sampledat1<-sampledat1[,-6]
save(sampledat1,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampledata_041315.robj')

data$MSI_STATUS_PRED<-treepred4[match(data$SAMPLE_NAME,sampledat1$sample_name)]
data1<-data
save(data1,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated_msi_values_converted_to_numeric_nans_corrected_msi_predicted.robj')

setkey(data,LOCUS_COORDINATES)
loc1<-data["10:100008321-100008335"]
sampleinfo3<-loc1

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_030115.robj')

alllocs<-rbind(sampleinfo1,sampleinfo3)
msifracs<-tapply(alllocs$MSI_STATUS,alllocs$TUMOR_TYPE,function(x)length(which(x=='MSI-H'))/length(x))
bardat<-data.frame(cancer=names(msifracs)[-c(1,7)],prop=msifracs[-c(1,7)])

#3 - core loci associated with MSI

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