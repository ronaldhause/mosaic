source('~/helper_functions.R')
require(qvalue)
require(caret)
require(doMC)
require(reshape2)
require(plyr)
require(scales)

load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/UCEC_1-7_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/READ_1-3_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/COAD_1-9_merged_full_analysis.robj')
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/UCEC_4-7_merged_full_analysis.robj')

#cleaning up duplicate sample entries in ucec
ucec2<-unique(ucec,by=c("SAMPLE_NAME","LOCUS_COORDINATES"))

#save exome results from cancers with MSI statuses for MOSAIC training
save(data,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated.robj')

#fix errors in MSI statuses from TCGA server
setkey(data,SAMPLE_NAME)
system.time(data['TCGA-B5-A1MU',TUMOR_TYPE := 'UCEC'])
system.time(data['TCGA-FI-A2EW',TUMOR_TYPE := 'UCEC'])
system.time(data['TCGA-B5-A1MU',MSI_STATUS := 'MSS'])
system.time(data['TCGA-FI-A2EW',MSI_STATUS := 'ND'])
data[,MSI_STATUS:=droplevels(MSI_STATUS)]
data[,TUMOR_TYPE:=droplevels(TUMOR_TYPE)]

#calculating average peak_diff etc. for each sample
system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #19s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=SAMPLE_NAME]) #8s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=SAMPLE_NAME]) #7s
system.time(num_called<-data[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==F)),by=SAMPLE_NAME]) #7s

#load info for these cancers
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_030115.robj')

#separating data into MSI-H and MSS
setkey(data,MSI_STATUS)
msidata<-data["MSI-H"]
mssdata<-data["MSS"]

#preliminarily examine of sites differentially unstable in MSI-H vs. MSS
setkey(msidata,LOCUS_COORDINATES)
setkey(mssdata,LOCUS_COORDINATES)
msi_mean_locus_peaks<-msidata[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=LOCUS_COORDINATES]
mss_mean_locus_peaks<-mssdata[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=LOCUS_COORDINATES]
msi_locus_nas<-msidata[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=LOCUS_COORDINATES]
mss_locus_nas<-mssdata[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=LOCUS_COORDINATES]
msi_sig_locs<-msidata[,length(intersect(which(KS_VALUE<0.05), which(PEAK_DIFFERENCE_VALUE>0)),by=LOCUS_COORDINATES]
mss_sig_locs<-mssdata[,length(intersect(which(KS_VALUE<0.05), which(PEAK_DIFFERENCE_VALUE>0)),by=LOCUS_COORDINATES]
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
#V4 = msi_pos_locs
#V5 = mss_pos_locs
#V6 = msi_locus_nas
#V7 = mss_locus_nas

#adjust for multiple hypotheses testing using Storey's q-value
fisher_pvals[,qvals:=qvalue(fisher_pvals$V1)$qvalues]

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
ucecsmall2$tumor_type<-loc1$TUMOR_TYPE[match(ucecsmall2$sample_name,loc1$SAMPLE_NAME)]
ucecsmall2$num_unstable_raw<-ucecsmall$num_unstable_raw[match(ucecsmall2$sample_name,ucecsmall$sample_name)]

ucecsmall_withsigsites<-data.frame(sample_name=peak_avg$SAMPLE_NAME,peak_avg=peak_avg$V1,peak_sd=peak_sd$V1,num_unstable=num_unstable$V1,msi=loc1$MSI_STATUS[match(peak_avg$SAMPLE_NAME,loc1$SAMPLE_NAME)],siglocs[,1:1000])
ucecsmall_withsigsites<-ucecsmall_withsigsites[which(ucecsmall_withsigsites$msi=='MSI-H' | ucecsmall_withsigsites$msi=='MSS'),]
ucecsmall_withsigsites$msi<-factor(ucecsmall_withsigsites$msi)
ucecsmall_withsigsites<-ucecsmall_withsigsites[,-1]
save(ucecsmall_withsigsites,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_training_data_sig_loc_incl_030315.robj')
save(ucecsmall2,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_training_data_030315.robj')

sampledat1<-data.frame(sample_name=ucecsmall2$sample_name,peak_avg=ucecsmall2$peak_avg,peak_sd=ucecsmall2$peak_sd,num_unstable_ks=ucecsmall2$num_unstable,num_unstable_raw=ucecsmall2$num_unstable_raw,msi_status=ucecsmall2$msi,tumor_type=ucecsmall2$tumor_type)
save(sampledat1,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampledata_030115.robj')

# define the control using a random forest selection function
#optimizing using accuracy
subsets<-c(1:25,50,100)
trainctrl <- trainControl(summaryFunction=twoClassSummary,classProb= TRUE)
rfFuncs$summary<-defaultSummary
ctrl <- rfeControl(functions=rfFuncs,method = "cv")
set.seed(101)
#variable selection
results <- rfe(ucecsmall2[,-4], ucecsmall2$msi, sizes=subsets,rfeControl=ctrl,trControl=trainctrl,metric="Accuracy")
longresults <- rfe(ucecsmall_withsigsites[,-5], ucecsmall_withsigsites$msi, sizes=subsets,rfeControl=ctrl,trControl=trainctrl,metric="Kappa"))
#7 variables achieve 97.7% accuracy; 3 variables achieves 97.4% accuracy, 92.6% kappa. 99.26% ROC, 91.25% sens, 98.52% spec

#Supplementary Figure 2B - accuracy as a function of the number of top features included in MOSAIC
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_accuracy_vs_features.pdf',width=7,height=7)
plotdatvar<-data.frame(y=longresults$results$Accuracy,ymin=longresults$results$Accuracy-2*longresults$results$AccuracySD,ymax=longresults$results$Accuracy+2*longresults$results$AccuracySD,variables=longresults$results$Variables,accuracy=longresults$results$Accuracy,kappa=longresults$results$Kappa)
plotdat<-data.frame(variables=longresults$results$Variables,accuracy=longresults$results$Accuracy,kappa=longresults$results$Kappa)
ggplot() + geom_line(aes(x=variables,y=accuracy),data=plotdat) + geom_smooth(aes(x=variables,y=accuracy,ymin=ymin,ymax=ymax),data=plotdatvar,stat="identity") + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylim(0.50,1.05) + scale_x_log10(breaks=c(1,2,10,100,1000),labels=comma(c(1,2,10,100,1000))) + geom_vline(xintercept=2,linetype="dashed")
dev.off()

#Figure 1F - boxplots of top variables separating MSI-H vs. MSS
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

#MOSAIC classifier

fitControl <- trainControl(## 10-fold CV
  method = "loocv",
  returnData = TRUE,
  returnResamp = "final",
  repeats =  1)

gbmGrid <-  expand.grid(mtry=2:3)

set.seed(101)

registerDoMC(cores = 4)
  system.time(msirffitt <- train(msi ~ ., data = ucecsmallt,
                   method = "rf",
                   trControl = fitControl,
                   verbose = FALSE,
                   prox=TRUE,
                   tuneGrid=gbmGrid,
                   n.trees=1000,
                   allowParallel=TRUE)) #3s, mtry=2, 97.1% accuracy, 91.7% kappa

ucecsmallt<-ucecsmall2
ucecsmallt$num_unstable_ks<-ucecsmall22$num_unstable_ks
ucecsmallt$loc1<-ucecsmall22$X698
ucecsmallt$loc2<-ucecsmall22$X333
ucecsmallt$loc3<-ucecsmall22$X64

#weighted classification to correct for unequal class balance
weight_test<-ucecsmall2$msi
levels(weight_test)<-rev(table(weight_test))
weight_test<-as.numeric(as.matrix(weight_test))

set.seed(101)
system.time(msirffitw1 <- train(msi ~ peak_avg + peak_sd + num_unstable, data = ucecsmall3,
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

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_classification_tree_1_041215.pdf',width=7*1.5,height=7*1.5)
plot(rpartfit2$finalModel, main="Classification Tree for MSI")
text(rpartfit2$finalModel, use.n=TRUE, all=TRUE, cex=.8)
dev.off()


fitControl <- trainControl(## 10-fold CV
method = "loocv",
returnData = TRUE,
returnResamp = "final",
repeats =  1)

fitControl <- trainControl(## 10-fold CV
method = "cv",
returnData = TRUE,
returnResamp = "final",
repeats =  1)

set.seed(101)
msirffitwscaled <- train(msi ~ peak_avg + peak_sd + num_unstable, data = ucecsmall3,
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

treepred <- predict(msirffit,newdata=ucecsmall2,type='raw')
treeconfmat<-table(pred = treepred, true = ucecsmall2$msi)