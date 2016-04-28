source('/net/shendure/vol1/home/hauser/scripts/shendure/miscellaneous/useful_functions.r')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_042515.robj')
#516,876 sites, 4224 samples
results_allloci$genomic_class<-factor(results_allloci$genomic_class)
results_allloci2<-results_allloci[which(results_allloci$num_present>(4224/2)),]
results_allloci2$prop<-((results_allloci2$simple_gain_site/results_allloci2$num_present)*100)
#222,357 sites have sufficient data in at least 2112 samples

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_training_data_030315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_tree_classifier_peak_avg_median_scaled_041315.robj')

require(caret)
colnames(ucecsmall2)[4]<-'num_unstable'
rpartPred<-predict(rpartfit2,ucecsmall2,type="raw")
confusionMatrix(rpartPred,ucecsmall2$msi)
#95.0% sensitivity, 97.8% specificity (2.22% FPR)

require(plyr)
ucecsmall2$msi<-revalue(ucecsmall2$msi,c("MSI-H"="MSI.H"))
weight_test<-ucecsmall2$msi
levels(weight_test)<-rev(table(weight_test))
weight_test<-as.numeric(as.matrix(weight_test))
set.seed(101)
system.time(rpartfit1 <- train(msi ~ ., data = ucecsmall2[,c(1,2,4,5)],
                   method = "rpart",
                   trControl = trainControl(method="loocv"),
                   control=rpart.control(minsplit=2),
                   tuneGrid=expand.grid(.cp=c(0,0.001,0.01,0.1,0.45,0.95)), 
                   weights=weight_test
                   ))

msiPred<-predict(msimodel,ucecsmall2,type="raw")
confusionMatrix(msiPred,ucecsmall2$msi)

require('ROCR')
rpartPred<-predict(rpartfit1,newdat=ucecsmall2[,c(1,2,4,5)],type="prob")
pred2<-prediction(rpartPred[,2],ucecsmall2$msi)
envision_perf2 <- performance(pred2,"tpr","fpr")
envision_auc <- round(performance(pred2,"auc")@y.values[[1]], digits=3) #0.683
minauct <- paste(c("AUC = "),envision_auc,sep="")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_roc.pdf',height=7,width=7)
plot(envision_perf2,main='ROC',colorize=T) #0.705 for Envision, #0.830 for CADD, 0.783 for SIFT, 0.802 for PolyPhen
text(0.8,0.2,minauct,cex=1.5)
dev.off()

library('ROCR')
pred2<-prediction(functplotdat$envision,functplotdat$functional_effect)
envision_perf2 <- performance(pred2,"tpr","fpr")
envision_auc <- round(performance(pred2,"auc")@y.values[[1]], digits=3) #0.683
minauct <- paste(c("AUC = "),envision_auc,sep="")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/dms_analysis/working_results/pmd_prelim_roc_dms_weighted_classes_012815.pdf',height=7,width=7)
plot(envision_perf2,main='ROC',colorize=T) #0.705 for Envision, #0.830 for CADD, 0.783 for SIFT, 0.802 for PolyPhen
text(0.8,0.2,minauct,cex=1.5)
dev.off()

#figure 2a
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
msifracs<-tapply(sampledat$msi_status_tree_scaled,sampledat$tumor_type,function(x)length(which(x=='MSI-H'))/length(x))
bardat<-data.frame(cancer=names(msifracs)[-c(1,11)],prop=msifracs[-c(1,11)]*100)
bardat$cancer<-factor(bardat$cancer,levels=bardat$cancer[order(bardat$prop,decreasing=TRUE)])
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_prop_tree_scaled_070115_final.pdf',width=7*2,height=7)
ggplot(data=bardat, aes(x=cancer, y=prop)) + geom_bar(aes(fill=cancer), stat="identity") + guides(fill=FALSE) + scale_fill_manual(values=getPalette(18)) + theme_bw(base_size = 20) + ylab('Proportion of MSI individuals (%)') + coord_cartesian(ylim = c(0,50)) + xlab('')
dev.off()

#figure 2b
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
require(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols=brewer.pal(9, "Set1")
msifracs<-tapply(sampledat$msi_status_tree_scaled,sampledat$tumor_type,function(x)length(which(x=='MSI-H'))/length(x))
bardat<-data.frame(cancer=names(msifracs)[-c(1,11)],prop=msifracs[-c(1,11)]*100)
bardat$cancer<-factor(bardat$cancer,levels=bardat$cancer[order(bardat$prop,decreasing=TRUE)])
sampledat2<-sampledat[-which(sampledat$tumor_type==''),]
sampledat3<-sampledat2
sampledat2$tumor_type<-factor(sampledat2$tumor_type,levels=levels(bardat$cancer))

require(ggplot2)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/frac_unstable_per_cancer_120115.pdf',width=7*2,height=7)
ggplot(sampledat2,aes(y=log10(num_unstable_raw+0.5),x=tumor_type)) + geom_boxplot(aes(fill=factor(tumor_type)),outlier.colour = NA) + geom_point(position = position_jitter(width = 0.2), size=1.25, aes(color=msi_status_tree_scaled)) +scale_fill_manual(values=getPalette(18)) + theme_bw(base_size = 20) + guides(fill=FALSE) + xlab('') + ylab(expression("log"[10]*" number of unstable loci")) + scale_color_manual(values=c(cols[1],'grey25'),guide=FALSE)
dev.off()

#figure 3c
r<-results_allloci2
r$type<-'stable'
r$type[which(results_allloci2$prop>0)]<-'unstable'
r$genomic_class<-revalue(results_allloci2$genomic_class, c("downstream"="noncoding","exonic"="coding","ncRNA_exonic"="noncoding","ncRNA_intronic"="noncoding","ncRNA_splicing"="noncoding","ncRNA_UTR3"="noncoding","ncRNA_UTR5"="noncoding","ncRNA_UTR5;ncRNA_UTR3"="noncoding","upstream"="noncoding","upstream;downstream"="noncoding","UTR5;UTR3"="noncoding","intergenic"="noncoding","intronic"="noncoding"))
require(reshape2)
a<-dcast(r,genomic_class~type,length)
a$stable<-a$stable/table(r$type)[1]
b<-melt(a,id="genomic_class")

require(scales)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_genomic_classes_070115.pdf',width=7,height=7)
ggplot(b,aes(x = variable, y = value,fill = genomic_class)) + 
    geom_bar(position = "fill",stat="identity") + 
    scale_y_continuous(labels = percent_format()) + theme_bw(base_size = 20) + scale_fill_brewer(palette="Set1") + ylab("% of microsatellites") + xlab("") + theme(legend.title=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
dev.off()

#figure 3d
r<-results_allloci2
r$type<-'stable'
r$type[which(results_allloci2$prop>0)]<-'unstable'
require(reshape2)
a<-dcast(r,repeat_type~type,length)
a$stable<-a$stable/table(r$type)[1]
b<-melt(a,id="repeat_type")

require(scales)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_repeat_lengths_070115.pdf',width=7,height=7)
ggplot(b,aes(x = variable, y = value,fill = repeat_type)) + 
    geom_bar(position = "fill",stat="identity") + 
    scale_y_continuous(labels = percent_format()) + theme_bw(base_size = 20) + scale_fill_brewer(palette="Set1") + ylab("% of microsatellites") + xlab("") + theme(legend.title=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
dev.off()

#2e - cosmic, two questions: enrichment for cancer genes for sites that are unstable across all cancer types (net proportion unstable > 0.05 or for each cancer) vs. sampled matched sizes from background of sites with sufficient reads?
#enrichment for cancer genes for sites that are significantly unstable in MSI-H cancers vs. MSS cancers vs. sampled matched sizes from background of sites with sufficient reads for MSI-H vs. MSS comparison
require(plyr)
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_042515.robj')
output3<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_050315.txt',header=T)

locinfo<-results_msi[which(results_msi$msi_status=='MSS'),]
locinfo$genomic_class<-as.character(locinfo$genomic_class)
locinfo$genomic_class<-revalue(locinfo$genomic_class, c("downstream"="noncoding","exonic"="coding","ncRNA_exonic"="noncoding","ncRNA_intronic"="noncoding","ncRNA_splicing"="noncoding","ncRNA_UTR3"="noncoding","ncRNA_UTR5"="noncoding","ncRNA_UTR5;ncRNA_UTR3"="noncoding","upstream"="noncoding","upstream;downstream"="noncoding","UTR5;UTR3"="noncoding","intergenic"="noncoding","intronic"="noncoding"))
locinfo<-locinfo[-which(locinfo$genomic_class=='noncoding'),]

backgroundgenes<-as.character(output3$gene[which(output3$msi_peak_)]) #21553 unique genes

output3$genomic_class<-revalue(output3$genomic_class, c("downstream"="noncoding","exonic"="coding","ncRNA_exonic"="noncoding","ncRNA_intronic"="noncoding","ncRNA_splicing"="noncoding","ncRNA_UTR3"="noncoding","ncRNA_UTR5"="noncoding","ncRNA_UTR5;ncRNA_UTR3"="noncoding","upstream"="noncoding","upstream;downstream"="noncoding","UTR5;UTR3"="noncoding","intergenic"="noncoding","intronic"="noncoding"))
siglocinfo<-output3[which(output3$qval<0.05),]
siglocinfo2<-siglocinfo[-which(siglocinfo$genomic_class=='noncoding'),]
siggenes<-as.character(siglocinfo2$gene) #5173 unique genes
siggenes<-unlist(lapply(siggenes,function(x)strsplit(x,',')))
cosmic<-read.csv('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cosmic_061515.csv',header=T) #573 unique genes
cosmicgenes<-as.character(unique(cosmic$Gene.Symbol))
#cosmicgenes<-c(as.character(unique(cosmic$Gene.Symbol)),'ACVR2A')

length(which(is.na(match(unique(siggenes),cosmicgenes))==FALSE)) #216/5173 = 4.2% of genes have at least one MS site associated with global MSI
length(which(is.na(match(unique(backgroundgenes),cosmicgenes))==FALSE)) #478/21553 = 2.2% of genes have at least one MS site associated with global MSI;

set.seed(101)
genes<-c()
permgenes<-c()
for (i in 1:1000){
	genes[i]<-length(unique(sample(backgroundgenes,length(siggenes))))
	#permgenes[i]<-length(which(is.na(match(unique(sample(backgroundgenes,length(siggenes))),cosmicgenes))==FALSE))
	permgenes[i]<-length(which(is.na(match(sample(unique(backgroundgenes),length(unique(siggenes))),cosmicgenes))==FALSE))
}

require(RColorBrewer)
cols<-brewer.pal(9, "Set1")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cosmic_enrichment_070115.pdf',width=7,height=7)
ggplot(data.frame(permgenes),aes(x=permgenes)) + geom_histogram(fill = cols[9],binwidth=5) + geom_point(x=216,y=0,size=3,color=cols[1]) +theme_bw(base_size = 20) + ylab("count") + xlab("# of overlapping cancer genes") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_cartesian(xlim=c(100,225))
dev.off()

require(RColorBrewer)
cols<-brewer.pal(9, "Set1")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cosmic_enrichment_011216.pdf',width=7,height=7)
ggplot(data.frame(permgenes),aes(x=permgenes)) + geom_histogram(fill = cols[9],binwidth=5) + geom_point(x=216,y=0,size=3,color=cols[1]) +theme_bw(base_size = 20) + ylab("count") + xlab("# of overlapping cancer genes") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_cartesian(xlim=c(100,275))
dev.off()

length(which(is.na(match(unique(siggenes),cosmicgenes))==FALSE)) #216/5173 = 4.2% of genes have at least one MS site associated with global MSI
length(which(is.na(match(unique(backgroundgenes),cosmicgenes))==FALSE)) #478/21553 = 2.2% of genes have at least one MS site associated with global MSI