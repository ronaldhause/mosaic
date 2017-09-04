source('/net/shendure/vol1/home/hauser/scripts/shendure/miscellaneous/useful_functions.r')

#load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
#load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
#load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_post_review_070816.robj')
#516,876 sites, 5938 samples

#figure 1a
require(ggplot2)
require(reshape)
require(plyr)
source('/net/shendure/vol1/home/hauser/scripts/shendure/miscellaneous/useful_functions.r')
dat<-read.csv('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sample_electropherogram.csv',header=T)
a<-melt(dat[,c(1:3,6:7)],id="position")
a$variable<-revalue(a$variable,c("tumor.unstable"="tumor unstable","tumor.stable.normalized"="tumor stable","normal.unstable"="normal unstable","normal.stable.normalized"="normal stable"))
a$variable<-factor(a$variable,levels=c("normal stable","normal unstable","tumor stable","tumor unstable"))
a$type<-rep('tumor',dim(a)[1])
a$type[grep('normal',a$variable)]<-'normal'
library(RColorBrewer)
cols<-brewer.pal(9, "Set1")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sample_virtual_electropherogram_test.pdf',width=7*1.25,height=7)
ggplot(data=a, aes(x=position,y=value,fill=type)) + geom_bar(stat="identity") + guides(fill=FALSE) + theme_bw(base_size = 20) + ylab('Fraction relative to major allele (%)') + xlab('Variant length relative to major allele') + coord_cartesian(ylim = c(0,1.1)) + facet_wrap(~variable) + theme(strip.background = element_blank()) + scale_fill_manual(values=cols[3:4]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#figure 1b
require(ggplot2)
require(plyr)
genomic_class_cats<-revalue(results_allloci$genomic_class, c("downstream"="noncoding","ncRNA_exonic"="noncoding","ncRNA_intronic"="noncoding","ncRNA_splicing"="noncoding","ncRNA_UTR3"="noncoding","ncRNA_UTR5"="noncoding","ncRNA_UTR5;ncRNA_UTR3"="noncoding","upstream"="noncoding","upstream;downstream"="noncoding","UTR5;UTR3"="noncoding","intergenic"="noncoding"))
dat<-as.data.frame(table(genomic_class_cats))
colnames(dat)<-c('genomic_class','totals')

alldata<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/full_genome_v2.annovar.variant_function.gz',header=F)
alldata$bp<-paste(alldata[,4],alldata[,5],sep='-')
alldata$locus<-paste(alldata[,3],alldata$bp,sep=':')
a1<-alldata
alldata$V1<-as.character(alldata$V1)
alldata$V1[match(results_allloci$locus,alldata$locus)]<-as.character(results_allloci$genomic_class)
alldata2<-alldata[-grep('[[:alpha::]]*.',as.character(alldata$V3)),]
#19075236 - 39634 from non-autosomes 
#19035602 loci in alldat, 516876 in exome

all_genomic_class_cats<-revalue(alldata2[,1], c("downstream"="noncoding","ncRNA_exonic"="noncoding","ncRNA_intronic"="noncoding","ncRNA_splicing"="noncoding","ncRNA_UTR3"="noncoding","ncRNA_UTR5"="noncoding","ncRNA_UTR5;ncRNA_UTR3"="noncoding","upstream"="noncoding","upstream;downstream"="noncoding","UTR5;UTR3"="noncoding","intergenic"="noncoding"))
alldat<-as.data.frame(table(all_genomic_class_cats))
colnames(alldat)<-c('genomic_class','totals')

require(reshape2)
require(plyr)
a<-merge(dat,alldat,by="genomic_class")
a[,2]<-a[,2]/sum(a[,2])
a[,3]<-a[,3]/sum(a[,3])
a<-melt(a,id.vars="genomic_class")
a$variable<-revalue(a$variable,c("totals.x"="exome","totals.y"="whole genome"))

require(scales)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_genomic_classes_020216.pdf',width=7,height=7)
ggplot(a,aes(x = variable, y = value,fill = genomic_class)) +
    geom_bar(position = "fill",stat="identity") + 
    scale_y_continuous(labels = percent_format()) + theme_bw(base_size = 20) + scale_fill_brewer(palette="Set1") + ylab("% of microsatellites") + xlab("") + theme(legend.title=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank())
dev.off()

#figure 1e - histograms of sample repeat lengths for MSI-H and MSS cancer
system.time(load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated_msi_values_converted_to_numeric_nans_corrected.robj')) #2 min, 21450540 x 20
require(data.table)
setkey(data,SAMPLE_NAME)
msih<-data["TCGA-NH-A5IV"]
mss<-data["TCGA-T9-A92H"]

msih<-data["TCGA-QG-A5Z2"]
mss<-data["TCGA-RU-A8FL"]

d<-rbind(msih,mss)
d$TYPE<-c(rep('MSI-H',dim(msih)[1]),rep('MSS',dim(msih)[1]))
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sample_peak_diff_histogram_020216.pdf',width=7,height=7)
ggplot(d,aes(x=PEAK_DIFFERENCE_VALUE,fill=TYPE,log10="y")) +
  geom_histogram(binwidth=1) +
  facet_wrap(~ TYPE, ncol=1) +
  scale_fill_brewer(palette="Set1") +
  scale_y_sqrt() +
  theme_bw()
dev.off()

require(plotrix)
e<-table(d$PEAK_DIFFERENCE_VALUE,d$TYPE)
e<-e[-c(1:4),]
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sample_peak_diff_histogram_gap_020216.pdf',width=7*0.5,height=7)
par(mfrow=c(2,1))
gap.barplot(e[,1], gap=c(7000,319900),xaxlab=rownames(e),col=rep(cols[1],10),ytics=c(seq(0,7000,1000),320000),xlab='allele difference in tumor relative to normal',ylab="number of microsatellites",main="TCGA-NH-A5IV (MSI-H)",cex.lab=0.5,cex.axis=0.5,cex.main=0.5)
gap.barplot(e[,2], gap=c(7000,339900),xaxlab=rownames(e),col=rep(cols[2],10),ytics=c(seq(0,7000,1000),340000),xlab='allele difference in tumor relative to normal',ylab="number of microsatellites",main="TCGA-T9-A92H (MSS)",cex.lab=0.5,cex.axis=0.5,cex.main=0.5)
dev.off()

#figure 1f and supplementary figure 1
load('/net/shendure/vol12/projects/msi_tcga/mosaic_training_data_062916.robj')
require(ggplot2)
require(scales)
mosaic_training$defbsite<-gsub("NA","ND",mosaic_training$defbsite)
ucecsmall2 <- mosaic_training[-(which(mosaic_training$msi=="ND")),]
ucecsmall2 <- ucecsmall2[-(which(ucecsmall2$msi=="MSI-L")),]
ucecsmall3<-mosaic_training[which(mosaic_training$msi=='MSI-L'),]
ucecsmall3$msi<-factor(ucecsmall3$msi)
ucecsmall3$lab<-rep('MSS',dim(ucecsmall3)[1])
ucecsmall3<-ucecsmall3[,-1]
defbdat<-as.data.frame(table(ucecsmall2$defbsite,ucecsmall2$msi))
defbdat$Var1<-factor(defbdat$Var1, levels = rev(c("unstable", "stable", "ND")))
defbdat<-defbdat[seq(6,1),]
num_unstable3<-ucecsmall3$num_unstable
a<-ggplot(ucecsmall2,aes(y=ucecsmall2$peak_avg,x=ucecsmall2$msi)) + geom_boxplot() + geom_jitter(data=ucecsmall3,aes(y=ucecsmall3$peak_avg,x=ucecsmall3$lab,color='a')) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('global average allele number difference') + xlab('') + theme(axis.title=element_text(size=16)) + theme(legend.position="none") + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
b<-ggplot(ucecsmall2,aes(y=ucecsmall2$peak_sd,x=ucecsmall2$msi)) + geom_boxplot() + geom_jitter(data=ucecsmall3,aes(y=ucecsmall3$peak_sd,x=ucecsmall3$lab,color='a')) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('global sd in allele number difference') + xlab('') + theme(axis.title=element_text(size=16)) + theme(legend.position="none") + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
c<-ggplot(ucecsmall2,aes(y=num_unstable+0.5,x=ucecsmall2$msi)) + geom_boxplot() + geom_jitter(data=ucecsmall3,aes(y=num_unstable3+0.5,x=ucecsmall3$lab,color='a')) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + scale_y_log10() + ylab('global number of significantly unstable loci') + xlab('') + theme(axis.title=element_text(size=16)) + theme(legend.position="none") + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
d<-ggplot(defbdat,aes(y=Freq, x=Var2, fill=factor(Var1))) + geom_bar(position = "fill", stat = "identity") + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('stability status at 8:7679723-7679741') + xlab('') + theme(axis.title=element_text(size=16)) + theme(legend.position="right") + scale_y_continuous(labels = percent_format()) +  theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + labs(fill="")

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_comparison_summary_with_msil_070716.pdf',width=7*1.5,height=7*1.5)
multiplot(a,c,d,b, cols=2)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/fig1f_only_mss_msi_comparison_summary_with_msil_070716.pdf',width=7*1.5,height=7*1.5)
multiplot(a,d, cols=2)
dev.off()

#MSI-L vs. MSS calculations
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_training_data_sig_loc_incl_030315.robj')
ucecsmall<-ucecsmall[,1:6]
ucecsmall$msi_pred<-sampledat$msi_status_tree_scaled[match(ucecsmall$sample_name,sampledat$sample_name)]
ucecsmall2<-ucecsmall[which(ucecsmall$msi=='MSI-H' | ucecsmall$msi=='MSS'),]
ucecsmall2$msi<-factor(ucecsmall2$msi)
wilcox.test(ucecsmall$num_unstable_raw[which(ucecsmall$msi=='MSS')],ucecsmall$num_unstable_raw[which(ucecsmall$msi=='MSI-L')])
wilcox.test(ucecsmall$peak_avg[which(ucecsmall$msi=='MSS')],ucecsmall$peak_avg[which(ucecsmall$msi=='MSI-L')])
wilcox.test(ucecsmall$peak_sd[which(ucecsmall$msi=='MSS')]^2,ucecsmall$peak_sd[which(ucecsmall$msi=='MSI-L')]^2)

#misclassification figure
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_training_data_062916.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_classifier_063016.robj')
require(ggplot2)
require(scales)
mosaic_training$defbsite<-gsub("ND","NA",mosaic_training$defbsite)
ucecsmall2 <- mosaic_training[-(which(mosaic_training$msi=="ND")),]
ucecsmall2$X8.7679723.7679741<-ucecsmall2$defbsite
ucecsmall2$mosaic<-predict(mosaic, ucecsmall2, type="raw")
ucecsmall2$misclass<-rep(0,dim(ucecsmall2)[1])
ucecsmall2$misclass[which(ucecsmall2$msi!=ucecsmall2$mosaic)]<-1
ucecsmall2 <- ucecsmall2[-(which(ucecsmall2$msi=="MSI-L")),]
ucecsmall3<-ucecsmall2[which(ucecsmall2$misclass==1),]
ucecsmall3$msi<-factor(ucecsmall3$msi)
ucecsmall3$lab<-ucecsmall3$mosaic
ucecsmall3<-ucecsmall3[,-1]
num_unstable3<-ucecsmall3$num_unstable
a<-ggplot(ucecsmall2,aes(y=ucecsmall2$peak_avg,x=ucecsmall2$msi)) + geom_boxplot() + geom_jitter(data=ucecsmall3, color = 'blue',aes(y=ucecsmall3$peak_avg,x=ucecsmall3$msi)) + theme_bw(base_size = 20) + ylab('global average allele number difference') + xlab('') + theme(axis.title=element_text(size=16)) + theme(legend.position="none") + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
c<-ggplot(ucecsmall2,aes(y=num_unstable+0.5,x=ucecsmall2$msi)) + geom_boxplot() + geom_jitter(data=ucecsmall3,color = 'blue',aes(y=num_unstable3+0.5,x=ucecsmall3$msi))  + theme_bw(base_size = 20) + scale_y_log10() + ylab('global number of significantly unstable loci') + xlab('') + theme(axis.title=element_text(size=16)) + theme(legend.position="none") + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_msi_misclassifications_overlaid_072116.pdf',width=7*1.75,height=7)
multiplot(a,c,cols=2)
dev.off()