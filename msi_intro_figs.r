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

#figure 1b - regenerate
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

#alldat<-data.frame(genomic_class=c('intergenic','intronic','ncRNA_intronic','UTR3','downstream','upstream','ncRNA_exonic','exonic','UTR5','upstream;downstream','splicing','UTR5;UTR3','ncRNA_splicing'),totals=c(1038315,809970,107657,17859,14171,10976,3223,150,1329,370,19,4,4))
#alldat$genomic_class<-revalue(alldat$genomic_class,c("downstream"="noncoding","ncRNA_exonic"="noncoding","ncRNA_intronic"="noncoding","ncRNA_splicing"="noncoding","ncRNA_UTR3"="noncoding","ncRNA_UTR5"="noncoding","ncRNA_UTR5;ncRNA_UTR3"="noncoding","upstream"="noncoding","upstream;downstream"="noncoding","UTR5;UTR3"="noncoding","intergenic"="noncoding"))
#colnames(alldat)<-c('genomic_class','totals')

ggpie <- function (dat, by, totals) {
ggplot(dat, aes_string(x=factor(1), y=totals, fill=by)) +
geom_bar(stat='identity', color='black') +
scale_fill_brewer(palette="Set1") +
guides(fill=guide_legend(override.aes=list(colour=NA))) + 
coord_polar(theta='y') +
theme(axis.ticks=element_blank(),
axis.text.y=element_blank(),
axis.text.x=element_text(colour='black'),
axis.title=element_blank(),
legend.position="none") +
scale_y_continuous(breaks=cumsum(dat[[totals]]) - dat[[totals]] / 2, labels=dat[[by]]) +
theme(panel.background = element_rect(fill = "white"))
}

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_genomic_classes_020216.pdf',width=7,height=7)
ggpie(dat,"genomic_class","totals")
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_genomic_classes_global_020216.pdf',width=7,height=7)
ggpie(alldat,"genomic_class","totals")
dev.off()

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

cosmic<-read.csv('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cosmic_061515.csv',header=T) #573 unique genes
cosmicgenes<-c(as.character(cosmic$Gene.Symbol),'ACVR2A')

msih2<-msih[which(msih$PEAK_DIFFERENCE_VALUE>0),]
msih3<-msih2[which(is.na(match(msih2$GENE,cosmicgenes))==F),]
output3<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_050315.txt',header=T)
siglocinfo<-output3[which(output3$qval<0.05),]
msih4<-msih3[which(is.na(match(msih3$GENE,siglocinfo$gene))==F),]
msih5<-msih4[which(msih4$GENOMIC_CLASS=='exonic'),]
m<-merge(msih5,mss,by="LOCUS_COORDINATES")

#MSI-L vs. MSS calculations
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_training_data_sig_loc_incl_030315.robj')
ucecsmall<-ucecsmall[,1:6]
ucecsmall$msi_pred<-sampledat$msi_status_tree_scaled[match(ucecsmall$sample_name,sampledat$sample_name)]
ucecsmall2<-ucecsmall[which(ucecsmall$msi=='MSI-H' | ucecsmall$msi=='MSS'),]
ucecsmall2$msi<-factor(ucecsmall2$msi)
wilcox.test(ucecsmall$num_unstable_raw[which(ucecsmall$msi=='MSS')],ucecsmall$num_unstable_raw[which(ucecsmall$msi=='MSI-L')])
wilcox.test(ucecsmall$peak_avg[which(ucecsmall$msi=='MSS')],ucecsmall$peak_avg[which(ucecsmall$msi=='MSI-L')])
wilcox.test(ucecsmall$peak_sd[which(ucecsmall$msi=='MSS')]^2,ucecsmall$peak_sd[which(ucecsmall$msi=='MSI-L')]^2)