source('/net/shendure/vol1/home/hauser/scripts/shendure/miscellaneous/useful_functions.r')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')

system.time(results_allloci<-results[,list(median_peak_diff=median(median_peak_diff,na.rm=T),num_unstable=sum(num_unstable,na.rm=T),num_missing=sum(num_missing,na.rm=T),num_present=sum(num_present,na.rm=T),total_num=sum(total_num,na.rm=T),genomic_class=genomic_class[1],gene=gene[1],repeat_type=repeat_type[1],repeat_dna_sequence=repeat_dna_sequence[1]),by="locus"])
save(results_allloci,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_post_review_070816.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_post_review_070816.robj')
#516,876 sites, 5930 samples
results_allloci$genomic_class<-factor(results_allloci$genomic_class)
results_allloci2<-results_allloci[which(results_allloci$num_present>(5930/2)),]
results_allloci2$prop<-((results_allloci2$num_unstable/results_allloci2$num_present)*100)
#223,082 sites have sufficient data in at least 2965 samples

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_small_msi_mss_training_data_sig_loc_incl_030315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_classifier_063016.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mosaic_training_data_062916.robj')

require(caret)
ucecsmall2<-mosaic_training
ucecsmall2$msi<-gsub('MSI-L','MSS',ucecsmall2$msi)
ucecsmall2<-ucecsmall2[-which(ucecsmall2$msi=='ND'),]
ucecsmall2$defbsite<-gsub('ND','stable',ucecsmall2$defbsite)
colnames(ucecsmall2)[13]<-'X8.7679723.7679741'
rpartPred<-predict(mosaic,ucecsmall2,type="raw")
confusionMatrix(rpartPred,ucecsmall2$msi)
#97.3% accuracy, 93.6% sensitivity, 98.5% specificity (1.54% FPR, 95.5% PPV)

mosaic_sub<-subset(mosaic_training, msi %in% c('MSS','MSI-H'))
output<-sampledat
output$msipcr<-mosaic_training$msi[match(sampledat$sample_name,mosaic_training$sample_name)]
output$training<-rep('no',dim(output)[1])
output$training[which(is.na(match(sampledat$sample_name,mosaic_sub$sample_name))==F)]<-'yes'
write.csv(output,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sample_information.csv',row.names=F)

#figure 2a
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
msifracs<-tapply(sampledat$msi,sampledat$tumor_type,function(x)length(which(x=='MSI-H'))/length(x))
bardat<-data.frame(cancer=names(msifracs),prop=msifracs*100)
bardat$cancer<-factor(bardat$cancer,levels=bardat$cancer[order(bardat$prop,decreasing=TRUE)])
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_prop_tree_scaled_070816_final.pdf',width=7*2,height=7)
ggplot(data=bardat, aes(x=cancer, y=prop)) + geom_bar(aes(fill=cancer), stat="identity") + guides(fill=FALSE) + scale_fill_manual(values=getPalette(18)) + theme_bw(base_size = 20) + ylab('Percentage of cancers classified as MSI-H') + coord_cartesian(ylim = c(0,30)) + xlab('') + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
dev.off()

#figure 2b
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')
require(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols=brewer.pal(9, "Set1")
msifracs<-tapply(sampledat$msi,sampledat$tumor_type,function(x)length(which(x=='MSI-H'))/length(x))
bardat<-data.frame(cancer=names(msifracs),prop=msifracs*100)
bardat$cancer<-factor(bardat$cancer,levels=bardat$cancer[order(bardat$prop,decreasing=TRUE)])
sampledat2<-sampledat
sampledat2$tumor_type<-factor(sampledat2$tumor_type,levels=levels(bardat$cancer))

require(ggplot2)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/frac_unstable_per_cancer_070816.pdf',width=7*2,height=7)
ggplot(sampledat2,aes(y=prop_unstable,x=tumor_type)) + geom_boxplot(aes(fill=factor(tumor_type)),outlier.colour = NA) + geom_point(position = position_jitter(width = 0.4), size=1.125, aes(color=msi, alpha=0.8)) + scale_fill_manual(values=getPalette(18)) + theme_bw(base_size = 20) + guides(fill=FALSE) + xlab('') + ylab("proportion of unstable loci") + scale_color_manual(values=c(cols[1],'grey25'),guide=FALSE) + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + theme(legend.position="none")
dev.off()

require(ggplot2)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/num_unstable_per_cancer_070816.pdf',width=7*2,height=7)
ggplot(sampledat2,aes(y=log10(num_unstable+0.5),x=tumor_type)) + geom_boxplot(aes(fill=factor(tumor_type)),outlier.colour = NA) + geom_point(position = position_jitter(width = 0.4), size=1.125, aes(color=msi, alpha=0.8)) +scale_fill_manual(values=getPalette(18)) + theme_bw(base_size = 20) + guides(fill=FALSE) + xlab('') + ylab(expression("log"[10]*" number of unstable loci")) + scale_color_manual(values=c(cols[1],'grey25'),guide=FALSE) + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + theme(legend.position="none")
dev.off()

#TO DO

#2e - cosmic, two questions: enrichment for cancer genes for sites that are unstable across all cancer types (net proportion unstable > 0.05 or for each cancer) vs. sampled matched sizes from background of sites with sufficient reads?
#enrichment for cancer genes for sites that are significantly unstable in MSI-H cancers vs. MSS cancers vs. sampled matched sizes from background of sites with sufficient reads for MSI-H vs. MSS comparison
require(plyr)
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_070816.robj')
output3<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_070816.txt', header=T)
#output3<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_050315.txt',header=T)

#testing for overlap with previously annotated recurrent mutational targets of colorectal cancer MSI
hits<-as.character(read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/recurrent_msi_targets.txt',header=F,stringsAsFactors=F)[,1])
hits<-unique(hits)
fisher.test(matrix(c(25,2,6882,18104-6882),byrow=T,nrow=2))
#25/27

locinfo<-results_msi[which(results_msi$msi_status=='MSS'),]
locinfo$genomic_class<-as.character(locinfo$genomic_class)
locinfo$genomic_class<-revalue(locinfo$genomic_class, c("downstream"="noncoding","exonic"="coding","ncRNA_exonic"="noncoding","ncRNA_intronic"="noncoding","ncRNA_splicing"="noncoding","ncRNA_UTR3"="noncoding","ncRNA_UTR5"="noncoding","ncRNA_UTR5;ncRNA_UTR3"="noncoding","upstream"="noncoding","upstream;downstream"="noncoding","UTR5;UTR3"="noncoding","intergenic"="noncoding","intronic"="noncoding"))
locinfo<-locinfo[-which(locinfo$genomic_class=='noncoding'),]

backgroundgenes<-unlist(lapply(as.character(locinfo$gene),function(x)strsplit(x,','))) #18104 unique genes

output3$genomic_class<-revalue(output3$genomic_class, c("downstream"="noncoding","exonic"="coding","ncRNA_exonic"="noncoding","ncRNA_intronic"="noncoding","ncRNA_splicing"="noncoding","ncRNA_UTR3"="noncoding","ncRNA_UTR5"="noncoding","ncRNA_UTR5;ncRNA_UTR3"="noncoding","upstream"="noncoding","upstream;downstream"="noncoding","UTR5;UTR3"="noncoding","intergenic"="noncoding","intronic"="noncoding"))
siglocinfo<-output3[which(output3$qval<0.05),]
siglocinfo2<-siglocinfo[-which(siglocinfo$genomic_class=='noncoding'),]
siggenes<-as.character(siglocinfo2$gene) #6821 unique genes
siggenes<-unlist(lapply(siggenes,function(x)strsplit(x,',')))
cosmic<-read.csv('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cosmic_061515.csv',header=T) #573 unique genes
cosmicgenes<-c(as.character(unique(cosmic$Gene.Symbol)),'ACVR2A')

length(which(is.na(match(unique(siggenes),cosmicgenes))==FALSE)) #272/6882 = 4.0% of genes have at least one MS site associated with global MSI
length(which(is.na(match(unique(backgroundgenes),cosmicgenes))==FALSE)) #475/18104 = 2.6% of genes have at least one MS site associated with global MSI

set.seed(101)
permgenes<-c()
test<-c()
for (i in 1:1000){
  permgenes[i]<-length(which(is.na(match(unique(sample(unique(backgroundgenes),length(unique(siggenes)))),cosmicgenes))==FALSE))
}

require(RColorBrewer)
cols<-brewer.pal(9, "Set1")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cosmic_enrichment_072716_2.pdf',width=7,height=7)
ggplot(data.frame(permgenes),aes(x=permgenes)) + geom_histogram(fill = cols[9],binwidth=5) + geom_point(x=272,y=0,size=3,color=cols[1]) +theme_bw(base_size = 20) + ylab("count") + xlab("Number of Genes") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_cartesian(xlim=c(150,275)) + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
dev.off()