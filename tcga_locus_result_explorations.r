a<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_tcga_loci_results_042215.txt',header=T)
b<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch2_tcga_loci_results_042215.txt',header=F,skip=1,fill=T)
b[which(b[,1]=='MSS' | b[,1]=='MSI-H'),1:13]<-c('Omitted',b[which(b[,1]=='MSS' | b[,1]=='MSI-H'),1:12])
c<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch3_tcga_loci_results_042215.txt',header=F,skip=1,fill=T)
d<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch4_tcga_loci_results_042215.txt',header=F,skip=1,fill=T)
e<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch5_tcga_loci_results_042215.txt',header=F,skip=1,fill=T)
f<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch6_tcga_loci_results_042215.txt',header=T,fill=T)
colnames(b)<-colnames(a)
colnames(c)<-colnames(a)
colnames(d)<-colnames(a)
colnames(e)<-colnames(a)
results<-rbind(a,b,c,d,e,f)
#results[which(results[,1]=='MSS' | results[,1]=='MSI-H' | results[,1]=='<NA>'),1:13]<-c('Omitted',results[which(b[,1]=='MSS' | results[,1]=='MSI-H' | results[,1]=='<NA>'),1:12])
results[which(results$repeat_dna_sequence==''),1:13]<-c('Omitted',results[which(results$repeat_dna_sequence==''),1:12])
#
save(results,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_042415.robj')
results2<-results[-which(is.na(results[,1])==TRUE),]
results2$msi_status<-droplevels(results2$msi_status)
results2$total_num<-as.numeric(as.matrix(results2$total_num))
results<-as.data.table(results2)
save(results,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')

#ANALYSIS
#1 - instability landscape across all cancer types
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')

system.time(results_allloci<-results[,list(median_peak_diff=median(median_peak_diff,na.rm=T),ks_gain_site=sum(ks_gain_site,na.rm=T),simple_gain_site=sum(simple_gain_site,na.rm=T),num_missing=sum(num_missing,na.rm=T),num_present=sum(num_present,na.rm=T),total_num=sum(total_num,na.rm=T),genomic_class=genomic_class[1],gene=gene[1],repeat_type=repeat_type[1],repeat_dna_sequence=repeat_dna_sequence[1]),by="locus"])
save(results_allloci,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_042515.robj')
#516,876 sites, 4224 samples
results_allloci2<-results_allloci[which(results_allloci$num_present>(4224/2)),]
results_allloci2$prop<-((results_allloci2$simple_gain_site/results_allloci2$num_present)*100)
#222,357 sites have sufficient data in at least 2112 samples
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_042515.robj')

library(ggplot2)
bardat<-data.frame(locus=results_allloci2$locus,prop=(results_allloci2$simple_gain_site/results_allloci2$num_present)*100)
bardat$locus<-factor(bardat$locus,levels=bardat$locus[order(bardat$prop,decreasing=TRUE)])
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_locus_msi.pdf',width=7*2,height=7)
ggplot(data=bardat[1:1000,], aes(x=locus,y=prop)) + geom_bar(stat="identity") + guides(fill=FALSE) + theme_bw(base_size = 20) + ylab('Proportion of samples (%)') + coord_cartesian(ylim = c(0,40)) + xlab('loci') + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_locus_msi_hist_050315.pdf',width=7*2,height=7)
ggplot(results_allloci2, aes(x=(results_allloci2$simple_gain_site/results_allloci2$num_present)*100)) + geom_histogram() + theme_bw(base_size = 20) + ylab('frequency') + xlab('proportion_of_samples')
dev.off()

r<-results_allloci2
r$type<-'stable'
r$type[which(results_allloci2$prop>0)]<-'unstable'
r$genomic_class<-droplevels(r$genomic_class)
table(r$genomic_class,r$type)
library(reshape2)
a<-dcast(r,genomic_class~type,length)
a$stable<-a$stable/table(r$type)[1]
a$stable<-a$stable/table(r$type)[1]
b<-melt(a,id="genomic_class")

library(scales)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_genomic_classes_042515.pdf',width=7,height=7)
colourCount = length(unique(b$genomic_class))
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
ggplot(b,aes(x = variable, y = value,fill = genomic_class)) + 
    geom_bar(position = "fill",stat="identity") + 
    scale_y_continuous(labels = percent_format()) + theme_bw(base_size = 20) + scale_fill_manual(values=getPalette(colourCount))
dev.off()

output1<-results_allloci2[order(results_allloci2$prop,decreasing=TRUE),]
write.table(output[,c(1,4,6,7,12,8:11),with=FALSE],'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_050315.txt',quote=F,row.names=F)

#2 - instability landscape of specific cancer types
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')
library(data.table)
results_per_cancer<-results[,list(median_peak_diff=median(median_peak_diff,na.rm=T),ks_gain_site=sum(ks_gain_site,na.rm=T),simple_gain_site=sum(simple_gain_site,na.rm=T),num_missing=sum(num_missing,na.rm=T),num_present=sum(num_present,na.rm=T),total_num=sum(total_num,na.rm=T),genomic_class=genomic_class[1],gene=gene[1],repeat_type=repeat_type[1],repeat_dna_sequence=repeat_dna_sequence[1]),by=list(locus,tumor_type)]
results_per_cancer$tumor_type<-droplevels(results_per_cancer$tumor_type)
save(results_per_cancer,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_per_cancer_042515.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_per_cancer_042515.robj')

require(reshape2)
results_per_cancer2<-results_per_cancer[-which(results_per_cancer$num_present<(results_per_cancer$total_num/2)),]
results_per_cancer3<-split(results_per_cancer2,results_per_cancer2$tumor_type)

results_per_cancer4 = Map(function(x, i) setNames(x, ifelse(names(x) %in% "locus",
      names(x), sprintf('%s.%d', names(x), i))), results_per_cancer3, seq_along(results_per_cancer3))
results_per_cancer5<-Reduce(function(...) merge(..., by="locus",all=FALSE),results_per_cancer4)

output3<-data.frame(locus=results_msi5$locus,msi_peak_unstable=results_msi5$simple_gain_site.1,msi_locus_calls=results_msi5$num_present.1,mss_peak_unstable=results_msi5$simple_gain_site.1.2,mss_locus_calls=results_msi5$num_present.1.2,pval=fisher_pvals$V1,qval=fisher_pvals$qvals,odds_ratio=fisher_odds$V1,genomic_class=results_msi5$genomic_class.1,gene=results_msi5$gene.1,repeat_type=results_msi5$repeat_type.1,repeat_dna_sequence=results_msi5$repeat_dna_sequence.1)
a<-c(5,7)
b<-rep(11,18*2)*rep(seq(0,17),each=2)
c<-a+b
output2<-results_per_cancer5[,c(1,c,196:199)]
d<-rep(names(results_per_cancer3),each=2)
e<-rep(c('num_unstable','num_called'),18)
f<-paste(d,e,sep='_')
colnames(output2)<-c('locus',f,'genomic_class','gene','repeat_type','repeat_dna_sequence')

dat3<-data.frame()
output2<-results_per_cancer[order(results_per_cancer$prop,decreasing=TRUE),]
write.table(output[,c(1,4,6,7,12,8:11),with=FALSE],'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_050315.txt',quote=F,row.names=F)

canceroutput<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_050315.txt',header=T)

#3 - core loci associated with MSI
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')
library(data.table)
results_msi<-results[,list(median_peak_diff=median(median_peak_diff,na.rm=T),ks_gain_site=sum(ks_gain_site,na.rm=T),simple_gain_site=sum(simple_gain_site,na.rm=T),num_missing=sum(num_missing,na.rm=T),num_present=sum(num_present,na.rm=T),total_num=sum(total_num,na.rm=T),genomic_class=genomic_class[1],gene=gene[1],repeat_type=repeat_type[1],repeat_dna_sequence=repeat_dna_sequence[1]),by=list(locus,msi_status)]
save(results_msi,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_042515.robj')
#4102 MSS samples and 122 MSI-H samples

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_042515.robj')

#'2:175246541-175246561'
library(reshape2)
results_msi2<-results_msi[-which(results_msi$num_present<(results_msi$total_num/2)),]
#DF$seq <- with(DF, ave(Value, ID, msi_status, FUN = seq_along))
results_msi3<-split(results_msi2,results_msi2$msi_status)

results_msi4 = Map(function(x, i) setNames(x, ifelse(names(x) %in% "locus",
      names(x), sprintf('%s.%d', names(x), i))), results_msi3, seq_along(results_msi3))
results_msi5<-Reduce(function(...) merge(..., by="locus",all=FALSE),results_msi4)
#209,459 sites have sufficient data for both groups

fisher.test(matrix(c(85,121-85,631,3879-631),nrow=2,byrow=T))

system.time(fisher_odds<-results_msi5[,fisher.test(matrix(c(simple_gain_site.1,num_present.1-simple_gain_site.1,simple_gain_site.1.2,num_present.1.2-simple_gain_site.1.2),nrow=2,byrow=T))$estimate,by=locus]) #246
setkey(results_msi5,locus)
system.time(fisher_pvals<-results_msi5[,fisher.test(matrix(c(simple_gain_site.1,num_present.1-simple_gain_site.1,simple_gain_site.1.2,num_present.1.2-simple_gain_site.1.2),nrow=2,byrow=T))$p.value,by=locus])
source('/net/shendure/vol1/home/hauser/Scripts/useful_functions.r')
fisher_pvals[,qvals:=qvalue2(fisher_pvals$V1)$qvalues]

output3<-data.frame(locus=results_msi5$locus,msi_peak_unstable=results_msi5$simple_gain_site.1,msi_locus_calls=results_msi5$num_present.1,mss_peak_unstable=results_msi5$simple_gain_site.1.2,mss_locus_calls=results_msi5$num_present.1.2,pval=fisher_pvals$V1,qval=fisher_pvals$qvals,odds_ratio=fisher_odds$V1,genomic_class=results_msi5$genomic_class.1,gene=results_msi5$gene.1,repeat_type=results_msi5$repeat_type.1,repeat_dna_sequence=results_msi5$repeat_dna_sequence.1)
output3<-output3[order(output3$qval,decreasing=FALSE),]
write.table(output3,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_050315.txt',quote=F,row.names=F)

output3<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_050315.txt',header=T)
write.csv(output3,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_050315.csv',row.names=F)

bitmap('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_locus_msi_vs_mss_qqplot_032316.tiff',width=7*1.25,height=7,units="in",type="tifflzw",res=300)
gg_qqplot(output3$pval)
dev.off()

#testing for significant enrichment in cancer genes
locinfo<-results_msi[which(results_msi$msi_status=='MSS'),]
locinfo<-locinfo[-which(locinfo$genomic_class=='intronic'),]
locinfo<-locinfo[-which(locinfo$genomic_class=='intergenic'),]
backgroundgenes<-as.character(locinfo$gene) #21553 unique genes
#backgroundgenes<-unique(unlist(lapply(backgroundgenes,function(x)strsplit(x,','))))

siglocinfo<-output3[which(output3$qval<0.05),]
siglocinfo2<-siglocinfo[-which(siglocinfo$genomic_class=='intronic'),]
siglocinfo2<-siglocinfo2[-which(siglocinfo2$genomic_class=='intergenic'),]
siggenes<-as.character(siglocinfo2$gene) #5173 unique genes
cosmic<-read.csv('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cosmic_061515.csv',header=T) #573 unique genes
cosmicgenes<-c(as.character(unique(cosmic$Gene.Symbol)),'ACVR2')

length(which(is.na(match(unique(siggenes),cosmicgenes))==FALSE)) #216/5173 = 4.2% of genes have at least one MS site associated with global MSI
length(which(is.na(match(unique(backgroundgenes),cosmicgenes))==FALSE)) #478/21553 = 2.2% of genes have at least one MS site associated with global MSI (2.3%)

#giving genes GO annotations
library(topGO)
library(biomaRt)
mart <- useMart( "ensembl", dataset="hsapiens_gene_ensembl")
get.go.biomart <- getGO(id=ipi.LL.sym.descrip.ug.unique_eg[,2],type="entrezgene",mart=mart)

library(org.Hs.eg.db)
genes=rep(0,length(unique(backgroundgenes)))
genes[match(siggenes,unique(backgroundgenes))]<-1
names(genes)<-unique(backgroundgenes)
godata <- new("topGOdata", ontology = "BP", allGenes=as.factor(genes), annot = annFUN.org, mapping="org.Hs.eg.db", ID = "symbol" )
resultFisher <- runTest(godata, algorithm = "classic", statistic = "fisher")

biocLite("mygene")
library(mygene)
res <- queryMany(siglocinfo, scopes='symbol', fields=c('entrezgene', 'go'), species='human')

fisher.test(matrix(c(216,5173-216,478,21553-478),nrow=2,byrow=T))

#4 - cancer-specific loci associated with MSI
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')
library(data.table)
results_msi_per_cancer<-results[,list(median_peak_diff=median(median_peak_diff,na.rm=T),ks_gain_site=sum(ks_gain_site,na.rm=T),simple_gain_site=sum(simple_gain_site,na.rm=T),num_missing=sum(num_missing,na.rm=T),num_present=sum(num_present,na.rm=T),total_num=sum(total_num,na.rm=T),genomic_class=genomic_class[1],gene=gene[1],repeat_type=repeat_type[1],repeat_dna_sequence=repeat_dna_sequence[1]),by=list(locus,msi_status,tumor_type)]
save(results_msi_per_cancer,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_per_cancer_042515.robj')

plotdat<-data.frame(locus=output2$locus,coad=output2$COAD_num_unstable/output2$COAD_num_called,ucec=output2$UCEC_num_unstable/output2$UCEC_num_called,read=output2$READ_num_unstable/output2$READ_num_called,stad=output2$STAD_num_unstable/output2$STAD_num_called)
plotdat2<-data.frame(locus=plotdat$locus,x=plotdat$read-plotdat$stad,y=plotdat$coad-plotdat$ucec,number=rowSums(output2[,c(2,4,6,30),with=FALSE]))
require(tpsmeta)
angle<-1.75*pi #rotating 45 degrees clockwise
plotdat3<-data.frame(locus=plotdat2$locus,x=cos(angle)*plotdat2$x-sin(angle)*plotdat2$y,y=sin(angle)*plotdat2$x+cos(angle)*plotdat2$y,number=rowSums(output2[,c(2,4,6,30),with=FALSE]))

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_unstable_sites_before_rotation.pdf',height=7,width=7)
ggplot(plotdat2,aes(x=x,y=y)) + geom_point(alpha=0.5) + xlab("") + ylab("") + theme_bw(base_size = 20) + coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_unstable_sites_after_rotation.pdf',height=7,width=7)
ggplot(plotdat3,aes(x=x,y=y)) + geom_point(alpha=0.5) + xlab("") + ylab("") + theme_bw(base_size = 20) + coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5)) + geom_hline(aes(yintercept=0),linetype='dashed')
dev.off()

library(grid)
my_grob = grobTree(textGrob("UCEC", x=-0.4,  y=0.4, hjust=0,
  gp=gpar(col="black", fontsize=12)))

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_unstable_sites_after_rotation_sized_060815.pdf',height=7,width=7)
ggplot(plotdat3,aes(x=x,y=y)) + geom_point(aes(size=number),alpha=0.75) + xlab("") + ylab("") + theme_bw(base_size = 20) + scale_size_continuous(range = c(0, 3)) + coord_cartesian(xlim=c(-0.5,0.5),ylim=c(-0.5,0.5)) + geom_hline(aes(yintercept=0),linetype='dashed') + geom_vline(aes(xintercept=0),linetype='dashed') + annotation_custom(my_grob)
dev.off()