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