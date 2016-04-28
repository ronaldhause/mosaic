#1 - instability landscape across all cancer types
load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_duplicates_eliminated_msi_values_converted_to_numeric_nans_corrected_msi_predicted.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampledata_041315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_sampleinfo_030115.robj')
#415 samples

data<-data1
setkey(data,LOCUS_COORDINATES)
system.time(median_locus_peaks<-data[,median(PEAK_DIFFERENCE_VALUE,na.rm=T),by=LOCUS_COORDINATES]) #60s
system.time(locus_nas<-data[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=LOCUS_COORDINATES]) #10s
system.time(sig_locs<-data[,length(which(KS_VALUE<0.05)),by=LOCUS_COORDINATES]) #9s
system.time(pos_locs<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=LOCUS_COORDINATES]) #9s

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_locusinfo.robj')
a<-data.table(locus=sig_locs$LOCUS_COORDINATES,median_locus_peak=median_locus_peaks$V1,sig_locs$V1,pos_locs$V1,locus_nas$V1)
locus_output<-a[order(V4,decreasing=T),]
locus_output$locus_calls<-dim(sampledat1)[1]-locus_output$V5
colnames(locus_output)<-c('locus','median_peak_diff','ks_sig_site','gain_sig_site','num_missing','num_present')
locus_output$genomic_class<-locusinfo$GENOMIC_CLASS[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$gene<-locusinfo$GENE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_type<-locusinfo$REPEAT_TYPE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output$repeat_dna_sequence<-locusinfo$REPEAT_DNA_SEQUENCE[match(locus_output$locus,locusinfo$LOCUS_COORDINATES)]
locus_output<-locus_output[order(locus_output$gain_sig_site/locus_output$num_present,decreasing=T),]
write.table(locus_output,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_tcga_loci_results_041415.txt',quote=F,row.names=F)

sub<-which(locus_output$num_present>dim(sampledat1)[1]/2)
sublocs<-as.character(locus_output$locus)[sub]
locus_output2<-locus_output[sub,]

library(ggplot2)
bardat<-data.frame(locus=locus_output2$locus,prop=(locus_output2$gain_sig_site/locus_output2$num_present)*100)
bardat$locus<-factor(bardat$locus,levels=bardat$locus[order(bardat$prop,decreasing=TRUE)])
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_tcga_locus_msi.pdf',width=7*2,height=7)
ggplot(data=bardat[1:1000,], aes(x=locus,y=prop)) + geom_bar(stat="identity") + guides(fill=FALSE) + theme_bw(base_size = 20) + ylab('Proportion of MSI individuals (%)') + coord_cartesian(ylim = c(0,50)) + xlab('loci') + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/batch1_tcga_locus_msi_hist_041515.pdf',width=7*2,height=7)
ggplot(locus_output2, aes(x=(locus_output2$gain_sig_site/locus_output2$num_present)*100)) + geom_histogram() + theme_bw(base_size = 20) + ylab('frequency') + xlab('proportion_of_samples')
dev.off()

