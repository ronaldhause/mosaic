results_allloci$repeat_type<-factor(results_allloci$repeat_type)

prop<-results_allloci$simple_gain_site/results_allloci$num_present

a<-tapply(prop,results_allloci$repeat_type,function(x)median(x,na.rm=T))
b<-tapply(results_allloci$median_peak_diff,results_allloci$repeat_type,function(x)median(x,na.rm=T))

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_per_cancer_042515.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_080115.robj')

results_msi_per_cancer$prop<-results_msi_per_cancer$simple_gain_site/results_msi_per_cancer$num_present
results_msi_per_cancer$repeat_type<-factor(results_msi_per_cancer$repeat_type)
a<-with(results_msi_per_cancer,tapply(prop,list(msi_status,repeat_type),function(x)mean(x,na.rm=T)))
b<-with(results_msi_per_cancer,tapply(median_peak_diff,list(msi_status,repeat_type),function(x)mean(x,na.rm=T)))

output3$msihprop<-output3$msi_peak_unstable/output3$msi_locus_calls
output3$mssprop<-output3$mss_peak_unstable/output3$mss_locus_calls

a2<-with(output3,tapply(msihprop,list(repeat_type),function(x)length(which(x>0.2))/length(x)))
a3<-with(output3,tapply(mssprop,list(repeat_type),function(x)length(which(x>0.2))/length(x)))

a22<-with(output3,tapply(msi_peak_unstable,list(repeat_type),function(x)length(x)))
a32<-with(output3,tapply(mss_peak_unstable,list(repeat_type),function(x)length(x)))

seq<-gsub('\\*','',results_msi_per_cancer$repeat_dna_sequence)
nums2<-as.numeric(gsub(".*([0-9]+)$","\\1",seq))
nums3<-as.numeric(gsub(".*([0-9][0-9])$","\\1",seq))
nums2[which(is.na(nums3)==F)]<-nums3[which(is.na(nums3)==F)]
results_msi_per_cancer$tract_repeat_number<-nums2

seq<-gsub('\\*','',output3$repeat_dna_sequence)
nums2<-as.numeric(gsub(".*([0-9]+)$","\\1",seq))
nums3<-as.numeric(gsub(".*([0-9][0-9])$","\\1",seq))
nums2[which(is.na(nums3)==F)]<-nums3[which(is.na(nums3)==F)]
output3$tract_repeat_number<-nums2

c<-with(results_msi_per_cancer,tapply(prop,list(msi_status,tract_repeat_number),function(x)mean(x,na.rm=T)))
d<-with(results_msi_per_cancer,tapply(median_peak_diff,list(msi_status,tract_repeat_number),function(x)mean(x,na.rm=T)))

c2<-with(output3,tapply(msihprop,list(tract_repeat_number),function(x)length(which(x>0.2))/length(x)))
c3<-with(output3,tapply(mssprop,list(tract_repeat_number),function(x)length(which(x>0.2))/length(x)))

require(reshape2)
plotdat<-melt(c)
colnames(plotdat)<-c('msi_status','num','prop')

library(RColorBrewer)
cols<-brewer.pal(9, "Set1")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tract_repeat_number_vs_instability.pdf',width=7*1.5,height=7)
ggplot(data=plotdat, aes(x=num,y=prop*100,fill=msi_status)) + geom_bar(stat="identity",position="dodge") + scale_fill_manual(values=cols[1:2]) + theme_bw(base_size = 20) + ylab('Proportion of unstable sites (%)') + coord_cartesian(ylim = c(0,55)) + xlab('Repeat tract length')
dev.off()

require(reshape2)
plotdat2<-melt(a)
colnames(plotdat2)<-c('msi_status','type','prop')

library(RColorBrewer)
cols<-brewer.pal(9, "Set1")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/repeat_type_vs_instability.pdf',width=7*1.5,height=7)
ggplot(data=plotdat2, aes(x=type,y=prop*100,fill=msi_status)) + geom_bar(stat="identity",position="dodge") + scale_fill_manual(values=cols[1:2]) + theme_bw(base_size = 20) + ylab('Proportion of unstable sites (%)') + coord_cartesian(ylim = c(0,12)) + xlab('Repeat tract length')
dev.off()
