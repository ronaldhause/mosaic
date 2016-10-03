#lsmat = load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')

results_msi_per_cancer<-results

system.time(results_msi_per_cancer<-results[,list(median_peak_diff=median(median_peak_diff,na.rm=T),num_unstable=sum(num_unstable,na.rm=T),num_missing=sum(num_missing,na.rm=T),num_present=sum(num_present,na.rm=T),total_num=sum(total_num,na.rm=T),genomic_class=genomic_class[1],gene=gene[1],repeat_type=repeat_type[1],repeat_dna_sequence=repeat_dna_sequence[1]),by=c("locus","msi_status")])

sites_with_insufficient_data<-unique(results_msi_per_cancer$locus[which(results_msi_per_cancer$num_present<(0.5*results_msi_per_cancer$total_num))])

results_msi_per_cancer<-results_msi_per_cancer[!(locus %in% sites_with_insufficient_data)]
#204,797 sites

results_msi_per_cancer$prop<-results_msi_per_cancer$num_unstable/results_msi_per_cancer$num_present
results_msi_per_cancer$repeat_type<-factor(results_msi_per_cancer$repeat_type)
results_msi_per_cancer$repeat_type<-revalue(results_msi_per_cancer$repeat_type, c("p5"="p4"))
a<-with(results_msi_per_cancer,tapply(prop,list(msi_status,repeat_type),function(x)length(which(x>0.2))/length(x)))

b<-with(results_msi_per_cancer,tapply(prop,list(msi_status,repeat_type),function(x)length(x)))
b2<-with(results_msi_per_cancer,tapply(prop,list(msi_status,repeat_type),function(x)length(which(x>0.2))))

seq<-gsub('\\*','',results_msi_per_cancer$repeat_dna_sequence)
nums2<-as.numeric(gsub(".*([0-9]+)$","\\1",seq))
nums3<-as.numeric(gsub(".*([0-9][0-9])$","\\1",seq))
nums2[which(is.na(nums3)==F)]<-nums3[which(is.na(nums3)==F)]
results_msi_per_cancer$tract_repeat_number<-nums2

c<-with(results_msi_per_cancer,tapply(prop,list(msi_status,tract_repeat_number),function(x)median(x,na.rm=T)))

require(reshape2)
plotdat<-melt(c)
colnames(plotdat)<-c('msi_status','num','prop')

library(RColorBrewer)
cols<-brewer.pal(9, "Set1")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tract_repeat_number_vs_instability_summed_071316.pdf',width=7*1.5,height=7)
ggplot(data=plotdat, aes(x=num,y=prop*100,fill=msi_status)) + geom_bar(stat="identity",position="dodge") + scale_fill_manual(values=cols[1:2]) + theme_bw(base_size = 20) + ylab('Proportion of sites unstable (%)') + coord_cartesian(ylim = c(0,62)) + xlab('Repeat tract length') + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
dev.off()

require(reshape2)
plotdat2<-melt(a)
colnames(plotdat2)<-c('msi_status','type','prop')

library(RColorBrewer)
cols<-brewer.pal(9, "Set1")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/repeat_type_vs_instability_summed_071316.pdf',width=7*1.5,height=7)
ggplot(data=plotdat2, aes(x=type,y=prop*100,fill=msi_status)) + geom_bar(stat="identity",position="dodge") + scale_fill_manual(values=cols[1:2]) + theme_bw(base_size = 20) + ylab('Proportion of sites unstable (%)') + coord_cartesian(ylim = c(0,12)) + xlab('Repeat tract length') + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
dev.off()