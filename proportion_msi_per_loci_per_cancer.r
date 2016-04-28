system.time(load('/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/TCGA_subset_with_MSI_merged_full_analysis_msi_values_converted_to_numeric.robj')) #13 min

loc1<-data[which(data$LOCUS_COORDINATES=='10:100008321-100008335'),]
table(loc1$TUMOR_TYPE,loc1$MSI_STATUS)
#58/209 ucec are MSI-H (27.8%)
#2/68 read (2.9%)
#30/167 coad (18.0%)

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_prop_batch1.pdf',width=7,height=7)
cancer_msi<-data.frame(type=c('UCEC','COAD','READ'),prop=c(27.8,18.0,2.9))
ggplot(data=cancer_msi, aes(x=type, y=prop)) + geom_bar(aes(fill=type), stat="identity") + guides(fill=FALSE) + scale_color_brewer(palette="Set1") + theme_bw(base_size = 20) + ylab('Proportion of MSI individuals (%)') + coord_cartesian(ylim = c(0,50)) + xlab('')
dev.off()

setkey(data,TUMOR_TYPE)
ucec<-data["UCEC"]
read<-data["READ"]
coad<-data["COAD"]

setkey(ucec,LOCUS_COORDINATES)
setkey(read,LOCUS_COORDINATES)
setkey(coad,LOCUS_COORDINATES)
system.time(ucec_locus_nas<-ucec[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=LOCUS_COORDINATES]) #12s
system.time(read_locus_nas<-read[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=LOCUS_COORDINATES]) #12s
system.time(coad_locus_nas<-coad[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=LOCUS_COORDINATES]) #12s
system.time(ucec_sig_locs<-ucec[,length(which(KS_VALUE<0.05)),by=LOCUS_COORDINATES]) #4s
system.time(read_sig_locs<-read[,length(which(KS_VALUE<0.05)),by=LOCUS_COORDINATES]) #10s
system.time(coad_sig_locs<-coad[,length(which(KS_VALUE<0.05)),by=LOCUS_COORDINATES]) #10s
system.time(ucec_pos_locs<-ucec[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=LOCUS_COORDINATES]) #4s
system.time(read_pos_locs<-read[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=LOCUS_COORDINATES]) #10s
system.time(coad_pos_locs<-coad[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=LOCUS_COORDINATES]) #4s
#209 UCEC, 68 READ, 167 COAD

subtype<-data.table(locus=ucec_sig_locs$LOCUS_COORDINATES,ucec_sig_locs$V1,read_sig_locs$V1,coad_sig_locs$V1,ucec_pos_locs$V1,read_pos_locs$V1,coad_pos_locs$V1,ucec_locus_nas$V1,read_locus_nas$V1,coad_locus_nas$V1)
colnames(subtype)<-c('locus','ucec_samples_ks_sig','read_samples_ks_sig','coad_samples_ks_sig','ucec_peak_unstable','read_peak_unstable','coad_peak_unstable','ucec_locus_nas','read_locus_nas','coad_locus_nas')
write.table(msi_mss_locus_output,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_subset_with_msi_calls_loci_results.txt',quote=F,row.names=F)
