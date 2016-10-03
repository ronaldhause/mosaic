#skcm

require(data.table)
dat<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/extended_exome_results/Metastatic_SKCM/SKCM_metastatic_tumor_normal_comparisons/all.txt', colClasses=c('character','character','integer'), skip=1, col.names=c('sample','locus','peak_diff'), sep='\t', showProgress=T)
sample_names<-gsub('.msi_diff.txt','',dat$sample)
dat2<-data.table(GENOMIC_CLASS=NA, GENE=NA, REPEAT_TYPE=NA, REPEAT_DNA_SEQUENCE=NA, LOCUS_COORDINATES=dat$locus, SAMPLE_NAME=sample_names, TUMOR_TYPE=rep('SKCM_met', dim(dat)[1]), MSI_STATUS=NA, EXO1=NA, MLH1=NA, MLH3=NA, MSH2=NA, MSH3=NA, MSH6=NA, PMS1=NA, PMS2=NA, POLD1=NA, POLE=NA, KS_VALUE=NA, PEAK_DIFFERENCE_VALUE=dat$peak_diff)
skcm_met_norm<-dat2
save(skcm_met_norm,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/skcm_metastatic_normal.robj')

dat<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/extended_exome_results/Metastatic_SKCM/SKCM_primary_tumor_normal_comparisons/all.txt', colClasses=c('character','character','integer'), skip=1, col.names=c('sample','locus','peak_diff'), sep='\t', showProgress=T)
sample_names<-gsub('.msi_diff.txt','',dat$sample)
dat2<-data.table(GENOMIC_CLASS=NA, GENE=NA, REPEAT_TYPE=NA, REPEAT_DNA_SEQUENCE=NA, LOCUS_COORDINATES=dat$locus, SAMPLE_NAME=sample_names, TUMOR_TYPE=rep('SKCM_primary', dim(dat)[1]), MSI_STATUS=NA, EXO1=NA, MLH1=NA, MLH3=NA, MSH2=NA, MSH3=NA, MSH6=NA, PMS1=NA, PMS2=NA, POLD1=NA, POLE=NA, KS_VALUE=NA, PEAK_DIFFERENCE_VALUE=dat$peak_diff)
skcm_primary_norm<-dat2

require(data.table)
dat<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/extended_exome_results/metastatic_analyses/THCA_metastatic_vs_normal/all.txt', colClasses=c('character','character','integer'), skip=1, col.names=c('sample','locus','peak_diff'), sep='\t', showProgress=T)
sample_names<-gsub('.msi_diff.txt','',dat$sample)
dat2<-data.table(GENOMIC_CLASS=NA, GENE=NA, REPEAT_TYPE=NA, REPEAT_DNA_SEQUENCE=NA, LOCUS_COORDINATES=dat$locus, SAMPLE_NAME=sample_names, TUMOR_TYPE=rep('THCA_met_norm', dim(dat)[1]), MSI_STATUS=NA, EXO1=NA, MLH1=NA, MLH3=NA, MSH2=NA, MSH3=NA, MSH6=NA, PMS1=NA, PMS2=NA, POLD1=NA, POLE=NA, KS_VALUE=NA, PEAK_DIFFERENCE_VALUE=dat$peak_diff)
thca_met_norm<-dat2
dat<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/extended_exome_results/metastatic_analyses/THCA_primary_vs_normal/all.txt', colClasses=c('character','character','integer'), skip=1, col.names=c('sample','locus','peak_diff'), sep='\t', showProgress=T)
sample_names<-gsub('.msi_diff.txt','',dat$sample)
dat2<-data.table(GENOMIC_CLASS=NA, GENE=NA, REPEAT_TYPE=NA, REPEAT_DNA_SEQUENCE=NA, LOCUS_COORDINATES=dat$locus, SAMPLE_NAME=sample_names, TUMOR_TYPE=rep('THCA_primary_norm', dim(dat)[1]), MSI_STATUS=NA, EXO1=NA, MLH1=NA, MLH3=NA, MSH2=NA, MSH3=NA, MSH6=NA, PMS1=NA, PMS2=NA, POLD1=NA, POLE=NA, KS_VALUE=NA, PEAK_DIFFERENCE_VALUE=dat$peak_diff)
thca_primary_norm<-dat2
dat<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/extended_exome_results/metastatic_analyses/BRCA_metastatic_vs_normal/all.txt', colClasses=c('character','character','integer'), skip=1, col.names=c('sample','locus','peak_diff'), sep='\t', showProgress=T)
sample_names<-gsub('.msi_diff.txt','',dat$sample)
dat2<-data.table(GENOMIC_CLASS=NA, GENE=NA, REPEAT_TYPE=NA, REPEAT_DNA_SEQUENCE=NA, LOCUS_COORDINATES=dat$locus, SAMPLE_NAME=sample_names, TUMOR_TYPE=rep('BRCA_met_norm', dim(dat)[1]), MSI_STATUS=NA, EXO1=NA, MLH1=NA, MLH3=NA, MSH2=NA, MSH3=NA, MSH6=NA, PMS1=NA, PMS2=NA, POLD1=NA, POLE=NA, KS_VALUE=NA, PEAK_DIFFERENCE_VALUE=dat$peak_diff)
brca_met_norm<-dat2
dat<-fread('/net/shendure/vol7/stevesal/MSI_TCGA_project/extended_exome_results/metastatic_analyses/BRCA_primary_vs_normal/all.txt', colClasses=c('character','character','integer'), skip=1, col.names=c('sample','locus','peak_diff'), sep='\t', showProgress=T)
sample_names<-gsub('.msi_diff.txt','',dat$sample)
dat2<-data.table(GENOMIC_CLASS=NA, GENE=NA, REPEAT_TYPE=NA, REPEAT_DNA_SEQUENCE=NA, LOCUS_COORDINATES=dat$locus, SAMPLE_NAME=sample_names, TUMOR_TYPE=rep('BRCA_primary_norm', dim(dat)[1]), MSI_STATUS=NA, EXO1=NA, MLH1=NA, MLH3=NA, MSH2=NA, MSH3=NA, MSH6=NA, PMS1=NA, PMS2=NA, POLD1=NA, POLE=NA, KS_VALUE=NA, PEAK_DIFFERENCE_VALUE=dat$peak_diff)
brca_primary_norm<-dat2
metdat<-do.call(rbind, list(thca_met_norm,thca_primary_norm,brca_met_norm,brca_primary_norm))
save(metdat,file='/net/shendure/vol7/stevesal/MSI_TCGA_project/robjs/thca_brca_metastatic_results.robj')

#skcm_met_norm and metdat
skcm_primary_norm[,PEAK_DIFFERENCE_VALUE:=as.numeric(as.matrix(PEAK_DIFFERENCE_VALUE))]
skcm_met_norm[,PEAK_DIFFERENCE_VALUE:=as.numeric(as.matrix(PEAK_DIFFERENCE_VALUE))]
metdat[,PEAK_DIFFERENCE_VALUE:=as.numeric(as.matrix(PEAK_DIFFERENCE_VALUE))]
data<-rbindlist(list(skcm_primary_norm,skcm_met_norm,metdat))

system.time(peak_avg<-data[,mean(PEAK_DIFFERENCE_VALUE,na.rm=T),by=c("SAMPLE_NAME","TUMOR_TYPE")]) #14s
system.time(peak_sd<-data[,sd(PEAK_DIFFERENCE_VALUE,na.rm=T),by=c("SAMPLE_NAME","TUMOR_TYPE")]) #15s
system.time(num_unstable<-data[,length(which(PEAK_DIFFERENCE_VALUE>0)),by=c("SAMPLE_NAME","TUMOR_TYPE")]) #10s
system.time(num_na<-data[,length(which(is.na(PEAK_DIFFERENCE_VALUE)==T)),by=c("SAMPLE_NAME","TUMOR_TYPE")]) #10s

metsampledat<-data.frame(sample_name=peak_avg$SAMPLE_NAME, tumor_type=peak_avg$TUMOR_TYPE, peak_avg=peak_avg$V1, peak_sd=peak_sd$V1, num_unstable=num_unstable$V1, num_called = 516876-num_na$V1)
metsampledat$prop_unstable<-metsampledat$num_unstable/metasampledat$num_called
metsampledat$met<-rep('primary',dim(metsampledat)[1])
metsampledat$met[grep('met',metsampledat$tumor_type)]<-'metastasis'
metsampledat$cancer_type<-substr(metsampledat$tumor_type,1,4)

plotdat<-metasampledat[-grep('SKCM',metasampledat$tumor_type),]
plotdat$tumor_type<-factor(plotdat$tumor_type)
plotdat$sample_name<-factor(plotdat$sample_name,levels(plotdat$sample_name)[c(1:4,6:8,5,9:13)])

library(RColorBrewer)
cols<-brewer.pal(9, "Paired")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/thca_brca_met_norm_instability.pdf',width=7*1.5,height=7)
ggplot(data=plotdat, aes(x=sample_name,y=prop_unstable*100,fill=tumor_type)) + geom_bar(stat="identity",position="dodge") + scale_fill_manual(values=cols[c(2,1,4,3)]) + theme_bw(base_size = 20) + ylab('Proportion of sites unstable (%)') + coord_cartesian(ylim = c(0,0.8)) + theme(axis.title.x=element_blank(), axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

mod <- lm(prop_unstable ~ cancer_type + cancer_type/met, data=metsampledat)
anova(mod)

plotdat<-metasampledat

library(RColorBrewer)
cols<-brewer.pal(9, "Paired")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/metastatic_samples_num_unstable_three_cancers.pdf',width=7*1.25,height=7)
ggplot(metsampledat,aes(y=prop_unstable*100,x=met)) + geom_point(position = position_jitter(width = 0.4), size=1.5, aes(color=cancer_type)) + theme_bw(base_size = 20) + guides(fill=FALSE) + xlab('') + ylab("proportion of unstable loci") + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(), panel.grid.major=element_blank(), panel.border = element_blank(), panel.background = element_blank())
dev.off()