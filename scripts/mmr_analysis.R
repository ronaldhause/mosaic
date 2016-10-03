load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mmr_mutational_data.robj')
mmrmutdat[mmrmutdat==0]<-'ND'
mmrmutdat[mmrmutdat==1]<-'nsSNV'
sampleinfo$EXO1[which(sampleinfo$EXO1 %in% c('','ND','nan','NaN'))]<-mmrmutdat$EXO1[match(sampleinfo$SAMPLE_NAME,rownames(mmrmutdat))][which(sampleinfo$EXO1 %in% c('','ND','nan','NaN'))]
sampleinfo$EXO1[is.na(sampleinfo$EXO1)]<-'ND'
for (i in 10:18){
	sampleinfo[which(sampleinfo[,i] %in% c('','ND','nan','NaN')),i]<-mmrmutdat[,i-8][match(sampleinfo$SAMPLE_NAME,rownames(mmrmutdat))][which(sampleinfo[,i] %in% c('','ND','nan','NaN'))]
sampleinfo[,i][is.na(sampleinfo[,i])]<-'ND'
}

empties<-c('ND','NaN','',NA)

mutbi<-matrix(rep(0,dim(sampleinfo)[1]*10),ncol=10)
for (i in 1:10){
	mutbi[which(is.na(match(sampleinfo[,i+8],empties))==T),i]<-1
}
colnames(mutbi)<-colnames(sampleinfo)[9:18]
rownames(mutbi)<-sampleinfo$SAMPLE_NAME
mutbi<-as.data.frame(mutbi)

library(reshape)
plotdat<-data.frame(sample=sampleinfo$SAMPLE_NAME,prop_unstable=sampledat$prop_unstable,mutbi, MLH1_silencing=sampleinfo$MLH1_silencing)
plotdat$MLH1_silencing[is.na(plotdat$MLH1_silencing)]<-0
plotdat2<-melt(plotdat,id=c('sample','prop_unstable'))

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_mmr_mutation_msi_barchart_072116.pdf',width=7*2,height=7)
ggplot(data=plotdat2, aes(x=variable,y=prop_unstable*100)) + geom_boxplot(aes(fill=factor(value))) + guides(fill=FALSE) + theme_bw(base_size = 20) + ylab('% of MSI events') + xlab('MMR gene') + scale_fill_brewer('Set1') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"), panel.background = element_blank(), panel.border = element_blank())
dev.off()

rowSums(mutbi)
cor.test(rowSums(mutbi),sampledat$prop_unstable) #r=0.16, p<2.2e-16, rho = 0.0058, p = 0.66

mutbi2<-data.frame(mutbi,MLH1_silencing=sampleinfo$MLH1_silencing)
mutbi2$MLH1_silencing[is.na(mutbi2$MLH1_silencing)]<-0
mmrmut<-rep('no_mut',dim(sampleinfo)[1])
mmrmut[which(rowSums(mutbi2)>0)]<-'mmr_mut'
sampleinfo$mmrmut<-factor(sampleinfo$mmrmut,levels=c('no_mut','mmr_mut'))
sampleinfo$num_unstable<-sampledat$num_unstable
sampleinfo$prop_unstable<-sampledat$prop_unstable

sig<-list()
for (i in 1:length(cancers)){
	subset<-sampleinfo[which(sampleinfo$TUMOR_TYPE==cancers[i]),]
	if (length(which(subset$mmrmut=='mmr_mut'))>1){
	sig[[i]]<-wilcox.test(subset$prop_unstable~subset$mmrmut,alternative="less")
	}
	else{sig[[i]]<-NA}
}
names(sig)<-cancers

require(ggplot2)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_num_unstable_mmr_mutation_072116.pdf',width=7*2,height=7)
ggplot(data=sampleinfo, aes(x=TUMOR_TYPE,y=log10(num_unstable+0.5))) + geom_boxplot(aes(fill=factor(mmrmut))) + guides(fill=FALSE) + theme_bw(base_size = 20) + ylab(expression(paste(log[10],' # of MSI events',sep=''))) + xlab('cancer type') + scale_fill_brewer('Set1')
dev.off()

require(ggplot2)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_num_unstable_mmr_mutation_prop_unstable_072116.pdf',width=7*2,height=7)
ggplot(data=sampleinfo, aes(x=TUMOR_TYPE,y=prop_unstable)) + geom_boxplot(aes(fill=factor(mmrmut))) + guides(fill=FALSE) + theme_bw(base_size = 20) + ylab("% of MSI events") + xlab('cancer type') + scale_fill_brewer('Set1')
dev.off()

output<-melt(data.frame(MSI_STATUS=sampleinfo[,8],mutbi2),id.vars='MSI_STATUS')

require(reshape2)
output2<-dcast(count(output),variable + value ~ MSI_STATUS, value.var="freq")

k<-1
res<-list()
for (i in seq(1,22,2)){
	res[[k]]<-fisher.test(matrix(c(output2[i+1,3],output2[i,3],output2[i+1,4],output2[i,4]),nrow=2,byrow=T))
	k<-k+1
}