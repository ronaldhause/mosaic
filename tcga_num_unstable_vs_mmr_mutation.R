load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')

sampleinfo$num_unstable<-sampledat$num_unstable_raw
empties<-c('ND','NaN','',NA)
mutbi<-matrix(rep(0,4478*10),ncol=10)
for (i in 1:10){
	mutbi[which(is.na(match(sampleinfo[,i+8],empties))==T),i]<-1
}
colnames(mutbi)<-colnames(sampleinfo)[9:18]
rownames(mutbi)<-sampleinfo$SAMPLE_NAME

rowSums(mutbi)
cor.test(rowSums(mutbi),sampledat$num_unstable_raw)

mmrmut<-rep('no_mut',4478)
mmrmut[which(rowSums(mutbi)>0)]<-'mmr_mut'

sampleinfo$mmrmut<-mmrmut
sampleinfo<-sampleinfo[-which(sampleinfo$num_unstable>50000),]
sampledat<-sampledat[-which(sampledat$num_unstable_raw>50000),]
sampleinfo$msi_status<-sampledat$msi_status_tree_scaled

cancers<-levels(sampleinfo$TUMOR_TYPE)[-c(1,11)]
sig<-list()
for (i in 1:length(cancers)){
	subset<-sampleinfo[which(sampleinfo$TUMOR_TYPE==cancers[i]),]
	if (length(which(subset$mmrmut=='mmr_mut'))>9){
	sig[[i]]<-wilcox.test(subset$num_unstable~subset$mmrmut,alternative="greater")
	}
	else{sig[[i]]<-NA}
}
names(sig)<-cancers

require(ggplot2)
sampleinfo$mmrmut<-mmrmut
sampleinfo<-sampleinfo[-which(sampledat$tumor_type==''),]
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_num_unstable_mmr_mutation.pdf',width=7*2,height=7)
ggplot(data=sampleinfo, aes(x=TUMOR_TYPE,y=log10(num_unstable+0.5))) + geom_boxplot(aes(fill=factor(mmrmut))) + guides(fill=FALSE) + theme_bw(base_size = 20) + ylab(expression(paste(log[10],' # of MSI events',sep=''))) + xlab('cancer type') + scale_fill_brewer('Set1')
dev.off()

empties<-c('ND','NaN','',NA)
mutbi<-matrix(rep(0,4478*10),ncol=10)
for (i in 1:10){
	mutbi[which(is.na(match(sampleinfo[,i+8],empties))==T),i]<-1
}
colnames(mutbi)<-colnames(sampleinfo)[9:18]
rownames(mutbi)<-sampleinfo$SAMPLE_NAME

library(reshape)
plotdat<-data.frame(sample=sampleinfo$SAMPLE_NAME,num_unstable=sampledat$num_unstable_raw,num_unstable_ks=sampledat$num_unstable_ks,mutbi,mmrmut=sampleinfo$mmrmut)
plotdat2<-melt(plotdat,id=c('sample','num_unstable'))
plotdat2<-plotdat2[-which(plotdat2$num_unstable>50000),]