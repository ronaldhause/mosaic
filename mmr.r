load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')

empties<-c('ND','NaN','',NA)
mutbi<-vector(rep(0,4478))
mutbi[which(is.na(match(sampleinfo$EXO1,empties))==T)]<-1

mutbi<-matrix(rep(0,4478*10),ncol=10)
for (i in 1:10){
	mutbi[which(is.na(match(sampleinfo[,i+8],empties))==T),i]<-1
}
colnames(mutbi)<-colnames(sampleinfo)[9:18]
rownames(mutbi)<-sampleinfo$SAMPLE_NAME

library(reshape)
plotdat<-data.frame(sample=sampleinfo$SAMPLE_NAME,num_unstable=sampledat$num_unstable_raw,num_unstable_ks=sampledat$num_unstable_ks,mutbi)
plotdat2<-melt(plotdat,id=c('sample','num_unstable','num_unstable_ks'))
plotdat2<-plotdat2[-which(plotdat2$num_unstable>50000),]
plotdat2$num_unstable<-log10(plotdat2$num_unstable)

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_mmr_mutation_msi_barchart_061515.pdf',width=7*2,height=7)
ggplot(data=plotdat2, aes(x=variable,y=num_unstable)) + geom_boxplot(aes(fill=factor(value))) + guides(fill=FALSE) + theme_bw(base_size = 20) + ylab(expression(paste(log[10],' # of MSI events',sep=''))) + xlab('MMR gene') + scale_fill_brewer('Set1')
dev.off()

rowSums(mutbi)
cor.test(rowSums(mutbi),sampledat$num_unstable_ks) #0.05, p=0.0006; 0.004, p=0.792 for raw

msih<-which(sampledat$msi_status_tree_scaled=='MSI-H')
cor.test(rowSums(mutbi)[-msih],sampledat$num_unstable_ks[-msih],method='s') #0.05, p=0.0006; 0.004, p=0.792 for raw
plot(rowSums(mutbi)[-msih],sampledat$num_unstable_ks[-msih])