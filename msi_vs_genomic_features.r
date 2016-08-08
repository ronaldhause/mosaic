load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_post_review_070816.robj')
results_allloci2<-results_allloci[which(results_allloci$num_present>(5930/2)),]
results_allloci2$prop<-((results_allloci2$num_unstable/results_allloci2$num_present)*100)

#making annotation, repeat type, and repeat seq identity log odds plots
require(plyr)
results_allloci2$genomic_class<-revalue(results_allloci2$genomic_class, c("downstream"="noncoding","ncRNA_exonic"="noncoding","ncRNA_intronic"="noncoding","ncRNA_splicing"="noncoding","ncRNA_UTR3"="noncoding","ncRNA_UTR5"="noncoding","ncRNA_UTR5;ncRNA_UTR3"="noncoding","upstream"="noncoding","upstream;downstream"="noncoding","UTR5;UTR3"="noncoding","intergenic"="noncoding"))

r<-as.data.frame(results_allloci2)
r$type<-'stable'
r$type[which(results_allloci2$prop>0)]<-'unstable'
r$genomic_class<-droplevels(r$genomic_class)
gclass<-table(r$genomic_class,r$type)

anno_test<-function(x,y){
	unstablehits<-length(intersect(which(r$type=='unstable'),which(r[,y]==x)))+1
	unstablemisses<-length(intersect(which(r$type=='unstable'),which(r[,y]!=x)))+1
	stablehits<-length(intersect(which(r$type=='stable'),which(r[,y]==x)))+1
	stablemisses<-length(intersect(which(r$type=='stable'),which(r[,y]!=x)))+1
	return(fisher.test(matrix(c(unstablehits,unstablemisses,stablehits,stablemisses),nrow=2,byrow=T)))
}

r$genomic_class<-droplevels(r$genomic_class)
gclass<-table(r$genomic_class,r$type)
results<-list()
features<-c()
mins<-c()
maxs<-c()
centers<-c()
for (i in 1:dim(gclass)[1]){
	dat<-anno_test(rownames(gclass)[i],'genomic_class')
	results[[i]]<-dat
	features[i]<-rownames(gclass)[i]
	mins[i]<-dat$conf.int[1]
	maxs[i]<-dat$conf.int[2]
	centers[i]<-dat$estimate
}
gclassdat<-data.frame(xlo=log2(mins),xhi=log2(maxs),x=log2(centers),y=features)

r$repeat_type<-droplevels(r$repeat_type)
rclass<-table(r$repeat_type,r$type)
results<-list()
features<-c()
mins<-c()
maxs<-c()
centers<-c()
for (i in 1:dim(rclass)[1]){
	dat<-anno_test(rownames(rclass)[i],'repeat_type')
	results[[i]]<-dat
	features[i]<-rownames(rclass)[i]
	mins[i]<-dat$conf.int[1]
	maxs[i]<-dat$conf.int[2]
	centers[i]<-dat$estimate
}
rclassdat<-data.frame(xlo=mins,xhi=maxs,x=centers,y=features)

p<-r[union(which(r$repeat_type=='p2'),which(r$repeat_type=='p1')),]
p$seq<-gsub(".*\\((.*)\\).*", "\\1", p$repeat_dna_sequence)
p$seq<-revalue(p$seq,c("GA"="GA","TC"="GA","A"="A/T","T"="A/T","C"="C/G","G"="C/G","TG"="CA","CA"="CA","GT"="CA","AC"="CA","CT"="GA","AG"="GA","TA"="AT","CG"="GC"))
sclass<-table(p$seq,p$type)
r<-p
results<-list()
features<-c()
mins<-c()
maxs<-c()
centers<-c()
for (i in 1:dim(sclass)[1]){
	dat<-anno_test(rownames(sclass)[i],'seq')
	results[[i]]<-dat
	features[i]<-rownames(sclass)[i]
	mins[i]<-dat$conf.int[1]
	maxs[i]<-dat$conf.int[2]
	centers[i]<-dat$estimate
}
sclassdat<-data.frame(xlo=log2(mins),xhi=log2(maxs),x=log2(centers),y=features)
sclassdat$y<-factor(sclassdat$y,levels=sclassdat$y[c(2,4,1,3,5,6)])

forestplot <- function(d, xlab="feature", ylab="log odds ratio"){
    require(ggplot2)
    p <- ggplot(d, aes(x=y, y=x, ymin=xlo, ymax=xhi)) + 
		geom_pointrange(size=1.25) + 
		geom_hline(yintercept=0, lty=2,size=1.25) +
		ylab(ylab) +
		xlab(xlab) +
		theme_bw(base_size = 20) +
		theme(axis.text.x=element_text(size=10)) +
		theme(panel.grid.minor = element_blank(), panel.border = element_blank())
    return(p)
}

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sequence_forestplot_071016.pdf',width=7*1.5,height=7)
forestplot(sclassdat,xlab="sequence")
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/genomic_class_forestplot_071016.pdf',width=7*1.5,height=7)
forestplot(gclassdat,xlab="genomic class")
dev.off()

#see msi_by_cancer_specific_figs.r
r<-msihresults
r$type<-'stable'
r$type[which(r$COAD>0)]<-'unstable'
p<-r[union(which(r$repeat_type=='p2'),which(r$repeat_type=='p1')),]
p$seq<-gsub(".*\\((.*)\\).*", "\\1", p$repeat_dna_sequence)
p$seq<-revalue(p$seq,c("GA"="GA","TC"="GA","A"="A/T","T"="A/T","C"="C/G","G"="C/G","TG"="CA","CA"="CA","GT"="CA","AC"="CA","CT"="GA","AG"="GA","TA"="AT","CG"="GC"))
sclass<-table(p$seq,p$type)
r<-p
results<-list()
features<-c()
mins<-c()
maxs<-c()
centers<-c()
for (i in 1:dim(sclass)[1]){
	dat<-anno_test(rownames(sclass)[i],'seq')
	results[[i]]<-dat
	features[i]<-rownames(sclass)[i]
	mins[i]<-dat$conf.int[1]
	maxs[i]<-dat$conf.int[2]
	centers[i]<-dat$estimate
}
coadclassdat<-data.frame(xlo=log2(mins),xhi=log2(maxs),x=log2(centers),y=features)
coadclassdat$y<-factor(coadclassdat$y,levels=coadclassdat$y[c(2,4,1,3,5,6)])

r<-msihresults
r$type<-'stable'
r$type[which(r$READ>0)]<-'unstable'
p<-r[union(which(r$repeat_type=='p2'),which(r$repeat_type=='p1')),]
p$seq<-gsub(".*\\((.*)\\).*", "\\1", p$repeat_dna_sequence)
p$seq<-revalue(p$seq,c("GA"="GA","TC"="GA","A"="A/T","T"="A/T","C"="C/G","G"="C/G","TG"="CA","CA"="CA","GT"="CA","AC"="CA","CT"="GA","AG"="GA","TA"="AT","CG"="GC"))
sclass<-table(p$seq,p$type)
r<-p
results<-list()
features<-c()
mins<-c()
maxs<-c()
centers<-c()
for (i in 1:dim(sclass)[1]){
	dat<-anno_test(rownames(sclass)[i],'seq')
	results[[i]]<-dat
	features[i]<-rownames(sclass)[i]
	mins[i]<-dat$conf.int[1]
	maxs[i]<-dat$conf.int[2]
	centers[i]<-dat$estimate
}
readclassddat<-data.frame(xlo=log2(mins),xhi=log2(maxs),x=log2(centers),y=features)
readclassddat$y<-factor(readclassddat$y,levels=readclassddat$y[c(2,4,1,3,5,6)])

r<-msihresults
r$type<-'stable'
r$type[which(r$UCEC>0)]<-'unstable'
p<-r[union(which(r$repeat_type=='p2'),which(r$repeat_type=='p1')),]
p$seq<-gsub(".*\\((.*)\\).*", "\\1", p$repeat_dna_sequence)
p$seq<-revalue(p$seq,c("GA"="GA","TC"="GA","A"="A/T","T"="A/T","C"="C/G","G"="C/G","TG"="CA","CA"="CA","GT"="CA","AC"="CA","CT"="GA","AG"="GA","TA"="AT","CG"="GC"))
sclass<-table(p$seq,p$type)
r<-p
results<-list()
features<-c()
mins<-c()
maxs<-c()
centers<-c()
for (i in 1:dim(sclass)[1]){
	dat<-anno_test(rownames(sclass)[i],'seq')
	results[[i]]<-dat
	features[i]<-rownames(sclass)[i]
	mins[i]<-dat$conf.int[1]
	maxs[i]<-dat$conf.int[2]
	centers[i]<-dat$estimate
}
ucecclassdat<-data.frame(xlo=log2(mins),xhi=log2(maxs),x=log2(centers),y=features)
ucecclassdat$y<-factor(ucecclassdat$y,levels=ucecclassdat$y[c(2,4,1,3,5,6)])

r<-msihresults
r$type<-'stable'
r$type[which(r$STAD>0)]<-'unstable'
p<-r[union(which(r$repeat_type=='p2'),which(r$repeat_type=='p1')),]
p$seq<-gsub(".*\\((.*)\\).*", "\\1", p$repeat_dna_sequence)
p$seq<-revalue(p$seq,c("GA"="GA","TC"="GA","A"="A/T","T"="A/T","C"="C/G","G"="C/G","TG"="CA","CA"="CA","GT"="CA","AC"="CA","CT"="GA","AG"="GA","TA"="AT","CG"="GC"))
sclass<-table(p$seq,p$type)
r<-p
results<-list()
features<-c()
mins<-c()
maxs<-c()
centers<-c()
for (i in 1:dim(sclass)[1]){
	dat<-anno_test(rownames(sclass)[i],'seq')
	results[[i]]<-dat
	features[i]<-rownames(sclass)[i]
	mins[i]<-dat$conf.int[1]
	maxs[i]<-dat$conf.int[2]
	centers[i]<-dat$estimate
}
stadclassdat<-data.frame(xlo=log2(mins),xhi=log2(maxs),x=log2(centers),y=features)
stadclassdat$y<-factor(stadclassdat$y,levels=stadclassdat$y[c(2,4,1,3,5,6)])

seqdat<-do.call(rbind,list(coadclassdat,readclassddat,ucecclassdat,stadclassdat))
seqdat$type<-c(rep('COAD',6),rep('READ',6),rep('UCEC',6),rep('STAD',6))

forestplot2 <- function(d, xlab="feature", ylab="log odds ratio"){
    require(ggplot2)
    p <- ggplot(d, aes(x=y, y=x, ymin=xlo, ymax=xhi)) + 
		geom_pointrange(aes(color=type),size=1.25,position = 
position_dodge(width = 0.65)) + 
		geom_hline(yintercept=0, lty=2,size=1.25) +
		ylab(ylab) +
		xlab(xlab) +
		theme_bw(base_size = 20) +
		theme(axis.text.x=element_text(size=10)) +
		scale_color_brewer(palette='Set1') +
		theme(panel.grid.minor = element_blank(),panel.border = element_blank())
    return(p)
}

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_sequence_forestplot_071216.pdf',width=7*1.75,height=7)
forestplot2(seqdat,xlab="sequence")
dev.off()

#MSI-H vs. all
output3<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_070816.txt', header=T)
r<-output3
r$type<-'stable'
r$type[which(r$STAD>0)]<-'unstable'
p<-r[union(which(r$repeat_type=='p2'),which(r$repeat_type=='p1')),]
p$seq<-gsub(".*\\((.*)\\).*", "\\1", p$repeat_dna_sequence)
p$seq<-revalue(p$seq,c("GA"="GA/TC","TC"="GA/TC","A"="A/T","T"="A/T","C"="C/G","G"="C/G","TG"="TG/CA","CA"="TG/CA","GT"="AC/GT","AC"="AC/GT","CT"="AG/CT","AG"="AG/CT"))
sclass<-table(p$seq,p$type)
r<-p
results<-list()
features<-c()
mins<-c()
maxs<-c()
centers<-c()
for (i in 1:dim(sclass)[1]){
	dat<-anno_test(rownames(sclass)[i],'seq')
	results[[i]]<-dat
	features[i]<-rownames(sclass)[i]
	mins[i]<-dat$conf.int[1]
	maxs[i]<-dat$conf.int[2]
	centers[i]<-dat$estimate
}
stadclassdat<-data.frame(xlo=log2(mins),xhi=log2(maxs),x=log2(centers),y=features)
stadclassdat$y<-factor(stadclassdat$y,levels=stadclassdat$y[c(4,6,1,2,7,10,3,9,5,8)])

r<-results_allloci2
r$type<-'stable'
r$type[which(results_allloci2$prop>0)]<-'unstable'
r$genomic_class<-droplevels(r$genomic_class)
table(r$genomic_class,r$type)
library(reshape2)
a<-dcast(r,genomic_class~type,length)
a$stable<-a$stable/table(r$type)[1]
a$unstable<-a$unstable/table(r$type)[2]
b<-melt(a,id="genomic_class")

anno<-read.csv('/net/shendure/vol3/data/nobackup/hauser/Shendure/databases/annovar/annovar/msi.hg19_multianno.csv')

anno$locus<-apply(anno,1,function(x)paste(x[1],paste(x[2],x[3],sep='-'),sep=':'))
anno$locus<-gsub(' ','',anno$locus)
anno2<-anno[match(r$locus,anno$locus),]

anno_test2<-function(x){
	unstablehits<-length(intersect(which(r$type=='unstable'),which(x!='.')))
unstablemisses<-length(intersect(which(r$type=='unstable'),which(x=='.')))
stablehits<-length(intersect(which(r$type=='stable'),which(x!='.')))
stablemisses<-length(intersect(which(r$type=='stable'),which(x=='.')))
return(fisher.test(matrix(c(unstablehits,unstablemisses,stablehits,stablemisses),nrow=2,byrow=T)))
}

a2<-apply(anno2[,6:14],2,anno_test2)

mins<-unlist(lapply(a2,function(x)x$conf.int[1]))
maxs<-unlist(lapply(a2,function(x)x$conf.int[2]))
centers<-unlist(lapply(a2,function(x)x$estimate))
features<-names(mins)
features<-c('phastCons','TFBS','miRNAs','miRNA targets','structural variants','DNase I HS sites','TFBS_bad','H3k27Ac','H3k4me3')
dat<-data.frame(xlo=mins,xhi=maxs,x=centers,y=features)
dat<-dat[-7,]

# d is a data frame with 4 columns
# d$x gives variable names
# d$y gives center point
# d$ylo gives lower limits
# d$yhi gives upper limits
forestplot <- function(d, xlab="Feature", ylab="Odds Ratio"){
    require(ggplot2)
    p <- ggplot(d, aes(x=y, y=x, ymin=xlo, ymax=xhi)) + 
		geom_pointrange(size=1.25) + 
		geom_hline(yintercept=1, lty=2,size=1.25) +
		ylab(ylab) +
		xlab(xlab) +
		theme_bw(base_size = 20) +
		theme(axis.text.x=element_text(size=10)) +
		theme(panel.grid.minor = element_blank(), panel.border = element_blank())
    return(p)
}

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/annovar_forestplot_071216.pdf',width=7*1.5,height=7)
forestplot(dat)
dev.off()