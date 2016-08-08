require(data.table)

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_070816.robj') #results_msi
dat<-results_msi
loci<-as.character(unique(dat$locus))
a<-strsplit(loci,'\\:')
b<-matrix(unlist(a),ncol=2,byrow=T)
c<-matrix(unlist(strsplit(b[,2],'-')),ncol=2,byrow=T)
d<-data.frame(locus=as.character(loci),chr=as.numeric(b[,1]),pos1=as.numeric(c[,1]),pos2=as.numeric(c[,2]))
d$mid<-apply(d[,c("pos1","pos2")],1,median)

convert.magic1 <- function(obj,types){
    out <- lapply(1:length(obj),FUN = function(i){switch(types[i],character = as.character,numeric = as.numeric,factor = as.factor); FUN(obj[,i])})
    names(out) <- colnames(obj)
    as.data.frame(out,stringsAsFactors = FALSE)
}

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')

results_msi_per_cancer<-results

require(reshape2)
require(dplyr)
results_msi_per_cancer$prop<-results_msi_per_cancer$num_unstable/results_msi_per_cancer$num_present
results_msi_per_cancer[,cancer_type:=paste(results_msi_per_cancer$tumor_type,results_msi_per_cancer$msi_status,sep='_')]
results_msi_per_cancer3<-results_msi_per_cancer[-which(results_msi_per_cancer$num_present<(0.5*results_msi_per_cancer$total_num)),]
dat<-results_msi_per_cancer3
dat$chr<-d$chr[match(dat$locus,d$locus)]
dat$pos<-d$mid[match(dat$locus,d$locus)]
dat2<-dcast(dat[,13:16,with=FALSE],chr + pos ~ cancer_type,value.var="prop")
save(dat2,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_msi_per_cancer_props_for_plotting_and_binning_072016.robj')
dat3<-dat2[complete.cases(dat2)]

#mss vs. msi-h
results_msi_per_cancer <- results_msi

results_msi_per_cancer$prop<-results_msi_per_cancer$num_unstable/results_msi_per_cancer$num_present
results_msi_per_cancer3<-results_msi_per_cancer[-which(results_msi_per_cancer$num_present<(0.5*results_msi_per_cancer$total_num)),]
dat<-results_msi_per_cancer3
dat$chr<-d$chr[match(dat$locus,d$locus)]
dat$pos<-d$mid[match(dat$locus,d$locus)]
msidat2<-dcast(dat[,c(2,12:14),with=FALSE],chr + pos ~ msi_status,value.var="prop")
msidat3<-msidat2[complete.cases(msidat2)]

#create 1 Mb bins across genome based on min and max coordinates per chr in coverage file
require(GenomicRanges)
require(ggbio)
require(BSgenome.Hsapiens.UCSC.hg19)
chr.lengths = seqlengths(Hsapiens)[1:22]
names(chr.lengths)<-seq(1,22)
g<-with(dat3,GRanges(chr,IRanges(pos,width=1)))
values(g)<-dat3[,-c(1:2), with=FALSE]
seqlengths(g)<-chr.lengths

averagePerBin <- function(x, binsize, mcolnames=NULL)
{
     if (!is(x, "GenomicRanges"))
         stop("'x' must be a GenomicRanges object")
     if (any(is.na(seqlengths(x))))
         stop("'seqlengths(x)' contains NAs")
     bins <- IRangesList(lapply(seqlengths(x),
                                function(seqlen)
                                  IRanges(breakInChunks(seqlen, binsize))))
     ans <- as(bins, "GRanges")
     seqinfo(ans) <- seqinfo(x)
     if (is.null(mcolnames))
         return(ans)
     averageMCol <- function(colname)
     {
         cvg <- coverage(x, weight=colname)
         views_list <- RleViewsList(
                           lapply(names(cvg),
                               function(seqname)
                                   Views(cvg[[seqname]], bins[[seqname]])))
         unlist(viewMeans(views_list), use.names=FALSE)
     }
     mcols(ans) <- DataFrame(lapply(mcols(x)[mcolnames], averageMCol))
     ans
}

sumPerBin <- function(x, binsize, mcolnames=NULL)
{
     if (!is(x, "GenomicRanges"))
         stop("'x' must be a GenomicRanges object")
     if (any(is.na(seqlengths(x))))
         stop("'seqlengths(x)' contains NAs")
     bins <- IRangesList(lapply(seqlengths(x),
                                function(seqlen)
                                  IRanges(breakInChunks(seqlen, binsize))))
     ans <- as(bins, "GRanges")
     seqinfo(ans) <- seqinfo(x)
     if (is.null(mcolnames))
         return(ans)
     averageMCol <- function(colname)
     {
         cvg <- coverage(x, weight=colname)
         views_list <- RleViewsList(
                           lapply(names(cvg),
                               function(seqname)
                                   Views(cvg[[seqname]], bins[[seqname]])))
         unlist(viewSums(views_list), use.names=FALSE)
     }
     mcols(ans) <- DataFrame(lapply(mcols(x)[mcolnames], averageMCol))
     ans
}

#sum coverage per 1Mb bins
bins <- IRangesList(lapply(seqlengths(g),
                                function(seqlen)
                                  IRanges(breakInChunks(seqlen, 1000000))))
ans <- as(bins, "GRanges")
seqinfo(ans) <- seqinfo(g)

sum_per_bin <- sumPerBin(g, 1000000, colnames(mcols(g)))
counts <- countOverlaps(ans,g)
res <- apply(data.frame(mcols(sum_per_bin)),2,function(x)x/counts)

g2<-with(msidat3,GRanges(chr,IRanges(pos,width=1)))
values(g2)<-msidat3[,-c(1:2), with=FALSE]
seqlengths(g2)<-chr.lengths
msi_sum_per_bin <- sumPerBin(g2, 1000000, colnames(mcols(g2)))
counts <- countOverlaps(ans,g2)
msi_res <- apply(data.frame(mcols(msi_sum_per_bin)),2,function(x)x/counts)

require(rtracklayer)
require(humarray)
mySession <- browserSession ()
genome(mySession) <- "hg19"
filelist <- tableNames(ucscTableQuery(mySession, track="wgEncodeUwRepliSeq"))
filelist2 <- filelist[grep('WaveSignal',filelist)]
filelist3 <- filelist2[c(1:10,16)]

gm12878 <- import.bw('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/repliseq/wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigwig')
seqlevelsStyle(gm12878)<-'NCBI'
gm12878_auto <- dropSeqlevels(select.autosomes(gm12878),"X")
gm12878_binned <- averagePerBin(gm12878_auto, 1000000, "score")

files = list.files(path='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/repliseq/', full.names=T)
test <- lapply(files,function(x){
	t <- import.bw(x)
	seqlevelsStyle(t)<-'NCBI'
	t2 <- dropSeqlevels(select.autosomes(t),"X")
	averagePerBin(t2, 1000000, "score")
})

#mcf7 for brca1, bj for skcm, imr90 for lusc/luad, hepg2 for lihc
repliseq <- do.call(cbind, lapply(test, function(x)mcols(x)))
names(repliseq)<-c('bg02','bj','gm12878','hela3','hepg2','huvec','imr90','k562','mcf7','nhek','sknsh')
repliseq2<-as.data.frame(repliseq)
repliseq2$median<-apply(repliseq2,1,median)
repliseq2$mean<-apply(repliseq2[,1:11],1,mean)

cor.test(msi_res[,1],repliseq2$median, method='s') #rho = 0.020, p = 0.33, MSI-H
cor.test(msi_res[,2],repliseq2$median, method='s') #rho = 0.016, p = 0.41, MSS

res<-as.data.frame(res)
cor.test(res$BRCA_MSS,repliseq2$median,method='s')
cor.test(res$BRCA_MSS,repliseq2$mcf7,method='s')

cor.test(res$LIHC_MSS,repliseq2$median,method='s')
cor.test(res$LIHC_MSS,repliseq2$hepg2,method='s')
cor.test(res$LIHC_MSI.H,repliseq2$median,method='s')
cor.test(res$LIHC_MSI.H,repliseq2$hepg2,method='s')

cor.test(res$SKCM_MSS,repliseq2$median,method='s')
cor.test(res$SKCM_MSS,repliseq2$bj,method='s')

apply(res,2,function(x)cor.test(x,repliseq2$median,method='s')$estimate)

library(ggplot2)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/repliseq_vs_instability_072716.pdf',height=7,width=7)
df<-data.frame(msih=msi_res[,1]*100,mss=msi_res[,2]*100,repliseq=repliseq2$median)
df2<-melt(df,id="repliseq")
ggplot(df2,aes(x=repliseq,y=value,color=variable)) + geom_point() + xlab("replication timing") + ylab("proportion unstable") + theme_bw(base_size = 20) + geom_smooth(method=lm) + theme(axis.line = element_line(colour = "black"), panel.grid.major=element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) + theme(legend.position="top",legend.justification="center",legend.key = element_rect(fill = NA, colou = NA), legend.title=element_blank()) + scale_color_brewer(palette="Set1")
dev.off()