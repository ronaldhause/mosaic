load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_070816.robj')
dat<-results_msi
loci<-as.character(unique(dat$locus))
a<-strsplit(loci,'\\:')
b<-matrix(unlist(a),ncol=2,byrow=T)
c<-matrix(unlist(strsplit(b[,2],'-')),ncol=2,byrow=T)
d<-data.frame(locus=as.character(loci),chr=as.numeric(b[,1]),pos1=as.numeric(c[,1]),pos2=as.numeric(c[,2]))
d$mid<-apply(d[,c("pos1","pos2")],1,median)

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')

results_msi_per_cancer<-results

require(reshape2)
require(dplyr)
require(data.table)
results_msi_per_cancer$prop<-results_msi_per_cancer$num_unstable/results_msi_per_cancer$num_present
results_msi_per_cancer[,cancer_type:=paste(results_msi_per_cancer$tumor_type,results_msi_per_cancer$msi_status,sep='_')]
results_msi_per_cancer3<-subset(results_msi_per_cancer, tumor_type %in% c('UCEC','READ','COAD','STAD'))
dat<-results_msi_per_cancer3
dat$chr<-d$chr[match(dat$locus,d$locus)]
dat$pos<-d$mid[match(dat$locus,d$locus)]
dat2<-dcast(dat[,c(3,13:16),with=FALSE],chr + pos + locus ~ cancer_type,value.var="prop")
colnames(dat2)<-c('chr','pos','locus','coad_msih','coad_mss','read_msih','read_mss','stad_msih','stad_mss','ucec_msih','ucec_mss')
dat2<-dat2[complete.cases(dat2)]

acvr2a <- results_msi_per_cancer3[which(results_msi_per_cancer3$locus=='2:148683681-148683698'),]
orc4 <- results_msi_per_cancer3[which(results_msi_per_cancer3$locus=='2:148701095-148701119'),]


#figure 4a avcr2a locus plot
require(GenomicRanges)
require(ggbio)
require(BSgenome.Hsapiens.UCSC.hg19)
chr.lengths = seqlengths(Hsapiens)[1:22]
names(chr.lengths)<-seq(1,22)
g<-with(dat2,GRanges(chr,IRanges(pos,width=1), coad_msih=dat2$coad_msih,coad_mss=-dat2$coad_mss,read_msih=dat2$read_msih,read_mss=-dat2$read_mss,stad_msih=dat2$stad_msih,stad_mss=-dat2$stad_mss,ucec_msih=dat2$ucec_msih,ucec_mss=-dat2$ucec_mss))
seqlengths(g)<-chr.lengths

require(RColorBrewer)
cols = brewer.pal(9, "Paired")
require(Gviz)
idxTrack <- IdeogramTrack(genome="hg19", chromosome="chr2")
data(geneModels)
axTrack <- GenomeAxisTrack()
genetrack <- BiomartGeneRegionTrack(genome="hg19", chromosome="chr2", start=148500000, end = 148830000, showId=T, collapseTranscripts="longest", shape="arrow", fill="#8282d2")
subdat <- g[seqnames(g)==2,]
subdat2 <- subsetByOverlaps(subdat,GRanges(2,IRanges(start=148600000,end=148730000)))
seqlevelsStyle(subdat2)<-"UCSC"
coad <- DataTrack(range=subdat2, data=as.data.frame(mcols(subdat2)[,1:2]), genome="hg19", name="COAD", type="h", groups=c("coad_msih","coad_mss"), col=cols[1:2], ylim=c(-1,1),lwd=c(6),grid=TRUE,col.grid='gray75',lty.grid=2)
read <- DataTrack(range=subdat2, data=as.data.frame(mcols(subdat2)[,3:4]), genome="hg19", name="READ", type="h", groups=c("read_msih","read_mss"), col=cols[3:4], ylim=c(-1,1),lwd=6,grid=TRUE,col.grid='gray75',lty.grid=2)
stad <- DataTrack(range=subdat2, data=as.data.frame(mcols(subdat2)[,5:6]), genome="hg19", name="STAD", type="h", groups=c("stad_msih","stad_mss"), col=cols[5:6], ylim=c(-1,1),lwd=6,grid=TRUE,col.grid='gray75',lty.grid=2)
ucec <- DataTrack(range=subdat2, data=as.data.frame(mcols(subdat2)[,7:8]), genome="hg19", name="UCEC", type="h", groups=c("ucec_msih","ucec_mss"), col=cols[7:8], ylim=c(-1,1),lwd=6,grid=TRUE,col.grid='gray75',lty.grid=2)
microsatellites<-AnnotationTrack(subdat2,genome="hg19",stacking="dense")

pdffile  <- "/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/gviz_test_072616_3.pdf"
pdf(pdffile, 8, 8)
plotTracks(list(idxTrack,axTrack,coad,read, stad, ucec, microsatellites,genetrack),from=148650000,to=148730000)
dev.off()

#figure 4b circos
dat2<-dat2[complete.cases(dat2),]

require(OmicCircos)
data(UCSC.hg19.chr)
options(stringsAsFactors=FALSE)

library(RColorBrewer)
colors = brewer.pal(9, "Paired")
colnames(dat2)<-c('chr','pos','coad_msih','coad_mss','read_msih','read_mss','stad_msih','stad_mss','ucec_msih','ucec_mss')

#create 1 Mb bins across genome based on min and max coordinates per chr in coverage file
require(GenomicRanges)
require(ggbio)
require(BSgenome.Hsapiens.UCSC.hg19)
chr.lengths = seqlengths(Hsapiens)[1:22]
names(chr.lengths)<-seq(1,22)
g<-with(dat2,GRanges(chr,IRanges(pos,width=1), coad_msih=dat2$coad_msih,coad_mss=dat2$coad_mss,read_msih=dat2$read_msih,read_mss=dat2$read_mss,stad_msih=dat2$stad_msih,stad_mss=dat2$stad_mss,ucec_msih=dat2$ucec_msih,ucec_mss=dat2$ucec_mss))
seqlengths(g)<-chr.lengths

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
sum_per_bin <- sumPerBin(g, 1000000, mcolnames=colnames(dat2)[3:10])
counts <- countOverlaps(ans,g)
res <- apply(data.frame(mcols(sum_per_bin)),2,function(x)x/counts)

max(abs(res),na.rm=T) #0.34 for all sites, 0.6 for just unstable sites

midpoints <- floor((start(avg_per_bin)+end(avg_per_bin))/2)
plotdat<-data.frame(chr=as.numeric(as.matrix(seqnames(avg_per_bin))),pos=midpoints,res)
plotdat[,c(4,6,8,10)]<-plotdat[,c(4,6,8,10)]*-1
dat4<-plotdat

source('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/modified_circos2.r')

pdffile  <- "/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/testcircos_doublesided_just_completely_ascertained_sites_072516_2.pdf"
pdf(pdffile, 8, 8)
par(mar=c(2, 2, 2, 2))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="")
legend(680,800,colnames(dat3)[3:10],lty=1,lwd=2,col=colors[1:8],cex=0.5,title='Proportion unstable',box.col="white")
circos2(R=300, type="chr", cir="hg19", print.chr.lab=T, W=4);
circos2(R=250, cir="hg19", W=50, mapping=dat4[,c(1:2,3:4)], type="ml", B=F, col=colors[1:2], cex=0.4,scale=T,lwd=1);
circos2(R=200, cir="hg19", W=50, mapping=dat4[,c(1:2,5:6)], type="ml", B=F, col=colors[3:4], cex=0.4,scale=T,lwd=1);
circos2(R=150, cir="hg19", W=50, mapping=dat4[,c(1:2,7:8)], type="ml", B=F, col=colors[5:6], cex=0.4,scale=T,lwd=1);
circos2(R=100, cir="hg19", W=50, mapping=dat4[,c(1:2,9:10)], type="ml", B=F, col=colors[7:8], cex=0.4,scale=T,lwd=1);
dev.off()