load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_042515.robj')
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

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_per_cancer_042515.robj')

require(reshape2)
require(dplyr)
results_msi_per_cancer$prop<-results_msi_per_cancer$simple_gain_site/results_msi_per_cancer$num_present
results_msi_per_cancer$cancer_type<-paste(results_msi_per_cancer$tumor_type,results_msi_per_cancer$msi_status,sep='_')
results_msi_per_cancer2<-results_msi_per_cancer[-which(results_msi_per_cancer$num_present<(0.5*results_msi_per_cancer$total_num)),]
results_msi_per_cancer3<-subset(results_msi_per_cancer2, tumor_type %in% c('UCEC','READ','COAD','STAD'))
dat<-results_msi_per_cancer3
dat$chr<-d$chr[match(dat$locus,d$locus)]
dat$pos<-d$mid[match(dat$locus,d$locus)]
dat2<-dat[,14:17]
require(data.table)
dat2<-dcast(dat[,14:17,with=FALSE],chr + pos ~ cancer_type,value.var="prop")

require(OmicCircos)
data(UCSC.hg19.chr)
options(stringsAsFactors=FALSE)

#colors <- rainbow(10, alpha=0.5);

library(RColorBrewer)
colors = brewer.pal(9, "Paired")

dat3<-dat2
dat3[is.na(dat3)]<-0

colnames(dat3)<-c('chr','pos','coad_msih','coad_mss','read_msih','read_mss','stad_msih','stad_mss','ucec_msih','ucec_mss')

#create 1 Mb bins across genome based on min and max coordinates per chr in coverage file
require(GenomicRanges)
require(ggbio)
require(BSgenome.Hsapiens.UCSC.hg19)
chr.lengths = seqlengths(Hsapiens)[1:22]
names(chr.lengths)<-seq(1,22)
g<-with(dat3,GRanges(chr,IRanges(pos,width=1)))
values(g)<-DataFrame(coad_msih=dat3$coad_msih,coad_mss=dat3$coad_mss,read_msih=dat3$read_msih,read_mss=dat3$read_mss,stad_msih=dat3$stad_msih,stad_mss=dat3$stad_mss,ucec_msih=dat3$ucec_msih,ucec_mss=dat3$ucec_mss)
g2<-with(dat3,GRanges(chr,IRanges(pos,width=1),coad_msih=dat3$coad_msih,coad_mss=dat3$coad_mss,read_msih=dat3$read_msih,read_mss=dat3$read_mss,stad_msih=dat3$stad_msih,stad_mss=dat3$stad_mss,ucec_msih=dat3$ucec_msih,ucec_mss=dat3$ucec_mss))
seqlengths(g)<-chr.lengths

#b<-tileGenome(seqlengths(g), tilewidth=1000000,cut.last.tile.in.chrom=T)
#score <- mcolAsRleList(g, "coad_msih")
#binres<-binnedAverage(b,score,"binned_coad_msih")

SumPerBin <- function(x, binsize, mcolnames=NULL)
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
avg_per_bin <- averagePerBin(g, 1000000, mcolnames=c("coad_msih","coad_mss"))

dat4<-as.data.frame(dat3)
dat4[,c(4,6,8,10)]<-dat4[,c(4,6,8,10)]*-1

source('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/modified_circos.r')

pdffile  <- "/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/testcircos.pdf"
pdf(pdffile, 8, 8)
par(mar=c(2, 2, 2, 2))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="")

legend(680,800,colnames(dat3)[3:10],lty=1,lwd=2,col=colors[1:8],cex=0.5,title='Proportion unstable',box.col="white")

circos(R=300, type="chr", cir="hg19", print.chr.lab=T, W=4);
circos(R=250, cir="hg19", W=50, mapping=dat4[,c(1:2,3:4)], type="ml", B=F, col=colors[1:2], cex=0.4,scale=T,lwd=0.1);
circos(R=200, cir="hg19", W=50, mapping=dat4[,c(1:2,5:6)], type="ml", B=F, col=colors[3:4], cex=0.4,scale=T,lwd=0.4);
circos(R=150, cir="hg19", W=50, mapping=dat4[,c(1:2,7:8)], type="ml", B=F, col=colors[5:6], cex=0.4,scale=T,lwd=0.4);
circos(R=100, cir="hg19", W=50, mapping=dat4[,c(1:2,9:10)], type="ml", B=F, col=colors[7:8], cex=0.4,scale=T,lwd=0.4);

dev.off()