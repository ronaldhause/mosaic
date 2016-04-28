load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_per_cancer_042515.robj')

output<-data.frame(tcga_id=sampledat$sample_name,num_unstable=sampledat$num_unstable_raw,predicted_msi_status=sampledat$msi_status_tree_scaled,tumor_type=sampledat$tumor_type)
#write.table(output,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sample_information_supplemental_table.txt',quote=F,row.names=F)

require(reshape2)
ucecmss<-results_msi_per_cancer$median_peak_diff[intersect(which(results_msi_per_cancer$tumor_type=='UCEC'),which(results_msi_per_cancer$msi_status=='MSS'))]
ucecmsih<-results_msi_per_cancer$median_peak_diff[intersect(which(results_msi_per_cancer$tumor_type=='UCEC'),which(results_msi_per_cancer$msi_status=='MSI-H'))]
coadmss<-results_msi_per_cancer$median_peak_diff[intersect(which(results_msi_per_cancer$tumor_type=='COAD'),which(results_msi_per_cancer$msi_status=='MSS'))]
coadmsih<-results_msi_per_cancer$median_peak_diff[intersect(which(results_msi_per_cancer$tumor_type=='COAD'),which(results_msi_per_cancer$msi_status=='MSI-H'))]
readmss<-results_msi_per_cancer$median_peak_diff[intersect(which(results_msi_per_cancer$tumor_type=='READ'),which(results_msi_per_cancer$msi_status=='MSS'))]
readmsih<-results_msi_per_cancer$median_peak_diff[intersect(which(results_msi_per_cancer$tumor_type=='READ'),which(results_msi_per_cancer$msi_status=='MSI-H'))]
stadmss<-results_msi_per_cancer$median_peak_diff[intersect(which(results_msi_per_cancer$tumor_type=='STAD'),which(results_msi_per_cancer$msi_status=='MSS'))]
stadmsih<-results_msi_per_cancer$median_peak_diff[intersect(which(results_msi_per_cancer$tumor_type=='STAD'),which(results_msi_per_cancer$msi_status=='MSI-H'))]
lsmat<-cbind(ucecmss,ucecmsih,coadmss,coadmsih,readmss,readmsih,stadmss,stadmsih)
rownames(lsmat)<-results_msi_per_cancer$locus[intersect(which(results_msi_per_cancer$tumor_type=='UCEC'),which(results_msi_per_cancer$msi_status=='MSS'))]
# save(lsmat,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_080115.robj')

lsmat2<-lsmat[complete.cases(lsmat),]

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_080115.robj')

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r = (cor(x, y,use="pairwise",method='s'))
         txt <- format(c(r, 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex)
     }

#pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_spearman_081315.pdf',height=7,width=7)
#pairs(lsmat,lower.panel=panel.smooth,upper.panel=panel.cor)
#dev.off()

require(reshape2)
require(dplyr)
results_msi_per_cancer$prop<-results_msi_per_cancer$simple_gain_site/results_msi_per_cancer$num_present
results_msi_per_cancer$cancer_type<-paste(results_msi_per_cancer$tumor_type,results_msi_per_cancer$msi_status,sep='_')
results_msi_per_cancer2<-results_msi_per_cancer[-which(results_msi_per_cancer$num_present<(0.5*results_msi_per_cancer$total_num)),]
results_msi_per_cancer3<-subset(results_msi_per_cancer2, tumor_type %in% c('UCEC','READ','COAD','STAD'))
dat<-results_msi_per_cancer3
adat2<-dcast(dat,locus ~ cancer_type, value.var="prop")
adat3<-adat2[complete.cases(adat2),]
#158,912 sites

require(RColorBrewer)
my.cols <- rev(brewer.pal(11, "RdYlBu"))
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_spearman_smoothscatter_021816.pdf',height=14,width=14)
pairs(lsmat,lower.panel=function(...) smoothScatter(..., nrpoints = 100, add = TRUE,colramp=colorRampPalette(my.cols)),upper.panel=panel.cor)
dev.off()

require(vegan)
data(varespec)
vare.dist <- vegdist(lsmat, method = "jaccard", binary=TRUE)

x1 <- matrix(sample(c(0,1), 8, rep = TRUE), ncol = 2)
a<-dist(lsmat,method="Jaccard")

stability_overlap<-function(x,y)
	{
	mss_unstable_sites<-which(x>0)
	msi_unstable_sites<-which(y>0)
	overlap<-length(intersect(which(x>0),which(y>0)))/min(c(length(which(x>0)),length(which(y>0))))
    jaccard<-length(intersect(which(x>0),which(y>0)))/length(union(which(x>0),which(y>0)))
	return(overlap)
    return(jaccard)
}

#jaccard
stability_overlap2<-function(x,y)
    {
    x<-lsmat[,x]
    y<-lsmat[,y]
    mss_unstable_sites<-which(x>0)
    msi_unstable_sites<-which(y>0)
    overlap<-length(intersect(which(x>0),which(y>0)))/length(union(which(x>0),which(y>0)))
    return(overlap)
}

#jaccard
stability_overlap2<-function(x){
    ncx <- ncol(x)
    r <- matrix(0, nrow = ncx, ncol = ncx)
    for (i in seq_len(ncx)) {
        for (j in seq_len(ncx)) {
            x2 <- x[, i]
            y2 <- x[, j]
            ok <- complete.cases(x2, y2)
            r[i,j]<-length(intersect(which(x2>0),which(y2>0)))/length(union(which(x2>0),which(y2>0)))
            rownames(r) <- colnames(x)
            colnames(r) <- colnames(x)
        }
    }
    r
}

#overlap coefficient
stability_overlap3<-function(x){
    ncx <- ncol(x)
    r <- matrix(0, nrow = ncx, ncol = ncx)
    for (i in seq_len(ncx)) {
        for (j in seq_len(ncx)) {
            x2 <- x[, i]
            y2 <- x[, j]
            ok <- complete.cases(x2, y2)
            r[i,j]<-length(intersect(which(x2>0),which(y2>0)))/min(c(length(which(x2>0)),length(which(y2>0))))
            rownames(r) <- colnames(x)
            colnames(r) <- colnames(x)
        }
    }
    r
}

#cosine
stability_overlap4<-function(x){
    ncx <- ncol(x)
    r <- matrix(0, nrow = ncx, ncol = ncx)
    for (i in seq_len(ncx)) {
        for (j in seq_len(ncx)) {
            x2 <- x[, i]
            y2 <- x[, j]
            ok <- complete.cases(x2, y2)
            x2 <- x2[which(ok==TRUE)]
            y2 <- y2[which(ok==TRUE)]
            x2[which(x2>0)]<-1
            x2[which(x2<1)]<-0            
            y2[which(y2>0)]<-1
            y2[which(y2<1)]<-0
            r[i,j]<-sum(x2*y2)/sqrt(sum(x2^2)*sum(y2^2))
            rownames(r) <- colnames(x)
            colnames(r) <- colnames(x)
        }
    }
    r
}




overlap1<-stability_overlap3(lsmat)
overlap2<-stability_overlap2(lsmat)
overlap3<-stability_overlap4(lsmat)


col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3","#92C5DE","#D1E5F0","#FFFFFF", "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_jiccard_final_overlap_021816.pdf',width=7,height=7)
corrplot(overlap2,is.corr = TRUE, method = "color", type = "lower", tl.col = rgb(0, 0, 0),col = c(rep('black',100),col2(200)[101:200]),cl.lim=c(0,1),addCoef.col='black',addCoefasPercent=FALSE)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_jiccard_final_overlap_021816.pdf',width=7,height=7)
overlap2[which(overlap2==1)]<-0.5
corrplot(overlap2,is.corr = TRUE, method = "color", type = "lower", tl.col = rgb(0, 0, 0),col = c(rep('black',100),col2(100)[50:100],rep('black',50)),cl.lim=c(0,0.5),addCoef.col='black',addCoefasPercent=FALSE,diag=FALSE)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_cosine_final_overlap_021816.pdf',width=7,height=7)
corrplot(overlap3,is.corr = TRUE, method = "color", type = "lower", tl.col = rgb(0, 0, 0),col = c(rep('black',100),col2(200)[101:200]),cl.lim=c(0,1),addCoef.col='black',addCoefasPercent=FALSE)
dev.off()


,tl.cex = 1/par("cex"),
    cl.cex = 1/par("cex"), addCoefasPercent = TRUE,

#(rep('black',100),col2(100))

panel.overlap <- function(x, y, digits=2, prefix="", cex.cor)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r = stability_overlap(x,y)
         txt <- format(c(r, 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex)
     }

require(RColorBrewer)
my.cols <- rev(brewer.pal(11, "RdYlBu"))
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_jaccard_overlap_smoothscatter_021816.pdf',height=14,width=14)
pairs(lsmat,lower.panel=function(...) smoothScatter(..., nrpoints = 100, add = TRUE,colramp=colorRampPalette(my.cols)),upper.panel=panel.overlap)
dev.off()