load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_post_review_070816.robj')

results_msi_per_cancer <- results

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
#save(lsmat,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_070816.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_070816.robj')

require(reshape2)
require(dplyr)
results_msi_per_cancer$prop<-results_msi_per_cancer$num_unstable/results_msi_per_cancer$num_present
results_msi_per_cancer$cancer_type<-paste(results_msi_per_cancer$tumor_type,results_msi_per_cancer$msi_status,sep='_')
results_msi_per_cancer2<-results_msi_per_cancer[-which(results_msi_per_cancer$num_present<(0.5*results_msi_per_cancer$total_num)),]
results_msi_per_cancer3<-subset(results_msi_per_cancer2, tumor_type %in% c('UCEC','READ','COAD','STAD'))
dat<-results_msi_per_cancer3
adat2<-dcast(dat,locus ~ cancer_type, value.var="prop")
adat3<-adat2[complete.cases(adat2),]
#121,154 sites

a<-subset(results_msi_per_cancer, tumor_type %in% c('UCEC','READ','COAD','STAD'))
orc4<-a[a$locus=='2:148701095-148701119',]

acvr2a<-a[a$locus=='2:148683681-148683698',]

jaccard_overlap<-function(x){
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

overlap_coefficient<-function(x){
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

cosine_similarity<-function(x){
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

jaccard_overlap2<-function(x){
    ncx <- ncol(x)
    r <- matrix(0, nrow = ncx, ncol = ncx)
    for (i in seq_len(ncx)) {
        for (j in seq_len(ncx)) {
            x2 <- x[, i]
            y2 <- x[, j]
            ok <- complete.cases(x2, y2)
            r[i,j]<-length(intersect(which(x2>=0.25),which(y2>=0.25)))/length(union(which(x2>=0.25),which(y2>=0.25)))
            rownames(r) <- colnames(x)
            colnames(r) <- colnames(x)
        }
    }
    r
}

overlap_coefficient2<-function(x){
    ncx <- ncol(x)
    r <- matrix(0, nrow = ncx, ncol = ncx)
    for (i in seq_len(ncx)) {
        for (j in seq_len(ncx)) {
            x2 <- x[, i]
            y2 <- x[, j]
            ok <- complete.cases(x2, y2)
            r[i,j]<-length(intersect(which(x2>=0.25),which(y2>=0.25)))/min(c(length(which(x2>=0.25)),length(which(y2>=0.25))))
            rownames(r) <- colnames(x)
            colnames(r) <- colnames(x)
        }
    }
    r
}

cosine_similarity2<-function(x){
    ncx <- ncol(x)
    r <- matrix(0, nrow = ncx, ncol = ncx)
    for (i in seq_len(ncx)) {
        for (j in seq_len(ncx)) {
            x2 <- x[, i]
            y2 <- x[, j]
            ok <- complete.cases(x2, y2)
            x2 <- x2[which(ok==TRUE)]
            y2 <- y2[which(ok==TRUE)]
            x2[which(x2>=0.25)]<-1
            x2[which(x2<0.25)]<-0            
            y2[which(y2>=0.25)]<-1
            y2[which(y2<0.25)]<-0
            r[i,j]<-sum(x2*y2)/sqrt(sum(x2^2)*sum(y2^2))
            rownames(r) <- colnames(x)
            colnames(r) <- colnames(x)
        }
    }
    r
}


overlap1<-overlap_coefficient(lsmat)
overlap2<-jaccard_overlap(lsmat)
overlap3<-cosine_similarity(lsmat)

adat22<-adat3[,-1]
overlap1n<-overlap_coefficient2(adat22)
overlap2n<-jaccard_overlap2(adat22)
overlap3n<-cosine_similarity2(adat22)

require(corrplot)
col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3","#92C5DE","#D1E5F0","#FFFFFF", "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_jiccard_final_overlap_070816_prop_0.25.pdf',width=7,height=7)
corrplot(overlap2n,is.corr = TRUE, method = "color", type = "lower", tl.col = rgb(0, 0, 0),col = c(rep('black',100),col2(200)[101:200]),cl.lim=c(0,1),addCoef.col='black',addCoefasPercent=FALSE)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_cosine_final_overlap_070816_prop_0.25.pdf',width=7,height=7)
corrplot(overlap3n,is.corr = TRUE, method = "color", type = "lower", tl.col = rgb(0, 0, 0),col = c(rep('black',100),col2(200)[101:200]),cl.lim=c(0,1),addCoef.col='black',addCoefasPercent=FALSE)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_mss_genomic_landscape_comparison_overlap_coefficient_070816_prop_0.25.pdf',width=7,height=7)
corrplot(overlap1n,is.corr = TRUE, method = "color", type = "lower", tl.col = rgb(0, 0, 0),col = c(rep('black',100),col2(200)[101:200]),cl.lim=c(0,1),addCoef.col='black',addCoefasPercent=FALSE)
dev.off()