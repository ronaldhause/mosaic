load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_042515.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_per_cancer_042515.robj')

canceroutput<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_global_loci_results_050315.txt',header=T)

# require(reshape2)
# results_per_cancer2<-results_per_cancer[-which(results_per_cancer$num_present<(results_per_cancer$total_num/2)),]
# results_per_cancer3<-split(results_per_cancer2,results_per_cancer2$tumor_type)

# results_per_cancer4 = Map(function(x, i) setNames(x, ifelse(names(x) %in% "locus",
#       names(x), sprintf('%s.%d', names(x), i))), results_per_cancer3, seq_along(results_per_cancer3))
# results_per_cancer5<-Reduce(function(...) merge(..., by="locus",all=FALSE),results_per_cancer4)
# save(results_per_cancer5,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_results_070115.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_results_070115.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_prop_test_results_070115.robj')

#117507 sites with sufficient data across all cancer types, testing for sites differentially unstable in one cancer independent of cancer status
unstable<-grep('simple_gain',colnames(results_per_cancer5))
all<-grep('num_present',colnames(results_per_cancer5))
#system.time(prop_results<-apply(results_per_cancer5,1,function(x)prop.test(as.numeric(x[unstable]),as.numeric(x[all])))) #2.5 min
unstablesites<-apply(results_per_cancer5, 1,function(x)sum(as.numeric(x[unstable])))
allcancers<-apply(results_per_cancer5,1,function(x)sum(as.numeric(x[all])))
propunstable<-unstablesites/allcancers
pvals<-unlist(lapply(prop_results,function(x)x$p.value))
props<-data.frame(matrix(unlist(lapply(prop_results,function(x)x$estimate)),nrow=dim(results_per_cancer5)[1],byrow=T))
#save(prop_results,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_prop_test_results_070115.robj')
prop_sums<-apply(props,1,sum)
#74418/117507 63.3% of sites are stable
prop.test(c(rep(0,17),5),as.numeric(x[all])) #5 or more total loci will give 95% power to identify cancer-specific loci
cancer_names<-as.character(as.matrix(results_per_cancer5[1,unstable-3]))
colnames(props)<-cancer_names
results<-data.frame(locus=results_per_cancer5$locus,genomic_class=results_per_cancer5$genomic_class.1,gene=results_per_cancer5$gene.1,repeat_type=results_per_cancer5$repeat_type.1,repeat_dna_sequence=results_per_cancer5$repeat_dna_sequence.1,props,pvals=pvals)
write.table(results,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_prop_test_results.txt',quote=F,row.names=F)
results2<-results[which(prop_sums>4),]
results2a<-results[which(propunstable>0.05),]
require(qvalue)
source('/net/shendure/vol1/home/hauser/scripts/shendure/miscellaneous/useful_functions.r')
#results3<-data.frame(results2,qvals=qvalue2(results2$pvals)$qvalues)
results3a<-data.frame(results2a,qvals=qvalue(results2a$pvals)$qvalues)
#hm<-results3[which(results3$qvals<0.05),]
#hm2<-as.matrix(hm[,-c(1:5,24:25)])
hm<-results3a[which(results3a$qvals<0.05),]
hm2<-as.matrix(hm[,-c(1:5,24:25)])

#require(gplots)
#pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_heatmap_fdr_5_bigger.pdf',width=14,height=14)
#heatmap.2(hm2,col=bluered,symkey=T,keysize=0.5,trace='none',labRow="",labCol="",scale='none',cexRow=0.4,cexCol=0.5)
#dev.off()

require(apcluster)
site.apclus <- apcluster(negDistMat(r=2), hm2,details=T,q=0)
cat("affinity propogation optimal number of clusters:", length(site.apclus@clusters), "\n") #5
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/affinity_propagation_sites.pdf',width=7,height=7)
heatmap(site.apclus)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/affinity_propagation_site_clustering.pdf',width=7,height=7)
plot(site.apclus,hm2)
dev.off()

 as.dendrogram2<-function(x,...)
 {
    .local <- function (object, base = 0.05, useNames = TRUE, 
        ...) 
    {
        if (all(dim(x@sim) <= 1)) 
            stop("similarity matrix not included in object")
        as.dendrogram(aggExCluster(x@sim, x, ...), 
            base = base, useNames = useNames)
    }
    .local(x, ...)
}

site.hclust<-as.hclust(site.apclus)
site.dendro<-as.dendrogram2(site.apclus)

hm2.1<-hm2[sample.hclust$order,site.hclust$order]

mydata <- hm2
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
  for (i in 2:20) wss[i] <- sum(kmeans(mydata,
                                       centers=i)$withinss)
mydata <- t(hm2)
wss2 <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
  for (i in 2:10) wss2[i] <- sum(kmeans(mydata,
                                       centers=i)$withinss)

hm2s<-scale(hm2)

require(NbClust)
nb <- NbClust(hm2, diss=NULL, distance = "euclidean", 
        min.nc=2, max.nc=15, method = "ward.D2", 
        index = "all", alphaBeale = 0.1) #2 or 4 for large

nb2 <- NbClust(t(hm2), diss=NULL, distance = "euclidean", 
        min.nc=2, max.nc=5, method = "ward.D2", 
        index = "all", alphaBeale = 0.1)

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/site_elbow_plot_large.pdf',width=7,height=7)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares", main="Sites")
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sample_elbow_plot_large.pdf',width=7,height=7)
plot(1:10, wss2, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares",main="Samples")
dev.off()

library(mclust)
d_clust <- Mclust(hm2,G=1:10)
m.best <- dim(d_clust$z)[2] #3 by 17; #3 or 4 x 9
cat("model-based optimal number of clusters:", m.best, "\n")

#for larger data
d_clust <- Mclust(hm2,G=1:10)
m.best <- dim(d_clust$z)[2] #6 site clusters

s_clust <- Mclust(t(hm2),G=1:5) #4 sample clusters

require("RColorBrewer")
blues<-colorRampPalette(brewer.pal(9,"Blues"))(256)
clusters<-cutree(hclust(dist(hm2),method='ward.D2'),k=4)
hm$clusters<-clusters
annotation<-data.frame(Cluster=as.character(hm$clusters))
rownames(annotation)<-rownames(hm)
c2<-cutree(hclust(dist(t(hm2)),method='ward.D2'),k=4)
annotationcol<-data.frame(Cancer=as.character(c2))
rownames(annotationcol)<-names(c2)
require(pheatmap)
cols<-brewer.pal(8, "Dark2")
cols2<-brewer.pal(8, "Set1")
anno_cols<-list(Cluster = c('1' = cols[1], '2' = cols[2], '3' = cols[3], '4' = cols[4]),Cancer= c('1' = cols2[1], '2' = cols2[2], '3' = cols2[6],'4' = cols2[7]))
pheatmap(hm2,col=blues,filename='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/figure3_cancer_specific_pheatmap_new_4clusters_blues_legend_030216_test.pdf',cutree_rows=4,cutree_cols=4,clustering_method="ward.D2",show_rownames=F,annotation_row=annotation,annotation_legend=T,cellwidth=10,annotation_colors=anno_cols,annotation_col=annotationcol)

hm$clusters2<-cutree(hclust(dist(hm2),method='ward.D2'),k=2)
hm$clusters3<-cutree(hclust(dist(hm2),method='ward.D2'),k=3)
hm$clusters4<-cutree(hclust(dist(hm2),method='ward.D2'),k=4)
hm$clusters6<-cutree(hclust(dist(hm2),method='ward.D2'),k=6)

require(devtools)
install_github(c("GuangchuangYu/DOSE", "GuangchuangYu/clusterProfiler"))

. ~hauser/scripts/shendure/miscellaneous/r_install_github.sh clusterProfiler GuangchuangYu

. ~hauser/scripts/shendure/miscellaneous/r_install_github.sh clusterProfiler GuangchuangYu

wget https://github.com/GuangchuangYu/clusterProfiler/archive/master.zip 

require(clusterProfiler)
require(dplyr)
require('biomaRt')
mart<-useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "jul2015.archive.ensembl.org")
ids <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","entrezgene","hgnc_symbol"),values=hm$gene,mart= mart)
#ids <- bitr(hm$gene, fromType="SYMBOL", toType=c("UNIPROT", "ENTREZID"), OrgDb="org.Hs.eg.db")
#head(ids)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
head(gene)
## [1] "4312"  "8318"  "10874" "55143" "55388" "991"
background<-results$gene
bids <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","entrezgene","hgnc_symbol"),values=background,mart= mart)
backgroundgenes<-as.factor(bids$entrezgene[match(background,bids$hgnc_symbol)])

hm$entrezgene<-as.factor(ids$entrezgene[match(hm$gene,ids$hgnc_symbol)])
hm$clusters<-as.factor(hm$clusters)

a<-groupGO(genes,ont="MF",level=3,readable=TRUE, OrgDb = 'org.Hs.eg.db')

geneclusterlist<-split(as.character(hm$entrezgene),hm$clusters4)

stablegenes<-as.factor(bids$entrezgene[match(results$gene[which(prop_sums==0)],bids$hgnc_symbol)])
unstablegenes<-as.factor(bids$entrezgene[match(results$gene[which(prop_sums>0)],bids$hgnc_symbol)])
stablegenes<-stablegenes[which(is.na(match(stablegenes,unstablegenes))==T)] #5772

geneclusterlist[[5]]<-as.character(stablegenes)

bpgoresults<-lapply(geneclusterlist,function(x)enrichGO(gene          = x,
                universe      = backgroundgenes,
                OrgDb = 'org.Hs.eg.db',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.10,
                readable      = TRUE,
                ))

bpgoresults2<-lapply(bpgoresults,function(x)simplify(x, cutoff=0.7, by="p.adjust", select_fun=min))

mfgoresults<-lapply(geneclusterlist,function(x)enrichGO(gene          = x,
                universe      = backgroundgenes,
                OrgDb = 'org.Hs.eg.db',
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.10,
                readable      = TRUE))

mfgoresults2<-lapply(mfgoresults,function(x)simplify(x, cutoff=0.7, by="p.adjust", select_fun=min))

ccgoresults<-lapply(geneclusterlist,function(x)enrichGO(gene          = x,
                universe      = backgroundgenes,
                OrgDb = 'org.Hs.eg.db',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.10,
                readable      = TRUE))

ccgoresults2<-lapply(ccgoresults,function(x)simplify(x, cutoff=0.7, by="p.adjust", select_fun=min))

keggresults <- lapply(geneclusterlist,function(x)enrichKEGG(gene = x, universe = backgroundgenes, organism="hsa", pvalueCutoff=0.10, pAdjustMethod="BH"))

require(data.table)
cluster1results<-rbind(lapply(bpgoresults2,summary)[[1]],lapply(mfgoresults2,summary)[[1]],lapply(ccgoresults2,summary)[[1]],lapply(keggresults,summary)[[1]])
cluster2results<-rbind(lapply(bpgoresults2,summary)[[2]],lapply(mfgoresults2,summary)[[2]],lapply(ccgoresults2,summary)[[2]],lapply(keggresults,summary)[[2]])
cluster3results<-rbind(lapply(bpgoresults2,summary)[[3]],lapply(mfgoresults2,summary)[[3]],lapply(ccgoresults2,summary)[[3]],lapply(keggresults,summary)[[3]])
cluster4results<-rbind(lapply(bpgoresults2,summary)[[4]],lapply(mfgoresults2,summary)[[4]],lapply(ccgoresults2,summary)[[4]],lapply(keggresults,summary)[[4]])
cluster5results<-rbind(lapply(bpgoresults2,summary)[[5]],lapply(mfgoresults2,summary)[[5]],lapply(ccgoresults2,summary)[[5]],lapply(keggresults,summary)[[5]])

cluster1results$fe<-sapply(strsplit(as.character(cluster1results$GeneRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))/sapply(strsplit(as.character(cluster1results$BgRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))
#cluster1results2<-cluster1results[-match(c('biological_process','molecular_function','cellular_component'),cluster1results$Description),]
cluster1results2<-cluster1results[rev(order(cluster1results$fe)),]
cluster1results2<-cluster1results2[-which(cluster1results2$Count<5),]
#cluster1results2<-cluster1results2[which(cluster1results2$fe>1.25),]

cluster2results$fe<-sapply(strsplit(as.character(cluster2results$GeneRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))/sapply(strsplit(as.character(cluster2results$BgRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))
#cluster2results2<-cluster2results[-match(c('biological_process','molecular_function','cellular_component'),cluster2results$Description),]
cluster2results2<-cluster2results[rev(order(cluster2results$fe)),]
cluster2results2<-cluster2results2[-which(cluster2results2$Count<5),]
#cluster2results2<-cluster2results2[which(cluster2results2$fe>1.25),]

cluster3results$fe<-sapply(strsplit(as.character(cluster3results$GeneRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))/sapply(strsplit(as.character(cluster3results$BgRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))
#cluster3results2<-cluster3results[-match(c('biological_process','molecular_function','cellular_component'),cluster3results$Description),]
cluster3results2<-cluster3results[rev(order(cluster3results$fe)),]
cluster3results2<-cluster3results2[-which(cluster3results2$Count<5),]

cluster4results$fe<-sapply(strsplit(as.character(cluster4results$GeneRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))/sapply(strsplit(as.character(cluster4results$BgRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))
#cluster4results2<-cluster4results[-match(c('biological_process','molecular_function'),cluster4results$Description),]
cluster4results2<-cluster4results[rev(order(cluster4results$fe)),]
cluster4results2<-cluster4results2[-which(cluster4results2$Count<5),]

cluster5results$fe<-sapply(strsplit(as.character(cluster5results$GeneRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))/sapply(strsplit(as.character(cluster5results$BgRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))
#cluster4results2<-cluster4results[-match(c('biological_process','molecular_function'),cluster4results$Description),]
cluster5results2<-cluster5results[rev(order(cluster5results$fe)),]


output<-rbind(cluster1results2,cluster2results2,cluster3results2,cluster4results2,cluster5results2)
gooutput<-data.frame(cluster=c(rep('1',dim(cluster1results2)[1]),rep('2',dim(cluster2results2)[1]),rep('3',dim(cluster3results2)[1]),rep('4',dim(cluster4results2)[1]),rep('stable',dim(cluster5results2)[1])),output)
#gooutput<-gooutput[-which(gooutput$Count<5),]
write.csv(gooutput,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/fig3_goenrichment_results_030216.csv',row.names=F)

#plotting fold changes
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/fig3_goenrichment_qvalues_030216.pdf',height=7*1.5,width=7*0.7)
par(mar=c(2,20,2,1),mfrow=c(5,1))
barplot(log2(cluster4results2$fe[1:10]+0.001), horiz=TRUE, names.arg=cluster4results2$Description[1:10], border="grey94", cex.names=.9, col=cols[4],las=1,main='Cluster 4',xlim=c(0,max(log2(gooutput$fe+0.001))))
barplot(log2(cluster1results2$fe[1:10]+0.001), horiz=TRUE, names.arg=cluster1results2$Description[1:10], border="grey94", cex.names=.9, col=cols[1],las=1,main='Cluster 1',xlim=c(0,max(log2(gooutput$fe+0.001))))
barplot(log2(cluster3results2$fe[1:10]+0.001), horiz=TRUE, names.arg=cluster3results2$Description[1:10], border="grey94", cex.names=.9, col=cols[3],las=1,main='Cluster 3',xlim=c(0,max(log2(gooutput$fe+0.001))))
barplot(log2(cluster2results2$fe[1:10]+0.001), horiz=TRUE, names.arg=cluster2results2$Description[1:10], border="grey94", cex.names=.9, col=cols[2],las=1,main='Cluster 2',xlim=c(0,max(log2(gooutput$fe+0.001))))
barplot(log2(cluster5results2$fe[1:10]+0.001), horiz=TRUE, names.arg=cluster5results2$Description[1:10], border="grey94", cex.names=.9, col="grey55",las=1,main='Stable genes',xlim=c(0,max(log2(gooutput$fe+0.001))))
dev.off()


pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/fig3_goenrichment_qvalues.pdf', height=7*3,width=7*2)
par(mar=c(0,6,0,.5),mfrow=c(2,2))
barplot()
plot(0, ylim=c(0, max(log2(cluster1results2$fe+0.0001))), type="n", yaxt="n", frame.plot=FALSE, xlab="", ylab="")
mypos<-barplot(log2(cluster1results2$fe[1:10]+0.001), horiz=FALSE, axes=FALSE, add=TRUE, names.arg=rep("", 10), border="grey94", space=0, cex.names=.9, col="grey55")
mtext("log2(fold change)", side=1, cex=.7, line=2.2, adj=1.1)
axis(2, at=mypos, labels=rev(cluster1results2$Description[1:10]), lwd=0, , cex.axis=1.1, line=-.2)
dev.off()

axis(2, at=mypos, labels=rev(hm$gene), lwd=0, , cex.axis=1.1, line=-.2)
mtext("-log10(q-value)", side=1, cex=.7, line=2.2, adj=1.1)
abline(v=-log10(.1), col="red", lwd=1)
abline(v=-log10(.25), col="purple", lty=2, lwd=1)

ego1 <- enrichGO(gene          = hm$entrezgene[which(hm$clusters4==1)],
                universe      = backgroundgenes,
                organism      = "human",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
simplify(ego1)
require(magrittr)
require(GOSemSim)
require(tidyr)
ego1simp <- simplify(ego1, cutoff=0.7, by="p.adjust", select_fun=min)

genes<-as.character(na.omit(hm$entrezgene))

ekegg1 <- enrichKEGG(gene          = hm$entrezgene[which(hm$clusters6==1)],
                universe      = backgroundgenes, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH",
                 qvalueCutoff=0.05, readable=TRUE)
require(magrittr)
require(GOSemSim)
require(tidyr)
ego1simp <- simplify(ego1, cutoff=0.7, by="p.adjust", select_fun=min)

geneclusterlist<-split(as.character(hm$entrezgene),hm$clusters4)
x=compareCluster(geneclusterlist, fun="enrichGO")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_all_go_enrichments.pdf',width=7*0.75,height=7*1.25)
plot(x)
dev.off()

#multiple groups
hm$entrezgene<-as.factor(ids$entrezgene[match(hm$gene,ids$hgnc_symbol)])
hm$clusters<-as.factor(hm$clusters)
a<-data.frame(entrez=hm$entrezgene,clusters=hm$clusters)
yy.formula <- compareCluster(entrez~clusters,data=a,fun='groupGO')
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/global_cluster_enrichments.pdf',width=7,height=7)
plot(yy.formula)
dev.off()

background<-results$gene
bids <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","entrezgene","hgnc_symbol"),values=background,mart= mart)
backgroundgenes<-as.factor(bids$entrezgene[match(background,bids$hgnc_symbol)])

stablegenes<-as.factor(bids$entrezgene[match(results$gene[which(prop_sums==0)],bids$hgnc_symbol)])

data(gcSample)
lapply(gcSample, head)
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG", use_internal_data = TRUE, universe=)
head(summary(ck))

#cluster1
ego1 <- enrichGO(gene          = hm$entrezgene,
                universe      = backgroundgenes,
                organism      = "human",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.20,
                readable      = TRUE)
require(magrittr)
require(GOSemSim)
require(tidyr)
ego1simp <- simplify(ego1, cutoff=0.7, by="p.adjust", select_fun=min)

# #DAVID
# geneclusterlist<-split(as.character(hm$entrezgene),hm$clusters)
# x=compareCluster(geneclusterlist, fun="enrichDAVID", annotation="KEGG_PATHWAY")
# pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_all_david_enrichments.pdf',width=7*0.75,height=7*1.25)
# plot(x)
# dev.off()

# #GSEA
# results_allloci$prop<-results_allloci$simple_gain_site/results_allloci$num_present
# results_allloci2<-results_allloci[-which(results_allloci$num_present<2112),]
# results_allloci2$gene<-factor(results_allloci2$gene)
# geneprops<-tapply(results_allloci2$prop,results_allloci2$gene,function(x)mean(x,na.rm=T))
# names(geneprops)<-bids$entrezgene[match(names(geneprops),bids$hgnc_symbol)]
# geneprops<-geneprops[na.omit(names(geneprops))]
# geneprops<-rev(sort(geneprops))
# geneprops<-as.numeric(geneprops)
# ego2 <- gseGO(geneList = geneprops, organism = "human", ont = "CC", nPerm = 100, minGSSize = 120, pvalueCutoff = 0.01,verbose = FALSE)
# head(summary(ego2))

output3<-read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_global_msi_050315.txt',header=T)
msisiglocs<-output3$locus[1:100]

hm3<-as.matrix(results[na.omit(match(msisiglocs,results$locus)),c(6:8,20)])
rownames(hm3)<-results$gene[na.omit(match(msisiglocs,results$locus))]
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_specific_loci_msi_cancer_distributions_quartiles_080115.pdf',width=7*0.75,height=7*1.25)
heatmap.2(hm3,col=bluered(256),key=TRUE,symkey=FALSE,symm=FALSE,keysize=0.5,trace='none',density.info="none",scale='none',cexRow=0.3,cexCol=0.5,Rowv=FALSE,lwid=c(3,4),lhei=c(2,4),breaks=c(seq(min(hm3),fivenum(hm3)[2],length.out=86),seq(fivenum(hm3)[2],fivenum(hm3)[4],length.out=85),seq(fivenum(hm3)[4],max(hm3),length.out=86)))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_specific_loci_msi_cancer_distributions_quartiles_row_ordered_080115.pdf',width=7*0.75,height=7*1.25)
heatmap.2(hm3,col=bluered(256),key=TRUE,symkey=FALSE,symm=FALSE,keysize=0.5,trace='none',density.info="none",scale='none',cexRow=0.3,cexCol=0.5,lwid=c(3,4),lhei=c(2,4),breaks=c(seq(min(hm3),fivenum(hm3)[2],length.out=86),seq(fivenum(hm3)[2],fivenum(hm3)[4],length.out=85),seq(fivenum(hm3)[4],max(hm3),length.out=86)))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_specific_loci_msi_cancer_distributions_raw_080115.pdf',width=7*0.75,height=7)
heatmap.2(hm3,col=bluered(256),symkey=FALSE,keysize=0.5,trace='none',scale='none',cexRow=0.3,cexCol=0.5,Rowv=FALSE,lwid=c(1.5,4),lhei=c(1.5,4))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_specific_loci_msi_cancer_distributions_scaled_080115.pdf',width=7*0.75,height=7)
heatmap.2(hm3,col=bluered(256),symkey=FALSE,keysize=0.5,trace='none',scale='row',cexRow=0.3,cexCol=0.5,Rowv=FALSE,lwid=c(1.5,4),lhei=c(1.5,4))
dev.off()

write.table(results2,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_instability_060715.txt',quote=F,row.names=F)

#FOUR WAY PLOT

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_msi_per_cancer_042515.robj')

plotdat<-data.frame(locus=output2$locus,coad=output2$COAD_num_unstable/output2$COAD_num_called,ucec=output2$UCEC_num_unstable/output2$UCEC_num_called,read=output2$READ_num_unstable/output2$READ_num_called,stad=output2$STAD_num_unstable/output2$STAD_num_called)
plotdat2<-data.frame(locus=plotdat$locus,x=plotdat$read-plotdat$coad,y=plotdat$stad-plotdat$ucec,number=rowSums(output2[,c(2,4,6,30)]))
#require(tpsmeta)
angle<-1.75*pi #rotating 45 degrees clockwise
plotdat3<-data.frame(locus=plotdat2$locus,x=cos(angle)*plotdat2$x-sin(angle)*plotdat2$y,y=sin(angle)*plotdat2$x+cos(angle)*plotdat2$y,number=rowSums(output2[,c(2,4,6,30)]))

library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols<-brewer.pal(9,'Set1')
steps<-c(cols[4],'black',cols[2]) #purple COAD to blue READ
pal<-color.palette(steps, c(96,96), space="rgb")
a<-hex2RGB(pal(117507))
steps2<-c(cols[1],'black',cols[3]) #red UCEC to green STAD
pal2<-color.palette(steps, c(96,96), space="rgb")
b<-hex2RGB(pal2(117507))
c<-mixcolor(0.5,a,b)
d<-hex(c)

library(grid)
my_grob = grobTree(textGrob("UCEC", x=-0.4,  y=0.4, hjust=0,
  gp=gpar(col="black", fontsize=12)))

plotdat2<-data.frame(locus=plotdat$locus,x=plotdat$read-plotdat$coad,y=plotdat$stad-plotdat$ucec,number=rowSums(output2[,c(2,4,6,30)]),color=d)

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_unstable_sites_sized_070615.pdf',height=7,width=7)
ggplot(plotdat2,aes(x=x,y=y)) + geom_point(aes(size=number),alpha=0.75) + xlab("") + ylab("") + theme_bw(base_size = 20) + scale_size_continuous(range = c(0, 3)) + coord_cartesian(xlim=c(-0.3,0.3),ylim=c(-0.3,0.3)) + geom_hline(aes(yintercept=0),linetype='dashed') + geom_vline(aes(xintercept=0),linetype='dashed') +
	xlab('proportion unstable in READ - COAD') +
	ylab('proportion unstable in STAD - UCEC')
dev.off()

#compareCluster for cancer-specific sites
canspec<-c()
for (i in 1:dim(hm2)[1]){
  site<-hm2[i,]
  sitecheck<-site[which(site<max(site))]<(0.75*max(site))
  if (sum(sitecheck)==17){
    canspec[i]<-colnames(hm2)[which(site==max(site))]
  }
  else{
    canspec[i]<-NA
  }
}

geneclusterlist<-split(as.character(hm$entrezgene),hm$clusters4)
x=compareCluster(geneclusterlist, fun="enrichGO")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_all_go_enrichments.pdf',width=7*0.75,height=7*1.25)
plot(x)
dev.off()
