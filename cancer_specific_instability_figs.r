load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_per_cancer_070816.robj')

# require(reshape2)
# results_per_cancer2<-results_per_cancer[-which(results_per_cancer$num_present<(results_per_cancer$total_num/2)),]
# results_per_cancer3<-split(results_per_cancer2,results_per_cancer2$tumor_type)

# results_per_cancer4 = Map(function(x, i) setNames(x, ifelse(names(x) %in% "locus", names(x), sprintf('%s.%d', names(x), i))), results_per_cancer3, seq_along(results_per_cancer3))
# results_per_cancer5<-Reduce(function(...) merge(..., by="locus",all=FALSE),results_per_cancer4)
# save(results_per_cancer5,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_results_071316.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_results_071316.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_prop_test_results_070115.robj')

#92,385 sites with sufficient data across all cancer types, testing for sites differentially unstable in one cancer independent of cancer status
unstable<-grep('num_unstable',colnames(results_per_cancer5))
all<-grep('num_present',colnames(results_per_cancer5))
#system.time(prop_results<-apply(results_per_cancer5,1,function(x)prop.test(as.numeric(x[unstable]),as.numeric(x[all])))) #2.5 min
unstablesites<-apply(results_per_cancer5, 1,function(x)sum(as.numeric(x[unstable])))
allcancers<-apply(results_per_cancer5,1,function(x)sum(as.numeric(x[all])))
propunstable<-unstablesites/allcancers
pvals<-unlist(lapply(prop_results,function(x)x$p.value))
props<-data.frame(matrix(unlist(lapply(prop_results,function(x)x$estimate)),nrow=dim(results_per_cancer5)[1],byrow=T))
#save(prop_results,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_prop_test_results_070115.robj')
prop_sums<-apply(props,1,sum)
#53004/92385 57.4% of sites are stable
#prop.test(c(rep(0,17),5),as.numeric(x[all])) #5 or more total loci will give 95% power to identify cancer-specific loci
cancer_names<-as.character(as.matrix(results_per_cancer5[1,])[unstable-2])
colnames(props)<-cancer_names
results<-data.frame(locus=results_per_cancer5$locus,genomic_class=results_per_cancer5$genomic_class.1,gene=results_per_cancer5$gene.1,repeat_type=results_per_cancer5$repeat_type.1,repeat_dna_sequence=results_per_cancer5$repeat_dna_sequence.1,props,pvals=pvals)
write.table(results,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_prop_test_results_071316.txt',quote=F,row.names=F)
write.csv(results,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_prop_test_results_071316.csv',row.names=F)
results2a<-results[which(propunstable>0.05),]
require(qvalue)
source('/net/shendure/vol1/home/hauser/scripts/shendure/miscellaneous/useful_functions.r')
results3a<-data.frame(results2a,qvals=qvalue(results2a$pvals)$qvalues)
hm<-results3a[which(results3a$qvals<0.05),]
hm2<-as.matrix(hm[,-c(1:5,24:25)])

require(apcluster)
site.apclus <- apcluster(negDistMat(r=2), hm2,details=T, q=0)
cat("affinity propogation optimal number of clusters:", length(site.apclus@clusters), "\n") #7
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/affinity_propagation_sites_071316.pdf',width=7,height=7)
heatmap(site.apclus)
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

sample.apclus <- apcluster(negDistMat(r=2), t(hm2),details=T) #5
sample.hclust<-as.hclust(sample.apclus)
sample.dendro<-as.dendrogram2(sample.apclus)

library(mclust)
d_clust <- Mclust(hm2,G=1:5)
m.best <- dim(d_clust$z)[2] #4 site clusters

s_clust <- Mclust(t(hm2),G=1:5) #4 sample clusters

require("RColorBrewer")
blues<-colorRampPalette(brewer.pal(9,"Blues"))(256)
clusters<-cutree(hclust(dist(hm2),method='ward.D2'),k=6)
hm$clusters<-clusters
annotation<-data.frame(Cluster=as.character(hm$clusters))
rownames(annotation)<-rownames(hm)
c2<-cutree(hclust(dist(t(hm2)),method='ward.D2'),k=5)
annotationcol<-data.frame(Cancer=as.character(c2))
rownames(annotationcol)<-names(c2)
require(pheatmap)
cols<-brewer.pal(8, "Dark2")
cols2<-brewer.pal(8, "Set1")
anno_cols<-list(Cluster = c('1' = cols[1], '2' = cols[2], '3' = cols[3], '4' = cols[4]),Cancer= c('1' = cols2[1], '2' = cols2[2], '3' = cols2[6],'4' = cols2[7]))
pheatmap(hm2,col=blues,filename='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/figure3_cancer_specific_pheatmap_new_4clusters_blues_legend_071316.pdf',cutree_rows=4,cutree_cols=4,clustering_method="ward.D2",show_rownames=F,annotation_row=annotation,annotation_legend=T,cellwidth=10,annotation_colors=anno_cols,annotation_col=annotationcol)

hm$clusters<-cutree(hclust(dist(hm2),method='ward.D2'),k=4)

require(clusterProfiler)
require(dplyr)
require('biomaRt')
mart<-useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "jul2015.archive.ensembl.org")
ids <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","entrezgene","hgnc_symbol"),values=hm$gene,mart= mart)
background<-results$gene
bids <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","entrezgene","hgnc_symbol"),values=background,mart= mart)
backgroundgenes<-as.factor(bids$entrezgene[match(background,bids$hgnc_symbol)])

hm$entrezgene<-as.factor(ids$entrezgene[match(hm$gene,ids$hgnc_symbol)])
hm$clusters<-as.factor(hm$clusters)

geneclusterlist<-split(as.character(hm$entrezgene),hm$clusters)

stablegenes<-as.factor(bids$entrezgene[match(results$gene[which(prop_sums==0)],bids$hgnc_symbol)])
unstablegenes<-as.factor(bids$entrezgene[match(results$gene[which(prop_sums>0)],bids$hgnc_symbol)])
stablegenes<-stablegenes[which(is.na(match(stablegenes,unstablegenes))==T)] #2285

geneclusterlist[[5]]<-as.character(stablegenes)

bpgoresults<-lapply(geneclusterlist,function(x)enrichGO(gene          = x,
                universe      = backgroundgenes,
                OrgDb = 'org.Hs.eg.db',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.20,
                readable      = TRUE,
                ))

bpgoresults2<-lapply(bpgoresults,function(x)simplify(x, cutoff=0.7, by="p.adjust", select_fun=min))

mfgoresults<-lapply(geneclusterlist,function(x)enrichGO(gene          = x,
                universe      = backgroundgenes,
                OrgDb = 'org.Hs.eg.db',
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.20,
                readable      = TRUE))

mfgoresults2<-lapply(mfgoresults,function(x)simplify(x, cutoff=0.7, by="p.adjust", select_fun=min))

ccgoresults<-lapply(geneclusterlist,function(x)enrichGO(gene          = x,
                universe      = backgroundgenes,
                OrgDb = 'org.Hs.eg.db',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.20,
                readable      = TRUE))

ccgoresults2<-lapply(ccgoresults,function(x)simplify(x, cutoff=0.7, by="p.adjust", select_fun=min))

keggresults <- lapply(geneclusterlist,function(x)enrichKEGG(gene = x, universe = backgroundgenes, organism="hsa", pvalueCutoff=0.20, pAdjustMethod="BH"))

require(data.table)
cluster1results<-rbind(lapply(bpgoresults2,summary)[[1]],lapply(mfgoresults2,summary)[[1]],lapply(ccgoresults2,summary)[[1]],lapply(keggresults,summary)[[1]])
cluster2results<-rbind(lapply(bpgoresults2,summary)[[2]],lapply(mfgoresults2,summary)[[2]],lapply(ccgoresults2,summary)[[2]],lapply(keggresults,summary)[[2]])
cluster3results<-rbind(lapply(bpgoresults2,summary)[[3]],lapply(mfgoresults2,summary)[[3]],lapply(ccgoresults2,summary)[[3]],lapply(keggresults,summary)[[3]])
cluster4results<-rbind(lapply(bpgoresults2,summary)[[4]],lapply(mfgoresults2,summary)[[4]],lapply(ccgoresults2,summary)[[4]],lapply(keggresults,summary)[[4]])
cluster5results<-rbind(lapply(bpgoresults2,summary)[[5]],lapply(mfgoresults2,summary)[[5]],lapply(ccgoresults2,summary)[[5]],lapply(keggresults,summary)[[5]])

cluster1results$fe<-sapply(strsplit(as.character(cluster1results$GeneRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))/sapply(strsplit(as.character(cluster1results$BgRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))
cluster1results2<-cluster1results[rev(order(cluster1results$fe)),]

cluster2results$fe<-sapply(strsplit(as.character(cluster2results$GeneRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))/sapply(strsplit(as.character(cluster2results$BgRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))
cluster2results2<-cluster2results[rev(order(cluster2results$fe)),]

cluster3results$fe<-sapply(strsplit(as.character(cluster3results$GeneRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))/sapply(strsplit(as.character(cluster3results$BgRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))
cluster3results2<-cluster3results[rev(order(cluster3results$fe)),]

cluster4results$fe<-sapply(strsplit(as.character(cluster4results$GeneRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))/sapply(strsplit(as.character(cluster4results$BgRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))
cluster4results2<-cluster4results[rev(order(cluster4results$fe)),]

cluster5results$fe<-sapply(strsplit(as.character(cluster5results$GeneRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))/sapply(strsplit(as.character(cluster5results$BgRatio), split = "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2]))
cluster5results2<-cluster5results[rev(order(cluster5results$fe)),]

output<-rbind(cluster1results2,cluster2results2,cluster3results2,cluster4results2,cluster5results2)
gooutput<-data.frame(cluster=c(rep('1',dim(cluster1results2)[1]),rep('2',dim(cluster2results2)[1]),rep('3',dim(cluster3results2)[1]),rep('4',dim(cluster4results2)[1]),rep('stable',dim(cluster5results2)[1])),output)
write.csv(gooutput,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/fig3_goenrichment_results_071316.csv',row.names=F)

#plotting fold changes
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/fig3_goenrichment_qvalues_071316.pdf',height=7*1.5,width=7*0.7)
par(mar=c(2,20,2,1),mfrow=c(5,1))
barplot(log2(cluster4results2$fe[1:10]+0.001), horiz=TRUE, names.arg=cluster4results2$Description[1:10], border="grey94", cex.names=.9, col=cols[4],las=1,main='Cluster 4',xlim=c(0,max(log2(gooutput$fe+0.001))))
barplot(log2(cluster1results2$fe[1:10]+0.001), horiz=TRUE, names.arg=cluster1results2$Description[1:10], border="grey94", cex.names=.9, col=cols[1],las=1,main='Cluster 1',xlim=c(0,max(log2(gooutput$fe+0.001))))
barplot(log2(cluster3results2$fe[1:10]+0.001), horiz=TRUE, names.arg=cluster3results2$Description[1:10], border="grey94", cex.names=.9, col=cols[3],las=1,main='Cluster 3',xlim=c(0,max(log2(gooutput$fe+0.001))))
barplot(log2(cluster2results2$fe[1:10]+0.001), horiz=TRUE, names.arg=cluster2results2$Description[1:10], border="grey94", cex.names=.9, col=cols[2],las=1,main='Cluster 2',xlim=c(0,max(log2(gooutput$fe+0.001))))
barplot(log2(cluster5results2$fe[1:10]+0.001), horiz=TRUE, names.arg=cluster5results2$Description[1:10], border="grey94", cex.names=.9, col="grey55",las=1,main='Stable genes',xlim=c(0,max(log2(gooutput$fe+0.001))))
dev.off()
