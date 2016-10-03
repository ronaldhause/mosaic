load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')

#just MSI-H cancers
results_msi_per_cancer<-results
require(reshape2)
require(dplyr)
results_msi_per_cancer$prop<-results_msi_per_cancer$num_unstable/results_msi_per_cancer$num_present
results_msi_per_cancer$cancer_type<-paste(results_msi_per_cancer$tumor_type,results_msi_per_cancer$msi_status,sep='_')
results_msi_per_cancer2<-results_msi_per_cancer[-which(results_msi_per_cancer$num_present<(0.5*results_msi_per_cancer$total_num)),]
results_msi_per_cancer3<-subset(results_msi_per_cancer2, tumor_type %in% c('UCEC','READ','COAD','STAD'))
dat<-results_msi_per_cancer3
dat<-subset(dat, msi_status %in% 'MSI-H')
dat2<-dcast(dat,locus ~ tumor_type, value.var="prop")
dat3<-dat2[complete.cases(dat2),]
#132,535 sites callable across all four MSI-H cancer samples

require(data.table)
fisherdat<-dcast(setDT(dat),locus ~ tumor_type, value.var=c("num_unstable","num_present"))
fisherdat2<-fisherdat[complete.cases(fisherdat),]
unstable<-grep('num_unstable',colnames(fisherdat2))
all<-grep('num_present',colnames(fisherdat2))
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msih_cancer_specific_prop_test_results_071216.robj')
#system.time(prop_results<-apply(fisherdat2,1,function(x)prop.test(as.numeric(x[unstable]),as.numeric(x[all])))) #2.5 min
pvals<-unlist(lapply(prop_results,function(x)x$p.value))
props<-data.frame(matrix(unlist(lapply(prop_results,function(x)x$estimate)),nrow=dim(fisherdat2)[1],byrow=T))
#save(prop_results,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msih_cancer_specific_prop_test_results_071216.robj')
prop_sums<-apply(props,1,sum)
#100663/132535 76% of sites are stable
colnames(props)<-c('COAD','READ','UCEC','STAD')
unstablesites<-apply(fisherdat2, 1,function(x)sum(as.numeric(x[unstable])))
allcancers<-apply(fisherdat2,1,function(x)sum(as.numeric(x[all])))
propunstable<-unstablesites/allcancers
msihresults<-data.frame(locus=fisherdat2$locus,genomic_class=dat$genomic_class[match(fisherdat2$locus,dat$locus)],gene=dat$gene[match(fisherdat2$locus,dat$locus)],repeat_type=dat$repeat_type[match(fisherdat2$locus,dat$locus)],repeat_dna_sequence=dat$repeat_dna_sequence[match(fisherdat2$locus,dat$locus)],props,pvals=pvals)
msihresults2<-msihresults[which(propunstable>0.05),]
require(qvalue)
msihresults3<-data.frame(msihresults2,qvals=qvalue(msihresults2$pvals)$qvalues)
msihresults3<-msihresults3[order(msihresults3$qvals),]
msimat<-msihresults3[which(msihresults3$qvals<0.05),]
msimat2<-as.data.frame(msimat[,-c(1:5,10:11)])

write.table(msimat,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msih_cancer_specific_instability_loci_071216.txt',quote=F,row.names=F)

msihresults4<-subset(msihresults3, !duplicated(msihresults3$gene))
msihresults4<-msihresults4[1:50,]
msimat3<-as.data.frame(msihresults4[,-c(1:5,10:11)])

#figure 4d - cancer-specific instability
require("RColorBrewer")
require(plyr)
blues<-colorRampPalette(brewer.pal(9,"Blues"))(256)
require(pheatmap)
rownames(msimat3)<-msihresults4$locus
anno<-msihresults4[,c(2,4)]
rownames(anno)<-rownames(msimat3)
anno$genomic_class<-revalue(anno$genomic_class, c("downstream"="noncoding","ncRNA_exonic"="noncoding","ncRNA_intronic"="noncoding","ncRNA_splicing"="noncoding","ncRNA_UTR3"="noncoding","ncRNA_UTR5"="noncoding","ncRNA_UTR5;ncRNA_UTR3"="noncoding","upstream"="noncoding","upstream;downstream"="noncoding","UTR5;UTR3"="noncoding","intergenic"="noncoding"))
rnames<-as.character(msihresults4$gene)
rnames[grep('dist',rnames)]<-paste0("chr",as.character(msihresults4$locus[grep('dist',rnames)]))

cols<-brewer.pal(8, "Set2")
cols2<-brewer.pal(9, "Set1")
pheatmap(msimat3,col=blues,filename='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_specific_loci_msi_cancer_distributions_072616.pdf',clustering_method="ward.D2",cellwidth=10,fontsize_row=10,treeheight_row=0,treeheight_col=0,annotation_row=anno,annotation_legend=T,labels_row=rnames,annotation_names_row=FALSE,annotation_colors=list(genomic_class = c(exonic = cols[4],intronic = cols[5], noncoding = cols[6],splicing = cols[7],UTR5 = cols[8],UTR3=cols2[1]),repeat_type=c(p1=cols[1],p2=cols[2],c=cols[3])))
dev.off()

#MSI-H and MSS cancers

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

require(corrplot)
scor <- cor(adat3[,-1,with=FALSE], method = "spearman")

fisherdat3<-fisherdat2[, 2:9, with = FALSE]
fisherdat3<-as.data.frame(fisherdat3)[which(propunstable>0.05),]

coadcol<-grep('COAD',colnames(fisherdat3))
coadsig<-apply(fisherdat3,1,function(x)fisher.test(matrix(c(x[coadcol[1]],x[coadcol[2]]-x[coadcol[1]],sum(x[-coadcol][1:3]),sum(x[-coadcol][4:6])-sum(x[-coadcol][1:3])),nrow=2,byrow=T),alternative="greater"))
coadpvals<-unlist(lapply(coadsig,function(x)x$p.value))
coadodds<-unlist(lapply(coadsig,function(x)x$estimate))
coadqvals<-qvalue(coadpvals)$qvalues
#716 at Q < 0.05

readcol<-grep('READ',colnames(fisherdat3))
readsig<-apply(fisherdat3,1,function(x)fisher.test(matrix(c(x[readcol[1]],x[readcol[2]]-x[readcol[1]],sum(x[-readcol][1:3]),sum(x[-readcol][4:6])-sum(x[-readcol][1:3])),nrow=2,byrow=T),alternative="greater"))
readpvals<-unlist(lapply(readsig,function(x)x$p.value))
readodds<-unlist(lapply(readsig,function(x)x$estimate))
readqvals<-qvalue(readpvals)$qvalues

uceccol<-grep('UCEC',colnames(fisherdat3))
ucecsig<-apply(fisherdat3,1,function(x)fisher.test(matrix(c(x[uceccol[1]],x[uceccol[2]]-x[uceccol[1]],sum(x[-uceccol][1:3]),sum(x[-uceccol][4:6])-sum(x[-uceccol][1:3])),nrow=2,byrow=T),alternative="greater"))
ucecpvals<-unlist(lapply(ucecsig,function(x)x$p.value))
ucecodds<-unlist(lapply(ucecsig,function(x)x$estimate))
ucecqvals<-qvalue(ucecpvals)$qvalues
#10 at Q < 0.10

stadcol<-grep('STAD',colnames(fisherdat3))
stadsig<-apply(fisherdat3,1,function(x)fisher.test(matrix(c(x[stadcol[1]],x[stadcol[2]]-x[stadcol[1]],sum(x[-stadcol][1:3]),sum(x[-stadcol][4:6])-sum(x[-stadcol][1:3])),nrow=2,byrow=T),alternative="greater"))
stadpvals<-unlist(lapply(stadsig,function(x)x$p.value))
stadodds<-unlist(lapply(stadsig,function(x)x$estimate))
stadqvals<-qvalue(stadpvals)$qvalues
#650 significant at Q < 0.05

#defining cancer-specific sites as having frequencies of instability at that site three-fold higher in that cancer relative to the other three
subdat<-fisherdat2[which(propunstable>0.05),]
subdat$gene<-dat$gene[match(subdat$locus,dat$locus)]
require(clusterProfiler)
require(dplyr)
require('biomaRt')
mart<-useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "jul2015.archive.ensembl.org")
ids <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","entrezgene","hgnc_symbol"),values=subdat$gene,mart= mart)
subdat$entrezgene<-as.factor(ids$entrezgene[match(subdat$gene,ids$hgnc_symbol)])

cluslist<-list(ucec=subdat$entrezgene[which(ucecodds>1.5)],coad=subdat$entrezgene[which(coadodds>3)],read=subdat$entrezgene[which(readodds>3)],stad=subdat$entrezgene[which(stadodds>3)])

x=compareCluster(cluslist, fun="enrichGO",pvalueCutoff=0.05,OrgDb='org.Hs.eg.db')
x2<-simplify(x,cutoff=0.7,by="p.adjust",select_fun=min)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_all_go_enrichments_msih_cancers_071216_simplified_oddsthresh.pdf',width=7*0.625,height=7)
plot(x2,showCategory=NULL,font.size=6)
dev.off()