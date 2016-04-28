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
dat<-subset(dat, msi_status %in% 'MSI-H')
dat2<-dcast(dat,locus ~ tumor_type, value.var="prop")
dat3<-dat2[complete.cases(dat2),]
#170,101 sites callable across all four MSI-H cancer samples

a<-dcast(dat,locus + gene ~ cancer_type, value.var="prop")
a2<-a[complete.cases(a),]
#158,912 sites callable across all four MSI-H and MSS samples

require(data.table)
fisherdat<-dcast(setDT(dat),locus ~ tumor_type, value.var=c("simple_gain_site","num_present"))
fisherdat2<-fisherdat[complete.cases(fisherdat),]
#results_per_cancer6<-results_per_cancer5[,c(1,2:12,13:23,24:34,156:166)]
unstable<-grep('simple_gain',colnames(fisherdat2))
all<-grep('num_present',colnames(fisherdat2))
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msih_cancer_specific_prop_test_results_021716.robj')
prop_results<-prop_results_msih
#system.time(prop_results<-apply(fisherdat2,1,function(x)prop.test(as.numeric(x[unstable]),as.numeric(x[all])))) #2.5 min
pvals<-unlist(lapply(prop_results,function(x)x$p.value))
props<-data.frame(matrix(unlist(lapply(prop_results,function(x)x$estimate)),nrow=dim(fisherdat2)[1],byrow=T))
#prop_results_msih<-prop_results
#save(prop_results_msih,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msih_cancer_specific_prop_test_results_021716.robj')
prop_sums<-apply(props,1,sum)
#74418/117507 63.3% of sites are stable
colnames(props)<-c('COAD','READ','UCEC','STAD')
unstablesites<-apply(fisherdat2, 1,function(x)sum(as.numeric(x[unstable])))
allcancers<-apply(fisherdat2,1,function(x)sum(as.numeric(x[all])))
propunstable<-unstablesites/allcancers
msihresults<-data.frame(locus=fisherdat2$locus,genomic_class=dat$genomic_class[match(fisherdat2$locus,dat$locus)],gene=dat$gene[match(fisherdat2$locus,dat$locus)],repeat_type=dat$repeat_type[match(fisherdat2$locus,dat$locus)],repeat_dna_sequence=dat$repeat_dna_sequence[match(fisherdat2$locus,dat$locus)],props,pvals=pvals)
msihresults2<-msihresults[which(propunstable>0.05),]
source('/net/shendure/vol1/home/hauser/scripts/shendure/miscellaneous/useful_functions.r')
require(qvalue)
msihresults3<-data.frame(msihresults2,qvals=qvalue(msihresults2$pvals)$qvalues)
msimat<-msihresults3[which(msihresults3$qvals<0.05),]
msimat2<-as.data.frame(msimat[,-c(1:5,10:11)])

write.table(msimat,'/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msih_cancer_specific_instability_loci.txt',quote=F,)

msihresults4<-msihresults3[which(msihresults3$qvals<1e-3),]

hits<-as.character(read.table('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/recurrent_msi_targets.txt',header=F,stringsAsFactors=F)[,1])

require("RColorBrewer")
blues<-colorRampPalette(brewer.pal(9,"Blues"))(256)
require(pheatmap)
rownames(msimat3)[50:51]<-paste0('chr',as.character(msihresults4$locus[50:51]))
rnames<-rownames(msimat3)
rownames(msimat3)<-msihresults4$locus
#pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_specific_loci_msi_cancer_distributions_021616.pdf')
anno<-msihresults4[,c(2,4)]
rownames(anno)<-rownames(msimat3)
anno$genomic_class<-revalue(anno$genomic_class,c("ncRNA_splicing"="splicing"))
cols<-brewer.pal(9, "Set2")
pheatmap(msimat3,col=blues,filename='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_specific_loci_msi_cancer_distributions_021816.pdf',clustering_method="ward.D2",cellwidth=10,fontsize_row=6,treeheight_row=0,treeheight_col=0,annotation_row=anno,annotation_legend=T,labels_row=rnames,annotation_names_row=FALSE,annotation_colors=list(genomic_class = c(exonic = cols[4],intronic = cols[5], intergenic = cols[6],splicing = cols[7],UTR5 = cols[8]),repeat_type=c(p1=cols[1],p2=cols[2],c=cols[3])))
dev.off()

# [1] "KIAA1598" "PICALM"   "TBX3"     "TSHR"     "DICER1"   "SPOP"
# [7] "ACVR2A"   "MAP3K13"  "UBR5"     "WRN" in cosmic

#require(gplots)
#msimat3<-as.matrix(msimat2)
#rownames(msimat3)<-as.character(msimat$gene)
#msimat3<-msimat3[which(msimat$qvals<1e-3),]
#pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_specific_loci_msi_cancer_distributions_021616_1e3.pdf',width=7*0.75,height=7*1.25)
#heatmap.2(msimat3,col=blues,key=TRUE,symkey=FALSE,symm=FALSE,keysize=0.5,trace='none',density.info="none",scale='none',cexRow=0.3,cexCol=0.5,lwid=c(3,4),lhei=c(2,4))
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

require(corrplot)
pcor <- cor(adat3[,-1,with=FALSE], method = "pearson")
scor <- cor(adat3[,-1,with=FALSE], method = "spearman")
summary(descrCor[upper.tri(descrCor)])

corrplot(descrCor, order = "alphabet", is.corr = TRUE, method = "color", type = "lower", tl.col = rgb(0, 0, 0),col = c(col2(200),col2(200)))

descrCor <- cor(msihresults2[,6:9], method = "spearman")

col2 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3","#92C5DE","#D1E5F0","#FFFFFF", "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msih_cancer_all_sufficient_site_spearman_correlations.pdf',width=7,height=7)
corrplot(scor, order = "hclust", is.corr = TRUE, method = "color", type = "lower", tl.col = rgb(0, 0, 0),col = c(rep('black',100),col2(100)),cl.lim=c(0,1))
dev.off()

fisherdat3<-fisherdat2[, 2:9, with = FALSE]
fisherdat3<-as.data.frame(fisherdat3)[which(propunstable>0.05),]

coadcol<-grep('COAD',colnames(fisherdat3))
coadsig<-apply(fisherdat3,1,function(x)fisher.test(matrix(c(x[coadcol[1]],x[coadcol[2]]-x[coadcol[1]],sum(x[-coadcol][1:3]),sum(x[-coadcol][4:6])-sum(x[-coadcol][1:3])),nrow=2,byrow=T),alternative="greater"))
coadpvals<-unlist(lapply(coadsig,function(x)x$p.value))
coadodds<-unlist(lapply(coadsig,function(x)x$estimate))
coadqvals<-qvalue(coadpvals)$qvalues
#1 < 0.05, 14 < 0.10; 

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
#3 at Q < 0.20

stadcol<-grep('STAD',colnames(fisherdat3))
stadsig<-apply(fisherdat3,1,function(x)fisher.test(matrix(c(x[stadcol[1]],x[stadcol[2]]-x[stadcol[1]],sum(x[-stadcol][1:3]),sum(x[-stadcol][4:6])-sum(x[-stadcol][1:3])),nrow=2,byrow=T),alternative="greater"))
stadpvals<-unlist(lapply(stadsig,function(x)x$p.value))
stadodds<-unlist(lapply(stadsig,function(x)x$estimate))
stadqvals<-qvalue(stadpvals)$qvalues
#194 significant at Q < 0.05

#defining cancer-specific sites as having frequencies of instability at that site three-fold higher in that cancer relative to the other three: 59 for UCEC, 758 for COAD, 1290 for READ, and 1230 for STAD
subdat<-fisherdat2[which(propunstable>0.05),]
subdat$gene<-dat$gene[match(subdat$locus,dat$locus)]
require(clusterProfiler)
require(dplyr)
require('biomaRt')
mart<-useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "jul2015.archive.ensembl.org")
ids <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","entrezgene","hgnc_symbol"),values=subdat$gene,mart= mart)
subdat$entrezgene<-as.factor(ids$entrezgene[match(subdat$gene,ids$hgnc_symbol)])
cluslist<-list(ucec=subdat$entrezgene[which(ucecodds>2)],coad=subdat$entrezgene[which(coadodds>3)],read=subdat$entrezgene[which(readodds>3)],stad=subdat$entrezgene[which(stadodds>3)])

cluslist<-list(ucec=subdat$entrezgene[which(ucecpvals<0.05)],coad=subdat$entrezgene[which(coadpvals<0.05)],read=subdat$entrezgene[which(readpvals<0.05)],stad=subdat$entrezgene[which(stadpvals<0.05)])

x=compareCluster(cluslist, fun="enrichGO",pvalueCutoff=0.05,OrgDb='org.Hs.eg.db')
x2<-simplify(x,cutoff=0.7,by="p.adjust",select_fun=min)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_all_go_enrichments_msih_cancers_021716_simplified_pvalthresh.pdf',width=7*0.625,height=7)
plot(x2,showCategory=NULL,font.size=6)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_all_go_enrichments_msih_cancers_021716_pvalthresh.pdf',width=7,height=7*1.5)
plot(x,showCategory=NULL,font.size=6)
dev.off()