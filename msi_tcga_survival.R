library('TCGA2STAT')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')

sampleinfo$num_unstable<-sampledat$num_unstable_raw
sampledat<-sampledat[complete.cases(sampledat),]

cancers<-levels(sampledat$tumor_type)[-c(1,11)]
cancers[5]<-'GBM'
cancers[14]<-'OV'
cancers[18]<-'LGG'

clindatlist<-list()
rnaseqlist<-list()
for (i in 1:length(cancers)){
	subset.rnaseq<-getTCGA(disease=cancers[i],data.type='RNASeq2',type='RPKM',clinical=TRUE)
	rnaseqlist[[i]]<-subset.rnaseq
	clindatlist[[i]]<-subset.rnaseq$clinical
}

#a<-getTCGA(disease='OV',data.type='RNASeq2',type='RPKM',clinical=TRUE)
#rnaseqlist[[14]]<-a
#clindatlist[[14]]<-a$clinical

for (i in 1:length(cancers)){
	clindatlist[[i]]<-cbind(clindatlist[[i]],tumor_type=rep(cancers[i],dim(clindatlist[[i]])[1]))
	clindatlist[[i]]<-cbind(clindatlist[[i]],sample_name=rownames(clindatlist[[i]]))
}

#ucecclin<-read.delim('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/ucec_tcga_clinical_data.tsv',header=T)
#clindatlist[[1]]<-cbind(clindatlist[[1]],yearstobirth=ucecclin[match(clindatlist[[1]][,16],ucecclin$Patient.ID),8],gender=rep('female',548))

save(rnaseqlist,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_rnaseq_data.robj')
save(clindatlist,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_clinical_data.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_rnaseq_data.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_clinical_data.robj')

clindat<-Reduce(function(...) merge(..., all=T), clindatlist)
clindat$tumor_type<-revalue(clindat$tumor_type,c('OVmu'='OV'))
clindat$num_unstable<-sampledat$num_unstable_raw[match(clindat$sample_name,sampledat$sample_name)]
clindat$msi_status<-sampledat$msi_status_tree_scaled[match(clindat$sample_name,sampledat$sample_name)]
clindat$yearstobirth<-as.numeric(as.matrix(clindat$yearstobirth))
clindat$peak_avg<-sampledat$peak_avg[match(clindat$sample_name,sampledat$sample_name)]
clindat<-clindat[-which(clindat$num_unstable>10000),]

#correlation between age and num_unstable = -0.019, p = 0.2188; no sig difference between gender (p = 0.91)
clindat<-clindat[-which(is.na(clindat$num_unstable)==T),]
clindat$g_quartile<-factor(cut(clindat$num_unstable, quantile(clindat$num_unstable), include.lowest = TRUE),labels = LETTERS[1:4])
clindat$s_quartile<-factor(NA,levels=LETTERS[1:4])
clindat$s_tertile<-factor(NA,levels=LETTERS[1:3])
clindat$half<-factor(NA,levels=c('low','high'))
clindat$os<-rowSums(cbind(as.numeric(as.matrix(clindat$daystodeath)),as.numeric(as.matrix(clindat$daystolastfollowup))),na.rm=T)

for (i in 1:length(cancers)){
	subset<-clindat[which(clindat$tumor_type==cancers[i]),]
	subset$s_quartile<-factor(cut(subset$num_unstable, quantile(subset$num_unstable), include.lowest = TRUE),labels = LETTERS[1:4])
	clindat$s_quartile[which(is.na(match(clindat$sample_name,subset$sample_name))==F)]<-na.omit(subset$s_quartile[match(clindat$sample_name,subset$sample_name)])
	subset$half<-factor(cut(subset$num_unstable, quantile(subset$num_unstable)[c(1,3,5)], include.lowest = TRUE),labels = c('low','high'))
	clindat$half[which(is.na(match(clindat$sample_name,subset$sample_name))==F)]<-na.omit(subset$half[match(clindat$sample_name,subset$sample_name)])
	subset$s_tertile<-factor(cut(subset$num_unstable, quantile(subset$num_unstable,probs=seq(0,1,1/3)), include.lowest = TRUE),labels = LETTERS[1:3])
	clindat$s_tertile[which(is.na(match(clindat$sample_name,subset$sample_name))==F)]<-na.omit(subset$s_tertile[match(clindat$sample_name,subset$sample_name)])
}

require(plyr)
clindat$vitalstatus<-as.numeric(as.matrix(clindat$vitalstatus))
clindat$os<-as.numeric(as.matrix(clindat$os))
clindat$pathologicstage<-revalue(clindat$pathologicstage,c('stage iia'='stage ii','stage iib'='stage ii','stage iic'='stage ii','stage iiia'='stage iii','stage iiib'='stage iii','stage iiic'='stage iii','stage iva'='stage iv','stage ia'='stage i','stage ivb'='stage iv','stage ib'='stage i','stage ivc'='stage iv','stage x' = NA,'i/ii nos'=NA,'stage 0'=NA))
clindat$residualtumor<-revalue(clindat$residualtumor,c('rx'=NA))

clindat2<-clindat[-intersect(which(clindat$msi_status=='MSI-H'),which(clindat$tumor_type %in% c('COAD','READ','UCEC','STAD'))),]
clindat2$g_quartile<-factor(cut(clindat2$num_unstable, quantile(clindat2$num_unstable), include.lowest = TRUE),labels = LETTERS[1:4])
clindat2$s_quartile<-factor(NA,levels=LETTERS[1:4])
clindat2$s_tertile<-factor(NA,levels=LETTERS[1:3])
clindat2$half<-factor(NA,levels=c('low','high'))
clindat2$os<-rowSums(cbind(as.numeric(as.matrix(clindat2$daystodeath)),as.numeric(as.matrix(clindat2$daystolastfollowup))),na.rm=T)

for (i in 1:length(cancers)){
	subset<-clindat2[which(clindat2$tumor_type==cancers[i]),]
	subset$s_quartile<-factor(cut(subset$num_unstable, quantile(subset$num_unstable), include.lowest = TRUE),labels = LETTERS[1:4])
	clindat2$s_quartile[which(is.na(match(clindat2$sample_name,subset$sample_name))==F)]<-na.omit(subset$s_quartile[match(clindat2$sample_name,subset$sample_name)])
	subset$half<-factor(cut(subset$num_unstable, quantile(subset$num_unstable)[c(1,3,5)], include.lowest = TRUE),labels = c('low','high'))
	clindat2$half[which(is.na(match(clindat2$sample_name,subset$sample_name))==F)]<-na.omit(subset$half[match(clindat2$sample_name,subset$sample_name)])
	subset$s_tertile<-factor(cut(subset$num_unstable, quantile(subset$num_unstable,probs=seq(0,1,1/3)), include.lowest = TRUE),labels = LETTERS[1:3])
	clindat2$s_tertile[which(is.na(match(clindat2$sample_name,subset$sample_name))==F)]<-na.omit(subset$s_tertile[match(clindat2$sample_name,subset$sample_name)])
}

a<-lm(clindat$num_unstable~clindat$yearstobirth+clindat$gender+clindat$msi_status+clindat$tumor_type)
b<-lm(clindat$daystodeath~clindat$yearstobirth+clindat$gender+clindat$msi_status+clindat$tumor_type)

require(party)
full_msi_model<-coxph(formula=Surv(clindat$os,clindat$vitalstatus)~clindat$msi_status+as.numeric(as.matrix(clindat$yearstobirth))+clindat$gender+clindat$tumor_type+clindat$radiationtherapy+clindat$pathologicstage) #P = 0.08
full_msi_model_zph<-cox.zph(full_msi_model)
full_msi_surv_tree <- ctree(formula=Surv(clindat$os,clindat$vitalstatus)~factor(clindat$msi_status)+as.numeric(as.matrix(clindat$yearstobirth))+clindat$gender+clindat$tumor_type+clindat$radiationtherapy+clindat$pathologicstage)

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/full_msi_survival_all_cancer_full_model_tree.pdf',width=7*4,height=7*1.5)
plot(full_msi_surv_tree)
dev.off()

advanced<-subset(clindat,pathologicstage %in% c('stage iii','stage iv'))
coxph(formula=Surv(advanced$os,advanced$vitalstatus)~advanced$msi_status+advanced$tumor_type+advanced$pathologicstage) #P = 0.070

coxph(formula=Surv(clindat$os,clindat$vitalstatus)~clindat$msi_status+clindat$tumor_type) #P = 0.047
coxph(formula=Surv(clindat$os,clindat$vitalstatus)~clindat$msi_status+strata(clindat$tumor_type)) #P = 0.03

coxph(formula=Surv(clindat$os,clindat$vitalstatus)~clindat$msi_status+as.numeric(as.matrix(clindat$yearstobirth))+clindat$gender+clindat$tumor_type+clindat$radiationtherapy+clindat$pathologicstage) #P = 0.082
coxph(formula=Surv(clindat$os,clindat$vitalstatus)~clindat$msi_status+as.numeric(as.matrix(clindat$yearstobirth))+clindat$gender+strata(clindat$tumor_type)+strata(clindat$radiationtherapy)+clindat$pathologicstage) #P = 0.083 after strata

coxph(formula=Surv(clindat2$os,clindat2$vitalstatus)~clindat2$num_unstable+clindat2$tumor_type) #P = 0.062

coxph(formula=Surv(clindat$os,clindat$vitalstatus)~clindat$num_unstable+as.numeric(as.matrix(clindat$yearstobirth))+clindat$gender+clindat$tumor_type+clindat$radiationtherapy+clindat$pathologicstage) #P = 0.49

fullmod<-coxph(formula=Surv(clindat$os,clindat$vitalstatus)~clindat$num_unstable+as.numeric(as.matrix(clindat$yearstobirth))+clindat$gender+strata(clindat$tumor_type)+strata(clindat$radiationtherapy)+clindat$pathologicstage) #P = 0.96

fullmod<-coxph(formula=Surv(clindat$os,clindat$vitalstatus)~clindat$peak_avg+as.numeric(as.matrix(clindat$yearstobirth))+clindat$gender+strata(clindat$tumor_type)+strata(clindat$radiationtherapy)+clindat$pathologicstage) #P = 0.93


cox.zph(fullmod)

coadread<-subset(clindat,tumor_type %in% c('COAD','READ'))
coadread2<-subset(clindat2,tumor_type %in% c('COAD','READ'))
msicancers<-subset(clindat,tumor_type %in% c('COAD','READ','UCEC','STAD'))
msicancers2<-subset(clindat2,tumor_type %in% c('COAD','READ','UCEC','STAD'))

lgg<-subset(clindat,tumor_type %in% 'LGG')
lgg<-subset(lgg,s_tertile %in% c('A','C'))

lgg$cluster<-kmeans(lgg$num_unstable,2)$cluster

coxph(formula=Surv(lgg$os,lgg$vitalstatus)~lgg$num_unstable+lgg$yearstobirth+lgg$gender+lgg$radiationtherapy) #P = 0.067

coxph(formula=Surv(lgg$os,lgg$vitalstatus)~lgg$s_tertile+lgg$yearstobirth+lgg$gender+lgg$radiationtherapy) #P = 0.067


lihc<-subset(clindat,tumor_type %in% 'LIHC')
lihc$cluster<-kmeans(lihc$num_unstable,2)$cluster
coxph(formula=Surv(lihc$os,lihc$vitalstatus)~lihc$cluster+lihc$yearstobirth+lihc$gender+lihc$radiationtherapy) #P = 0.052

ucec<-subset(clindat,tumor_type %in% 'UCEC')
#ucec<-subset(ucec,s_quartile %in% c('A','D'))
ucec$cluster<-kmeans(ucec$num_unstable,2)$cluster
coxph(formula=Surv(ucec$os,ucec$vitalstatus)~ucec$s_quartile+ucec$yearstobirth+ucec$radiationtherapy) #P = 0.052
coxph(formula=Surv(ucec$os,ucec$vitalstatus)~ucec$cluster+ucec$yearstobirth+ucec$radiationtherapy) #P = 0.052

gbm<-subset(clindat,tumor_type %in% 'GBM')
gbm2<-subset(gbm,s_quartile %in% c('A','D'))
gbm$cluster<-kmeans(gbm$num_unstable,2)$cluster
coxph(formula=Surv(gbm2$os,gbm2$vitalstatus)~gbm2$s_quartile+gbm2$yearstobirth+gbm2$radiationtherapy) #P = 0.052
coxph(formula=Surv(gbm$os,gbm$vitalstatus)~gbm$cluster+gbm$yearstobirth+gbm$radiationtherapy) #P = 0.052

coxph(formula=Surv(coadread$os,coadread$vitalstatus)~as.numeric(coadread$num_unstable)+coadread$tumor_type) #P = 0.23
coxph(formula=Surv(coadread2$os,coadread2$vitalstatus)~coadread2$num_unstable+coadread2$tumor_type) #P = 0.44
coxph(formula=Surv(msicancers$os,msicancers$vitalstatus)~msicancers$num_unstable+msicancers$tumor_type) #P = 0.0368
coxph(formula=Surv(msicancers$os,msicancers$vitalstatus)~msicancers$msi_status+msicancers$tumor_type) #P = 0.029
coxph(formula=Surv(msicancers$os,msicancers$vitalstatus)~msicancers$msi_status+msicancers$tumor_type+msicancers$yearstobirth+msicancers$gender+msicancers$radiationtherapy+msicancers$pathologicstage) #P = 0.021

coxph(formula=Surv(msicancers2$os,msicancers2$vitalstatus)~msicancers2$num_unstable+msicancers2$tumor_type+msicancers2$yearstobirth+msicancers2$gender+msicancers2$radiationtherapy+msicancers2$pathologicstage) #P = 0.209

coxph(formula=Surv(coadread$os,coadread$vitalstatus)~coadread$msi_status+coadread$tumor_type + as.numeric(as.matrix(coadread$yearstobirth))+coadread$gender+coadread$tumor_type+coadread$radiationtherapy+coadread$pathologicstage) #P = 0.32

coxph(formula=Surv(coadread$os,coadread$vitalstatus)~as.numeric(coadread$num_unstable)+coadread$tumor_type + as.numeric(as.matrix(coadread$yearstobirth))+coadread$gender+coadread$tumor_type+coadread$radiationtherapy+coadread$pathologicstage) #P = 0.32
coxph(formula=Surv(coadread2$os,coadread2$vitalstatus)~coadread2$num_unstable+coadread2$tumor_type+as.numeric(as.matrix(coadread2$yearstobirth))+coadread2$gender+coadread2$tumor_type+coadread2$radiationtherapy+coadread2$pathologicstage) #P = 0.52

b<-coxph(formula=Surv(msicancers$os,msicancers$vitalstatus)~msicancers$num_unstable+msicancers$tumor_type+as.numeric(as.matrix(msicancers$yearstobirth))+msicancers$gender+msicancers$tumor_type+msicancers$radiationtherapy+msicancers$pathologicstage) #P = 0.0377

cox.zph(b)

#Figure 6A
a<-coxph(formula=Surv(msicancers$os/365,msicancers$vitalstatus)~msicancers$msi_status+msicancers$tumor_type+as.numeric(as.matrix(msicancers$yearstobirth))+msicancers$gender+msicancers$tumor_type+msicancers$radiationtherapy+msicancers$pathologicstage) #P = 0.021

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/fig6a_msi_cancers_mss_vs_msi.pdf',width=7*1.5,height=7)
plot(survfit(formula=Surv(msicancers$os/365,msicancers$vitalstatus)~msicancers$msi_status), col=cols[1:2], lty=1, lwd=c(2,2), xlab="Years", ylab="Probability of survival",bty='n')
legend(0.01 * max(msicancers$os/365), 0.2, col = cols[1:2], lty = 1, lwd=2, legend = c('MSI-H','MSS'), bty='n')
text(0.1 * max(msicancers$os/365), 0.25, paste("p =", as.character(round(summary(a)$coef[1,5], 4))), cex=0.8)
dev.off()

#Figure 6B
b<-coxph(formula=Surv(msicancers$os,msicancers$vitalstatus)~msicancers$num_unstable+msicancers$tumor_type+as.numeric(as.matrix(msicancers$yearstobirth))+msicancers$gender+msicancers$tumor_type+msicancers$radiationtherapy+msicancers$pathologicstage) #P = 0.0377

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/fig6b_msi_cancers_num_unstable.pdf',width=7*1.5,height=7)
plot(survfit(formula=Surv(msicancers$os/365,msicancers$vitalstatus)~msicancers$s_quartile), col=cols[1:4], lty=1, lwd=2, xlab="Years", ylab="Probability of survival", bty='n')
legend(0.01 * max(msicancers$os/365), 0.2, col=cols[1:4], lty = 1, lwd=2, legend=c('1st','2nd','3rd','4th'), bty='n')
text(0.1 * max(msicancers$os/365), 0.25, paste("p =", as.character(round(summary(b)$coef[1,5], 4))), cex=0.8)
dev.off()

mypamr.plotsurvival <- function (group, survival.time, censoring.status, cols, lty, lwd, main) 
{
  require(survival)
  n.class <- length(unique(group))
  junk <- survfit(Surv(survival.time, censoring.status) ~ as.factor(group))
  junk2 <- coxph(Surv(survival.time, censoring.status) ~ as.factor(group))
  pv <- 1 - pchisq(2 * (junk2$loglik[2] - junk2$loglik[1]), 
                   df = n.class - 1)
  plot(junk, col = cols, xlab = "Years", ylab = "Probability of survival", lty=lty, lwd=lwd,main=main)
  legend(0.01 * max(survival.time), 0.2, col = cols, lty = lty, lwd=lwd, legend = as.character(1:n.class), box.lwd=1)
  text(0.1 * max(survival.time), 0.25, paste("p =", as.character(round(pv, 4))), cex=0.8)
  return()
}

coxph(formula=Surv(msicancers2$os,msicancers2$vitalstatus)~msicancers2$num_unstable+msicancers2$tumor_type+as.numeric(as.matrix(msicancers2$yearstobirth))+msicancers2$gender+msicancers2$tumor_type+msicancers2$radiationtherapy+msicancers2$pathologicstage) #P = 0.362

coxph(formula=Surv(msicancers$os,msicancers$vitalstatus)~msicancers$num_unstable+msicancers$tumor_type) #P = 0.0368
coxph(formula=Surv(msicancers$os,msicancers$vitalstatus)~msicancers$msi_status+msicancers$tumor_type) #P = 0.029



coxph(formula=Surv(coadread$os,coadread$vitalstatus)~coadread$num_unstable+as.numeric(as.matrix(clindat$yearstobirth))+clindat$gender)


survival3<-list()
for (i in 1:length(cancers)){
	subset<-clindat2[which(clindat2$tumor_type==cancers[i]),]
	#subset2<-subset(subset,s_tertile %in% c('A','C'), select=c('os','vitalstatus','s_tertile','yearstobirth','radiationtherapy','pathologicstage','gender'))
	#subset2$s_tertile<-droplevels(subset2$s_tertile)
	subset2<-subset(subset,select=c('os','vitalstatus','half','yearstobirth','radiationtherapy','pathologicstage','gender'))
	subset2<-transform(subset2,os=as.numeric(os),vitalstatus=as.numeric(vitalstatus),yearstobirth=as.numeric(yearstobirth),radiationtherapy=as.factor(radiationtherapy),pathologicstage=as.factor(pathologicstage),gender=as.factor(gender))
	covarsel<-apply(subset2,2,function(x)length(unique(x)))
	subset3<-subset2[,which(covarsel>1)]
	survival3[[i]]<-coxph(formula=Surv(os,vitalstatus)~.,data=subset3)
}

names(survival3)<-cancers

survival<-list()
for (i in 1:length(cancers)){
	subset<-clindat2[which(clindat2$tumor_type==cancers[i]),]
	subset2<-subset(subset,select=c('os','vitalstatus','num_unstable','yearstobirth','radiationtherapy','pathologicstage','gender'))
	subset2<-transform(subset2,os=as.numeric(os),vitalstatus=as.numeric(vitalstatus),yearstobirth=as.numeric(yearstobirth),radiationtherapy=as.factor(radiationtherapy),pathologicstage=as.factor(pathologicstage),gender=as.factor(gender))
	covarsel<-apply(subset2,2,function(x)length(unique(x)))
	subset3<-subset2[,which(covarsel>1)]
	survival[[i]]<-coxph(formula=Surv(os,vitalstatus)~.,data=subset3)
}

names(survival)<-cancers

binarysurvival<-list()
k<-1
for (i in c(1,2,3,13)){
	subset<-clindat2[which(clinda2$tumor_type==cancers[i]),]
	subset2<-subset(subset,select=c('os','vitalstatus','msi_status','yearstobirth','radiationtherapy','pathologicstage','gender'))
	subset2<-transform(subset2,os=as.numeric(os),vitalstatus=as.numeric(vitalstatus),msi_status=as.factor(msi_status),yearstobirth=as.numeric(yearstobirth),radiationtherapy=as.factor(radiationtherapy),pathologicstage=as.factor(pathologicstage),gender=as.factor(gender))
	covarsel<-apply(subset2,2,function(x)length(unique(x)))
	subset3<-subset2[,which(covarsel>1)]
	binarysurvival[[k]]<-coxph(formula=Surv(os,vitalstatus)~.,data=subset3)
	k<-k+1
}

names(binarysurvival)<-cancers[c(1,2,3,13)] #0.20 for UCEC, 1 for READ, 0.409 for COAD, 0.054 for STAD

binarysurvivalraw<-list()
k<-1
for (i in c(1,2,3,13)){
	subset<-clindat[which(clindat$tumor_type==cancers[i]),]
	subset2<-subset(subset,select=c('os','vitalstatus','msi_status','yearstobirth','radiationtherapy','pathologicstage','gender'))
	binarysurvivalraw[[k]]<-coxph(formula=Surv(os,vitalstatus)~msi_status,data=subset2)
	k<-k+1
}

names(binarysurvivalraw)<-cancers[c(1,2,3,13)] #0.20 for UCEC, 1 for READ, 0.409 for COAD, 0.054 for STAD

msisurvival<-list()
for (i in 1:length(cancers)){
	subset<-clindat2[which(clindat2$tumor_type==cancers[i]),]
	subset2<-subset(subset,select=c('os','vitalstatus','num_unstable','yearstobirth','radiationtherapy','pathologicstage','gender'))
	subset2<-transform(subset2,os=as.numeric(os),vitalstatus=as.numeric(vitalstatus),yearstobirth=as.numeric(yearstobirth),radiationtherapy=as.factor(radiationtherapy),pathologicstage=as.factor(pathologicstage),gender=as.factor(gender))
	covarsel<-apply(subset2,2,function(x)length(unique(x)))
	subset3<-subset2[,which(covarsel>1)]
	msisurvival[[i]]<-coxph(formula=Surv(os,vitalstatus)~.,data=subset3)
}

names(msisurvival)<-cancers

survivalraw<-list()
for (i in 1:length(cancers)){
	subset<-clindat[which(clindat$tumor_type==cancers[i]),]
	subset2<-subset(subset,select=c('os','vitalstatus','num_unstable','yearstobirth','radiationtherapy','pathologicstage','gender'))
	survivalraw[[i]]<-coxph(formula=Surv(os,vitalstatus)~num_unstable,data=subset2)
}

names(survivalraw)<-cancers

msisurvivalraw<-list()
for (i in 1:length(cancers)){
	subset<-clindat2[which(clindat2$tumor_type==cancers[i]),]
	subset2<-subset(subset,select=c('os','vitalstatus','num_unstable','yearstobirth','radiationtherapy','pathologicstage','gender'))
	msisurvivalraw[[i]]<-coxph(formula=Surv(os,vitalstatus)~num_unstable,data=subset2)
}

names(msisurvivalraw)<-cancers

classes<-apply(clindat,2,class)

test<-coxph(Surv(clindat$os, as.numeric(as.matrix(clindat$vitalstatus))) ~ clindat$gender + clindat$yearstobirth + clindat$tumor_type + clindat$s_tertile)

plotdat <- data.frame(msi=clindat$msi_status,status=clindat$vitalstatus,os=clindat$os,g_quartiles=clindat$g_quartile,s_quartiles=clindat$s_quartile,s_tertiles=clindat$s_tertile,tumor_type=clindat$tumor_type,half=clindat$half,num_unstable=clindat$num_unstable,age=clindat$yearstobirth,gender=clindat$gender,sample_name=clindat$sample_name)
plotdat <- plotdat[order(plotdat[,1]),]
plotdat <- plotdat[complete.cases(plotdat),]
plotdat$status <- as.numeric(as.matrix(plotdat$status))
plotdat$os <- as.numeric(as.matrix(plotdat$os))

plotdat2 <- data.frame(msi=clindat2$msi_status,status=clindat2$vitalstatus,os=clindat2$os,g_quartiles=clindat2$g_quartile,s_quartiles=clindat2$s_quartile,s_tertiles=clindat2$s_tertile,tumor_type=clindat2$tumor_type,half=clindat2$half,num_unstable=clindat2$num_unstable,age=clindat2$yearstobirth,gender=clindat2$gender,sample_name=clindat2$sample_name)
plotdat2 <- plotdat2[order(plotdat2[,1]),]
plotdat2 <- plotdat2[complete.cases(plotdat2),]
plotdat2$status <- as.numeric(as.matrix(plotdat2$status))
plotdat2$os <- as.numeric(as.matrix(plotdat2$os))

junk <- survfit(Surv(plotdat$os, plotdat$status) ~ plotdat$tumor_type+plotdat$s_quartile+plotdat$gender)

mypamr.plotsurvival <- function (group, survival.time, censoring.status, cols, lty, lwd, main) 
{
  require(survival)
  n.class <- length(unique(group))
  junk <- survfit(Surv(survival.time, censoring.status) ~ as.factor(group))
  junk2 <- coxph(Surv(survival.time, censoring.status) ~ as.factor(group))
  #pv <- 1 - pchisq(2 * (junk2$loglik[2] - junk2$loglik[1]), df = n.class - 1)
  pv <- summary(survival[[i]])$coef[1,5]
  plot(junk, col = cols, xlab = "Years", ylab = "Probability of survival", lty=lty, lwd=lwd,main=main)
  legend(0.01 * max(survival.time), 0.2, col = cols[1:2], lty = lty, lwd=c(2,2), legend = c('top','bottom'), bty='n')
  text(0.1 * max(survival.time), 0.25, paste("p =", as.character(round(pv, 4))), cex=0.8)
  return()
}

plot(survfit(formula=Surv(msicancers$os/365,msicancers$vitalstatus)~msicancers$s_quartile), col=cols[1:4], lty=1, lwd=2, xlab="Years", ylab="Probability of survival", bty='n')
legend(0.01 * max(msicancers$os/365), 0.2, col=cols[1:4], lty = 1, lwd=2, legend=c('1st','2nd','3rd','4th'), bty='n')
text(0.1 * max(msicancers$os/365), 0.25, paste("p =", as.character(round(summary(b)$coef[1,5], 4))), cex=0.8)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/full_msi_survival_all_cancers.pdf',width=7*2,height=7)
mypamr.plotsurvival(plotdat[,1], plotdat[,3]/365, plotdat[,2], c(colors()[c(153,220)]), lty=c(2,1), lwd=c(2,2))
dev.off()

require(RColorBrewer)
cols<-brewer.pal(9, "Set1")

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/full_msi_survival_continuous_global.pdf',width=7*2,height=7)
mypamr.plotsurvival(plotdat[,4], plotdat[,3]/365, plotdat[,2], cols[1:4], lty=c(2,1), lwd=c(2,2))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/full_msi_survival_continuous_tumor_specific.pdf',width=7*2,height=7)
mypamr.plotsurvival(plotdat[,5], plotdat[,3]/365, plotdat[,2], cols[1:4], lty=c(2,1), lwd=c(2,2))
dev.off()

mss.plotdat<-plotdat[which(plotdat$msi=='MSS'),]

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_only_msi_survival_continuous_all_cancers.pdf',width=7*2,height=7)
mypamr.plotsurvival(mss.plotdat[,5], mss.plotdat[,3]/365, mss.plotdat[,2], cols[1:4], lty=c(2,1), lwd=c(2,2))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_survival_continuous_all_cancers_split_by_cancer_042516.pdf',width=20,height=12)
par(mfrow=c(3,5),mar=c(4,4,2,1))
for (i in 1:18){
subset<-plotdat[which(plotdat$tumor_type==cancers[i]),]
if (survival[[i]]$n>0 & survival[[i]]$nevent>10){
mypamr.plotsurvival(subset[,8], subset[,3]/365, subset[,2], cols[1:4], lty=c(2,1), lwd=c(2,2),main=cancers[i])
}
}
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_survival_continuous_all_cancers_split_by_cancer_2_042516.pdf',width=12,height=12)
par(mfrow=c(3,3),mar=c(4,4,2,1))
for (i in 10:18){
subset<-plotdat[which(plotdat$tumor_type==cancers[i]),]
mypamr.plotsurvival(subset[,8], subset[,3]/365, subset[,2], cols[1:4], lty=c(2,1), lwd=c(2,2),main=cancers[i])
}
dev.off()

plotdatsub2<-subset(plotdat2,s_quartiles %in% c('A','D'))
oldplotdat<-plotdat
oldplotdat2<-plotdat2

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_only_msi_survival_continuous_top_vs_bottom_all_cancers_split_by_cancer.pdf',width=12,height=12)
par(mfrow=c(3,3),mar=c(4,4,2,1))
for (i in 1:9){
subset<-plotdat[which(plotdat$tumor_type==cancers[i]),]
mypamr.plotsurvival(subset[,5], subset[,3]/365, subset[,2], cols[1:2], lty=c(2,1), lwd=c(2,2),main=cancers[i])
}
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_only_msi_survival_continuous_top_vs_bottom_all_cancers_split_by_cancer_2.pdf',width=12,height=12)
par(mfrow=c(3,3),mar=c(4,4,2,1))
for (i in 10:18){
subset<-plotdat[which(plotdat$tumor_type==cancers[i]),]
mypamr.plotsurvival(subset[,5], subset[,3]/365, subset[,2], cols[1:2], lty=c(2,1), lwd=c(2,2),main=cancers[i])
}
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_only_msi_survival_continuous_all_cancers_split_by_cancer.pdf',width=12,height=12)
par(mfrow=c(3,3),mar=c(4,4,2,1))
for (i in 1:9){
subset<-mss.plotdat[which(mss.plotdat$tumor_type==cancers[i]),]
mypamr.plotsurvival(subset[,5], subset[,3]/365, subset[,2], cols[1:4], lty=c(2,1), lwd=c(2,2),main=cancers[i])
}
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_only_msi_survival_continuous_all_cancers_split_by_cancer_2.pdf',width=12,height=12)
par(mfrow=c(3,3),mar=c(4,4,2,1))
for (i in 10:18){
subset<-mss.plotdat[which(mss.plotdat$tumor_type==cancers[i]),]
mypamr.plotsurvival(subset[,5], subset[,3]/365, subset[,2], cols[1:4], lty=c(2,1), lwd=c(2,2),main=cancers[i])
}
dev.off()

mypamr.msi.plotsurvival <- function (group, survival.time, censoring.status, cols, lty, lwd, main) 
{
  require(survival)
  n.class <- length(unique(group))
  junk <- survfit(Surv(survival.time, censoring.status) ~ as.factor(group))
  junk2 <- coxph(Surv(survival.time, censoring.status) ~ as.factor(group))
  pv <- 1 - pchisq(2 * (junk2$loglik[2] - junk2$loglik[1]), 
                   df = n.class - 1)
  plot(junk, col = cols, xlab = "Years", ylab = "Probability of survival", lty=lty, lwd=lwd,main=main)
  legend(0.01 * max(survival.time), 0.2, col = cols, lty = lty, lwd=lwd, legend = as.character(c('MSI-H','MSS')), bty='n')
  text(0.1 * max(survival.time), 0.25, paste("p =", as.character(round(pv, 4))), cex=0.8)
  return()
}

mypamr.msi.plotsurvival <- function (group, survival.time, censoring.status, cols, lty, lwd, main, legend) 
{
  require(survival)
  n.class <- length(unique(group))
  junk <- survfit(Surv(survival.time, censoring.status) ~ as.factor(group))
  junk2 <- coxph(Surv(survival.time, censoring.status) ~ as.factor(group))
  pv <- 1 - pchisq(2 * (junk2$loglik[2] - junk2$loglik[1]), 
                   df = n.class - 1)
  plot(junk, col = cols, xlab = "Years", ylab = "Probability of survival", lty=lty, lwd=lwd,main=main)
  legend(0.01 * max(survival.time), 0.2, col = cols, lty = lty, lwd=lwd, legend = legend, bty='n')
  text(0.1 * max(survival.time), 0.25, paste("p =", as.character(round(pv, 4))), cex=0.8)
  return()
}

#COADREAD
coadread<-plotdat[union(which(plotdat$tumor_type=='COAD'),which(plotdat$tumor_type=='READ')),]
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/coadread_msih_vs_mss.pdf',width=7*2,height=7)
mypamr.msi.plotsurvival(coadread[,1], coadread[,3]/365, coadread[,2], cols[1:2], lty=c(2,1), lwd=c(2,2),main='')
dev.off()

coadread<-plotdat[union(which(plotdat$tumor_type=='COAD'),which(plotdat$tumor_type=='READ')),]
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/coadread_msih_vs_mss_continuous.pdf',width=7*2,height=7)
mypamr.msi.plotsurvival(coadread[,5], coadread[,3]/365, coadread[,2], cols[1:4], lty=c(2,1), lwd=c(2,2),main='',legend=c('first','second','third','fourth'))
dev.off()

coadread2<-plotdat2[union(which(plotdat2$tumor_type=='COAD'),which(plotdat2$tumor_type=='READ')),]
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/coadread_mss_continuous.pdf',width=7*2,height=7)
mypamr.msi.plotsurvival(coadread2[,5], coadread2[,3]/365, coadread2[,2], cols[1:4], lty=c(2,1), lwd=c(2,2),main='',legend=c('first','second','third','fourth'))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/survival_trees_split_by_cancer.pdf',width=12,height=12)
par(mfrow=c(3,3),mar=c(4,4,2,1))
for (i in 2:9){
subset<-clindat[which(clindat$tumor_type==cancers[i]),]
plot(ctree(Surv(os,vitalstatus)~num_unstable+as.numeric(as.matrix(yearstobirth))+gender+radiationtherapy+pathologicstage,data=subset))
}
dev.off()

full_msi_surv_tree <- ctree(formula=Surv(clindat$os,clindat$vitalstatus)~factor(clindat$msi_status)+as.numeric(as.matrix(clindat$yearstobirth))+clindat$gender+clindat$tumor_type+clindat$radiationtherapy+clindat$pathologicstage)

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/survival_trees_split_by_cancer_2.pdf',width=12,height=12)
par(mfrow=c(3,3),mar=c(4,4,2,1))
for (i in 10:18){
subset<-mss.plotdat[which(mss.plotdat$tumor_type==cancers[i]),]
mypamr.plotsurvival(subset[,5], subset[,3]/365, subset[,2], cols[1:4], lty=c(2,1), lwd=c(2,2),main=cancers[i])
}
dev.off()