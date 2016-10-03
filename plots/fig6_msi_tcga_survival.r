library('TCGA2STAT')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_post_review_070816.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_post_review_070816.robj')

cancers<-unique(sampledat$tumor_type)

clindatlist<-list()
rnaseqlist<-list()
for (i in 1:length(cancers)){
	subset.rnaseq<-getTCGA(disease=cancers[i],data.type='RNASeq2',type='RPKM',clinical=TRUE)
	rnaseqlist[[i]]<-subset.rnaseq
	clindatlist[[i]]<-subset.rnaseq$clinical
}

for (i in 1:length(cancers)){
	clindatlist[[i]]<-cbind(clindatlist[[i]],tumor_type=rep(cancers[i],dim(clindatlist[[i]])[1]))
	clindatlist[[i]]<-cbind(clindatlist[[i]],sample_name=rownames(clindatlist[[i]]))
}

save(rnaseqlist,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_rnaseq_data.robj')
save(clindatlist,file='/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_clinical_data.robj')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_rnaseq_data.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_clinical_data.robj')

clindat<-Reduce(function(...) merge(..., all=T), clindatlist)
clindat$tumor_type<-revalue(clindat$tumor_type,c('OVmu'='OV'))
clindat$tumor_type<-as.character(clindat$tumor_type)
clindat$num_unstable<-sampledat$num_unstable[match(clindat$sample_name,sampledat$sample_name)]
clindat$msi_status<-sampledat$msi[match(clindat$sample_name,sampledat$sample_name)]
clindat$yearstobirth<-as.numeric(as.matrix(clindat$yearstobirth))
clindat$peak_avg<-sampledat$peak_avg[match(clindat$sample_name,sampledat$sample_name)]
clindat<-clindat[-which(is.na(clindat$num_unstable)==T),]
#clinical data on 5921 samples

#correlation between age and prop_unstable rho = 0.02, p = 0.07; sig difference between gender (p = 0.018, 0.48% in female, 0.46% in male); p = 0.0005 female = 0.0005, male = 0.0002 for peak avg
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

ucec_msi_status<-read.delim('/net/shendure/vol7/stevesal/MSI_TCGA_project/ucec_tcga_pub_clinical_data.tsv', header=T, stringsAsFactors=F)
ucec_msi_status$Tumor.Stage.2009<-tolower(ucec_msi_status$Tumor.Stage.2009)
ucec_msi_status$os<-clindat$os[match(ucec_msi_status$Patient.ID, clindat$sample_name)]
require(plyr)
clindat$vitalstatus<-as.numeric(as.matrix(clindat$vitalstatus))
clindat$os<-as.numeric(as.matrix(clindat$os))
clindat$pathologicstage<-as.character(clindat$pathologicstage)
clindat$pathologicstage[which(is.na(match(clindat$sample_name, ucec_msi_status$Patient.ID))==F)]<-ucec_msi_status$Tumor.Stage.2009[na.omit(match(clindat$sample_name, ucec_msi_status$Patient.ID))]
clindat$pathologicstage<-revalue(clindat$pathologicstage,c('stage iia'='stage ii','stage iib'='stage ii','stage iic'='stage ii','stage iiia'='stage iii','stage iiib'='stage iii','stage iiic'='stage iii','stage iva'='stage iv','stage ia'='stage i','stage ivb'='stage iv','stage ib'='stage i','stage ivc'='stage iv','stage x' = NA,'i/ii nos'=NA,'stage 0'=NA,'unknown'=NA))

#eliminating MSI-H samples in the major four cancer types
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

require(party)
require(survival)

#just MSI cancers
msicancers<-subset(clindat,tumor_type %in% c('COAD','READ','UCEC','STAD'))
#just MSS cancers
msicancers2<-subset(clindat2,tumor_type %in% c('COAD','READ','UCEC','STAD'))

#Figure 6A
a<-coxph(formula=Surv(msicancers$os/365,msicancers$vitalstatus)~msicancers$msi_status+msicancers$yearstobirth+msicancers$gender+msicancers$tumor_type+msicancers$radiationtherapy+msicancers$pathologicstage) #P = 0.234

sub<-msicancers[which(msicancers$msi_status=='MSI-H'),]
sub$msi_status<-factor(sub$msi_status)
msih<-survfit(coxph(formula=Surv(sub$os/365,sub$vitalstatus)~sub$yearstobirth+sub$gender+sub$tumor_type+sub$radiationtherapy+sub$pathologicstage))
sub<-msicancers[which(msicancers$msi_status=='MSS'),]
sub$msi_status<-factor(sub$msi_status)
mss<-survfit(coxph(formula=Surv(sub$os/365,sub$vitalstatus)~sub$yearstobirth+sub$gender+sub$tumor_type+sub$radiationtherapy+sub$pathologicstage))

cols<-brewer.pal(9, "Set1")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/fig6a_msi_cancers_mss_vs_msi_071316_adjusted.pdf',width=7*1.5,height=7)
plot(mss, col=cols[2], lty=1, lwd=c(2,2), xlab="Years", ylab="Probability of survival",bty='n', mark.time=T, conf.int=F, xlim=c(0,20))
lines(msih, col=cols[1], lty=1, lwd=2, mark.time=T, conf.int=F)
legend(0.01 * max(msicancers$os/365), 0.2, col = cols[1:2], lty = 1, lwd=2, legend = c('MSI-H','MSS'), bty='n')
text(0.1 * max(msicancers$os/365), 0.25, paste("p =", as.character(round(summary(a)$coef[1,5], 4))), cex=0.8)
dev.off()

#Figure 6B
b<-coxph(formula=Surv(msicancers$os,msicancers$vitalstatus)~msicancers$num_unstable+msicancers$yearstobirth+msicancers$gender+msicancers$tumor_type+msicancers$radiationtherapy+msicancers$pathologicstage) #P = 0.016

sub<-msicancers[which(msicancers$s_quartile=='A'),]
sub$s_quartile<-factor(sub$s_quartile)
first<-survfit(coxph(formula=Surv(sub$os/365,sub$vitalstatus)~sub$yearstobirth+sub$gender+sub$tumor_type+sub$radiationtherapy+sub$pathologicstage))
sub<-msicancers[which(msicancers$s_quartile=='B'),]
sub$s_quartile<-factor(sub$s_quartile)
second<-survfit(coxph(formula=Surv(sub$os/365,sub$vitalstatus)~sub$yearstobirth+sub$gender+sub$tumor_type+sub$radiationtherapy+sub$pathologicstage))
sub<-msicancers[which(msicancers$s_quartile=='C'),]
sub$s_quartile<-factor(sub$s_quartile)
third<-survfit(coxph(formula=Surv(sub$os/365,sub$vitalstatus)~sub$yearstobirth+sub$gender+sub$tumor_type+sub$radiationtherapy+sub$pathologicstage))
sub<-msicancers[which(msicancers$s_quartile=='D'),]
sub$s_quartile<-factor(sub$s_quartile)
fourth<-survfit(coxph(formula=Surv(sub$os/365,sub$vitalstatus)~sub$yearstobirth+sub$gender+sub$tumor_type+sub$radiationtherapy+sub$pathologicstage))

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/fig6b_msi_cancers_num_unstable_071316_adjusted.pdf',width=7*1.5,height=7)
plot(first, col=cols[1], lty=1, lwd=2, xlab="Years", ylab="Probability of survival", bty='n', mark.time=T, conf.int=F, xlim=c(0,20))
lines(second, col=cols[2], lty=1, lwd=2, mark.time=T, conf.int=F)
lines(third, col=cols[3], lty=1, lwd=2, mark.time=T, conf.int=F)
lines(fourth, col=cols[4], lty=1, lwd=2, mark.time=T, conf.in=F)
legend(0.01 * max(msicancers$os/365), 0.2, col=cols[1:4], lty = 1, lwd=2, legend=c('1st','2nd','3rd','4th'), bty='n')
text(0.1 * max(msicancers$os/365), 0.25, paste("p =", as.character(round(summary(b)$coef[1,5], 4))), cex=0.8)
dev.off()

#Figure 6C
msicancers<-msicancers2

c<-coxph(formula=Surv(msicancers$os,msicancers$vitalstatus)~msicancers$num_unstable+msicancers$yearstobirth+msicancers$gender+msicancers$tumor_type+msicancers$radiationtherapy+msicancers$pathologicstage) #0.0039

sub<-msicancers[which(msicancers$s_quartile=='A'),]
sub$s_quartile<-factor(sub$s_quartile)
first<-survfit(coxph(formula=Surv(sub$os/365,sub$vitalstatus)~sub$yearstobirth+sub$gender+sub$tumor_type+sub$radiationtherapy+sub$pathologicstage))
sub<-msicancers[which(msicancers$s_quartile=='B'),]
sub$s_quartile<-factor(sub$s_quartile)
second<-survfit(coxph(formula=Surv(sub$os/365,sub$vitalstatus)~sub$yearstobirth+sub$gender+sub$tumor_type+sub$radiationtherapy+sub$pathologicstage))
sub<-msicancers[which(msicancers$s_quartile=='C'),]
sub$s_quartile<-factor(sub$s_quartile)
third<-survfit(coxph(formula=Surv(sub$os/365,sub$vitalstatus)~sub$yearstobirth+sub$gender+sub$tumor_type+sub$radiationtherapy+sub$pathologicstage))
sub<-msicancers[which(msicancers$s_quartile=='D'),]
sub$s_quartile<-factor(sub$s_quartile)
fourth<-survfit(coxph(formula=Surv(sub$os/365,sub$vitalstatus)~sub$yearstobirth+sub$gender+sub$tumor_type+sub$radiationtherapy+sub$pathologicstage))

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/fig6b_msi_cancers_num_unstable_071316_adjusted_mss_only.pdf',width=7*1.5,height=7)
plot(first, col=cols[1], lty=1, lwd=2, xlab="Years", ylab="Probability of survival", bty='n', mark.time=T, conf.int=F, xlim=c(0,20))
lines(second, col=cols[2], lty=1, lwd=2, mark.time=T, conf.int=F)
lines(third, col=cols[3], lty=1, lwd=2, mark.time=T, conf.int=F)
lines(fourth, col=cols[4], lty=1, lwd=2, mark.time=T, conf.in=F)
legend(0.01 * max(msicancers$os/365), 0.2, col=cols[1:4], lty = 1, lwd=2, legend=c('1st','2nd','3rd','4th'), bty='n')
text(0.1 * max(msicancers$os/365), 0.25, paste("p =", as.character(round(summary(c)$coef[1,5], 4))), cex=0.8)
dev.off()

#Figure 6A
a<-coxph(formula=Surv(msicancers$os/365,msicancers$vitalstatus)~msicancers$msi_status+msicancers$yearstobirth+msicancers$gender+msicancers$tumor_type+msicancers$radiationtherapy+msicancers$pathologicstage) #P = 0.234

aplot<-coxph(formula=Surv(msicancers$os/365,msicancers$vitalstatus)~strata(msicancers$msi_status)+msicancers$yearstobirth+msicancers$gender+msicancers$tumor_type+msicancers$radiationtherapy+msicancers$pathologicstage) #P = 0.12

sub<-msicancers[which(msicancers$msi_status=='MSI-H'),]
sub$msi_status<-factor(sub$msi_status)
msih<-survfit(coxph(formula=Surv(sub$os/365,sub$vitalstatus)~sub$yearstobirth+sub$gender+sub$tumor_type+sub$radiationtherapy+sub$pathologicstage))
sub<-msicancers[which(msicancers$msi_status=='MSS'),]
sub$msi_status<-factor(sub$msi_status)
mss<-survfit(coxph(formula=Surv(sub$os/365,sub$vitalstatus)~sub$yearstobirth+sub$gender+sub$tumor_type+sub$radiationtherapy+sub$pathologicstage))

#evaluating across all cancers
survival<-list()
for (i in 1:length(cancers)){
	subset<-clindat[which(clindat$tumor_type==cancers[i]),]
	subset2<-subset(subset,select=c('os','vitalstatus','num_unstable','yearstobirth','radiationtherapy','pathologicstage','gender'))
	subset2<-transform(subset2,os=as.numeric(os),vitalstatus=as.numeric(vitalstatus),yearstobirth=as.numeric(yearstobirth),radiationtherapy=as.factor(radiationtherapy),pathologicstage=as.factor(pathologicstage),gender=as.factor(gender))
	covarsel<-apply(subset2,2,function(x)length(unique(x)))
	subset3<-subset2[,which(covarsel>1)]
	survival[[i]]<-coxph(formula=Surv(os,vitalstatus)~.,data=subset3)
}

names(survival)<-cancers

#plotting cancer-specific survival outcomes
plotdat <- data.frame(msi=clindat$msi_status,status=clindat$vitalstatus,os=clindat$os,g_quartiles=clindat$g_quartile,s_quartiles=clindat$s_quartile,s_tertiles=clindat$s_tertile,tumor_type=clindat$tumor_type,half=clindat$half,num_unstable=clindat$num_unstable,age=clindat$yearstobirth,gender=clindat$gender,sample_name=clindat$sample_name)
plotdat <- plotdat[order(plotdat[,1]),]
plotdat$status <- as.numeric(as.matrix(plotdat$status))
plotdat$os <- as.numeric(as.matrix(plotdat$os))
plotdat$half<-factor(plotdat$half,levels=c("high","low"))

mypamr.plotsurvival2 <- function (group, survival.time, censoring.status, cols, lty, lwd, main) 
{
  require(survival)
  n.class <- length(unique(group))
  junk <- survfit(Surv(survival.time, censoring.status) ~ as.factor(group))
  pv <- summary(survival[[i]])$coef[1,5]
  plot(junk, col = cols, xlab = "Years", ylab = "Probability of survival", lty=lty, lwd=lwd,main=main, bty='n', mark.time=T, conf.int=F)
  legend(0.01 * max(survival.time), 0.2, col = cols[1:3], lty = lty, lwd=c(2,2,2), legend = c('1st', '2nd', '3rd'), bty='n')
  text(0.1 * max(survival.time), 0.25, paste("p =", as.character(round(pv, 4))), cex=0.8)
  return()
}

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_survival_continuous_all_cancers_split_by_cancer_tertiles_1_072116.pdf',width=12,height=12)
par(mfrow=c(3,3),mar=c(4,4,2,1))
for (i in 1:9){
subset<-plotdat[which(plotdat$tumor_type==cancers[i]),]
mypamr.plotsurvival2(subset[,6], subset[,3]/365, subset[,2], cols[1:4], lty=c(1,1,1), lwd=c(2,2,2),main=cancers[i])
}
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msi_survival_continuous_all_cancers_split_by_cancer_tertiles_2_072116.pdf',width=12,height=12)
par(mfrow=c(3,3),mar=c(4,4,2,1))
for (i in 10:18){
subset<-plotdat[which(plotdat$tumor_type==cancers[i]),]
mypamr.plotsurvival2(subset[,6], subset[,3]/365, subset[,2], cols[1:4], lty=c(1,1,1), lwd=c(2,2,2),main=cancers[i])
}
dev.off()