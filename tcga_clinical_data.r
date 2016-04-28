require(cgdsr)
cbiop <- CGDS("http://www.cbioportal.org/public-portal/")
studies<-getCancerStudies(cbiop)[,1]
genetics<-getGeneticProfiles(cbiop,'ucec_tcga')[,c(1:2)]
mutations<-getGeneticProfiles(cbiop,'ucec_tcga')[8,1]
cases<-getCaseLists(cbiop,"ucec_tcga")[,c(1:3)]
cases2<-getCaseLists(cbiop,"ucec_tcga_pub")[,c(1:3)]
ucec_clin_data1<-getClinicalData(cbiop,cases[2,1])
ucec_clin_data2<-getClinicalData(cbiop,cases2[1,1])
b<-merge(ucec_clin_data1,ucec_clin_data2,by='row.names')
drops<-c('DFS_MONTHS','DFS_STATUS','AGE','OS_MONTHS','OS_STATUS')
ucec_clin_data<-merge(ucec_clin_data1[,!(names(ucec_clin_data1) %in% drops)],ucec_clin_data2,by='row.names')

#ucec_samples<-gsub('\\.','-',rownames(ucec_clin_data))
#ucec_samples<-substr(ucec_samples,start=1,stop=12)
#match(sampleinfo1$SAMPLE_NAME[which(sampleinfo1$TUMOR_TYPE=='UCEC')],ucec_samples)
#ucec_mut_data<-getProfileData(cbiop,'MLH1',genetics[,1],cases[2,1])

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')

#questions to address
#1 differences in pt outcomes between MSI-H and MSS cancers
#2 associations between rnaseq and unstable sites
#3 associations between # of unstable sites per sample and another quantitative trait

#COAD/READ new
#install.packages('TCGA2STAT',repo = 'http://cran.fhcrc.org/')
mypamr.plotsurvival <- function (group, survival.time, censoring.status, cols, lty, lwd) 
{
  require(survival)
  n.class <- length(unique(group))
  junk <- survfit(Surv(survival.time, censoring.status) ~ as.factor(group))
  junk2 <- coxph(Surv(survival.time, censoring.status) ~ as.factor(group))
  pv <- 1 - pchisq(2 * (junk2$loglik[2] - junk2$loglik[1]), 
                   df = n.class - 1)
  plot(junk, col = cols, xlab = "Years", ylab = "Probability of survival", lty=lty, lwd=lwd)
  legend(0.01 * max(survival.time), 0.2, col = cols, lty = lty, lwd=lwd, legend = as.character(1:n.class), box.lwd=1)
  text(0.1 * max(survival.time), 0.25, paste("p =", as.character(round(pv, 4))), cex=0.8)
  return()
}

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')

sampleinfo$num_unstable<-sampledat$num_unstable_raw

library('TCGA2STAT')
coad.rnaseq<-getTCGA(disease='COAD',data.type='RNASeq2',type='RPKM',clinical=TRUE)
coad.mut<-getTCGA(disease='COAD',data.type='Mutation',type='RPKM',clinical=TRUE)
coad.rnaseq.tum.norm<-TumorNormalMatch(coad.rnaseq$dat)

coad.os<-coad.rnaseq$merged.dat[,1:3]
coad.os$num_unstable<-sampledat$num_unstable_raw[match(coad.os$bcr,sampledat$sample_name)]
coad.os$msi_status<-sampledat$msi_status_tree_scaled[match(coad.os$bcr,sampledat$sample_name)]
coad.os2<-coad.os[complete.cases(coad.os),]
plotdat <- data.frame(msi=coad.os2$msi_status,status=coad.os2$status,os=coad.os2$OS,quartiles=factor(cut(coad.os2$num_unstable, quantile(coad.os2$num_unstable), include.lowest = TRUE),labels = LETTERS[1:4]))
plotdat <- plotdat[order(plotdat[,1]),]
coad.plotdat <- plotdat

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/coad_msi_survival.pdf',width=7*2,height=7)
mypamr.plotsurvival(plotdat[,1], plotdat[,3]/365, plotdat[,2], c(colors()[c(153,220)]), lty=c(2,1), lwd=c(2,2))
dev.off()

require(RColorBrewer)
cols<-brewer.pal(9, "Set1")

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/coad_msi_survival_continuous.pdf',width=7*2,height=7)
mypamr.plotsurvival(plotdat[,4], plotdat[,3]/365, plotdat[,2], cols[1:4], lty=c(2,1), lwd=c(2,2))
dev.off()

#read
read.rnaseq<-getTCGA(disease='READ',data.type='RNASeq2',type='RPKM',clinical=TRUE)

read.os<-read.rnaseq$merged.dat[,1:3]
read.os$num_unstable<-sampledat$num_unstable_raw[match(read.os$bcr,sampledat$sample_name)]
read.os$msi_status<-sampledat$msi_status_tree_scaled[match(read.os$bcr,sampledat$sample_name)]
read.os2<-read.os[complete.cases(read.os),]
plotdat <- data.frame(msi=read.os2$msi_status,status=read.os2$status,os=read.os2$OS,quartiles=factor(cut(read.os2$num_unstable, quantile(read.os2$num_unstable), include.lowest = TRUE),labels = LETTERS[1:4]))
plotdat <- plotdat[order(plotdat[,1]),]
read.plotdat <- plotdat

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/read_msi_survival.pdf',width=7*2,height=7)
mypamr.plotsurvival(plotdat[,1], plotdat[,3]/365, plotdat[,2], c(colors()[c(153,220)]), lty=c(2,1), lwd=c(2,2))
dev.off()

require(RColorBrewer)
cols<-brewer.pal(9, "Set1")

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/read_msi_survival_continuous.pdf',width=7*2,height=7)
mypamr.plotsurvival(plotdat[,4], plotdat[,3]/365, plotdat[,2], cols[1:4], lty=c(2,1), lwd=c(2,2))
dev.off()

#UCEC
ucec.rnaseq<-getTCGA(disease='UCEC',data.type='RNASeq2',type='RPKM',clinical=TRUE)
ucec.rnaseq.old<-getTCGA(disease='UCEC',data.type='RNASeq',type='RPKM',clinical=TRUE)

ucec.os<-ucec.rnaseq$merged.dat[,1:3]
ucec.os$num_unstable<-sampledat$num_unstable_raw[match(ucec.os$bcr,sampledat$sample_name)]
ucec.os$msi_status<-sampledat$msi_status_tree_scaled[match(ucec.os$bcr,sampledat$sample_name)]
ucec.os2<-ucec.os[complete.cases(ucec.os),]
plotdat <- data.frame(msi=ucec.os2$msi_status,status=ucec.os2$status,os=ucec.os2$OS,quartiles=factor(cut(ucec.os2$num_unstable, quantile(ucec.os2$num_unstable), include.lowest = TRUE),labels = LETTERS[1:4]))
plotdat <- plotdat[order(plotdat[,1]),]
ucec.plotdat <- plotdat

#STAD
stad.rnaseq<-getTCGA(disease='STAD',data.type='RNASeq2',type='RPKM',clinical=TRUE)

stad.os<-stad.rnaseq$merged.dat[,1:3]
stad.os$num_unstable<-sampledat$num_unstable_raw[match(stad.os$bcr,sampledat$sample_name)]
stad.os$msi_status<-sampledat$msi_status_tree_scaled[match(stad.os$bcr,sampledat$sample_name)]
stad.os2<-stad.os[complete.cases(stad.os),]
plotdat <- data.frame(msi=stad.os2$msi_status,status=stad.os2$status,os=stad.os2$OS,quartiles=factor(cut(stad.os2$num_unstable, quantile(stad.os2$num_unstable), include.lowest = TRUE),labels = LETTERS[1:4]))
plotdat <- plotdat[order(plotdat[,1]),]
stad.plotdat <- plotdat

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/stad_msi_survival.pdf',width=7*2,height=7)
mypamr.plotsurvival(plotdat[,1], plotdat[,3]/365, plotdat[,2], c(colors()[c(153,220)]), lty=c(2,1), lwd=c(2,2))
dev.off()

require(RColorBrewer)
cols<-brewer.pal(9, "Set1")

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/stad_msi_survival_continuous.pdf',width=7*2,height=7)
mypamr.plotsurvival(plotdat[,4], plotdat[,3]/365, plotdat[,2], cols[1:4], lty=c(2,1), lwd=c(2,2))
dev.off()

#META

meta.os2<-rbind(ucec.os2,coad.os2,read.os2,stad.os2)
plotdat <- data.frame(msi=meta.os2$msi_status,status=meta.os2$status,os=meta.os2$OS,quartiles=factor(cut(meta.os2$num_unstable, quantile(meta.os2$num_unstable), include.lowest = TRUE),labels = LETTERS[1:4]))
plotdat <- plotdat[order(plotdat[,1]),]

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/full_msi_survival.pdf',width=7*2,height=7)
mypamr.plotsurvival(plotdat[,1], plotdat[,3]/365, plotdat[,2], c(colors()[c(153,220)]), lty=c(2,1), lwd=c(2,2))
dev.off()

require(RColorBrewer)
cols<-brewer.pal(9, "Set1")

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/full_msi_survival_continuous.pdf',width=7*2,height=7)
mypamr.plotsurvival(plotdat[,4], plotdat[,3]/365, plotdat[,2], cols[1:4], lty=c(2,1), lwd=c(2,2))
dev.off()

mss.plotdat<-plotdat[which(plotdat$msi=='MSS'),]

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/mss_only_msi_survival_continuous.pdf',width=7*2,height=7)
mypamr.plotsurvival(mss.plotdat[,4], mss.plotdat[,3]/365, mss.plotdat[,2], cols[1:4], lty=c(2,1), lwd=c(2,2))
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/full_continuous_survival.pdf',width=7*2,height=7)
plot(survfit(junk2))
dev.off()

#all cancers

#COAD/READ
require(cgdsr)
cbiop <- CGDS("http://www.cbioportal.org/public-portal/")
studies<-getCancerStudies(cbiop)[,1]
cases<-getCaseLists(cbiop,"stad_tcga")[,c(1:3)]
cases2<-getCaseLists(cbiop,"coadread_tcga_pub")[,c(1:3)]
coadread_clin_data1<-getClinicalData(cbiop,cases[2,1])
coadread_clin_data2<-getClinicalData(cbiop,cases2[1,1])
a<-merge(coadread_clin_data1,coadread_clin_data2,by='row.names')
drops<-c('GENDER','PRIMARY_SITE','OS_MONTHS','OS_STATUS')
coadread_clin_data<-merge(coadread_clin_data1,coadread_clin_data2[,!(names(coadread_clin_data2) %in% drops)],by='row.names')
coadread_clin_data$SAMPLE_ID<-substr(gsub('\\.','-',coadread_clin_data$Row.names),start=1,stop=12) #57

#source("http://bioconductor.org/biocLite.R")
#biocLite("RTCGA.clinical")
library(devtools)
library(dplyr)
library(RTCGA.clinical)
#library(devtools);biocLite("mi2-warsaw/RTCGA.tools") 
#library(RTCGA.tools)
library(survival)
library(survMisc)

LUAD.clinical %>%
   mutate(
      patient.vital_status = ifelse(LUAD.clinical$patient.vital_status %>% as.character() =="dead",1,0),
      barcode = patient.bcr_patient_barcode %>% as.character(),
      times = ifelse( !is.na(patient.days_to_last_followup),
                 patient.days_to_last_followup %>% as.character() %>% as.numeric(),
                 patient.days_to_death %>% as.character() %>% as.numeric() ),
      stage = patient.stage_event.pathologic_stage
   ) %>%
   rename(
      therapy = patient.drugs.drug.therapy_types.therapy_type
   ) %>%
   filter( !is.na(times) ) -> LUAD.clinical.selected 

LUAD.clinical.selected$stage<-revalue(LUAD.clinical.selected$stage,c("stage i"="I","stage ia"="I","stage ib"="I","stage ii"="II","stage iia"="II","stage iib"="II","stage iiia"="III","stage iiib"="III","stage iv"="IV"))

LUAD.clinical.selected %>%
   survfit( Surv(times, patient.vital_status) ~ stage, data = .) %>%
   survMisc:::autoplot.survfit( titleSize=12, type="CI" ) %>%
   .[[2]] -> km_plot_luad

luadres<-LUAD.clinical.selected %>%
   survfit( Surv(times, patient.vital_status) ~ stage, data = .)

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/luad_surv_plot_by_stage.pdf',height=7,width=7)
km_plot_luad
dev.off()

