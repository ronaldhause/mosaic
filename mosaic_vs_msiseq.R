require(MSIseq)
data(test.mutationNum)
data(train.mutationNum)
data(NGStrainclass)
msiseqdata<-rbind(test.mutationNum,train.mutationNum)
msiseqdata<-train.mutationNum
msiseqdata$msiseq<-rep('MSS',dim(msiseqdata)[1])
msiseqdata$msiseq[which(msiseqdata$S.ind>0.395)]<-'MSI-H'
msiseqdata$msi_status<-sampleinfo1n$MSI_STATUS[match(rownames(msiseqdata),sampleinfo1n$SAMPLE_NAME)]
msiseqdata$mosaic<-sampledat1n$old_msi_pred[match(rownames(msiseqdata),sampledat1n$sample_name)]
msiseqdata$mosaicnew<-sampledat1n$new_msi_pred[match(rownames(msiseqdata),sampledat1n$sample_name)]
msiseqdata$mosaicneww<-sampledat1n$new_msi_predw[match(rownames(msiseqdata),sampledat1n$sample_name)]
msiseqdata$prop_unstable<-sampledat1n$prop_unstable[match(rownames(msiseqdata),sampledat1n$sample_name)]
msiseqdata$peak_avg<-sampledat1n$peak_avg[match(rownames(msiseqdata),sampledat1n$sample_name)]
msiseqdata$corrected_msi<-sampledat1n$corrected_msi[match(rownames(msiseqdata),sampledat1n$sample_name)]
msiseqdata$tumor_type<-sampledat1n$tumor_type[match(rownames(msiseqdata),sampledat1n$sample_name)]

NGStrainclass$msi_status<-sampleinfo1n$MSI_STATUS[match(NGStrainclass$Tumor_Sample_Barcode,sampleinfo1n$SAMPLE_NAME)]
ngs_msih_samples<-NGStrainclass$Tumor_Sample_Barcode[which(NGStrainclass$MSI_status=='MSI-H')]
sampleinfo1n$MSI_STATUS[na.omit(match(ngs_msih_samples,sampleinfo1n$SAMPLE_NAME))]<-'MSI-H'

#correlation of 0.8 between S.ind and prop.unstable or peak.avg
msiseqdata$corrected_msi[msiseqdata$corrected_msi=='MSI-L']<-'MSS'
msiseqdata$msi_status[msiseqdata$msi_status=='MSI-L']<-'MSS'
msiseqdata<-msiseqdata[-which(msiseqdata$msi_status=='ND'),]
msiseqdata<-msiseqdata[complete.cases(msiseqdata),]

confusionMatrix(msiseqdata$msiseq,msiseqdata$msi_status) #98.9% accuracy, 94.4% sensitivity, 100% specificity
confusionMatrix(msiseqdata$mosaic,msiseqdata$msi_status) #97.4% accuracy, 93.1% sensitivity, 98.6% specificity
confusionMatrix(msiseqdata$mosaicnew,msiseqdata$msi_status) #98.0% accuracy, 94.4% sensitivity, 98.9% specificity
confusionMatrix(msiseqdata$mosaicneww,msiseqdata$msi_status) #98.0% accuracy, 93.1% sensitivity, 99.3% specificity

#just test set
confusionMatrix(msiseqdata$msiseq,msiseqdata$msi_status) #99.1% accuracy, 96.0% sensitivity, 100% specificity
confusionMatrix(msiseqdata$mosaic,msiseqdata$msi_status) #97.4% accuracy, 92% sensitivity, 98.9% specificity
confusionMatrix(msiseqdata$mosaicnew,msiseqdata$msi_status) #99.1% accuracy, 96.0% sensitivity, 100% specificity
confusionMatrix(msiseqdata$mosaicneww,msiseqdata$msi_status) #98.2% accuracy, 92% sensitivity, 100% specificity

#just training set
confusionMatrix(msiseqdata$msiseq,msiseqdata$msi_status) #98.7% accuracy, 93.6% sensitivity, 100% specificity
confusionMatrix(msiseqdata$mosaic,msiseqdata$msi_status) #97.5% accuracy, 93.6% sensitivity, 98.4% specificity
confusionMatrix(msiseqdata$mosaicnew,msiseqdata$msi_status) #97.5% accuracy, 93.6% sensitivity, 98.4% specificity
confusionMatrix(msiseqdata$mosaicneww,msiseqdata$msi_status) #97.9% accuracy, 93.6% sensitivity, 99.0% specificity

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msiseq_sind_vs_prop_unstable.pdf',width=7*1.5,height=7*1.5)
plot(msiseqdata$S.ind,msiseqdata$prop_unstable)
abline(v=0.395,lty=2,lwd=2)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msiseq_sind_vs_peak_avg.pdf',width=7,height=7)
plot(msiseqdata$S.ind,msiseqdata$peak_avg, pch=16, col=as.factor(msiseqdata$msi_status))
legend('bottomright', bty='n', legend = levels(as.factor(msiseqdata$msi_status)), col = 1:4, cex = 1, pch = 16)
abline(v=0.395,lty=2,lwd=2)
abline(h=0.0053,lty=2,lwd=2)
dev.off()

require(DTComPair)
compar<-tab.paired(as.numeric(revalue(msiseqdata$msi_status, c("MSI-H"=1,"MSS"=0))),as.numeric(revalue(msiseqdata$msiseq, c("MSI-H"=1,"MSS"=0))),as.numeric(revalue(msiseqdata$mosaic, c("MSI-H"=1,"MSS"=0))),testnames=c('msiseq','mosaic'))
sesp.exactbinom(compar) #exact binomial test P > 0.13

cols<-

require(scatterplot3d)
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msiseq_sind_vs_peak_avg.pdf',width=7,height=7)
scatterplot3d(msiseqdata$S.ind,msiseqdata$peak_avg, msiseqdata$prop_unstable, pch=16, col=as.factor(msiseqdata$msi_status))
legend('bottomright', bty='n', legend = levels(as.factor(msiseqdata$msi_status)), col = 1:4, cex = 1, pch = 16)
abline(v=0.395,lty=2,lwd=2)
abline(h=0.0053,lty=2,lwd=2)
dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/msiseq_sind_vs_peak_avg_3d.pdf',width=7,height=7)
msiseqdata$pcolor[msiseqdata$msi_status=='MSI-H'] <- "red"
msiseqdata$pcolor[msiseqdata$msi_status=='MSS'] <- "blue"
with(msiseqdata, {
    s3d <- scatterplot3d(S.ind, peak_avg, prop_unstable,        # x y and z axis
                  color=pcolor, pch=19,        # circle color indicates no. of cylinders
                  type="h", lty.hplot=2,       # lines to the horizontal plane
                  scale.y=.75,                 # scale y axis (reduce by 25%)
                  xlab="S.ind",
                  ylab="peak_avg",
                  zlab="prop_unstable")
# add the legend
legend("topleft", inset=.05,      # location and inset
    bty="n", cex=.5,              # suppress legend box, shrink text 50%
    c("MSI-H", "MSS"), fill=c("red", "blue"))
})
dev.off()