source('/net/shendure/vol1/home/hauser/scripts/shendure/miscellaneous/useful_functions.r')

load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_loci_results_050315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampleinfo_042315.robj')

#MSI-H
msihtumor<-read.delim('/net/shendure/vol7/hauser/msi_tcga/individual_patient_calls/UCEC_4-7/UCECnew_3/TCGA-D1-A17U.tumor.msi.txt',header=T)
msihnormal<-read.delim('/net/shendure/vol7/hauser/msi_tcga/individual_patient_calls/UCEC_4-7/UCECnew_3/TCGA-D1-A17U.normal.msi.txt',header=T)
msihresult<-read.delim('/net/shendure/vol7/hauser/msi_tcga/individual_patient_calls/UCEC_4-7/UCECnew_3/TCGA-D1-A17U.MSI_full_analysis.txt',header=F)

cor.test(as.numeric(as.matrix(msihtumor$Avg_read_depth)),as.numeric(as.matrix((msihtumor$Number_Peaks))),method='s') #rho = -0.03
cor.test(as.numeric(as.matrix(msihnormal$Avg_read_depth)),as.numeric(as.matrix((msihnormal$Number_Peaks))),method='s') #rho = -0.02
header<-read.delim('/net/shendure/vol7/hauser/msi_tcga/TCGA_Analysis_header.txt',header=F)

msihresult$peak_diff<-as.numeric(as.matrix(msihresult$V20))
msihresult$msihtumor_read_depth<-as.numeric(as.matrix(msihtumor$Avg_read_depth[match(msihresult$V5,msihtumor$Position)]))
msihresult$msihnormal_read_depth<-as.numeric(as.matrix(msihnormal$Avg_read_depth[match(msihresult$V5,msihnormal$Position)]))
msihresult$joint_read_depth<-msihresult$msihtumor_read_depth+msihresult$msihnormal_read_depth
msihresult$diff<-msihresult$msihtumor_read_depth-msihresult$msihnormal_read_depth

cor.test(msihresult$joint_read_depth,msihresult$peak_diff,method='s') #rho = -0.06
cor.test(msihresult$msihtumor_read_depth,msihresult$peak_diff,method='s') #rho = -0.06
cor.test(msihresult$msihnormal_read_depth,msihresult$peak_diff,method='s') #rho = -0.06
cor.test(msihresult$diff,msihresult$peak_diff,method='s') #rho = -0.04

plotdat<-data.frame(summed_read_depth=msihresult$joint_read_depth,tumor_read_depth=msihresult$msihtumor_read_depth,normal_read_depth=msihresult$msihnormal_read_depth,tumor_minus_normal_read_depth=msihresult$diff,peak_diff=msihresult$peak_diff)

#bottom, left, top, right
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/read_depth_vs_instability.pdf',width=5.5,height=5.5,bg="transparent")
par(mfrow=c(2,2),mar=c(4,4,1,1))
smoothScatter(plotdat$peak_diff,plotdat$normal_read_depth,xlab='peak difference',ylab='normal read depth',lwd=2,bty='n',las=2,cex=1.6)
smoothScatter(plotdat$peak_diff,plotdat$tumor_read_depth,xlab='peak difference',ylab='tumor read depth',lwd=2,bty='n',las=2,cex=1.6)
smoothScatter(plotdat$peak_diff,plotdat$tumor_minus_normal_read_depth,xlab='peak difference',ylab='tumor + normal read depth',lwd=2,bty='n',las=2,cex=1.6)
smoothScatter(plotdat$peak_diff,plotdat$tumor_minus_normal_read_depth,xlab='peak difference',ylab='tumor - normal read depth',lwd=2,bty='n',las=2,cex=1.6)
dev.off()

ggplot(df) + stat_density2d(aes(x = x, y = y, fill = ..level..), geom = "polygon") + facet_grid(.~group) + theme_bw() + labs(x = expression(B73~~deviation~~(mu)))

ggplot(mydata) + aes(x=x, y=y) + scale_x_log10() + scale_y_log10() + 
  stat_density2d(geom="tile", aes(fill=..density..^0.25, alpha=1), contour=FALSE) + 
  geom_point(size=0.5) +
  stat_density2d(geom="tile", aes(fill=..density..^0.25,     alpha=ifelse(..density..^0.25<0.4,0,1)), contour=FALSE) + 
  scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256))


a<-msihresult$peak_diff
a[which(msihresult$peak_diff>0)]<-'unstable'
a[which(msihresult$peak_diff<=0)]<-'stable'
msihresult$stable<-a

x <- rnorm(20000)
y <- x + rnorm(20000, 0.05)
df <- data.frame(x = x, y = y,
  d = densCols(x, y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
p <- ggplot(df) +
    geom_point(aes(x, y, col = d), size = 1) +
    scale_color_identity() +
    theme_bw()
print(p)

require(reshape2)
plotdat2<-melt(a)
colnames(plotdat2)<-c('msi_status','type','prop')

library(RColorBrewer)
cols<-brewer.pal(9, "Set1")
pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/repeat_type_vs_instability.pdf',width=7*1.5,height=7)
ggplot(data=plotdat2, aes(x=type,y=prop*100,fill=msi_status)) + geom_bar(stat="identity",position="dodge") + scale_fill_manual(values=cols[1:2]) + theme_bw(base_size = 20) + ylab('Proportion of unstable sites (%)') + coord_cartesian(ylim = c(0,12)) + xlab('Repeat tract length')
dev.off()

#extracting read depth for each individual
load('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/sampledata_042315.robj')
require(XML)
xml = xmlTreeParse("/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/everything.xml", useInternalNodes=TRUE)
filename<-xpathApply(xml,"////filename",xmlValue)
filesize<-xpathApply(xml,"////filesize",xmlValue)
dat<-data.frame(filename=unlist(filename),filesize=as.numeric(as.matrix(unlist(filesize))))
dat<-dat[-grep('bai',dat$filename),]
samplenames<-sampledat$sample_name
a<-c()
for (i in 1:dim(dat)[1]){
	print(i*100/dim(dat)[1])
	loc<-regexpr('TCGA-',dat$filename[i])
	if (loc>0){
		a[i]<-substr(dat$filename[i],loc,loc+14)
	}
	else{
		a[i]<-substr(dat$filename[i],1,15)
	}
}
dat$sample_name2<-a
dat$sample_name<-substr(dat$sample_name2,1,12)
dat$num_unstable<-sampledat$num_unstable_raw[match(dat$sample_name,samplenames)]
b<-substr(dat$sample_name2,14,16)
dat$sample_type<-rep(NA,dim(dat)[1])
dat$sample_type[which(b=='01')]<-'tumor'
dat$sample_type[which(b=='10')]<-'normal'
require(dplyr)
dat %>% group_by(sample_name) %>% summarize(diff=filesize[match('tumor',sample_type)]-filesize[match('normal',sample_type)]) -> dat2
dat$msi_status<-sampledat$msi_status_tree_scaled[match(dat$sample_name,samplenames)]
dat$tumor_type<-sampledat$tumor_type[match(dat$sample_name,samplenames)]
dat$diff<-dat2$diff[match(dat$sample_name,dat2$sample_name)]
tumors<-dat[which(dat$sample_type=='tumor'),]
normals<-dat[which(dat$sample_type=='normal'),]
normals2<-normals[which(is.na(normals$msi_status)==F),]
normals2<-normals2[grep('Illumina',normals2$filename),]
normals3<-normals2[grep('hg19',normals2$filename),]

tumors2<-tumors[which(is.na(tumors$msi_status)==F),]
tumors2<-tumors2[grep('Illumina',tumors2$filename),]
tumors3<-tumors2[grep('hg19',tumors2$filename),]

#plotting estimated read depth vs. number of unstable satellites identified
filenames <- list.files('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/tcga_readdepth', pattern="*.stats", full.names=TRUE)
ldf <- lapply(filenames, function(x)read.delim(x,header=F))
mappedreads <- unlist(lapply(ldf, function(x)sum(x$V3)))
unmappedreads <- unlist(lapply(ldf, function(x)sum(x$V4)))
names(mappedreads) <- substr(filenames,69,80)
names(unmappedreads) <- substr(filenames,69,80)
totalreads<-mappedreads+unmappedreads
filesizes<-tumors3$filesize[match(names(mappedreads),tumors3$sample_name)]
a<-lm(totalreads~filesizes)
new.df<-data.frame(filesizes=tumors3$filesize)
tumors3$estimatedreads<-predict(a,new.df)
new.df<-data.frame(filesizes=normals3$filesize)
normals3$estimatedreads<-predict(a,new.df)

merged<-merge(normals3,tumors3,by='sample_name') #y = tumor, x = normal
merged$diff<-merged$estimatedreads.x-merged$estimatedreads.y

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/global_read_depth_vs_instability.pdf',width=5.5,height=5.5,bg="transparent")
par(mar=c(4,4,1,1))
smoothScatter(merged$num_unstable.x,merged$diff,xlab='# of unstable microsatellites',ylab='tumor - normal read depth',lwd=2,bty='n',las=2,cex=1.6,mgp=c(3,0.5,0))
dev.off()

one<-tumors2$sample_name[which(abs(tumors2$filesize-quantile(tumors2$filesize)[2])==min(abs(tumors2$filesize-quantile(tumors2$filesize)[2])))]
two<-tumors2$sample_name[which(abs(tumors2$filesize-quantile(tumors2$filesize)[3])==min(abs(tumors2$filesize-quantile(tumors2$filesize)[3])))]
three<-tumors2$sample_name[which(abs(tumors2$filesize-quantile(tumors2$filesize)[4])==min(abs(tumors2$filesize-quantile(tumors2$filesize)[4])))]
samps<-c(one,two[1],three)
a<-tumors2[intersect(which(tumors2$filesize>11950000000),which(tumors2$filesize<12000000000)),][2,]
b<-tumors2[intersect(which(tumors2$filesize>17325000000),which(tumors2$filesize<17375000000)),]
tumors2[match(samps,tumors2$sample_name),]

