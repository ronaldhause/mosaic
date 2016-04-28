vsize <- dim(hm2)[1]/10+3
extra.vsize <- 0

genematrix<-hm2
rownames(genematrix)<-hm$gene
hcr <- hclust(dist(hm2))
ddr <- as.dendrogram(hcr)
row.ord <- order.dendrogram(ddr)
hcc <- hclust(dist(t(hm2)))
ddc <- as.dendrogram(hcc)
col.ord <- order.dendrogram(ddc)
hm2 <- (hm2)[row.ord,col.ord]

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_advanced_heatmap_2.pdf', 8.5, vsize+1+extra.vsize)
#plot.coMut()
layout(matrix(c(6, 4, 7, 2, 1, 3, 0, 9, 0, 8, 5, 0), ncol=3, byrow=T), widths=c(1,2,1), heights=c(1, vsize-2, extra.vsize, 2))
par(mar=c(0,0,0,0), las=1)
image(1:ncol(genematrix), 1:nrow(genematrix),t(hm2), col=bluered(256), xlab="", ylab="", axes=FALSE)
segments(.5 + 1:(ncol(genematrix)-1), .5, .5 + 1:(ncol(genematrix)-1), .5 + nrow(genematrix), col="grey96", lwd=ifelse(ncol(genematrix)>200, .2, .5))
segments(.5, .5 + 1:(nrow(genematrix)-1), .5 + ncol(genematrix), .5 + 1:(nrow(genematrix)-1), col="grey96", lwd=ifelse(ncol(genematrix)>200, .2, .5))

plot(0,0)

## q-values
par(mar=c(0,6,0,.5))
plot(0, xlim=c(min(.4, min(-log10(hm$qvals+0.0001))), max(-log10(hm$qvals+0.0001))), type="n", yaxt="n", frame.plot=FALSE, xlab="", ylab="")
par(usr=c(par("usr")[1:2], 0, nrow(hm)), lwd=.8)
mypos <- barplot(rev(-log10(hm$qvals+0.001)), horiz=TRUE, axes=FALSE, add=TRUE, names.arg=rep("", nrow(hm)), border="grey94", space=0, cex.names=.9, col=c("grey70", "grey55")[factor((rev(hm$qvals)<=.1)+1, levels=c(1,2))])
axis(2, at=mypos, labels=rev(hm$gene), lwd=0, , cex.axis=1.1, line=-.2)
mtext("-log10(q-value)", side=1, cex=.7, line=2.2, adj=1.1)
abline(v=-log10(.1), col="red", lwd=1)
abline(v=-log10(.25), col="purple", lty=2, lwd=1)

dev.off()

pdf('/net/shendure/vol3/data/nobackup/hauser/Shendure/msi/cancer_specific_advanced_heatmap_trial_2.pdf', height=7*3,width=7*2)
heatmap.2(genematrix,col=bluered(256),trace='none',scale='none',cexCol=1,cexRow=0.1,keysize=0.75,density.info='none')
dev.off()
