library(gdsfmt)
library(SNPRelate)
gdsfile = "plots/snps_selected.gds"
vcf.fn <- "vcf/FIXME.All.SNP.combined_selected.vcf.gz"

if(!file.exists(gdsfile)){
	snpgdsVCF2GDS_R(vcf.fn, gdsfile,method="biallelic.only")
	                #option=snpgdsOption(CM002236=1,CM002237=2,CM002238=3,CM002239=4,CM002240=5,CM002241=6,CM002242=7))
}

snpgdsSummary(gdsfile)
genofile <- snpgdsOpen(gdsfile)
chroms <- read.gdsn(index.gdsn(genofile,"snp.chromosome"))
chr <- strtoi(sub("SCAF_([0-9]+)","\\1",chroms,perl=TRUE))

pca <- snpgdsPCA(genofile,num.thread=2,autosome.only=FALSE)

pc.percent <- pca$varprop*100

#pca$sample.id

head(round(pc.percent, 2))
pdf("plots/PCA_snp_plots.pdf")
tab <- data.frame(sample.id = pca$sample.id,
                 # pop = pheno$MinimalMediaGrowth,
                  EV1=pca$eigenvect[,1], # PCA vector 1
                  EV2=pca$eigenvect[,2], # PCA vector 2
		  stringsAsFactors=FALSE)

plot(tab$EV2, tab$EV1,
     #, col=as.integer(tab$pop),
xlab="eigenvector 2", ylab="eigenvector 1", main="PCA SNP plot")
text(x = pca$eigenvect[,2], y = pca$eigenvect[,1], labels = tab$sample.id, pos = 1 ,cex =0.8, offset = 0.5)

set.seed(100)
# recode the snp.gds to support chromosomes?
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2,autosome.only=FALSE))
rv <- snpgdsCutTree(ibs.hc)
plot(rv$dendrogram, leaflab="none", main="A. socialis CCGP Strains")

snpgdsDrawTree(rv, type="z-score", main="A. socialis CCGP Strains")
snpgdsDrawTree(rv, main="A. socialis CCGP Strains",
               edgePar=list(col=rgb(0.5,0.5,0.5, 0.75), t.col="black"))

table(rv$samp.group)
df = data.frame(group = rv$samp.group)
rownames(df) = pca$sample.id
write.csv(df,"plots/popset_inferred.csv")
tab <- data.frame(sample.id = pca$sample.id,
                  pop = rv$samp.group,
                  EV1=pca$eigenvect[,1], # PCA vector 1
                  EV2=pca$eigenvect[,2], # PCA vector 2
                  stringsAsFactors=FALSE)
plot(tab$EV2, tab$EV1,
     col=as.integer(tab$pop),
     xlab="eigenvector 2", ylab="eigenvector 1", main="PCA SNP plot")

CORRSNP <- snpgdsPCACorr(pca, genofile, eig.which=1:4,num.thread=2)

#savepar <- par(mfrow=c(3,1), mai=c(0.3, 0.55, 0.1, 0.25))
#for (i in 1:3)
#{
#  plot(abs(CORRSNP$snpcorr[i,]), ylim=c(0,1), xlab="", ylab=paste("PC", i),
#       col=factor(chr), pch="+")
#}

