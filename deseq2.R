######################################
## Differential Expression Analysis ##
######################################

# generate read ount table from htseq-count results
setwd("/home/hugj2006/ML/htseq_count")
files<-grep(".namesort.bam.txt",list.files(),value=TRUE)
for(f in files)
{
    x <- read.table(f,sep="\t")
    names(x) <- c("gene", gsub(".namesort.bam.txt","",f))
    if(!exists("allcount")) {allcount = x}
    else {allcount <- merge(allcount, x, by="gene")}
}
write.table(allcount, file="count.raw.txt", row.names=FALSE, sep="\t")


allcount<-read.table(file="count.raw.txt", header=TRUE, sep="\t")
row.names(allcount) <- allcount$gene
count13<-allcount[grep("Gorai.0",allcount$gene,value=TRUE),-1]
dim(count13) # 37223    47

# prepare column information table
coldata<- as.data.frame(t(allcount[grep("__",allcount$gene,value=TRUE),-1]))
coldata$sample<-rownames(coldata)
coldata$flowering <-"F"
coldata$flowering[grep("NF",coldata$sample)] <-"NF"
coldata$condition<-gsub("[.].*","",rownames(coldata))
coldata$sample<-gsub("[.]S.*|[.]L.*","",coldata$sample)
coldata$genome<-gsub(".*[.]","",coldata$sample)

# prepare for exploratory plots
library(genefilter)
sumPCA<-
function (x, genes, intgroup = "condition", ntop = 500)
{
    rv = rowVars(assay(x)[genes,])
    select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(x)[genes,][select, ]))
    tt = pca$x[,c("PC1","PC2")]
    tt = cbind(tt, coldata[rownames(tt),] )
    attr(tt,"percentVar")<-summary(pca)$importance[2,1:2]
    return(tt)
}
# consider flowering time genes from Corrinne
ftgenes<-as.character(read.table("~/ML/floweringtime.genes")$V1 )
library(DESeq2)
library(ggplot2)
dds <- DESeqDataSetFromMatrix( countData = count13, colData = coldata, design = ~ sample)
rld <- rlog(dds, blind=FALSE)
save(coldata, count13, rld, ftgenes, file="MLcounts.RData")

# now ploting
pdf("MLplots.pdf")
### 1: plotPCA with top 500, default
plotPCA(rld, intgroup=c("condition", "genome"))
data <- plotPCA(rld, intgroup=c("condition", "genome"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
### 2. same as plotPCA with shape
ggplot(data, aes(PC1, PC2, color=condition, shape=genome)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
ggtitle("Top 500 genes - plotPCA default")
### 3. consider only flowing genes
data <- sumPCA(rld, genes= ftgenes, intgroup=c("condition", "genome"), ntop = length(ftgenes) )
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data=data, aes(PC1, PC2, color=condition, shape=genome)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
ggtitle("Flowering genes")
### 4. Heatmap of the count matrix for flower genes
library("pheatmap")
pheatmap(assay(rld)[ftgenes,], cluster_rows=TRUE, show_rownames=TRUE,cluster_cols=TRUE, annotation_col=coldata[,c("condition","genome")], fontsize_row = 4, fontsize_col = 8)
###
dev.off()

pdf("ML.ftgenes.pdf")
perSum<- colSums(count13[ftgenes,])/colSums(count13) * 100
boxplot(perSum~coldata$condition, main = "Percentage Count Sum, by condition" )
boxplot(perSum~coldata$genome, main = "Percentage Count Sum, by genome"  )
boxplot(perSum~coldata$condition*coldata$genome , las=2, , main = "Percentage Count Sum, by genome and condition")
dev.off()

#####  Some conclusions
# Expression profiles of flowering genes are dominated by difference between sunrise(LD7 and SD7) and sunset(LD9 and SD5)
# Total expression of flowering genes (percentage of total counts) are higher in sunrise vs sunset samples.



# multifactor design
dds <- DESeqDataSetFromMatrix( countData = count13, colData = coldata, design = ~ condition + genome)
dds<-DESeq(dds)
batch <- rbind(
c("condition","LD7","LD9"),
c("condition","SD7","SD5"),
c("condition","LD7","SD7"),
c("condition","LD9","SD5"),
c("genome", "A2","D5"),
c("genome", "Maxxa","A2xD5"),
c("genome", "Maxxa","A2"),
c("genome", "Maxxa","D5"),
c("genome", "A2xD5","A2"),
c("genome", "A2xD5","D5")
)
apply(batch,1,function(x) { res <- results(dds, x); print(x); print( summary(res,alpha=.05) ); write.table(res, file=paste("DE/",paste(x[2:3], collapse="vs"),".txt", sep=""), sep="\t")})

# "condition" "LD7"       "LD9"      
# out of 35814 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 5818, 16% 
# LFC < 0 (down)   : 6467, 18% 

# "condition" "SD7"       "SD5"      
# LFC > 0 (up)     : 7270, 20% 
# LFC < 0 (down)   : 7287, 20% 

#  "condition" "LD7"       "SD7"      
# LFC > 0 (up)     : 1236, 3.5% 
# LFC < 0 (down)   : 1012, 2.8% 

# "condition" "LD9"       "SD5"      
# LFC > 0 (up)     : 1877, 5.2% 
# LFC < 0 (down)   : 2407, 6.7% 

# So, morning and night are quite different, and long day and short day are more different at sunset than sunrise.


