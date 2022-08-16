#Analises de DEGs de ARVC do GSE164490
library(DESeq2)
getwd()

cts = read.delim("AllReadCount.tab", sep="", header = T)
dim(cts)


Name=colnames(cts)
Time=c(rep("Control", 6), rep("Case",12))
coldata=data.frame(Name,Time)
head(coldata,18)
dim(coldata)

log_cts<- log(cts+1, 10) # This is log base 10 + 1 for "0"
table(rowSums(log_cts>log10(3))>=18)

keep.exprs<- rowSums(log_cts>log10(3))>=18
cts_filt<-cts[keep.exprs,]
dim(cts_filt)
head(cts_filt,10)

dds <- DESeqDataSetFromMatrix(countData = round(cts_filt), colData = coldata, design = ~ Time)
dds$Time <- factor(dds$Time, levels = c("Control","Case"))
head(dds)

dds <- DESeq(dds)
res <- results(dds)
head(res)
dim(res)
data.res <- as.data.frame(res)
write.csv(data.res,"res_arvc.csv", row.names = T)
