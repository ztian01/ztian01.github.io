# RNASeq Analysis - Discovery of Differentially Expressed Genes
**Languages:** R <br>
**Softwares / packages:** RStudio, DESeq2 <br>
**Data:**  Himes, Blanca E., et al. "RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive
 gene that modulates cytokine function in airway smooth muscle cells." PloS one 9.6 (2014): e99625.
## Normalization
RNASeq data of untreated and dexamethasone-treated samples was downloaded and processed. 
Reads quantification was conducted as described in [Reads Visualization and Quantification](./RNASeq-IGV.md)<br>
![rawCounts](./RNASeq%20Implementation/20220623/Records/rawCounts.png)<br>
```<R>
options(stringsAsFactors = FALSE)
raw_count <- read.table("readCount-edited.txt", sep = "\t", header = T)
row.names(raw_count)<-make.names(raw_count[,1],TRUE)
raw_count<-raw_count[,-1]
# Delete 0 rows
raw_count= raw_count[which(rowSums(raw_count==0)==0),]

condition <- factor(c("control","treat","control","treat","control","treat","control","treat"),levels = c("control","treat"))
colData <- data.frame(row.names=colnames(raw_count), condition)

colnames(raw_count) <-c("X1-UT", "X1-T", "X2-UT", "X2-T", "X3-UT", "X3-T", "X4-UT", "X4-T")
library(DESeq2)
dds <- DESeqDataSetFromMatrix(raw_count, colData, design= ~ condition)
dds = dds[rowSums(counts(dds)) > 10 ,]
dim(dds)

# normalization
rld <- rlogTransformation(dds)
# Output
write.csv(assay(rld), file="rld.csv")

```
![rld](./RNASeq%20Implementation/20220623/Records/rld.png)<br>

## DEG analysis
```<R>
# dds from normalization.R
norm_dds <- DESeq(dds)
res <- results(norm_dds)
res05 <- results(norm_dds, alpha=0.05)
res = res[order(res$padj),]
head(res)

# Delete NA
which(is.na(res$padj)) -> na
res = res[-na,]
dim(res)
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

# Converting gene id to gene names
library(biomaRt)
listMarts()
plant <- useMart("ensembl")
listDatasets(plant)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listFilters(mart)
ensembl_gene_id_version<-row.names(res)
hg_symbols<- getBM(attributes=c('ensembl_gene_id_version','external_gene_name',"description"), filters = 'ensembl_gene_id_version', values = ensembl_gene_id_version, mart = mart)
res <- cbind(ensembl_gene_id_version, res)

res <- merge(res,hg_symbols, by="ensembl_gene_id_version")
write.csv(res,file="All_results.csv")

diff_gene <-subset(res, padj < 0.05 & abs(log2FoldChange) > 2)
write.csv(diff_gene,file= "DEG_treat_vs_control.csv")
```
![AllResults](./RNASeq%20Implementation/20220623/Records/AllResults.png)<br>

## Visualization in volcano plot
```<R>
library(ggplot2)
library(ggpubr)
library(ggthemes)

res$logp <- -log10(res$padj)
res = data.frame(res)
ggscatter(res,x= "log2FoldChange", y="logp") +theme_base()
res$type <- ifelse(res$padj < 0.05,ifelse(abs(res$log2FoldChange) > 2, ifelse(res$log2FoldChange < -2,'down','up'),'noSig'),'noSig')
ggscatter(res,x= "log2FoldChange", y="logp",color = "type", palette=c("green", "gray", "red"),size=1, xlab="log2FoldChange", ylab="-log10(padj)") +theme_base()+geom_hline(yintercept=-log(0.05,10), linetype="dashed")+geom_vline(xintercept=c(-2,2), linetype="dashed")
res$label = ""

res = res[order(res$padj), ]

up.gene = head(res$external_gene_name[which(res$type=="up")],10)
down.gene = head(res$external_gene_name[which(res$type=="down")],10)

res.top10.gene =c(as.character(up.gene),as.character(down.gene))
res$label[match(res.top10.gene,res$external_gene_name)]=res.top10.gene
ggscatter(res,x= "log2FoldChange", y="logp",color="type", palette=c("#2f5688", "#BBBBBB", "#CC0000"),size=1, label=res$label, font.label=8, repel=T, xlab="log2FoldChange", ylab="-log10(p-value)") +theme_base()+geom_hline(yintercept=-log(0.05,10), linetype="dashed")+geom_vline(xintercept=c(-2,2), linetype="dashed")
```
![volcano-with-label](./RNASeq%20Implementation/20220623/Records/volcano-with-label.png)