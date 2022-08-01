### 1. PCA analysis(StatQuest). From: https://github.com/StatQuest/pca_demo/blob/master/pca_demo.R

#open as csv so that gene names can be row names
#always check if gene names are being read as dates by excel
TPMc <- read.csv("C:/Users/zVishwesh2/Desktop/TPMc.csv", row.names=1) #used TPM from edgeR by company analysis. can used normCounts from DESeq2 also
#transposing so that gene names are columns
TPMt = t(TPMc)
#removing 0 variance columns
TPM1 = TPMt[, colSums(TPMt != 0) > 0] 

pca <- prcomp(TPM1, scale=TRUE)

#scree plot
pca.var <- pca$sdev^2 #variation per PC
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) #percentages
barplot(pca.var.per, main="Scree plot", xlab= "Principal Component", ylab="Percent Variation")

#PCA1-2 T1D Organoids
library(ggplot2)
pca.data_12 <- data.frame(Sample=rownames(pca$x), 
                          X=pca$x[,1],
                          Y=pca$x[,2])
pca.data_12

ggplot(data=pca.data_12, aes(x=X, y=Y, label=Sample)) + 
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("T1D organoids PCA1-2")

#PCA1-3 T1D Organoids
pca.data_13 <- data.frame(Sample=rownames(pca$x), 
                          X=pca$x[,1],
                          Y=pca$x[,3])
pca.data_13

ggplot(data=pca.data_13, aes(x=X, y=Y, label=Sample)) + 
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC3 - ", pca.var.per[3], "%", sep="")) +
  theme_bw() +
  ggtitle("T1D organoids PCA1-3")

#The top 50 genes that contribute most to pc1
loading_scores_1 <- pca$rotation[,1]
gene_scores_1 <- abs(loading_scores_1) #Magnitudes
gene_score_ranked_1 <- sort(gene_scores_1, decreasing=TRUE)
top_50_genes_PC1 <- names(gene_score_ranked_1[1:50])

top_50_genes_PC1 #Names of the top 50 genes

pca$rotation[top_50_genes_PC1,1] #The scores (and +/- sign)

#The top 50 genes that contribute most to pc2
loading_scores_2 <- pca$rotation[,2]
gene_scores_2 <- abs(loading_scores_2) #Magnitudes
gene_score_ranked_2 <- sort(gene_scores_2, decreasing=TRUE)
top_50_genes_PC2 <- names(gene_score_ranked_2[1:50])

top_50_genes_PC2 #Names of the top 20 genes 

pca$rotation[top_50_genes_PC2,1] #The scores (and +/- sign)

### 2.DESeq2 analysis (Mike Vandewege YT video) From: https://youtu.be/wPzeea1Do18

library(DESeq2)
library(apeglm)
library(ggplot2)
library(pheatmap)
library(cowplot)

#length (Counts) file
#Analyse UD and Dif separately, as toooo many differences to analyse together.
Lengthc <- read.csv("C:/Users/zVishwesh2/Desktop/NumReadUD.csv", row.names=1)
#loading column data
coldata <- read.csv("C:/Users/zVishwesh2/Desktop/coldataUD.csv", row.names=1)

#deseq2 using Diabetes as differentiator
dds <- DESeqDataSetFromMatrix(round(Lengthc), coldata, ~Diabetes) #round to remove decimals

#removing lowly expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#main DESeq
ddsDE <- DESeq(dds)

#export normalized read counts
normCounts <- counts(ddsDE, normalized = T)
write.csv(normCounts, "normal.T1D(UD).csv")

#DESeq results
res <- results(ddsDE, alpha =0.05) #p value<0.05, results
#summary(res)
resOrdered <- res[order(res$padj),] #ordering by padj
write.csv(resOrdered, "deSeq.T1D(UD).csv")

plotMA(ddsDE)

#load saved csv, results of DESeq2
deSeq.T1D <- read.csv("C:/Users/zVishwesh2/Desktop/DESeq2 (UD)/deSeq.T1D(UD).csv", row.names=1)
#add column w yes if padj<0.05
deSeq.T1D$sig <- ifelse(deSeq.T1D$padj <= 0.05, "yes", "no")
#reomving NAs
deSeq.T1D <- na.omit(deSeq.T1D)

maPlot <- ggplot(deSeq.T1D, aes(x = log10(baseMean), y = log2FoldChange, color = sig)) + geom_point()

#volcano plot
volcano <- ggplot(deSeq.T1D, aes(x = log2FoldChange, y = -log10(padj), color = sig)) + geom_point()

#pheatmap
signi <- subset(deSeq.T1D, padj <= 0.05)
#merging the significant genes w normalized counts, by 0 ie rows
##need to merge as pheatmap needs normalized read counts
allSig <- merge(normCounts, signi, by = 0) 
#need only the normalized read counts
sigCounts <- allSig[,2:4]
row.names(sigCounts) <- allSig$Row.names

pheatmap(sigCounts)
pheat <- pheatmap(log2(sigCounts + 1), scale = 'row', show_rownames = F, treeheight_row = 0) #log to correct for outlier 
##can add pheatmap(log2( +1),scale = , show_rownames = F, treeheight_row = 0, treeheight_col = 0)

plot_grid(maPlot, volcano)
plot_grid(volcano,pheat[[4]], labels = c("A","B"))

### 3. GSEA analysis (Tutorial). From: https://yulab-smu.top/biomedical-knowledge-mining-book/index.html
#Need output of DESeq2 for GSEA (ranked list, normCounts)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggridges)
library(ggplot2)
library(enrichplot)
library(ReactomePA)

#getting entrezgene_id
library(biomaRt)
gene_id <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
                 mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl")))
df <- read.csv("C:/Users/zVishwesh2/Desktop/DESeq2 (UD)/deSeq.T1D(UD).csv")

#Only differentially expressed genes (Up and Down) for ORA
up_df <- na.omit(df[(df$log2FoldChange > 0) & (df$padj <= 0.05), ])
down_df <- na.omit(df[(df$log2FoldChange < 0) & (df$padj <= 0.05), ])
up_df <- merge(up_df, gene_id[,c(2,3)], by.x = "X", by.y = "external_gene_name")
down_df <- merge(down_df, gene_id[,c(2,3)], by.x = "X", by.y = "external_gene_name")


#All genes for GSEA
df <- merge(df, gene_id[,c(2,3)], by.x = "X", by.y = "external_gene_name") #all genes GSEA
geneList = df[,3] ## feature 1: numeric vector (log2FoldChange)
names(geneList) = as.character(df[,8]) ## feature 2: named vector (entrezgeneid)
geneList = sort(geneList, decreasing = TRUE) ## feature 3: decreasing order

##ORA:
#GO
ego_up  <- enrichGO(gene       = up_df$entrezgene_id,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENTREZID',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
barplot(ego_up, showCategory = 20) + ggtitle("barplot for Up(UD) GO_ORA")
dotplot(ego_up, showCategory=20) + ggtitle("dotplot for Up(UD) GO_ORA")

ego_down  <- enrichGO(gene        = down_df$entrezgene_id,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENTREZID',
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
barplot(ego_down, showCategory = 20) + ggtitle("barplot for Down(UD) GO_ORA")
dotplot(ego_down, showCategory=20) + ggtitle("dotplot for Down(UD) GO_ORA")
#KEGG
ekk_up <- enrichKEGG(gene     = up_df$entrezgene_id,
                     organism     = 'hsa',
                     keyType       = 'kegg',
                     pvalueCutoff = 0.05)
barplot(ekk_up, showCategory = 20) + ggtitle("barplot for Up(UD) KEGG_ORA")
dotplot(ekk_up, showCategory=20) + ggtitle("dotplot for Up(UD) KEGG_ORA")

ekk_down <- enrichKEGG(gene     = down_df$entrezgene_id,
                       organism     = 'hsa',
                       keyType       = 'kegg',
                       pvalueCutoff = 0.05)
barplot(ekk_down, showCategory = 20) + ggtitle("barplot for Down(UD) KEGG_ORA")
dotplot(ekk_down, showCategory=20) + ggtitle("dotplot for Down(UD) KEGG_ORA")
#Reactome
reac <- enrichPathway(gene=up_df$entrezgene_id, pvalueCutoff = 0.05, readable=TRUE)
barplot(reac, width = 0.5, xlim= 2, space = 10) + ggtitle("Up(UD) Reactome_ORA")
dotplot(reac) + ggtitle("Dotplot: Up(UD) Reactome_ORA")

reac_down <- enrichPathway(gene=down_df$entrezgene_id, pvalueCutoff = 0.05, readable=TRUE)
barplot(reac_down)
#WikiPathways
WP <- enrichWP(up_df$entrezgene_id, organism = "Homo sapiens")
barplot(WP)

WP_down <- enrichWP(down_df$entrezgene_id, organism = "Homo sapiens")
barplot(WP_down)

##GSEA
#GO
gsGO  <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               ont          = "ALL",
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
ridgeplot(gsGO, showCategory = 25)
#KEGG
gsKG <- gseKEGG(geneList     = geneList,
                organism     = 'hsa',
                minGSSize    = 120,
                pvalueCutoff = 0.05,
                verbose      = FALSE)
ridgeplot(gsKG, showCategory = 25)
#Reactome
gsRC <- gsePathway(geneList, 
                   pvalueCutoff = 0.2,
                   pAdjustMethod = "BH", 
                   verbose = FALSE)
ridgeplot(gsRC)
#best GSEA:Reactome. treeplot
p1 <- pairwise_termsim(gsRC)
treeplot(p1, hclust_method = "average", fontsize = 2.2, margin =0.2) + ggtitle("GSEA_Reactome(UD)")
#WP
gsWP <- gseWP(geneList, organism = "Homo sapiens")
ridgeplot(gsWP)

### 4. Custom GSEA (Using self gene set .gmt format)
# Here: using cell type markers from Lu J Sci Rep

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(dplyr)

# reading in data from deseq2
df = read.csv("C:/Users/zVishwesh2/Desktop/DESeq2 (UD)/deSeq.T1D(UD).csv")
#getting entrezgene_id
library(biomaRt)
gene_id <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
                 mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl")))

df <- merge(df, gene_id[,c(2,3)], by.x = "X", by.y = "external_gene_name") #all genes GSEA

## feature 1: numeric vector
geneList = df[,3]
## feature 2: named vector
names(geneList) = as.character(df[,8])
## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

gse <- gseGO(geneList=geneList, 
             ont ="ALL",
             OrgDb = organism)

ridgeplot(gse)

#loading cell type data
celltypes <- read.gmt("C:/Users/zVishwesh2/Desktop/celltypesF.gmt.txt")
#getting entrez id
my.symbols <- celltypes$gene
entrz <- mapIds(org.Hs.eg.db, my.symbols, 'ENTREZID', 'SYMBOL')
celltypes$entrz = na.omit(entrz)

t2g <- celltypes %>% select(-gene) #t2g - term:cell type & g:entrez

y <- GSEA(geneList, TERM2GENE = t2g)

p1 <- gseaplot2(y, geneSetID = 1, title = "Early_Enterocytes_Progenitors")
p2 <- gseaplot2(y, geneSetID = 2, title = "Enteroendocrine cells")
p3 <- gseaplot2(y, geneSetID = 3, title = "Goblet cells")

cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
