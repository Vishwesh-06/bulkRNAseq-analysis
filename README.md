# RNAseq-analysis
RNA sequencing analysis in R from scratch (only need NumReads)

RNA sequencing gives you raw reads as NumReads. I've done further analysis using the same:
(Project -  I had RNA sequencing data from 6 samples (2*2 diabetic and 2*1 control))

# 1. DESeq2
Commonly used tool for Differential Expression Gene Analysis.
Gives us an output as normalised read counts (used in PCA) and log2FoldChange and adjusted p value for each gene (used in ORA and GSEA)

# 2. PCA
Exploratory analysis using normalised counts from DESeq2/ EdgeR.
Tells you how your data is clustered, you can then analyse it accordingly.
(Project - I had 2 separate clusters forming on PC1 unrelated to the diabetic status, so I analysed the 2 clusters (2 diabetic and 1 control) separately, as I was interested in the changes in diabetes)

# 3. GSEA
Gives you biological insights about the differentially expressed genes. 
Using various libraries (Gene Ontology, Reactome, WikiPathways etc) and custom gene set.
