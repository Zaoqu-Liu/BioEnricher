## :bar_chart: Overview
## Integrate analysis and visualization for bioinformatic enrichment analyzer
BioEnricher lies in addressing two issues: firstly, it facilitates the seamless integration for enrichment analysis, encompassing diverse functionalities such as GO, KEGG, WikiPathways, Reactome, MsigDB, Disease Ontology, Cancer Gene Network, DisGeNET, CellMarker, and CMAP (drugs); secondly, it encapsulates advanced visualization functions, streamlining the process for faster and more convenient data presentation.

## :arrow_double_down: Installation
**You can install the released version of BioEnricher from Github with:**
```R
devtools::install_github("Zaoqu-Liu/BioEnricher")
```
## :beginner: Examples
### Get an interested gene list (for ORA) or an order-ranked geneList (for GSEA)
**You should identify an interested gene list or an order-ranked geneList by employing differential analysis or other methods.**
```R
library(airway)
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
data(airway)
se <- airway
se$dex <- relevel(se$dex, "untrt") 
res <- DESeqDataSet(se, design = ~ cell + dex)%>%
  estimateSizeFactors()%>%DESeq()%>%
  results()%>%as.data.frame()%>%na.omit()
ann <- bitr(rownames(res),'ENSEMBL','SYMBOL',org.Hs.eg.db)
res <- merge(ann,res,by.x=1,by.y=0)%>%distinct(SYMBOL,.keep_all = T) # Very crude, just as an example
```
```R
# Define an up-regulated gene list
up.genes <- res$SYMBOL[res$log2FoldChange > 2 & res$padj < 0.05]
# Define a down-regulated gene list
down.genes <- res$SYMBOL[res$log2FoldChange < -2 & res$padj < 0.05]
```
## :paperclip: ORA
**This function will perform over-representative analysis including GO, KEGG, WikiPathways, Reactome, MsigDB, Disease Ontoloty, Cancer Gene Network, DisGeNET, CellMarker, and CMAP.**
```R
# Set enrich.type using an enrichment analysis method mentioned above.
kegg <- lzq_ORA(
  genes = res$SYMBOL[res$log2FoldChange > 0 & res$padj < 0.05],
  enrich.type = 'KEGG'
)

# This function will output its calculation process.
+++ Updating gene symbols...
Maps last updated on: Thu Oct 24 12:31:05 2019
+++ Transforming SYMBOL to ENTREZID...
'select()' returned 1:1 mapping between keys and columns
+++ Performing KEGG enrichment...
+++ 109 significant terms were detected...
+++ Done!
```
## :page_facing_up: Simple visualization of KEGG pathway based on the pathview package
```R
res2 <- res[res$log2FoldChange > 0 & res$padj < 0.05,c(2,4)]
res2 <- data.frame(row.names = res2$SYMBOL,R=res2$log2FoldChange)

lzq_KEGGview(gene.data = res2,pathway.id = 'hsa04218')
```
<img src="man/hsa04510.png" width="100%" />

## :paperclip: ORA.integrated
**This function will perform an integration for ORA enrichment analysis, including GO, KEGG, WikiPathways, Reactome, MsigDB, Disease Ontology, Cancer Gene Network, DisGeNET, CellMarker, and CMAP (drugs).**
```R
library(BioEnricher)
# Integrative enrichment analysis of the up-regulated gene list
up.enrich <- lzq_ORA.integrated(
  genes = up.genes,
  background.genes = NULL,
  GO.ont = 'ALL',
  perform.WikiPathways = T,
  perform.Reactome = T,
  perform.MsigDB = T,
  MsigDB.category = 'ALL',
  perform.Cancer.Gene.Network = T,
  perform.disease.ontoloty = T,
  perform.DisGeNET = T,
  perform.CellMarker = T,
  perform.CMAP = T,
  min.Geneset.Size = 3
)

# This function will output its calculation process.
+++ Updating gene symbols...
Maps last updated on: Thu Oct 24 12:31:05 2019
+++ Transforming SYMBOL to ENTREZID...
'select()' returned 1:1 mapping between keys and columns
+++ Performing GO-ALL enrichment...
+++ Symplifying GO results...
+++ Performing KEGG enrichment...
+++ Performing Module KEGG enrichment...
+++ Performing WikiPathways enrichment...
+++ Performing Reactome pathways enrichment...
+++ Performing Disease Ontoloty enrichment...
+++ Performing Cancer Gene Network enrichment...
+++ Performing DisGeNET enrichment...
+++ Performing CellMarker enrichment...
+++ Performing MsigDB-ALL enrichment...                                               
+++ Performing CMAP enrichment...
+++ 1765 significant terms were detected...
+++ Done!
```
```R
# Integrative enrichment analysis of the down-regulated gene list
down.enrich <- lzq_ORA.integrated(
  genes = down.genes,
  background.genes = NULL,
  GO.ont = 'ALL',
  perform.WikiPathways = T,
  perform.Reactome = T,
  perform.MsigDB = T,
  MsigDB.category = 'ALL',
  perform.Cancer.Gene.Network = T,
  perform.disease.ontoloty = T,
  perform.DisGeNET = T,
  perform.CellMarker = T,
  perform.CMAP = T,
  min.Geneset.Size = 3
)

# This function will output its calculation process.
+++ Updating gene symbols...
Maps last updated on: Thu Oct 24 12:31:05 2019
+++ Transforming SYMBOL to ENTREZID...
'select()' returned 1:1 mapping between keys and columns
+++ Performing GO-ALL enrichment...
+++ Symplifying GO results...
+++ Performing KEGG enrichment...
+++ Performing Module KEGG enrichment...
+++ Performing WikiPathways enrichment...
+++ Performing Reactome pathways enrichment...
+++ Performing Disease Ontoloty enrichment...
+++ Performing Cancer Gene Network enrichment...
+++ Performing DisGeNET enrichment...
+++ Performing CellMarker enrichment...
+++ Performing MsigDB-ALL enrichment...                                               
+++ Performing CMAP enrichment...
+++ 1426 significant terms were detected...
+++ Done!
```
## :page_facing_up: Visualization for one ORA enrichment object
**barplot**
```R
lzq_ORA.barplot1(enrich.obj = up.enrich$simplyGO)
```
<img src="man/GO1.jpg" width="60%" />

**dotplot**
```R
lzq_ORA.dotplot1(enrich.obj = up.enrich$simplyGO)
```
<img src="man/GO2.jpg" width="60%" />

## :page_facing_up: Visualization for two types of ORA enrichment objects
```R
lzq_ORA.barplot2(
  enrich.obj1 = up.enrich$simplyGO,
  enrich.obj2 = down.enrich$simplyGO,
  obj.types = c('Up','Down')
)
```
<img src="man/Two-types-GO.png" width="60%" />

### You can translate the terms in the graph into Chinese using use.Chinese = T
```R
lzq_ORA.barplot2(
  enrich.obj1 = up.enrich$simplyGO,
  enrich.obj2 = down.enrich$simplyGO,
  obj.types = c('Up','Down'),
  use.Chinese = T
)
```
<img src="man/Two-types-GO-Chinese.jpg" width="35%" />

**Note: use.Chinese exists all the plot functions.**

## :paperclip: GSEA
**This function will perform gene-set enrichment analysis including GO, KEGG, WikiPathways, Reactome, MsigDB, Disease Ontoloty, Cancer Gene Network, DisGeNET, CellMarker, and CMAP.**
```R
# Obtain an order ranked geneList.
grlist <- res$log2FoldChange; names(grlist) <- res$SYMBOL
grlist <- sort(grlist,decreasing = T)

# Set enrich.type using an enrichment analysis method mentioned above.
fit <- lzq_GSEA(grlist,enrich.type = 'KEGG')

# This function will output its calculation process.
+++ Updating gene symbols...
Maps last updated on: Thu Oct 24 12:31:05 2019
+++ Transforming SYMBOL to ENTREZID...
'select()' returned 1:many mapping between keys and columns
+++ Performing KEGG enrichment...
+++ 8 significant terms were detected...
+++ Done!
```
## :paperclip: GSEA.integrated
**This function will perform an integration for GSEA enrichment analysis, including GO, KEGG, WikiPathways, Reactome, MsigDB, Disease Ontology, Cancer Gene Network, DisGeNET, CellMarker, and CMAP (drugs).**
```R
# Integrative enrichment analysis of the up-regulated gene list
fit2 <- lzq_GSEA.integrated(
  genes = grlist,
  gene.type = 'SYMBOL',
  GO.ont = 'ALL',
  perform.WikiPathways = T,
  perform.Reactome = T,
  perform.MsigDB = T,
  MsigDB.category = 'ALL',
  perform.Cancer.Gene.Network = T,
  perform.disease.ontoloty = T,
  perform.DisGeNET = T,
  perform.CellMarker = T,
  perform.CMAP = T,
  min.Geneset.Size = 3
)

# This function will output its calculation process.
+++ Updating gene symbols...
Maps last updated on: Thu Oct 24 12:31:05 2019
+++ Transforming SYMBOL to ENTREZID...
'select()' returned 1:many mapping between keys and columns
+++ Performing GO-ALL enrichment...
+++ Symplifying GO results...
+++ Performing KEGG enrichment...
+++ Performing Module KEGG enrichment...
+++ Performing WikiPathways enrichment...
+++ Performing Reactome pathways enrichment...
+++ Performing Disease Ontoloty enrichment...
+++ Performing Cancer Gene Network enrichment...
no term enriched under specific pvalueCutoff...
+++ Performing DisGeNET enrichment...
+++ Performing CellMarker enrichment...
+++ Performing MsigDB-ALL enrichment...                                               
+++ Performing CMAP enrichment...
no term enriched under specific pvalueCutoff...
+++ 311 significant terms were detected...
+++ Done!
```
## :page_facing_up: Visualization for positive or negative GSEA enrichment results
**Visualize analyzing result of GSEA**
```R
lzq_gseaplot(
  fit2$simplyGO,
  Pathway.ID = 'GO:0030016',
  rank = F,
  statistic.position = c(0.71,0.85),
  rel.heights = c(1, 0.4)
)
```
<img src="man/GSEA.jpg" width="60%" />

**Enrichment barplot for positive or negative GSEA results**
```R
lzq_GSEA.barplot1(enrich.obj = fit2$simplyGO,type = 'pos')
```
<img src="man/GSEA1.jpg" width="60%" />

**Enrichment dotplot for positive or negative GSEA results**
```R
lzq_GSEA.dotplot1(enrich.obj = fit2$simplyGO,type = 'pos')
```
<img src="man/GSEA2.jpg" width="60%" />

## :page_facing_up: Enrichment barplot for positive and negative GSEA results
```R
lzq_GSEA.barplot2(enrich.obj = fit2$simplyGO)
```
<img src="man/GSEA4.jpg" width="60%" />

### You can translate the terms in the graph into Chinese using use.Chinese = T
```R
lzq_GSEA.barplot2(enrich.obj = fit2$simplyGO,use.Chinese = T)
```
<img src="man/GSEA3.jpg" width="45%" />

**Note: use.Chinese exists all the plot functions.**


