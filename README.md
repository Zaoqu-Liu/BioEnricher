# BioEnricher
## Integrate analysis and visualization for bioinformatic enrichment analyzer
The primary strength of BioEnricher lies in addressing two issues: firstly, it facilitates the seamless integration for enrichment analysis, encompassing diverse functionalities such as GO, KEGG, WikiPathways, Reactome, MsigDB, Disease Ontology, Cancer Gene Network, DisGeNET, CellMarker, and CMAP (drugs); secondly, it encapsulates advanced visualization functions, streamlining the process for faster and more convenient data presentation.

## Installation
You can install the released version of BioEnricher from Github with:
```R
devtools::install_github("Zaoqu-Liu/BioEnricher")

# Examples
## Get an interested gene list (for ORA) or an order-ranked geneList (for GSEA)
You should identify an interested gene list or an order-ranked geneList by employing differential analysis or other methods.
```R
library(airway)
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
data(airway)
se <- airway
se$dex <- relevel(se$dex, "untrt") 
dds <- DESeqDataSet(se, design = ~ cell + dex)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)                   
res <- as.data.frame(results(dds))%>%na.omit()
ann <- bitr(rownames(res),'ENSEMBL','SYMBOL',org.Hs.eg.db)
res <- merge(ann,res,by.x=1,by.y=0)%>%distinct(SYMBOL,.keep_all = T) # Very crude, just as an example
