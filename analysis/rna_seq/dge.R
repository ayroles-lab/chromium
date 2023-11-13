library(data.table)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

sampleData = read_csv("/Genomics/ayroleslab2/scott/git/chromium/metadata/acute_metdata.csv")
rownames(sampleData) = sampleData$ENA_RUN
sampleData$individual = as.factor(sampleData$individual)
sampleData$paris_classification = as.factor(sampleData$paris_classification)

# Use relevel() to set adjacent normal samples as reference
sampleData$paris_classification = relevel(sampleData$paris_classification, ref = "normal")

files = list.files(paste0(dir, "star"), "*ReadsPerGene.out.tab$", full.names = T)