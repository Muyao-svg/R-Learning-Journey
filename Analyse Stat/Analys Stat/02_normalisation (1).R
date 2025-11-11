# <muyao.guo@etu-upsaclay.fr>
# 19/05/2024
# réaliser un travail de “normalisation” des valeurs de comptage dans le fichier de données. 
# L’objectif est de rendre comparables les mesures d’expression obtenues dans des expériences différentes. 
# Normalisation by "DESeq2". Ce script suppose que 01_import_data.R a déjà été exécuté via main.R
#===============================

# Package installation (if needed)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

library(DESeq2)
countData <- matrix_RNA_seq_sort

# Setting the condition
group_annotation_all = data.frame(
  Sample = c(group_a,group_b,group_c),
  Group = c(
    rep("A", length(group_a)),
    rep("B", length(group_b)),
    rep("C", length(group_c))
  )
)
group_annotation = group_annotation_all %>% 
  filter(Sample %in% colnames(matrix_RNA_seq_sort)) %>% 
  column_to_rownames("Sample")
group_annotation = group_annotation[colnames(matrix_RNA_seq_sort), ,drop = FALSE]

condition = group_annotation$Group
length(condition)  # 36

# colData
colData = data.frame(row.names = colnames(countData),
                     condition = factor(condition))
#  DESeq2
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~condition)

# Affichage aux données de comptage
View(counts(dds))    
# Réalisation de la normalisation des comptages
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
# Table des comptages normalisés
normCounts <- counts(dds, normalized = TRUE)
# Ecriture des résultats
write.table(normCounts, file = "Résultats/matrix_RNA_seq_norm",
            quote = F, sep = "\t", col.names = T, row.names = T)
