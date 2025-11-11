## Load packages ---------------------------------------------------------------

## Installation:
## install.packages("BiocManager")
## BiocManager::install(c("FactoMineR", "ggplot2","edgeR", "limma", "csaw", "GenomicRanges","tidyr", "dplyr"))
# 在时间轴上同时追踪“表达层（转录）变化”和“调控层（染色质可及性/TF 调控）变化”，从而重建促进视网膜神经节细胞
#（RGC）轴突再生的基因调控程序。

## Graphics plots
library(FactoMineR)
library(ggplot2)

## Differential analysis
library(edgeR)
library(limma)

## ATAC-seq packages
library(csaw)
library(GenomicRanges)

## General data manipulation
library(tidyr)
library(dplyr)

## Load data and create experimental design ------------------------------------

atacseq <- read.table("TP3Zebrafish/ZebraFish/Dhara_atacseq.txt", header=TRUE, sep = "\t")
rnaseq <- read.table("TP3Zebrafish/ZebraFish/Dhara_rnaseq.txt", header=TRUE, sep = "\t")
expdes <-
    data.frame(sample = colnames(atacseq),
               time = unlist(strsplit(substr(colnames(atacseq),start=5, stop=6),
                                      split = "_")))

# substr(x, start=n1, stop=n2)  Extract or replace substrings in a character vector.
head(expdes)
## => Check the dimension and characteristics of each dataset

## Exploratory analyses --------------------------------------------------------

pseudo_rnaseq <- log2(rnaseq + 0.5)
pseudo_atacseq <- log2(atacseq + 0.5)
pca_rnaseq <- PCA(t(pseudo_rnaseq), scale.unit = TRUE, graph = FALSE)
pca_atacseq <- PCA(t(pseudo_atacseq), scale.unit = TRUE, graph = FALSE)
plot(pca_rnaseq, title = "RNA-seq")
plot(pca_atacseq, title = "ATAC-seq")

## RNA-seq  --------------------------------------------------------------------

### Normalization --------------------------------------------------------------

dge_rnaseq <- DGEList(rnaseq, group = rep(1, ncol(rnaseq)),
               remove.zeros = TRUE, samples = expdes) #给所有样本分配同一个组（group 1）。
# 这通常表示你现在还没有分组（例如只是想先过滤或质控，还没做差异分析）。
# 这一步会为每个样本计算一个 normalization factor（标准化因子）;
# calcNormFactors() 是在计算“校准系数”，让不同样本的 read counts 变得在统计上可比。
dge_tmm <- calcNormFactors(dge_rnaseq, method = "TMM")
# 
dge_tmm$samples  # 查看normalisation size factor
# RNA-seq 数据标准化之后的数据转换步骤，常用于可视化（比如PCA、热图）或聚类分析
pseudo_tmm <- log2(cpm(dge_tmm) + 0.5) 
# CPM:count per millions,这样消除了测序深度的影响，让不同样本之间可比较
#RNA-seq 数据中很多基因在某些样本中读数为 0。
#而 log(0) = −∞，数学上不定义。
#所以我们加一个很小的值（比如 0.5 或 1）来避免报错，这个小值就叫 伪计数 (pseudo-count)。
normc_rnaseq <- cpm(dge_tmm)

## => Compare boxplots before and after normalization
## => Check PCA after normalization

### Differential analysis ------------------------------------------------------

design_matrix <- model.matrix(~factor(time), data = dge_tmm$samples)
# 它会把分组变量转换为哑变量（dummy variables），即 0/1 矩阵;
# 用于比较不同condition下的拟合比如wide type 是beta0，mutant是beta0+beta1（相对wild type）
# 描述每个样本的实验条件（组别、时间、批次等）的数学结构expression = intercept + X*beta; X-->design_matrix
keep <- filterByExpr(dge_tmm, design_matrix)
head(keep)  # 会给出true和false的结果
# filterByExpr():Determine which genes have sufficiently large counts to be retained in a statistical analysis.
dge_tmm <- dge_tmm[keep,]
dge_tmm <- estimateDisp(dge_tmm, design_matrix)
# edgeR 假设每个基因的 read count 服从 负二项分布：Y_ij~NB(u_ij,dispersion_i);
#RNA-seq 差异表达分析的核心是：
#“比较不同条件下基因的表达差异是否显著，而不被技术噪音或生物变异误导。”
#离散度估计是这一步的关键，因为它决定了统计检验的置信度：
#离散度大 → 数据波动大 → 检验显著性降低
#离散度小 → 数据更一致 → 更容易检出差异表达
fit_tmm <- glmFit(dge_tmm, design_matrix)  # Fit a negative binomial generalized log-linear model to the read counts for each gene.
lrt_tmm <- glmLRT(fit_tmm, coef=2:5)  # 采取2:5位的factor作为coef，见matrix design，这里代表factor time(12,2,4,7),time(0)作为reference
# 分别比较time(12)与reference的logFC差别，每个time都和time(0)比
# Rappel：u_ij = cj*q_ij
#         log2(q_ij) = xi*beta_j  xi此处就是X就是matrix design，样本原本的数据，即样本-时间矩阵
# glmLRT()用于比较difference entre time是否显著;对模型中的一个或多个系数（coefficient）执行似然比检验，以判断这些系数是否显著偏离 0。
# H0： beta0=beta1=beta2...=0; H1:至少有一个beta不是0。
# 只能比较time(12,2,4,7)相对于time(0)是否显著，不能比较time(12)相对于time(4)是否显著。
res_tmm <- topTags(lrt_tmm, n = nrow(dge_tmm))  # 查看检验结果 res_tmm$table
de_genes <- rownames(res_tmm$table)[which(res_tmm$table$FDR < 0.05)]

## => Check p-value distributions
hist(res_tmm$table$PValue)
## => How many differentially expressed genes?
length(de_genes)  # 7475
head(de_genes)

## ATAC-seq --------------------------------------------------------------------

### Normalization --------------------------------------------------------------

dge_atacseq <- DGEList(atacseq, group = rep(1, ncol(atacseq)),
                      remove.zeros = TRUE, samples = expdes)
# DGEList():它是一个 专门容纳“计数矩阵 + 样本信息 + 归一化参数 + 结果”的容器
# sample= 样本元信息就是描述每个样本“背景信息”的表格。它不是实验数据本身（不是测序读数），而是告诉分析软件：
# “这些列（样本）代表什么实验条件、来自哪里、属于哪个批次、是哪个个体的重复…
dge_loess <- normOffsets(dge_atacseq)
# normOffsets(): Normalize ChIP-Seq read counts to input control values, then test for significant enrichment relative to the control.
# normaliser le biais entre diffenrent condition echan.
# calcNormFactors() in rna-seq normalisation.ATAC不适用于TMM，因为TMM适用于大多数基因表达量不变，少数基因差异表达；
# ATAC-seq / ChIP-seq 这类“富集型”测序中，峰的数量和强度可能在不同程度上整体上升或下降；
# 它通过对 logCPM 值（对数化的 counts per million）执行 Loess 局部回归（LOcally Estimated Scatterplot Smoothing） 来估计 样本间系统偏移（offsets）
# Loess 归一化。在rna-seq中是：norm.factors
# offset：biais de composition. Chromatine ouvrir pas de meme maniere
# ATAC-seq 中，每个样本的 测序深度、文库复杂度、背景噪声 都不同。
#如果直接比较原始计数，显然会被这些差异误导。
#offset 的作用就是：
#在建模时告诉统计模型“每个样本的有效测序量是多少”，
#让模型在比较条件效应（β）时排除这些技术性偏差。

# TMM只适用于rnaseq是因为它们的GC%都是很稳定的cindition一样，但是TAXA-seq在不同condition'下染色体打开程度不一样
pseudo_loess <- cpm(dge_loess, offset = dge_loess$offset, log = TRUE,
                    prior.count = 0.5)
# offset=偏移量
# prior.count = 0.5: log(data+0.5)
normc_atacseq <- cpm(dge_loess, offset = dge_loess$offset)

## => Compare boxplots before and after normalization
## => Check PCA after normalization

limma::plotMA(atacseq[,c(1,3)], xlab="M", ylab="A", main = "")
abline(h = 0, col = "red")
# 草图,我们的目的是让散开的点，在归一化后可以尽量围绕在0附近
limma::plotMA(normc_atacseq[,c(1,3)], xlab="M", ylab="A", main = "")
abline(h = 0, col = "red")
# normaliser后的数据
# MA图:moyenne logFC-comptage

### Differential analysis ------------------------------------------------------

design_matrix <- model.matrix(~factor(time), data = dge_loess$samples)
keep <- filterByExpr(dge_loess, design_matrix)
dge_loess <- dge_loess[keep,]
dge_loess <- estimateDisp(dge_loess, design_matrix)
# plotBCV(dge_loess) 看一眼图
fit_loess <- glmFit(dge_loess, design_matrix)
lrt_loess <- glmLRT(fit_loess, coef=2:5)
res_loess <- topTags(lrt_loess, n = nrow(dge_loess))
de_peaks <- rownames(res_loess$table)[which(res_loess$table$FDR < 0.05)]
# 这一步找出了显著差异开放的染色质区域（differentially accessible peaks）
#在 ATAC-seq 里，“peak” 并不是“峰值”本身，而是：
#一个统计显著的 DNA 区域，代表该区域染色质开放、Tn5 酶可进入、被测序检测到的富集信号。
#这些区域通常是：启动子（promoter）区域,增强子（enhancer）转录因子结合位点（TF binding sites）
# 生物学解释：如果某个 peak 在某一组 组信号更强 → 说明该区域变得更开放（可能被激活）；如果信号更弱 → 表示染色质变紧（可能被抑制）。
# 为什么在 ATAC-seq 中更容易出现假阳性？
#因为：没有对照样本（不像 ChIP-seq 有 input DNA）；信号强度变化剧烈；染色质结构有局部噪音（mapping bias、GC 含量高区域）；
#PCR bias 或 duplication。所以如果不进行 FDR 控制，就容易把随机信号误判成“真正的开放区域”。



## => Check p-value distributions
hist(res_loess$table$PValue) 
## => How many differentially accessible peaks?
length(de_peaks)  # 492
## Genomic annotation of ATAC-seq ----------------------------------------------

### All peaks ------------------------------------------------------------------

atac_pos <- data.frame(peak_IDs = rownames(dge_loess$counts)) %>%
    separate(peak_IDs, into = c("chr", "start", "end")) %>%
    filter(chr %in% c(1:25))
#peak 的坐标信息转换成基因组区间对象（GRanges），
#从而让后续的分析（比如注释、可视化、重叠计算）都能在“基因组空间”上进行
atac_ranges <- GRanges(seqnames=paste0("chr", atac_pos[,1]),
                       IRanges(start=as.numeric(atac_pos[,2]),
                               end=as.numeric(atac_pos[,3])))

# library(ChIPseeker)
# library(TxDb.Drerio.UCSC.danRer10.refGene)
# anno <-
#   annotatePeak(atac_ranges, TxDb = TxDb.Drerio.UCSC.danRer10.refGene)
# plotAnnoPie(anno)
# gr_anno <- as.GRanges(anno)
# saveRDS(gr_anno, file="gr_anno.rds")
gr_anno <- readRDS("TP3Zebrafish/ZebraFish/gr_anno.rds")

## Focus on TSS regions
gr_tss<- gr_anno[abs(gr_anno$distanceToTSS) < 500]
gr_tss[1, ]

### Differential ATAC-seq peaks ------------------------------------------------
diffatac_pos <- data.frame(peak_IDs = rownames(dge_loess$counts),
                       FDR = res_loess$table$FDR) %>%
    separate(peak_IDs, into = c("chr", "start", "end")) %>%
    filter(chr %in% c(1:25), FDR < 0.05)
diffatac_ranges <- GRanges(seqnames=paste0("chr", diffatac_pos[,1]),
                       IRanges(start=as.numeric(diffatac_pos[,2]),
                               end=as.numeric(diffatac_pos[,3])))
diffanno <- subsetByOverlaps(gr_anno, diffatac_ranges)

## Integration of RNA-seq + ATAC-seq -------------------------------------------

### Qualitative -----------------------------------------------------------------

ensembl_to_refseq <- read.table("TP3Zebrafish/ZebraFish/ensembl_to_refseq.txt", sep = "\t",
                                header=TRUE)
diffanno_df <- left_join(data.frame(diffanno), ensembl_to_refseq,
                         by = "transcriptId")
de_atac_genes <- na.omit(unique(diffanno_df$ensemblId))
interesting_genes <- intersect(de_atac_genes, de_genes)
interesting_genes

### Quantitative ----------------------------------------------------------------

gene_choice <- "ENSDARG00000063518" # interesting_genes[1]
peak_choice <- diffanno_df %>% dplyr::filter(ensemblId == gene_choice) %>%
    dplyr::select(seqnames, start, end) %>%
    mutate(seqnames = paste0(seqnames, ":", start, "-", end)) %>%
    dplyr::select(seqnames) %>%
    mutate(seqnames = substr(seqnames, start=4, 20)) %>%
    slice(2) %>%
    unlist()
df <-
    data.frame(
        rnaseq = pseudo_tmm[which(rownames(pseudo_tmm) == gene_choice),],
        atacseq = pseudo_loess[which(rownames(pseudo_loess) == peak_choice),])
df <- cbind(df, do.call("rbind", strsplit(rownames(df), split = "_")))
colnames(df)[3:4] <- c("time","rep")
ggplot(df, aes(x=rnaseq, y = atacseq)) +
    geom_smooth(method = "lm", color = "black") +
    geom_point(aes(color=time)) +
    ggtitle(gene_choice)
cor(df$rnaseq, df$atacseq, method="spearman")


## Bonus: GO enrichment --------------------------------------------------------
#BiocManager::install(c("enrichplot","org.Dr.eg.db","clusterProfiler"))

#BiocManager::install(c("org.Dr.eg.db","clusterProfiler"))

library(enrichplot)
library(org.Dr.eg.db)
library(clusterProfiler)
go <- enrichGO(de_genes, OrgDb = "org.Dr.eg.db", ont = "BP", keyType="ENSEMBL")
dotplot(go)
#Voire sth gene surexpression: enrichissenemnt；
#  在你的差异基因（或差异开放峰关联基因）中，哪些生物学功能（GO terms）被显著富集。
#换句话说：它告诉你：这些基因更集中地出现在哪些生物过程（Biological Process, BP）里。
# test exact fisher,loi hyper什么什么