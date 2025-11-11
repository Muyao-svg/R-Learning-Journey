# Hypothese: 不同的echantillo之间基因的表达量有没有显著差别？
# https://www.cnblogs.com/swindler/p/16909670.html
# Condition A: tissu cancereux; Condition B: tissu normal
# H0: mu_A = Mu_B;
# Donnees EMT
data=read.table("TP1 DESeq/countsEMT.txt")
boxplot(data)
boxplot(log(data+1))
barplot(colSums(data),main = "Row data")
# 去除0表达量数据
data = data[-which(rowSums(data)==0),]
# 或者去除表达量<10的数据
data = data[-which(rowSums(data)<=10),]
# 验证是否符合Poisson 分布
vectMean = apply(data,1, mean)
vectVar = apply(data,1,var)
plot(log(vectMean+1), log(vectVar+1))
abline(0,1,col="red")
# E(Yij) != var(Eij) 因此数据不符合泊松分布
# Nornalisation des donnees
# CPM: depend la profondeur de sequancage
library(edgeR)
cpm(data)
cpm = data/colSums(data)*10^6  # 每个echantillons是一列
boxplot(cpm,main="cpm")
boxplot(log(cpm+1), main="log cpm")
# 但这样做不太好，忽略了基因的长度

# La normalisation Upper quantile
# quantile在75%作为reference. 只用每个样本中 非零 counts 的第 75 百分位数 
# 来代表这个样本的测序深度。
calcNormFactors(data)

# Relative Log Expression


# TMM: Trimmed mean of m-value


# 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
help(DESeq)
# Condition day0,day7
cond = factor(c(rep("day0",3),rep("day7",3)))
# object construction
dds <- DESeqDataSetFromMatrix(data, DataFrame(cond), ~ cond)
# standard analysis
dds <- DESeq(dds)
res <- results(dds)
variable.names(res)  #有一个是pvalue
# 在 DESeq2 中，每个基因都会被检验一个 统计假设：
# 原假设（H₀）：
# 该基因在两组条件下的平均表达量相同。
# 即 log₂FoldChange = 0。
# 备择假设（H₁）：
# 该基因在两组间表达不同（log₂FoldChange ≠ 0）。
# DESeq2 对每个基因进行 Wald 检验（Wald test），
# 计算出一个统计量 stat = log2FoldChange / lfcSE，
# 然后根据标准正态分布求出对应的 p-value
table(res$pvalue<0.05)
# FALSE  TRUE 
# 9445 12870
table(res$padj<0.05)  # BH调整后的pvalue,代表通过了 FDR（假发现率）控制: 如果但看pvalue,在数据里面他是uniform 分布的，有5%的可能性会出现假阳性
# FALSE  TRUE 
# 8665 12363 
table(res$pvalue<0.05/dim(data)[1])
# 数量总和并不一样，说明有paval null
table(complete.cases(res$pvalue))
cnts = data[!table(complete.cases(res$pvalue)),]

#!! 在之前就要先去除表达量为0的数据

# moderated log2 fold changes
# Yij = log2(Yij+1)
# 在之前的mean-var图里面，虽然他们不相等，但是有一定的relation，trouver la fonction de la variance depandant de la moyenne
# stabiliser la variance 在均值范围内具有恒定方差）的数值矩阵。该转换同时根据样本量进行标准化处理。对数转换（rlog）对样本量因子敏感度较低，当样本量因子波动较大时尤为适用
help("varianceStabilizingTransformation")

# CLustering 
vsd <- vst(dds,nsub = 500)  # variance stabilizing transformation (VST)
plotPCA(vsd)

help(plotPCA)




# Standard analysis
help("DESeq")
estimateSizeFactors(dds)
dds$sizeFactor

dds <- estimateSizeFactors(dds)
estimateDispersions(dds)    # alpha_i. dispersion离散趋势
# dispersion 是一个可以将均值和var联系起来的统计量。 alpha越大，var越大，越远离loi Boisson
# 是否tres variance 
# dispersion 将一组数据的标准差除以其均值，用来测度数据离散程度的相对数
plotDispEsts(dds)
# 数据量越大，dispersion越小
# 粉色：estimation dispersion ajuster; bleu: estimation dispersion final
which.max(dispersions(dds))  # 4893

# Wald检验是一种在统计模型中测试特定解释变量的重要性的方法。
# Wald检验是测试与一组解释变量相关的参数是否为零的多种方法之一
# wald检验用于检测log2 fold changes: 表示两样品（组）间表达量的比值，对其取以2为底的对数之后即为log2FC。 一般默认取log2FC绝对值大于1为差异基因的筛选标准
# 以及用于计算pvalue
res["CDH1",]
hist(res$pvalue)






