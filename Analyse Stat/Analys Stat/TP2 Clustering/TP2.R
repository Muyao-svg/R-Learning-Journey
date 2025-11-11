# Data input
data = read.table("TP2 Clustering/GSE94396_iTreg_raw_counts.txt",header = TRUE)
summary(data)
boxplot(log(data[,-1]+1),main="Boxplot of row data")
rownames(data)=data[,1]
data = data[,-1]

# BiocManager::install("coseq")
??coseq
library(coseq)
set.seed(9170)
data = data[-which(rowSums(data)<=10),]
## === Clusting probabilist ===
# We assume the data follows GMM distribution.
runGMM <- coseq(data,K=2:13,model="Normal", transformation="voom", GaussianModel = "Gaussian_pk_Lk_I",set.seed=9170)
summary(runGMM)  # Selected number of clusters via ICL: 13
#Why GaussianModel = "Gaussian_pk_Lk_I" :
#Occasionally, the default form of Gaussian mixture model (Gaussian_pk_Lk_CK,
#which is the most general form available with per-cluster proportions, volumes, 
#shapes, and orientations for covariance matrices) used in coseq is not estimable 
#(via the Rmixmod package used by coseq) due to non-invertible per-cluster covariance
#matrices. The error message thus suggests using a slightly more restrictive form of 
#Gaussian mixture model, such as the Gaussian_pk_Lk_Bk (which imposes spherical covariance matrices) 
#or Gaussian_pk_Lk_I (which imposes diagonal covariance matrices). See ?Rmixmod::mixmodGaussianModel 
#for more details about the nomenclature and possible forms of Gaussian mixture models.

plot(runGMM)$ICL
# Le curve continue de baisser mais maisn a un rythm moins soutenu, ce qui indique que le "gain" lie a
# l'ajoute de clusters sont limetee.(avec une risque de overfitting)


## === Clusting non probabilist:kmeans ===
runKmeans <- kmeans(data, centers = 10)
cl_KM <- runKmeans$cluster

plot_profiles_km <- function(data, cl, main="k-means profiles (K=10)") {
  data <- as.matrix(data)
  rng <- range(data, na.rm = TRUE)
  K <- length(unique(cl))
  op <- par(mfrow = c(ceiling(K/3), 3), mar=c(3,3,2,1))
  on.exit(par(op), add = TRUE)
  
  for (k in sort(unique(cl))) {
    mat <- data[cl == k, , drop = FALSE]
    mu  <- colMeans(mat, na.rm = TRUE)

    graphics::plot(1:ncol(data), mu, type="n", ylim=rng,
                   xlab="", ylab="", main=paste("Cluster", k))
    apply(mat, 1, function(z) graphics::lines(z, col=adjustcolor("grey", 0.4)))
    graphics::lines(mu, lwd=3)
  }
  mtext(main, outer=TRUE, line=-1.5, cex=1.1)
}

plot_profiles_km(data, cl_KM, main="k-means profiles (K=10)")

## === Comparaison 2 methodes===
# 1.Time cost
#The Gaussian mixture model is more time-intensive but enables the estimation of
#per-cluster correlation structures among samples. The K-means algorithm has a 
#much smaller computational cost, at the expense of assuming a diagonal per-cluster 
#covariance structures. Generally speaking, we feel that the K-means approach is a 
#good first approach to use, although this can be somewhat context-dependent. Finally, 
#although we feel the Poisson mixture model was a good first approach for RNA-seq co-expression, 
#we now recommend either the Gaussian mixture model or K-means approach instead.

# 2.Stability of result
plot(runGMM)$profile
plot_profiles_km(data, cl_KM, main="k-means profiles (K=10)")
# GMM a une profile plus flexible et stable.


# 3. Sensibility de la methode a diverses parametres
#Modèle de mélange gaussien (GMM, méthode probabiliste):
#Le GMM est sensible au choix du nombre de composantes K et au type de covariance utilisé. 
#Ces paramètres déterminent la forme des clusters dans l’espace.
#Il peut aussi être influencé par le choix de la transformation des données
#(ex. voom, logCLR) qui affecte la normalité supposée des distributions.
#Comme le modèle est probabiliste, il est souvent moins sensible à l’initialisation :
#les probabilités d’appartenance permettent une adaptation progressive des clusters,
#contrairement au k-means qui fait une assignation “tout ou rien”.


# 4. Facilite interpretation
#Le k-means est plus simple et visuellement intuitif, idéal pour une première exploration,
#tandis que le modèle gaussien offre une interprétation plus complète mais nécessite une 
#lecture plus analytique.






