# Q1
#Le problème d’apprentissage supervisé associé à ce jeu de données est un problème de 
#classification binaire.
#Il s’agit de prédire la variable Diag (diagnostic positif ou négatif du diabète) à 
#partir des autres variables physiologiques telles que glucose, pressure, insulin, etc.

# Q2 OOB
# model$err.rate 用来显示error OOB的数值
# L’erreur OOB est un estimateur de l’erreur de généralisation du modèle，
# 也就是说它估计的是模型在从未见过的新数据上的平均错误率（test error）。
# 换句话说：
# Elle fournit une estimation de la performance prédictive du modèle sur des données indépendantes.
# L’erreur OOB (Out-Of-Bag) correspond à l’erreur moyenne commise sur les observations non utilisées 
# pour construire chaque arbre de la forêt.
# Elle est un estimateur de l’erreur de généralisation du modèle, c’est-à-dire de 
# la performance prédictive sur de nouvelles données.
# Elle est considérée comme un estimateur sans biais car, pour chaque observation, la prédiction OOB 
# est réalisée à partir d’arbres qui n’ont pas été entraînés sur cette observation. 
# Ainsi, l’évaluation est effectuée sur des données indépendantes du processus d’apprentissage, 
# comme dans un vrai test externe.
# 每个点的OOB预测都是“像测试集那样独立”的预测结果，因此OOB误差不依赖于训练集内部信息，是无偏的估计量
#（estimateur sans biais）de l’erreur de test.



# Q3 比较OOB error and Bagging
# Bagging 是随机森林的特例（当 mtry = p 时，即每次分裂节点时使用全部变量）
# En recalculant l’erreur OOB à l’aide de la fonction predict, on obtient une valeur proche de celle
# fournie par le modèle.
# En comparaison avec la méthode de Bagging (où mtry = p), la forêt aléatoire présente une erreur 
# OOB plus faible.
# Cela montre que la sélection aléatoire d’un sous-ensemble de variables à chaque division réduit 
# la corrélation entre arbres et améliore la capacité de généralisation du modèle.
# Conclusion : la forêt aléatoire surpasse le Bagging en termes de performance prédictive.