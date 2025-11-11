###################### MAUREL DYLANE TCHUISSE II ###############################


################################################################################

########################### DEVOIR À RENDRE ####################################

################################################################################



#### OBJECTIF : Estimation d'un paramètre de bioéquivalence et de sa 
###variabilité (par bootstrap)



# 1) Implémentation du jeu de données des patients traités


n = 8 #nombre de patients

subject = 1:n #numéro des patients

# Concentration hormonale après prise du placebo
placebo = c(9243,9671,11792,13357,9055,6290,12412,18806)
# Concentration hormonale après prise de l'ancien traitement
old = c(17649,12013,19979,21816,13850,9806,17208,29044)
# Concentration hormonale après prise du nouveau traitement
new = c(16449,14614,17274,23798,12560,10157,16570,26325)

donnees_patients = data.frame(subject, placebo, old, new, 
                              old_minus_placebo = old - placebo, 
                              new_minus_old = new - old)

donnees_patients
# subject placebo   old   new old_minus_placebo new_minus_old
#       1    9243 17649 16449              8406         -1200
#       2    9671 12013 14614              2342          2601
#       3   11792 19979 17274              8187         -2705
#       4   13357 21816 23798              8459          1982
#       5    9055 13850 12560              4795         -1290
#       6    6290  9806 10157              3516           351
#       7   12412 17208 16570              4796          -638
#       8   18806 29044 26325             10238         -2719


# 2) Estimations

## 2-1) Paramètre de bioéquivalence theta

# Étant donné : 

### Z: la variable aléatoire de loi inconnue, désignant la différence, 
### pour chaque sujet, entre la concentration hormonale après administration
### de l'ancien traitement (old) et celle après prise du placebo (placebo), 
### soit la colonne "old_minus_placebo" dans le jeu de données précédent;
Z = donnees_patients$old_minus_placebo

### Y: la variable aléatoire de loi inconnue, désignant la différence, 
###pour chaque sujet, entre la concentration hormonale après administration
### du nouveau traitement (new) et celle après prise de l'ancien traitement(old),
###  soit la colonne "new_minus_old" dans le jeu de données précédent;
Y = donnees_patients$new_minus_old

# La bioéquivalence theta est définie comme étant le rapport entre l'espérance de
# Z et celle de Y (theta = E(Y)/E(Z)).  Ces deux espérances étant inconnues,
# theta peut être estimé par la méthode plug-in, qui consiste dans ce cas à remplacer
# le rapport des espérances par un rapport de moyennes empiriques.

# Estimateur plug-in de theta
theta_estim = mean(Y)/mean(Z)
theta_estim
#-0.0713061

## 2-2) Estimation par bootstrap de l'écart-type de l'estimateur theta_estim

# Une appréciation de la variabilité (calcul de l'écart-type par exemple) 
# de l'estimateur theta_estim du paramètre theta, nécessite que l'on ait accès à 
# la loi de cet estimateur (qui elle même dépend de la loi des variables aléa-
# toires Y et Z dont est fonction l'estimateur), c'est-à-dire accès à une distri- 
# bution de l'ensemble de ses valeurs possibles. 

# Ceci suppose dans ce cas précis d'étude, de disposer de plusieurs 
# n-échantillons du type [(y1,z1), (y2,z2), ...., (yn,zn)] à partir desquels on 
# calculera, pour chaque échantillon, une valeur de l'estimateur. Ces valeurs, 
# une fois obtenues pourront servir à représenter la distribution (et donc la loi)
# de l'estimateur et à calculer des paramètres de dispersion (écart-type) qui 
# renseignent sur sa variabilité.

# Cependant, il est en pratique impossible de disposer d'une multitude de n-échan-
# tillons pour le même couple de variables aléatoires (Y,Z) et qui pourraient être 
# exploités pour le calcul des paramètres de variabilité de l'estimateur.

# C'est une problématique à laquelle la méthode bootstrap non-paramétrique répond. 
# Elle permet d'estimer la distribution de l'estimateur theta_estim ainsi que ses 
# paramètres de variabilité :  

# ======> On génère par ré-échantillonnage d'un n-échantillon de départ 
#         [(y1,z1), (y2,z2), ...., (yn,zn)] de loi F inconnue, plusieurs 
#         n-échantillons [(y1,z1)*, (y2,z2)*, ...., (yn,zn)*], appelés
#         échantillons bootstrap; le tirage des éléments du n-échantillon 
#         [(y1,z1), (y2,z2), ...., (yn,zn)] étant réalisé selon une loi discrète
#         empirique de fonction de répartition F_hat, qui représente une estimation 
#         de la fonction de répartition F associée à la loi de 
#         [(y1,z1), (y2,z2), ...., (yn,zn)]. 

# ======> Pour chacun des échantillons bootstrap, on calcule une estimation 
#         "theta_estim_bootstrap" du paramètre theta. 

#           - La distribution des estimations "theta_estim_bootstrap" représente  
#         la distribution d'échantillonnage bootstrap de l'estimateur theta_estim;
#       
#           - L'estimation de l'écart-type de theta_estim (estimateur du paramètre  
#         de bioéquivalence theta), est donné par la racine carée de la variance
#         empirique de la précédente distribution d'échantillonnage, c'est-à-dire
#         la racine carée de la variance des estimations theta_estim_bootstrap de   
#         theta, calculées sur chacun des échantillons bootstrap.




#Graine aléatoire fixée pour la reproductibilité des résultats
set.seed(23) 

#Nombre d'échantillons bootstrap à générer 
B = 1000

# Génération des échantillons bootstrap

# L'idée est que chaque échantillon bootstrap contienne les numéros des patients
# tirés selon une loi discrète uniforme (tirage avec remise) à partir de 
# l'ensemble des numéros de départ {1,2,3,4,5,6,7,8}.
#On échantillonne selon les indices ici afin de tenir de la dépendance
#entre les variables Y et Z prises sur un même individu. 

echantillons_bootstrap_patients = replicate(B, sample(1:n, n, replace = TRUE))

echantillons_bootstrap_patients #Matrice de 1000 échantillons de numéros de patients, chaque 
#échantillon étant disposé en colonne.

# Fonction calculant une valeur de theta_estim à partir du vecteur des numéros des 
#patients formant un échantillon
theta_hat_star = function(num_patient) mean(Y[num_patient]) / mean(Z[num_patient])

# Calcul des estimations de theta pour chacun des échantillons bootstrap 
theta_estim_bootstrap = apply(echantillons_bootstrap_patients, 2, theta_hat_star)
theta_estim_bootstrap

#Loi d'échantillonnage bootstrap de l'estimateur 
hist(theta_estim_bootstrap, breaks = 30, xlab = "Estimations de theta sur les
     échantillons bootstrap", ylab = "fréquences absolues", 
     main= "Distribution d'échantillonnage boostrap \n de l'estimateur theta_estim
     de theta")
######### On obtient pour cette simulation, une distribution d'échantillonnage
######### asymétrique, que l'on aurait pas pu capturer si l'on avait fait 
######### l'hypothèse d'une distribution normale

# Calcul de l'estimation bootstrap de l'écart-type de l'estimateur theta_estim
std_err_bootstrap_estimator = sqrt((B-1)*var(theta_estim_bootstrap)/B)
std_err_bootstrap_estimator
#0.1009134 valeur proche de celle trouvée par Wasserman. La différence étant liée
# aux différences d'échantillonnage. La graine aléatoire assurant la reproductibilité
# n'étant pas forcément la même que celle prise par Wasserman


# 3) Calcul de l'intervalle de confiance de niveau 1-alpha = 95% du paramètre de bio-
# équivalence theta. 

########################## Méthode pivotale #################################### 

# Soit "theta_estim" l'estimateur de "theta", considérons la statistique 
# (theta_estim - theta). 

# Si q(alpha/2) et q(1 - alpha/2) sont respectivement les quantiles d'ordre 
# alpha/2 et 1-alpha/2  de la statistique (theta_estim - theta) alors 
# par définition : 

#                  P[(theta_estim - theta) <= q(alpha/2)] = alpha/2 
#                                     et 
#                 P[(theta_estim - theta) >= q(1-alpha/2)] = alpha/2


# On peut donc écrire que : 

#       P[q(alpha/2) <= (theta_estim - theta) <= q(1 - alpha/2)] = 1 - alpha

# On en déduit que : 

#      P[q(alpha/2) <= (theta_estim - theta) <= q(1 - alpha/2)] = 1 - alpha

# ==> P[theta_estim - q(1 - alpha/2) <= theta <= theta_estim - q(alpha/2)] = 1 - alpha (i)

# Or, la loi de (theta_estim - theta) étant inconnue, il est impossible d'accéder
# aux valeurs de q(1 - alpha/2) et q(alpha/2). 

# On approche donc la loi de (theta_estim - theta) par celle de 
# (theta_estim_bootstrap - theta_estim) où "theta_estim_bootstrap" représente 
# la variable des estimations de theta sur les échantillons bootstrap. 

# On a donc que :

#                     q(alpha/2) = a(alpha/2) - theta_estim  (ii)
#                                    et 
#                 q(1 - alpha/2) = a(1 - alpha/2) - theta_estim (iii)

# avec a(alpha/2) et a(1 - alpha/2), les quantiles respectivement d'ordre alpha/2
# et 1-alpha/2 de la loi de theta_estim_bootstrap

# (ii) et (iii) dans (i) permettent d'écrire que : 

# P[2*theta_estim - a(1 - alpha/2) <= theta <= 2*theta_estim - a(alpha/2)] = 1 - alpha

# D'où IC(theta) = [2*theta_estim - a(1 - alpha/2), 2*theta_estim - a(alpha/2)]
# est un intervalle de confiance de niveau 1-alpha de theta


########################### Application numérique ##############################
alpha = 0.05 

# Borne inférieure de l'intervalle de confiance
ICinf = as.numeric(2*theta_estim - quantile(theta_estim_bootstrap,probs=1-alpha/2))

# Borne supérieure de l'intervalle de confiance
ICsup = as.numeric(2*theta_estim - quantile(theta_estim_bootstrap,probs=alpha/2))

# Intervalle de confiance
IC = c(ICinf, ICsup)
IC
#-0.30082947  0.09137424


# Interprétation : La bioéquivalence est admise, si et seulement si 
# -0.2 <= theta <= 0.2. Or notre intervalle de confiance de niveau 95% est 
# [-0.30, 0.091] et n'est pas entièrement inclus dans [-0.2, 0.2]. On n'a donc pas
# démontré la bioéquivalence au seuil 95%. 

# Il se pourrait donc qu'il y ait une différence entre les effets du nouveau trai- 
# tement hormonal et ceux de l'ancien traitement. 