# Muyao Guo

rm(list=objects())     # supprime les objets existant en session 
graphics.off()         # supprime les graphiques existant en session

## ------------------------------------------------------------
## 1) Données et statistiques d’intérêt
##    - Contexte bioéquivalence : chaque sujet reçoit placebo, old patch, new patch
##    - Donnée appariée : les trois mesures d’un même sujet sont corrélées
## ------------------------------------------------------------
# Create data 
df = data.frame(
  "placebo" = c(9243,9671,11792,13357,9055,6290,12412,18806),
  "old" = c(17649,12013,19979,21816,13850,9806,17208,29044),
  "new" = c(16449,14614,17274,23798,12560,10157,16570,26325),
  "old_placebo" = c(8406,2342,8187,8459,4795,3516,4796,10238),
  "new_old" = c(-1200,2601,-2705,1982,-1290,351,-638,-2719)
)
size = nrow(df)

## Estimateur plug-in (theta chapeau)
origine_theta = mean(df$new_old)/mean(df$old_placebo)
origine_theta   # -0.0713061


## ------------------------------------------------------------
## 2) Bootstrap non paramétrique
##    Méthodo : ré-échantillonnage PAR LIGNES (appariement conservé)
##    — à chaque réplication b : tirer des indices i1..in avec remise
##      puis calculer theta* = mean(Y[idx]) / mean(Z[idx])
## ------------------------------------------------------------
# Bootstrapt simulation
set.seed(1907)
B = 1000
bootstrapt_theta = replicate(B, {
  idx <- sample.int(size, size, replace = TRUE)      # ré-échantillonner les SUJETS (lignes)
  mean(df$new_old[idx]) / mean(df$old_placebo[idx])
})
length(bootstrapt_theta)   # 1000


## ------------------------------------------------------------
## 3) Écart-type bootstrap et figure de la loi d’échantillonnage
## ------------------------------------------------------------ 
# Car le vrai theta est inconnu, on utilise "origine_theta" a la place de "vrai_theta" et "bootstrapt_theta"
# a la place de "origine_theta"
# Autrement dire, "T_boot = theta_hat - theta", mais theta est inconnu. Donc on utilise "T_boot = bootstrapt_theta - origine_theta"
# Loi d'echantillonage bootstrap centrée
T_boot = bootstrapt_theta - origine_theta 
hist(bootstrapt_theta,prob=TRUE, breaks = 25)
hist(T_boot, prob=TRUE, breaks = 25)
se = sd(T_boot)
se   # 0.09825901

## ------------------------------------------------------------
## 4) Intervalles de confiance bootstrap
##    Deux méthodes usuelles :
##    (a) Percentile : [q_0.025(theta*), q_0.975(theta*)]
##    (b) Basic (centré) : T* = theta* - theta_hat ; IC = [theta_hat - q_0.975(T*),
##                                                         theta_hat - q_0.025(T*)]
##    Remarque : basic et percentile sont proches ici, mais pas identiques.
## ------------------------------------------------------------
# Pour repondre a la question bioequivalence, on compute IC
# IC basic
IC_b = data.frame(binf = as.numeric(origine_theta-quantile(T_boot, 0.975)),
                bsup = as.numeric(origine_theta - quantile(T_boot, 0.025)))

IC_b
#   binf      bsup
#1 -0.3064328 0.08632264
# Interprétation :
# 95% intervalle bootstrap basic ≈ (-0.31, 0.09).
# Cet intervalle n'est pas entièrement inclus dans l'intervalle d'équivalence (-0.20, 0.20).
# Conclusion : on ne peut pas conclure à une bioéquivalence au niveau 95%.

# IC pourcentage
CI_perc <- quantile(bootstrapt_theta, c(.025, .975))
CI_perc
#   2.5%      97.5% 
#  -0.2289348  0.1638206
# Interprétation :
# 95% intervalle bootstrap percentile ≈ (-0.23, 0.16).
# Cet intervalle n'est pas entièrement inclus dans l'intervalle d'équivalence (-0.20, 0.20).
# Conclusion : on ne peut pas conclure à une bioéquivalence au niveau 95%.

