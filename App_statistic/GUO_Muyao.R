# GUO Muyao
college = read.table("College.dat",header=TRUE)
dim(college)
str(college)
college$Private=factor(college$Private)

# Enlever la 1er colonne qui n'est pas une covariable continue
college = college[,-1]
head(college)
set.seed(1907)
# 1
# (a)
B=1000
dim(college_boot)
n=sample(1:nrow(college), B, replace = TRUE)
college_boot = college[n,]
dim(college_boot)
# 1000 18
head(college_boot)

# (b)
# On estime par regression lineaire 
# lm sur jeu de donnee College. Parce que c'est regression lineaire.
reg1 = lm(Apps~.,data=college)
summary(reg1)
# New data
xnew= college[115,-1]
pred1 = predict(reg1, newdata = xnew)
pred1    # 2918.669
# lm sur jeu de donnee College_boostrap
reg2 = lm(Apps~.,data=college_boot)
summary(reg2)
pred2 =  predict(reg2, newdata = xnew)
pred2   # 2688.634


## (b) T_boot(Predire Apps sur newdata)
B <- 1000
T_boot <- numeric(B)
n0 <- nrow(college)

for (b in 1:B) {
  idx_b <- sample.int(n0, n0, replace = TRUE)      
  dat_b <- college[idx_b, ]
  fit_b <- lm(Apps ~ ., data = dat_b)
  T_boot[b] <- predict(fit_b, newdata = xnew)     
}
hist(T_boot)
# ici l'histogramme présente une allure gaussienne symétrique
boot_mean <- mean(T_boot)   # 2883.307
boot_se   <- sd(T_boot)   # 244.7917


# (c)
boot_ci   <- quantile(T_boot, c(0.025, 0.975))
#     2.5%    97.5% 
# 2420.497   3367.741
  
# 2. Comparaison des perfprmance
# separer les donnee test et train
n = 330
test = sample(1:nrow(college_boot),n)
train = -test
train=college_boot[train,]
test = college_boot[test,]

# (a) Lasso
## -------- LASSO（cv.glmnet） --------
library(glmnet)

## model.matrix 
xtrain <- model.matrix(Apps ~ ., data = train)[, -1]
ytrain <- train$Apps
xtest  <- model.matrix(Apps ~ ., data = test)[, -1]
ytest  <- test$Apps

## Validation croisee pour choisir lambda, family="gaussian"
set.seed(1907)
cv.out <- cv.glmnet(xtrain, ytrain, alpha = 1, family = "gaussian")
plot(cv.out)  
# On choisir le lambda d la valeur minimum 
## Prediction avec lambda_min
pred_lasso <- predict(cv.out, newx = xtest, s = "lambda.min")
rmse_lasso <- sqrt(mean((pred_lasso - ytest)^2))
# 1028.076
predict(cv.out,type="coefficients",s="lambda.min")
length(predict(cv.out,type="nonzero",s="lambda.min")[,1])
# 13 covariable dans le modele ajuste final.


# (b) random forest
library(randomForest)
set.seed(1907)
fit.rf  <- randomForest(Apps ~ ., data = train, ntree = 2000)
fit.rf

## Presente courbe mse
plot(1:fit.rf$ntree, fit.rf$mse, type = "l", xlab = "Number of trees", ylab = "MSE")

## Essayer mtry different. Par defaut mtry = p/3
p <- ncol(train) - 1 
set.seed(1907)
fit.rf2 <- randomForest(Apps ~ ., data = train, ntree = 2000, mtry = 5)
fit.rf3 <- randomForest(Apps ~ ., data = train, ntree = 2000, mtry = 8)

lines(1:fit.rf2$ntree, fit.rf2$mse,col="red") 
lines(1:fit.rf3$ntree, fit.rf3$mse,col="blue")

## D'apred le plot, fit.rf present la meilleur MSE. Donc on choit le mtry par defaut
## Prediction
pred_rf   <- predict(fit.rf, newdata = test)
rmse_rf   <- sqrt(mean((pred_rf - ytest)^2))
# 932.2591

# Regle apprentissage bagging
fit_bag <- randomForest(Apps ~ ., data = train,ntree = 2000,mtry  = p)
pred_rf_bag   <- predict(fit_bag, newdata = test)
rmse_rf   <- sqrt(mean((pred_rf_bag - ytest)^2))

# (c) Comparaison avec validation croisee 10 blocs (LASSO,RF,bagging)
K=10
set.seed(0719)
ind_fold=sample(1:K,nrow(college),replace=TRUE)

error_reg=numeric()
error_rf=numeric()
error_lasso=numeric()

mse=function(y,ypred) round((mean((y-ypred)**2)),digits=2)
rmse=function(y,ypred) sqrt(mse(y,ypred))

for (j in 1:K)
{
  fit.lm=lm(Apps~.,data=college[ind_fold!=j,])
  y.lm=predict(fit.lm,newdata=college[ind_fold==j,])
  y.test=college[ind_fold==j,"Apps"]
  error_reg[j]=rmse(y.test,y.lm)
  
  # lasso
  xtrain <- model.matrix(Apps ~ ., data = college[ind_fold!=j,])[, -1]
  ytrain <- train$Apps
  xtest  <- model.matrix(Apps ~ ., data = college[ind_fold==j,])[, -1]
  ytest  <- test$Apps
  cv.out <- cv.glmnet(xtrain, ytrain, alpha = 1, family = "gaussian")
  pred_lasso <- predict(cv.out, newx = xtest, s = "lambda.min")
  error_lasso[j]=rmse(y.test,y.lm)
}
cv.error_reg=mean(error_reg)
cv.error_reg  # 1071.305



## =========================
## 10-fold CV: lm / LASSO / RF / Bagging
## =========================
set.seed(0719)
K <- 10
n <- nrow(college)
p <- ncol(college) 

fold_id <- sample(rep(1:K, length.out = n))  # 每个样本一个折号，1..K

mse  <- function(y, yhat) mean((y - yhat)^2)
rmse <- function(y, yhat) sqrt(mse(y, yhat))

## Resultat
err_lm     <- numeric(K)
err_lasso  <- numeric(K)
err_rf     <- numeric(K)
err_bag    <- numeric(K)

for (j in 1:K) {

  tr_idx <- fold_id != j
  te_idx <- !tr_idx
  
  train_dat <- college[tr_idx, ]
  test_dat  <- college[te_idx, ]
  
  y_test <- test_dat$Apps
  
  ## ============ 1) Regression lm ============
  fit_lm <- lm(Apps ~ ., data = train_dat)
  pred_lm <- predict(fit_lm, newdata = test_dat)
  err_lm[j] <- rmse(y_test, pred_lm)
  
  ## ============ 2) LASSO (glmnet) ============
  ## 
  x_tr <- model.matrix(Apps ~ ., data = train_dat)[, -1, drop = FALSE]
  y_tr <- train_dat$Apps
  x_te <- model.matrix(Apps ~ ., data = test_dat)[, -1, drop = FALSE]
  
  cv_out <- cv.glmnet(
    x_tr, y_tr,
    alpha = 1,  
    family = "gaussian"
  )
  pred_lasso <- predict(cv_out, newx = x_te, s = "lambda.min")
  err_lasso[j] <- rmse(y_test, as.numeric(pred_lasso))
  
  ## ============ 3) Random Forest ============

  fit_rf <- randomForest(
    Apps ~ ., data = train_dat,
    ntree = 1000
  )
  pred_rf <- predict(fit_rf, newdata = test_dat)
  err_rf[j] <- rmse(y_test, pred_rf)
  
  ## ============ 4) Bagging ============
  fit_bag <- randomForest(
    Apps ~ ., data = train_dat,
    ntree = 1000,
    mtry  = p
  )
  pred_bag <- predict(fit_bag, newdata = test_dat)
  err_bag[j] <- rmse(y_test, pred_bag)
}

res <- data.frame(
  Fold = 1:K,
  RMSE_lm = round(err_lm, 3),
  RMSE_lasso = round(err_lasso, 3),
  RMSE_rf = round(err_rf, 3),
  RMSE_bag = round(err_bag, 3)
)
print(res)





# (d) Conclusion
# rmse_lasso: 1028.076
# rmse_rf : 932.2591

#     Fold  RMSE_lm RMSE_lasso  RMSE_rf RMSE_bag
#1     1  917.343    881.424  872.098  710.304
#2     2 1422.435   1417.346 1236.499 1238.126
#3     3 1077.933   1068.962  896.670  820.254
#4     4  813.859    801.881  615.146  590.816
#5     5 1904.661   1975.950 3819.038 3741.438
#6     6 1084.482   1080.416  828.756  808.705
#7     7  755.756    765.594  779.844  805.675
#8     8  864.614    860.867  734.192  882.279
#9     9 1074.651   1083.428 1194.812 1149.988
#10   10  947.627    929.236  949.350  923.356

mean(res$RMSE_lm)  # 1086.336
mean(res$RMSE_lasso)  # 1086.51
mean(res$RMSE_rf)  # 1192.64
mean(res$RMSE_bag)  # 1167.094

# D'apres les resultats k fold, lm et lasso ont la meilleur performance






