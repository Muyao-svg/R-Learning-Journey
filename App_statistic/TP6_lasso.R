## TP 6
data = read.table("wdbc.data", sep = ",")
data=data[,-1] # 去掉第一列
dim(data)   # 569 31
names(data)[1]="Y"
summary(data)
# Que variable continue dans ce probleme
# 要把良性和恶性转化为factor 0，1
data[Y] = ifelse(data[1]=="M",1,0)
data$Y = as.numeric(data$Y)
data$Y=factor(data$Y)
str(data)

set.seed(9170)
test = sample(1:length(data$Y),200)
train = -test
test=data[test,]
train=data[train,]

x=model.matrix(Y~.,train)[,-1]
dim(x)
y=train$Y
xtest=model.matrix(Y~.,test)[,-1]
ytest=test$Y

library(glmnet)
fit.lasso=glmnet(x,y,family = "binomial",alpha = 1) # alpha = 1 表示 lasso penalty
names(fit.lasso)
summary(fit.lasso)
head(coef(fit.lasso))



# 用交叉验证自动选最优 λ
cv.out = cv.glmnet(x,y,family="binomial") # 会产生erreur de generalisation de lambda estimee par validation croisee
plot(cv.out)
# 选虚线内选最小的(和选cp的原理一模一样)。因为这里是 -log(lambda),按照1-SE原则，modele plus parcimonieux (moins variable)所以选靠近-4那条线的
cv.out$lambda.min  #  0.00171973
cv.out$lambda.1se  # 0.004360142  更保守的λ（误差在最小值1个标准差内）

fit.lasso1=glmnet(x,y,family = "binomial",alpha = 1,lambda = 0.00171973)  
pred_lasso1 = predict(fit.lasso1,xtest,type="link")  
## link:combinaision lineaire, 会返回x*beta_hat
## reponse: Prob(Y=1)
## class: Y_hat (0 or 1)

pred_lasso2 = predict(fit.lasso1,s=cv.out$lambda.min ,type="coefficients") 
pred_lasso2  # logit(P(Y=1))= -47.1607839 + 9.3812741*V9 + ....  
## 列出在给定 λ 下，哪些变量（特征）系数 ≠ 0, type = "nonzero"：告诉 predict() 返回非零系数的变量索引
predict(fit.lasso1,s=cv.out$lambda.min ,type="nonzero")
# 10个variables使用了. lasso在30个参数中选择了10个作为modele的参数


pred_lasso3 = predict(fit.lasso1,s=cv.out$lambda.min ,type="class",newx=xtest) 
table(pred_lasso3)
table(ytest)
# Taux erreur
mean(pred_lasso3 != ytest)  #0.04


########################################
##### Comparer avec Random forest
#######################################
library(nnet)
library(randomForest)
fit.rf = randomForest(Y~., data=train,ntree=2000)
fit.rf
# OOB estimate of  error rate: 3.79%

# Calibration de mtry et ntree
plot(1:fit.rf$ntree, fit.rf$err.rate[,1],type="l")
fit.rf1 = randomForest(Y~., data=train,ntree=2000, mtry = 5)
lines(1:fit.rf1$ntree, fit.rf1$err.rate[,1],col="red")
fit.rf2 = randomForest(Y~., data=train,ntree=2000, mtry = 4)
lines(1:fit.rf2$ntree, fit.rf2$err.rate[,1],col="blue")

