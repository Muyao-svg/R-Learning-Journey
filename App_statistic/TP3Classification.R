library(MASS)
data(birthwt)
?birthwt
head(birthwt)
dim(birthwt)

summary(birthwt)  # R导入的文件默认数字就是数字，要手动把numeric 转换成factor
str(birthwt)
unique(birthwt$ftv)   # variable discrete 6 modalite (5 variables)
unique(ui)    # variable discrete
summary(factor(birthwt$ftv))

## === 1. 重新构建data ===
birthwt = within(birthwt,{
  race = factor(race, labels=c("white", "black", "other"))
  smoke = factor(smoke,labels=c("No", "Yes"))
  ptl = factor(ptl > 0)
  ht= factor(ht>0)
  ui = factor(ui,labels=c("No", "Yes"))
  ftv = factor(ftv)
  levels(ftv)[-(1:2)] = "2+"
})

birthwt = birthwt[,1:ncol(birthwt)-1]
summary(birthwt)   # 现在factor 都是factor了
attach(birthwt)


## === Q1 ===
fit = glm(low~., data = birthwt, family = binomial)
summary(fit)
# Null deviance: 234.67  on 188  degrees of freedom   这里写的就是deviance。
# NUll: modele M0: logit_P(Y=1) = beta0; M1 :logit_P(Y=1) = beta0 + sum(beta*X)
# Trv=deviance(M0)-deviance(M1) = -2log(M0) + 2log(M1)
# Trv ~ X2(df = 两个模型的参数差-1)
# pvalue = P(Trv >= Trv_obs)
# Residual deviance: 195.48  on 178  degrees of freedom
# AIC: 217.48
pval = 1-pchisq(234.67-195.48, 10)
pval   # 2.351524e-05
pred1 = predict(fit, newdata = data.frame(age = 20, lwt=105,race="white", smoke="Yes", ptl="FALSE",ht="FALSE", ui="No", ftv="1"),type = "response")
pred1    # 0.2233787 

# 在fit结果里面，race 旁边有一颗星，我们想知道它是不是真的significative?
# Trv = Deviance(H0) - Deviance(H1)
# H0:modele sans race; H1:modele complete
# 注意这里df=2，因为race包含了“black”和“other”
fit.sans_race = glm(low~.-race,data = birthwt, family = binomial)
summary(fit.sans_race)
Trv = deviance(fit.sans_race) - deviance(fit)
Trv   # 5.751273
pval2 = 1-pchisq(Trv,2)
pval2   # 0.05638025   race n'est pas significative
anova(fit.sans_race, fit, test = "LRT")
# anova 算出来的pvalue的结果是一样的

# IC
newdata1 = data.frame(age = 20, lwt=105,race="white", smoke="Yes", ptl="FALSE",ht="FALSE", ui="No", ftv="1")
prob1 = pred1   # 0.2233787
x1 = c(1,20,105,0,0,1,0,0,0,1,0)
se = sqrt(t(x)%*%vcov(fit)%*%x)
se   # 0.5352285
IC = cbind(prob1-qnorm(0.975)*se, prob1+qnorm(0.975)*se)
IC = exp(IC)/(1+exp(IC))
IC
#0.3045657 0.7811546 包含了0.5，这不是一个好的预测

# 或者直接用
p=predict(fit, newdata = newdata1,type = "link", se.fit = TRUE)
p1 = p$fit
IC = cbind(p1-qnorm(0.975)*p$se.fit, p1+qnorm(0.975)*p$se.fit)
IC = exp(IC)/(1+exp(IC))
IC
# !!!! 注意这里一定要用 type = "link"
#想直接得到概率 ⇒ 用 type="response".
#想做标准误/置信区间或数值更稳定的计算 ⇒ 用 type="link" 拿到 𝜂
#η 和 se.fit，再 plogis() 映射回概率更稳妥。
# "link" 返回线性预测值,log(P)； “response”返回响应尺度上的E（Y），直接给概率P,[0,1]


newdata2 = data.frame(age = 20, lwt=105,race="white", smoke="Yes", ptl="FALSE",ht="FALSE", ui="No", ftv="1")




# === Q2 ===
##  Erreur test ===
# === 1.Decoupage train/test ===
test = sample(1:length(birthwt$low),60)
train = -test
train = birthwt[train, ]
test = birthwt[test,]

model_train = glm(low~.,data = train,family = binomial)
test_pred = predict(model_train, newdata = test, type = "response")
head(test_pred)
length(test_pred)
pred = (0,60)
pred = as.numeric(test_pred>0.5)
table(test$low, pred)
#pred
#   0  1
# 0 37  4
# 1 14  5

accuracy =  (37+5)/60
accuracy    # 0.7
mean(test$low != pred)  # 0.3



# Erreur de validation croisee
## === Validation coisee 交叉验证===
# 把数据多次拆分成“训练集 + 验证集”，轮流训练和测试，再取平均误差 → 更可靠地估计模型的泛化能力。
K=5
ind_fold=sample(1:K,nrow(birthwt),replace=TRUE)
ind_fold
table(ind_fold)
error=numeric(K)
error
for (j in 1:K)
{
  fit.glm=glm(low~.,data=birthwt[ind_fold!=j,],family = binomial)       # 不属于ind_fold 第j个格子的数据集，作为data_train
  y.glm=predict(fit.glm,newdata=birthwt[ind_fold==j,], type = "response")   # 属于ind_fold 第j个格子的数据集，作为data_test
  y.test=birthwt[ind_fold==j,"low"]
  pred = rep(0,length(y.test))
  pred[y.glm>0.5] = 1   # 把 y.glm 中大于 0.5 的那些位置，在 pred 向量里对应的位置赋值为 1
  #y.glm > 0.5 会生成一个逻辑向量（TRUE/FALSE），长度与 y.glm 相同。
  #用这个逻辑向量作为下标：pred[ y.glm > 0.5 ]，会选出 pred 里对应为 TRUE 的那些元素。
  #= 1 则把这些元素改成 1。其余位置保持原值（通常你先把 pred 初始化为 0）。
  error[j]= mean(y.test != pred)   # 在这里不急着开根号
}

CVerror= mean(error) 
CVerror   # 0.321661



# === Q3 ===



