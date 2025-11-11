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
#  0   1   2   3   4   6 
#100  47  30   7   4   1   这个访问6次的只有一组家庭，是否要保留这一组是个问题
# 我们可以选择一个seuil，<2的为一组，>2的为一组
  
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
birthwt=birthwt[,-1]
summary(birthwt)   # 现在factor 都是factor了
attach(birthwt)
# df = 189 - 11


# Representation Graphique
hist(bwt)
plot(lwt,bwt,main="relation entre lwt et bwt")
boxplot(bwt~race)


## === 2. Previson par regression lineaire ===
model = lm(bwt~.,data = birthwt)
summary(model)
# 如何测试某个variable 是否significative? --> anova( model_sans_ui,model_complet)  这是test Ficher
# P_value < 0.05, on rejette H0 model_sans_age, ui是variable sinificative

# Modele bwt~age+ui+race
model_1 = lm(bwt~age+ui+race,data=birthwt)
summary(model_1)

predict(model,newdata = data.frame(age = 23,lwt = 94,race="other",smoke ="Yes",
                                   ptl="FALSE",ht="FALSE",ui="No",ftv="0"))   # 2596.384 
# Erreur d'apprentissage (Training error): RSS/n
err_app = mean((birthwt$bwt - predict(model))**2)
err_app    # 393096.3
err_app = sqrt(err_app)   # 这个也叫“RMSE”
err_app    # 626.9739
# 结果并不好，因为是biaise 的estimateur



## === 3. Erreur test ===
# === 1.Decoupage train/test ===
test = sample(1:length(birthwt$bwt),60)
train = -test
train = birthwt[train, ]
test = birthwt[test,]

model_train = lm(bwt~.,data = train)
MSE = mean((test$bwt - predict(model_train,newdata = test))**2)
MSE   # 423787  Erreur test
RMSE_test = sqrt(MSE)
RMSE_test   # 650.9892  这个结果通常比RMSE_train 高
# Comment rendre le resultat reproductible? 同样的代码大家得出来的结果并不一样,而且非常fluctuation!!!


## === 4. Validation coisee 交叉验证===
# 把数据多次拆分成“训练集 + 验证集”，轮流训练和测试，再取平均误差 → 更可靠地估计模型的泛化能力。
K=5
ind_fold=sample(1:K,nrow(birthwt),replace=TRUE)
table(ind_fold)
error=numeric(K)
for (j in 1:K)
{
  fit.lm=lm(bwt~.,data=birthwt[ind_fold!=j,])       # 不属于ind_fold 第j个格子的数据集，作为data_train
  y.lm=predict(fit.lm,newdata=birthwt[ind_fold==j,])   # 属于ind_fold 第j个格子的数据集，作为data_test
  y.test=birthwt[ind_fold==j,"bwt"]
  error[j]= mean((y.test - y.lm)**2)    # 在这里不急着开根号
}
CVerror= mean(error) # à compléter
CVerror   # 451655.3
CVerror = sqrt(CVerror)   # 在全部算完后再开根号
CVerror   # 672.053   对比别的同学的结果，variability 明显比RSS/n小，大家的结果都在672附近

