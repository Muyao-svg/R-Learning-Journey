library(VGAM)

x = c(3.1,2.4,2.6,2.2,1.9,2.8,1.1,0.7,2.3,4.3)
n = length(x)
hist(x,prob=T)
# Histogram estimateur non paramerique

# Q1
# 用optimise 返回EMV
fdensite=function(x,a) x/(a**2)*exp(-x**2/2/(a**2))
logL = function(a) sum(log(fdensite(x=x,a=a)))
a_hat = optimize(logL,c(0.01,10),maximum = T)

vlogL = Vectorize(logL,"a")
curve(vlogL, from=0.01, to=10)
points(a_hat$maximum,a_hat$objective, col="red",pch = 4)
a_hat$maximum   #1.787453


# Q2
x1 = sample(x,10,replace = TRUE)
hist(x1)
estimateur_1 = mean(x1)   # estimateur bootstrapt 1.82

# Q3
set.seed(1907)
B=1000
#help(replicate)

a_i_hat = function(){
  # (1) 重新抽样
  x_i = sample(x,n,replace = TRUE)
  # (2) 定义方程式
  fdensite=function(x,a) x/(a**2)*exp(-x**2/2/(a**2))
  logL = function(a) sum(log(fdensite(x=x_i,a=a)))
  # (3) 计算该样本的最大似然估计值
  a_hat = optimize(logL,c(0.01,10),maximum = T)
  return(a_hat$maximum)
}

# 1000个a_i_*_hat
bootstrap_a = replicate(B,a_i_hat())
hist(bootstrap_a,probability = T)
lines(density(bootstrap_a), col = "blue", lwd = 2)

 


# Q4
a_origine = a_hat$maximum   #1.787453
MSE = mean((bootstrap_a - a_origine)**2)
MSE  #0.04637799
#或者用MSE = Biais^2(a_hat) - var(a_hat)
var_a_hat = var(bootstrap_a)
biais = mean(bootstrap_a) - a_origine
MSE_2 = biais**2 + var_a_hat
MSE_2  # 0.04642429



# Q5 
help("quantile")
IC = quantile(bootstrap_a, c(0.025,0.975))
IC   # 1.370749 2.193389





