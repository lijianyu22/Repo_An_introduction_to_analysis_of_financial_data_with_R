# read data, use table, param header
da=read.table("w-petroprice.txt",header = T)
da1=read.table("w-gasoline.txt")
# set as log func
pus=log(da$US)
pgs=log(da1[,1])
# transfer the discrete data to series.
tdx=c(1:717)/52+1997
# plot as some mode.
par(mfcol=c(2,1))
plot(tdx,pgs,xlab='year',ylab='ln(price)',type='l')
title(main='(a) Gasline')
plot(tdx,pus,xlab='year',ylab='ln(price)',type='l')
title(main='(b) Crude oil')
# diff the ts data, get stable data.
dpgs=diff(pgs)
acf(dpgs,lag=20)
pacf(dpgs,lag=20)

# 1.从acf和pacf数据来看，acf具备五阶相关，且偏指数下降 -> AR模型
# 2.pacf前五阶明显，且第一阶最为明显，说明与第一阶强烈相关，相关程度比其他阶次都高
# s:考虑AR(5)
m1=ar(diff(pgs),method = 'mle')
# 使用ar()直接拟合，主要参数1.系数协方差矩阵va.coef -> 每个系数的标准误差 -> t值（>=2 ？可以估计出该参数对整个模型影响程度
#  2.AIC （R2） 
t.test(dpgs)

# 使用arima（）来拟合，c(5,0,0)为(AR,DIF,MA)，比如其与ar()区别，发现几无区别；
m2=arima(dpgs,order=c(5,0,0),include.mean = F)
# 使用fixed函数来去除影响小系数
m1=arima(dpgs,order=c(5,0,0),include.mean = F,fixed = c(NA,NA,NA,0,NA))
# 残差计算，三张图：1.标准化残差的分布 2. 残差的ACF PACF，观察残差是否有相关性；
tsdiag(m1,gof=20)

#  加入pus，即外部变量
# 
# 
dpus=diff(pus)
# lm( ~1 + x) 含义为添加两者的线性关系并包含常数项
m3=lm(dpgs~1+dpus)
# 观察残差，发现残差仍具备相关性，且R2为33%，不是很接近1
# 注：因变量，自变量的关系用R2，而其他模型一般用AIC BIC
# 再注：m3输出没有特定的r，但是有r.squared
acf(m3$residuals,lag=20)
pacf(m3$residuals,lag=20)
summary(m3)$r.squared

# 使用arima中的xreg=dpus自动拟合，同时以相同方式去因子
m4=ar(m3$residuals,method='mle')
m4=arima(dpgs,order=c(6,0,0),include.mean = F,xreg=dpus)
m4=arima(dpgs,order=c(5,0,0),include.mean = F,xreg=dpus)
m4=arima(dpgs,order=c(5,0,0),include.mean = F,xreg=dpus,fixed = c(NA,NA,NA,0,NA,NA))
tsdiag(m4,gof=20)

# 回测函数，从序列的第316个开始预测
c1=c(NA,NA,NA,0,NA)
source('backtest.R')
pm1=backtest(m1,dpgs,316,1,fixed = c1,inc.mean = F)
c4=c(NA,NA,NA,0,NA,NA)
pm4=backtest(m4,dpgs,316,1,xre=dpus,inc.mean = F,fixed = c4)

# 预测结果绘图
tdx=tdx[2:717]
pm4fit=dpgs[317:716]-pm4$error
pm1fit=dpgs[317:716]-pm1$error
plot(tdx[317:716],dpgs[317:716],xlab='year',ylab='growth',type='l')
points(tdx[317:716],pm1fit,pch='*')
plot(tdx[317:716],dpgs[317:716],xlab='year',ylab='growth',type='l')
points(tdx[317:716],pm4fit,pch='*')

# 滞后一个周期，lm -> ar -> arima -> bt
m6=lm(dpgs[2:716]~-1+dpus[1:715])
acf(m6$residuals,lag=20)
pacf(m6$residuals,lag=20)

m7=ar(m6$residuals,method = 'mle')
m7=arima(dpgs[2:716],order = c(9,0,0),xreg=dpus[1:715],include.mean = F)
m7=arima(dpgs[2:716],order=c(9,0,0),xreg=dpus[1:715],include.mean = F,fixed=c(NA,NA,NA,0,NA,0,0,0,NA,NA))
tsdiag(m7,gof=20)
c7=c(NA,NA,NA,0,NA,0,0,0,NA,NA)
pm7=backtest(m7,dpgs[2:716],315,1,xre=dpus[1:715],inc.mean = F,fixed = c7)

# 预测结果展示：
# 经典模型
# [1] "RMSE of out-of-sample forecasts"
# [1] 0.02171235
# [1] "Mean absolute error of out-of-sample forecasts"
# [1] 0.01537881
# 
# # 带其他参数的模型
# [1] "RMSE of out-of-sample forecasts"
# [1] 0.01925732
# [1] "Mean absolute error of out-of-sample forecasts"
# [1] 0.01285104
# 
# # 带其他参数（滞后一个周期）的模型
# [1] "RMSE of out-of-sample forecasts"
# [1] 0.0216638
# [1] "Mean absolute error of out-of-sample forecasts"
# [1] 0.01548401
# 可以看出，经典模型效果 - 带其他参数


# 全球温度变化数据研究
# 
# 观察数据的一阶差分特点，基于此给出ARIMA模型的预估
Gt=scan(file='m-GLBTs.txt')
Gtemp=ts(Gt,frequency = 12,start = c(1880,1))
plot(Gtemp,xlab='year',ylab='temperature',type='l')
par(mfcol=c(2,1))
acf(diff(Gt),lag=36)
pacf(diff(Gt),lag=36)

# 1. AR 一阶主 2. MA 两阶（非指数） 3. 可以看出季节性，在24，12处
m1=arima(Gt,order=c(1,1,2))
tsdiag(m1,gof=20)
acf(m1$residuals,lag=36)

# 添加季节参数
season = list(order=c(0,0,1),period=24)
m2=arima(Gt,order = c(1,1,2),seasonal = list(order=c(0,0,1),period=24))
tsdiag(m1,gof=36)

# 比较两个模型
source('backtest.R')
pm1=backtest(m1,Gt[2:716],315,1,inc.mean = F)
pm2=backtest(m2,Gt[2:716],315,1,inc.mean = F)
#  AIC基本无差距，回测也是，说明24阶季节性影响很小

# [1] "RMSE of out-of-sample forecasts"
# [1] 15.74282
# [1] "Mean absolute error of out-of-sample forecasts"
# [1] 12.3788
# > pm2=backtest(m2,Gt[2:716],315,1,inc.mean = F)
# [1] "RMSE of out-of-sample forecasts"
# [1] 15.66913
# [1] "Mean absolute error of out-of-sample forecasts"
# [1] 12.29691

#  通过和时间趋势做解释，看到残差在8,29阶次有显著。
time=c(1:1568)
m2=lm(Gt~time)
par(mfcol=c(2,1))
acf(m2$residuals,lag=36)
pacf(m2$residuals,lag=36)
# 为什么非差分，使用AR2 MA1加时间可以消除。
# 由于时间吗，那为什么是2和1
m2=arima(Gt,order = c(2,0,1),xreg=time)
tsdiag(m2,gof=36)

m3=arima(Gt,order=c(2,0,1),seasonal = list(order=c(0,0,1),period=24),xreg=time)
tsdiag(m3,gof=36)


pm1=backtest(m1,Gt,1368,1)
time=as.matrix(time)
pm2=backtest(m2,Gt,1368,1,xre=time)


time=c(1:1568)
time1=c(rep(0,1212),time[1213:1568])
mm1=lm(Gt~time+time1)
x1=cbind(time,time1)
mm1=arima(Gt,order=c(2,0,1),seasonal = list(order=c(0,0,1),period=24),xreg=x1)
tsdiag(mm1,gof=36)
Box.test(mm1$residuals,lag=8,type='Ljung')

da=read.table('m-unrate.txt',header=T)
unemp=da$rate
unrate=ts(unemp,frequency = 12,start = c(1948,1))
plot(unrate,xlab='year',ylab='unrate',type='l')
par(mfcol=c(2,2))
acf(unemp,lag=36)
pacf(unemp,lag=36)
acf(diff(unemp),lag=36)
pacf(diff(unemp),lag=36)
m1=arima(unemp,order=c(1,1,5),seasonal = list(order=c(1,0,1),period=12))
c1=c(NA,NA,NA,0,0,NA,NA,NA)
m1=arima(unemp,order=c(1,1,5),seasonal = list(order=c(1,0,1),period=12),fixed=c1)
tsdiag(m1,gof=36)
Box.test(m1$residuals,lag=24,type='Ljung')
Box.test(m1$residuals,lag=36,type='Ljung')
mm=arima(unemp,order=c(0,1,0),seasonal = list(order=c(1,0,1),period=12))
par(mfcol=c(2,1))
acf(mm$residuals,lag=24)
pacf(mm$residuals,lag=24)
mm1=arima(unemp,order = c(5,1,0),seasonal = list(order=c(1,0,1),period=12))
cc1=c(0,NA,NA,NA,NA,NA,NA)
mm1=arima(unemp,order=c(5,1,0),seasonal = list(order=c(1,0,1),period=12),fixed = cc1)
tsdiag(mm1,gof=36)
source('backtest.R')
pm1=backtest(m1,unemp,700,1,fixed = c1,inc.mean = F)
pmm1=backtest(mm1,unemp,700,1,fixed = cc1,inc.mean = F)

da=read.table('m-unrateic.txt',header = T)
unrate=da$rate
x=da[,5:9]/1000
nm1=lm(unrate~icm1,data=x)
