# Project2 STA108
# read in data
diabetesoriginal=read.table("~/Desktop/Fall Quarter 2017/STA 108/project2/diabetes.txt", header=T)
dim(diabetesoriginal)
head(diabetesoriginal)

# arrange the data set a little bit, assign variables, value the factors
diabetes=cbind(diabetesoriginal$glyhb, diabetesoriginal[,-5])
head(diabetes)

# assign variables
#   Y=diabetes$glyhb,  X1=diabetes$chol,   X2=diabetes$stab.glu, X3=diabetes$hdl,     X4=diabetes$ratio,  X5=diabetes$location
#  X6=diabetes$age,    X7=diabetes$gender, X8=diabetes$height,   X9=diabetes$weight, X10=diabetes$frame, X11=diabetes$bp.1s,
# X12=diabetes$bp.1d, X13=diabetes$waist, X14=diabetes$hip,      X15=diabetes$time.ppn

length(colnames(diabetes))
colnames(diabetes)=c("Y","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15")
head(diabetes)

# 1
# str(filename): to see which are quantitative variables and qualitative variables
str(diabetes)

# group, view, and take out the factors, then view the data without factors
factor=c(6,8,11)                    # to group: 6,8,11 are column number
head(diabetes[,factor])             # to view the factors
head(diabetes[,-factor])            # to view the data without factors

# draw histogram
#   Y=diabetes$glyhb,  X1=diabetes$chol,   X2=diabetes$stab.glu, X3=diabetes$hdl,     X4=diabetes$ratio,  X5=diabetes$location
#  X6=diabetes$age,    X7=diabetes$gender, X8=diabetes$height,   X9=diabetes$weight, X10=diabetes$frame, X11=diabetes$bp.1s,
# X12=diabetes$bp.1d, X13=diabetes$waist, X14=diabetes$hip,      X15=diabetes$time.ppn
hist(diabetes$X1, main="Histogram of diabetes$chol", breaks=40)
hist(diabetes$X2, main="Histogram of diabetes$stab.glu", breaks=40)
hist(diabetes$X3, main="Histogram of diabetes$hdl", breaks=40)
hist(diabetes$X4, main="Histogram of diabetes$ratio", breaks=40)
hist(diabetes$X6, main="Histogram of diabetes$age", breaks=40)
hist(diabetes$X8, main="Histogram of diabetes$height", breaks=40)
hist(diabetes$X9, main="Histogram of diabetes$weight", breaks=40)
hist(diabetes$X11, main="Histogram of diabetes$bp.1s", breaks=40)
hist(diabetes$X12, main="Histogram of diabetes$bp.1d", breaks=40)
hist(diabetes$X13, main="Histogram of diabetes$waist", breaks=40)
hist(diabetes$X14, main="Histogram of diabetes$hip", breaks=40)
hist(diabetes$X15, main="Histogram of diabetes$time.ppn", breaks=40)

# draw pie chart: table(variable.name)-show the weight of each frame
pie(table(diabetesoriginal$location), main="Pie Chart of diabetes$location", col=c("pink","white"))
pie(table(diabetesoriginal$gender), main="Pie Chart of diabetes$gender", col=c("yellow","white"))
pie(table(diabetesoriginal$frame), main="Pie Chart of diabetes$frame", col=c("black","white","grey"))

# draw scatterplot matrix
pairs(~X1+X2+X3+X4+X6+X8+X9+X11+X12+X13+X14+X15, data=diabetes, main="Scatterplot Matrix for Diabetes Data", cex.main=0.8)

# obtain pairwise correlation matrix for all quantitative variables and view them
View(cor(diabetes[, -factor]))





# 2
model1=lm(Y~., data=diabetes)
par(mfrow = c(2,2))
plot(model1)
# when encountering error "figure margins too large": type graphics.off() in console, and try again





# 3
library(MASS)                      # for boxcox()
boxcox(model1)
diabetes$Ytrans=1/diabetes$Y       # transform Y to 1/Y and attach it to the data
diabetes=diabetes[, names(diabetes)!="Y"]    # extract all variables not named Y
View(diabetes)
head(diabetes)

model2=lm(Ytrans~., data=diabetes)
par(mfrow = c(2,2))
plot(model2)
boxcox(model2)



# 4
set.seed(10)                       # set seed for random variable generator so everyone gets the same split of the data
N=nrow(diabetes)                   # number of cases in diabetes
N
# randomly sample N/2 cases to form the training data
index=sample(1:N, size=N/2, replace=FALSE)
data.t=diabetes[index,]            # get the training data set
data.v=diabetes[-index,]           # use the remaining cases to form the validation set




# 5
model3=lm(Ytrans~., data=data.t)
summary(model3)
countbeta=length(model3$coefficients)
countbeta
MSEmodel3=summary(model3)$sigma^2
MSEmodel3



# 6
install.packages("leaps")
library(leaps)
# nbest-number of best model to be selected, nvmax-number of maximum models in pool
best=regsubsets(Ytrans~., data=data.t, nbest=1, nvmax=16)
summary.subset=summary(best)
summary.subset                     # return the top models in their subset
summary.subset$which               # clearly see which X to be included in top models

# SSEp, R2p, Ra2p, Cp, AICp, BICp for top models and none model
n=nrow(data.t)                     # number of cases in data.t
n
p.m=as.integer(as.numeric(rownames(summary.subset$which))+1)
p.m                                # coding "p.m=2:17" to get p.m is fine too

sse=summary.subset$rss
sse
r2=summary.subset$rsq
r2
ra2=summary.subset$adjr2
ra2
c=summary.subset$cp
c
aic=n*log(sse)+2*p.m-n*log(n)
aic
bic=n*log(sse)+log(n)*p.m-n*log(n)
bic


list.result=cbind(summary.subset$which,sse,r2,ra2,c,aic,bic)
list.result

modelnone=lm(Ytrans~1,data=data.t) # fit the none model
ssenone=sum(modelnone$residuals^2)
ssenone
p=1
r2none=0
ra2none=0
cnone=ssenone/MSEmodel3-(n-2*p)
cnone
aicnone=n*log(ssenone)+2*p-n*log(n)
aicnone
bicnone=n*log(ssenone)+log(n)*p-n*log(n)
bicnone

list.resultnone=c(1,rep(0,16),ssenone,r2none, ra2none,cnone,aicnone,bicnone)
list.resultnone

# combine list.resultnone to list.result, set columm names
list.result=rbind(list.resultnone,list.result)
colnames(list.result)=c(colnames(summary.subset$which), "SSE", "R2", "Ra2", "C","AIC","BIC")
colnames(list.result)
list.result

Frametype=model.matrix(~data.t$X10-1)
Frametype                          # show the columns of Frametype in terms of different frames


# base on the result, best models:
model3.1=lm(Ytrans~X2+X4+X6+Frametype[,3]+X13, data=data.t)            # in terms of AIC
model3.2=lm(Ytrans~X2+X6+X13, data=data.t)                             # in terms of BIC
model3.3=lm(Ytrans~X2+X4+X6+Frametype[,3]+X13+X15, data=data.t)        # in terms of Ra2



# 7
model4=lm(Ytrans~.^2, data=data.t)
countcof=length(model4$coefficients)
countcof
MSEmodel4=summary(model4)$sigma^2
MSEmodel4


# 8
fsp1=stepAIC(modelnone, scope=list(lower=modelnone,upper=model4), direction="both", k=2, data=data.t)
summary(fsp1)   # In the results, "X2:X4" means X2*X4.
sse.fsp1=sum(fsp1$residuals^2)
sse.fsp1
p.fsp1=length(fsp1$coefficients)
p.fsp1


# 9
fsp2=stepAIC(model3, scope=list(lower=modelnone,upper=model4), direction="both", k=2, data=data.t)
summary(fsp2)
sse.fsp2=sum(fsp2$residuals^2)
sse.fsp2
p.fsp2=length(fsp2$coefficients)
p.fsp2



# 10
bic.fsp1=n*log(sse.fsp1)+log(n)*p.fsp1
bic.fsp1
bic.fsp2=n*log(sse.fsp2)+log(n)*p.fsp2
bic.fsp2


model4.1=fsp2        # in terms of AIC
model4.2=fsp1        # in terms of BIC


# 11
# PRESS=sum( ( ei/(1-hii) )^2) where e=(yi-yiminusihat)

press.model3.1=sum(model3.1$residuals^2/(1-lm.influence(model3.1)$hat)^2)
press.model3.1

press.model3.2=sum(model3.2$residuals^2/(1-lm.influence(model3.2)$hat)^2)
press.model3.2

press.model3.3=sum(model3.3$residuals^2/(1-lm.influence(model3.3)$hat)^2)
press.model3.3

press.model4.1=sum(model4.1$residuals^2/(1-lm.influence(model4.1)$hat)^2)
press.model4.1

press.model4.2=sum(model4.2$residuals^2/(1-lm.influence(model4.2)$hat)^2)
press.model4.2
# Since model4.1 has the smallest PRESS, we choose it as the best model among the five.


# 12
size=nrow(data.t)
size

mspr=function(model){
    yhat=predict(model, data.v)
    mspr = sum((data.v$Ytrans-yhat)^2) / size
    print(mspr)
    return(mspr)
}

mspr.model3.1=mspr(model3.1)
mspr.model3.2=mspr(model3.2)
mspr.model3.3=mspr(model3.3)
mspr.model4.1=mspr(model4.1)
mspr.model4.2=mspr(model4.2)
# choose model3.3 because it has the smallest mspr.

press.model3.1/size
press.model3.2/size
press.model3.3/size
press.model4.1/size
press.model4.2/size


# 13
model3.3
Frametypenew=model.matrix(~diabetes$X10-1)
model5=lm(Ytrans ~ X2 + X4 + X6 + Frametypenew[, 3] + X13 + X15, data=diabetes)
summary(model5)
anova(model5)
