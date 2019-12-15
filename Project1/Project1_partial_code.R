# third attempt

# open the data file
dataUN = read.table("~/Desktop/Fall Quarter 2017/STA 108/project1/UN.txt", header=T)

# view all data
head(dataUN)

# assign values to variables
y = dataUN$Fertility
x = dataUN$PPgdp

# plot the data with a title and variable names 
plot(x, y, xlab = "PPgdp", ylab = "Fertility", main = "Relationship Between Fertility and PPgdp")

# fix the data to a model
model=lm(y~x)
par(mfrow=c(2,2))
plot(model)

# linearity doesn't hold, transform x first
# data points distribution: smaller values that are together spread further out
# take log on x
xtrans1=log(x)
plot(xtrans1, y, xlab = "log(PPgdp)", ylab = "Fertility", main = "Relationship Between Fertility and log(PPgdp)")
modeltrans1=lm(y~xtrans1)
par(mfrow=c(2,2))
plot(modeltrans1)


# linearity still doesn't hold but looks better
# equal variance doesn't hold
# try transformation on y
boxcox(modeltrans1)
ytrans1 = log(y)
plot(xtrans1, ytrans1, xlab = "log(PPgdp)", ylab = "log(Fertility)", main = "Relationship Between log(Fertility) and log(PPgdp)")
modeltrans2=lm(ytrans1~xtrans1)
par(mfrow=c(2,2))
plot(modeltrans2)


summary(modeltrans2)
summary(modeltrans2)$r.squared





# **********Model fitting and diagnostic**********
# Fit the simple linear model on the transformed data through three ways. 
# Report the least square estimates for the coefficients and R2. 
# Add the fitted line to the scatter plot on the transformed data and comment on the fit.
# **********              3a            **********
# Plain coding (not using the ‘lm’ function or matrix manipulation).

# plain coding
# xbar, ybar
xbar=mean(xtrans1)
ybar=mean(ytrans1)
# SSx, SSxy, b1hat, b0hat, ytrans1hat
SSx=sum((xtrans1-xbar)^2)
SSxy=sum((xtrans1-xbar)*(ytrans1-ybar))
b1hat=SSxy/SSx
b1hat
b0hat=ybar-b1hat*xbar
b0hat
ytrans1hat=b0hat+b1hat*xtrans1
# SSR, SSTO, rsq
SSR=sum((ytrans1hat-ybar)^2)
SSTO=sum((ytrans1-ybar)^2)
rsq=SSR/SSTO
rsq

# add the fitted line
plot(xtrans1, ytrans1, xlab = "Log(PPgdp)", ylab = "Log(Fertility)", main = "Relationship Between Log(Fertility) and Log(PPgdp)")
abline(a=b0hat, b=b1hat)


# **********Model fitting and diagnostic**********
# **********              3b            **********
# Using the ‘lm’ function.


#using lm()
modeltrans2=lm(ytrans1~xtrans1)
# get b1hat, b0hat
modeltrans2$coefficients
summary(modeltrans2)$r.squared
# add the fitted line
plot(xtrans1, ytrans1, xlab = "Log(PPgdp)", ylab = "Log(Fertility)", main = "Relationship Between Log(Fertility) and Log(PPgdp)")
abline(modeltrans2$coefficients)


# **********Model fitting and diagnostic**********
# **********              3c            **********
# Through matrix manipulation.

# xtm: xtrans1 matrix // ytm: ytrans1 matrix 
# cbind: combine the columns // rep(1,n): repeat 1 for n times
xtm = cbind(rep(1,n), xtrans1)
ytm = cbind(ytrans1)

# if a%*%X=b, solve(a,b): solve for X
# if only solve(a): gives the inverse of a
# t(a): gives the transpose of a
# In least square equation, Y=AX+E, A=[solve(Xt%*%X)]%*%[Xt%*%Y], where Xt is the transpose of X
# betam: beta matrix // xttm: xtrans1 transpose
xttm = t(xtm)
betam = (solve(xttm%*%xtm))%*%(xttm%*%ytm)
betam

# ythatm: yhat matrix // em: epsilon matrix 
ythatm = xtm%*%betam

# find rsq: SSR=sum((ytrans1hat-ybar)^2), SSTO=sum((ytrans1-ybar)^2), rsq=SSR/SSTO
# rowMeans(a):returns vector of row means
# colMeans(a):returns vector of column means
# rowSums(a): returns vector of row sums
# colSums(a): returns vector of column means
# yhyb: ytrans1hat-ybar // yyb: ytrans1-ybar // ybarm: ybar matrix
ybarm = colSums(ytm)/n
yhyb = ythatm-ybarm
yyb = ytm-ybarm
SSR = colSums(yhyb*yhyb)
SSTO = colSums(yyb*yyb)
rsq = SSR/SSTO
rsq

# get b0, b1 by betam(row, column) and draw the line
betam[1,1]
betam[2,1]
plot(xtrans1, ytrans1, xlab = "Log(PPgdp)", ylab = "Log(Fertility)", main = "Relationship Between Log(Fertility) and Log(PPgdp)")
abline(a=betam[1,1],b=betam[2,1])




# **********Model fitting and diagnostic**********
# **********              4             **********
# Draw the diagnostic plots and comment.




# **********Making inferences based on the model**********
# **********                   5                **********
# Test whether there is a linear relationship between the transformed variables at 0.05 significance level.

# calculate T*=abs(b1hat/SEb1hat)
# known b1hat, b0hat from previous calculation
b1hat
b0hat
# n, e, SSE, MSE, sigmahat
n=dim(dataUN)[1]
n
e=ytrans1-ytrans1hat
SSE=sum(e^2)
MSE=SSE/(n-2)
sigmahat=sqrt(MSE)

# standard error for b1hat
SEb1hat=sigmahat/sqrt(SSx)
SEb1hat
# standard error for b0hat
SEb0hat=sqrt(MSE*((1/n)+((xbar^2)/SSx)))
SEb0hat
b1hat/SEb1hat
abs(b1hat/SEb1hat)
# check with summary()
summary(modeltrans2)





# **********Making inferences based on the model**********
# **********                   6                **********
# Provide a 99% confidence interval on the expected Fertility for a region with PPgdp 20,000 US dollars in 2001.
a=1-0.99
1-a/2
xnew=20000
xtrans1new=log(xnew)
ytrans1new=b0hat+b1hat*xtrans1new



# SEytrans1new
SEytrans1new = sqrt(MSE*((1/n)+((xtrans1new-xbar)^2/SSx)))
SEytrans1new
# lbytrans1new, ubytrans1new
lbytrans1new=ytrans1new-qt(1-a/2,n-2)*SEytrans1new
lbytrans1new
ubytrans1new=ytrans1new+qt(1-a/2,n-2)*SEytrans1new
ubytrans1new

# calculate the lb, ub for y
lb=exp(lbytrans1new)
lb
ub=exp(ubytrans1new)
ub


# **********Making inferences based on the model**********
# **********                   7                **********
# Provide a 95% confidence band for the relation between the expected Fertility and PPgdp. 
# Add the bands to the scatter plot of the original data.


W=sqrt(2*qf(.95, df1=2, df2=n-2))

yhat=function(x){ b0hat+b1hat*x }
sigmahat2=(summary(modeltrans2)$sigma)^2
xbar=mean(xtrans1)
ybar=mean(ytrans1)
SSx=sum((xtrans1-xbar)^2)
seyhat=function(xh){ sqrt(sigmahat2*(1/n + (xh-xbar)^2/SSx)) }

cband=function(xvec){
  d=length(xvec)
  CIs=matrix(0,d,2)
  colnames(CIs)=c("Lower", "Upper")
  
  for(i in 1:d){  CIs[i,]=c(yhat(x=xvec[i])+c(-1,1)*W*seyhat(xh=xvec[i]))  }
  as.data.frame(CIs)
}

# get the range of x, create the bands
maxx=max(x)
minx=min(x)
range = seq(from=minx, to=maxx, length.out=1000)
bands = cband(log(range))

# draw bands on the original plot
plot(x, y, xlab = "PPgdp", ylab = "Fertility", main = "Relationship Between Fertility and PPgdp", pch=20, cex=.5)
points(range, exp(bands$Lower), col = "red", type = "l")
points(range, exp(bands$Upper), col = "red", type = "l")










# **********Making inferences based on the model**********
# **********                   8                **********
# Assuming that the same relationship between Fertility and PPgdp holds, 
# give a 99% prediction interval on Fertility for a region with PPgdp 25,000 US dollars in 2018.

a=1-0.99
1-a/2
xpred=25000
xtrans1pred=log(xpred)
ytrans1pred=b0hat+b1hat*xtrans1pred

# SEytrans1pred
SEytrans1pred = sqrt(MSE*(1+(1/n)+((xtrans1pred-xbar)^2/SSx)))
# lbytrans1pred, ubytrans1pred
lbytrans1pred=ytrans1pred-qt(1-a/2,n-2)*SEytrans1pred
lbytrans1pred
ubytrans1pred=ytrans1pred+qt(1-a/2,n-2)*SEytrans1pred
ubytrans1pred

# calculate the lbpred, ubpred for y
lbpred=exp(lbytrans1pred)
lbpred
ubpred=exp(ubytrans1pred)
ubpred


# **********Making inferences based on the model**********
# **********                   9                **********
# Based on the diagnostic plots in Part 4, do you have any concern on the above hypothesis testing and inferences? 
# If so, what are the concerns?