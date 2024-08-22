library(MASS)
library(lattice)
library(olsrr)
library(car)
library(L1pack)
setwd("D:\\PG files\\Projects\\Regression-Analysis-Project-main")
X = read.csv(file = "population_drinking1.txt",header = TRUE,sep = "\t")
names(X) <- c("Ind","Ind_1","Urban population", "Late births", "Wine consumption per capita", "Liquor consumption per capita","Cirrhosis death rate")
head(X)
#### Type of the covariates
str(X)
names(X) <- c("I","1","A1","A2","A3","A4","Y")
summary(X[,-c(1,2)])
pairs(X[,-c(1,2)])
cor(X[,-c(1,2,7)])
crPlots(lm(Y~A1+A2+A3+A4,data = X))
boxplot(X[,-c(1,2)])
colnames(X)=c("I","1","A1","A2","A3","A4","Y")
attach(X)
reg <- lm(Y~A1+A2+A3+A4)
summary(reg)
plot(1:nrow(X),X$Y,type = "o",pch = 20,ylab = "observed & predicted values",xlab = "Index")
lines(1:nrow(X),reg$fitted.values,type = "o",pch = 20,col = "red",lty = 2)
abline(v = 1:nrow(X),lty = 2,col = rgb(0,1,0,alpha = 0.3))
legend("topleft",legend = c("Obs","Pred"),fill = c("black","red"))
resi<-residuals(reg)
qqnorm(resi,ylab="Residuals",main="")
qqline(resi)
shapiro.test(resi)
plot(fitted(reg),residuals(reg),xlab="Fitted",ylab="Residuals")
abline(h=0)
A = as.matrix(X[,-1])
H = A%*%solve(t(A)%*%A)%*%t(A)
H_i = diag(H)
e_i = residuals(reg)
b_i = e_i^2/(1-H_i)
plot(fitted(reg),b_i,pch = 20,main = bquote("Plot of" ~ b[i]^2 ~ "vs fitted values"),xlab = "",ylab = "")

library(lmtest)
bptest(reg)

acf(resi,ylab = "",xlab = "",main = "")
title(xlab="Lag", ylab="ACF", line=2)

pacf(resi,ylab = "",xlab = "",main = "")
title(xlab="Lag", ylab="PACF", line=2)

require(lmtest)
dwtest(Y~A1+A2+A3+A4,data=X)

require(lmtest)
bgtest(reg,order = 20)

hat_d <- hatvalues(reg)
head(sort(hat_d,dec=TRUE))

n = nrow(X);p = 5
hat_d[hat_d > 2*(p/n)]

stud <- rstudent(reg)
plot(stud,ylim = c(-3,3),pch=20,col = "blue")
abline(h=c(0))

ols_plot_resid_stud_fit(reg)

par(mfrow = c(2,3))
DFBETAS = dfbetas(reg)
for(i in 1:5)
{
  plot(DFBETAS[,i],main=bquote("DFBETAS for" ~ beta[.(i-1)]),ylab="",ylim=c(-1.5,1.5),xlab="",pch=20)
  abline(h=c(-2/sqrt(n),2/sqrt(n)))
  ind = which(abs(DFBETAS[,i]) > 2/sqrt(n)) # beta_0
  text = text(ind,DFBETAS[ind,i],pos = 3,labels = ind)
}

ols_plot_dffits(reg) 

COVRATIO = covratio(reg)
plot(abs(COVRATIO-1),ylab=expression(abs(COVRATIO-1)),pch = 20)
abline(h = 3*p/n)
ind = which(abs(COVRATIO-1) >= 3*p/n)
text(ind,abs(COVRATIO[ind]-1),pos = 1,labels = ind)

COOKSD = cooks.distance(reg)
plot(COOKSD,pch=20)
abline(h = 3*mean(COOKSD))
ind = which(COOKSD > 3*mean(COOKSD)) # beta_0
text = text(ind,COOKSD[ind],pos = 1,labels = ind)

summary(lm(Y~A1+A2+A3+A4,data=X[-c(12,20,30,36),]))

reg_out = lm(Y~A1+A2+A3+A4,data = X[-c(12,20,30,36),])
hist(rstandard(reg_out),breaks = 5,main = "Histogram Of Residual Values") 

reg_out = lm(Y~A1+A3+A4,data = X[-c(12,20,30,36),])
qqnorm(rstandard(reg_out)) 
qqline(rstandard(reg_out)) 

X_mdl = model.matrix(reg)[,-1]
kappa(scale(X_mdl))

require("faraway")
vif(lm(Y~A1+A2+A3+A4,data = X[-c(12,20,30,36),]))#high implies collinearity 

require("faraway")
vif(lm(Y~A1+A2+A3+A4,data = X[-c(12,20,30,36),]))#high implies collinearity 

require("faraway")
vif(lm(Y~A1+A2+A3+A4,data = X[]))#high implies collinearity 

set.seed(2124)
lmod_per = lm(Y+10*rnorm(nrow(X)-4,s = 5) ~ A1+A2+A3+A4,data = X[-c(12,20,30,36),])
reg$coefficients
lmod_per$coefficients

vif(reg_out)

avPlots(lm(Y~A1+A2+A3+A4,data = X[-c(12,20,30,36),])) 

mod_wt_out = lm(Y ~ A1+A2+A3+A4, data = X[-c(12,20,30,36),])
ols_step_both_aic(mod_wt_out,details = TRUE)

X_m = X[,-c(1,2)]
X1 <- X_m[-c(12,20,30,36),]
names(X1) <- c("V1","V2","V3","V4","V5")
models<-list()
models[["V1"]]<-lm(V5~V1,X1)
models[["V2"]]<-lm(V5~V2,X1)
models[["V3"]]<-lm(V5~V3,X1)
models[["V4"]]<-lm(V5~V4,X1)
models[["V12"]]<-lm(V5~V1+V2,X1)
models[["V13"]]<-lm(V5~V1+V3,X1)
models[["V14"]]<-lm(V5~V1+V4,X1)
models[["V23"]]<-lm(V5~V2+V3,X1)
models[["V24"]]<-lm(V5~V2+V4,X1)
models[["V34"]]<-lm(V5~V3+V4,X1)
models[["V123"]]<-lm(V5~V1+V2+V3,X1)
models[["V124"]]<-lm(V5~V1+V2+V4,X1)
models[["V134"]]<-lm(V5~V1+V3+V4,X1)
models[["V234"]]<-lm(V5~V2+V3+V4,X1)
models[["V1234"]]<-lm(V5~V1+V2+V3+V4,X1)

mnames<-factor(names(models),levels = names(models))

R2 <- sapply(models, function(fit) summary(fit)$r.squared)
dotplot(R2 ~ mnames, type = "o", pch = 16,auto.key=list(space="right"),xlab="Models",ylab=expression(R^2))


adj.R2 <- sapply(models, function(fit) summary(fit)$adj.r.squared)
dotplot(R2 + adj.R2 ~ mnames, type = "o", pch = 16,auto.key=list(space="right"),xlab="Models")

sigma.sq <- summary(models[["V1234"]])$sigma^2
Cp <- sapply(models, function(fit) extractAIC(fit, scale = sigma.sq)[2])
dotplot(Cp ~ mnames, type = "o", pch = 16)

AIC <- sapply(models, function(fit) AIC(fit))
dotplot(AIC ~ mnames, type = "o", pch = 16,xlab="Models")

BIC <- sapply(models, function(fit) BIC(fit))
dotplot(BIC ~ mnames, type = "o", pch = 16,xlab="Models")

CV = NULL
for(i in 1:15)
{
  X_mdl_mat = model.matrix(models[[i]])
  head(X_mdl_mat)
  Y_vec = X1$V5
  H = X_mdl_mat%*%solve(t(X_mdl_mat)%*%X_mdl_mat)%*%t(X_mdl_mat)
  h = diag(H)
  n_h = nrow(X1)
  CV[i] = (1/n_h)*sum((Y_vec-H%*%Y_vec)^2/(1-h)^2)
}
dotplot(CV ~ mnames, type = "o", pch = 16,xlab="Models")

data.frame(R2,adj.R2,Cp,AIC,BIC,CV)

X1 = X[-c(12,20,30,36),]
mod_opt = lm(Y ~ A1+A3, data = X[-c(12,20,30,36),])
plot(1:nrow(X1),X1$Y,type = "o",pch = 20,ylab = "observed & predicted values",xlab = "Index")
lines(1:nrow(X1),mod_opt$fitted.values,type = "o",pch = 20,col = "red",lty = 2)
abline(v = 1:nrow(X1),lty = 2,col = rgb(0,1,0,alpha = 0.3))
legend("topleft",legend = c("Obs","Pred"),fill = c("black","red"))

mod_wt_out = lm(Y ~ A1+A2+A3+A4, data = X[-c(12,20,30,36),])
STEP_REG = stepAIC(mod_wt_out,scope = list(upper = ~(A1+A2+A3+A4)^2, lower = ~1),trace = TRUE)
STEP_REG

k_grid = 10^seq(-2,1/2,length.out = 100)
PE = NULL
X_R = as.matrix(X[-c(12,20,30,36),c(2,3,4,5,6)])
Y_R = as.matrix(X[-c(12,20,30,36),c(7)])
n = nrow(X_R)
for(i in 1:length(k_grid))
{
  k = k_grid[i]
  beta_k = solve(t(X_R)%*%X_R + k*diag(rep(1,5)))%*%t(X_R)%*%Y_R
  A_k = X_R%*%solve(t(X_R)%*%X_R + k*diag(rep(1,5)))%*%t(X_R)
  Y_ft_R = X_R%*%beta_k
  A_K_diag = diag(A_k)
  PE[i] = (1/n)*sum((Y_R-Y_ft_R)^2/(1-A_K_diag)^2)
}

plot(k_grid,PE,type = "l",ylab = "CV(1)",xBlab = "k")
k_opt = k_grid[which(PE == min(PE))]
abline(v = k_opt)

beta_opt = solve(t(X_R)%*%X_R + k_opt*diag(rep(1,5)))%*%t(X_R)%*%Y_R
plot(1:nrow(X_R),X_R%*%beta_opt,type = "o",pch = 20,ylab = "observed & predicted values",xlab = "Index",ylim = c(30,110))
lines(1:nrow(X_R),Y_R,type = "o",pch = 20,col = "red",lty = 2)
abline(v = 1:nrow(X),lty = 2,col = rgb(0,1,0,alpha = 0.3))
legend("topleft",legend = c("Obs","Pred"),fill = c("black","red"))

library(glmnet)
x <- as.matrix(X[-c(12,20,30,36),c(3,4,5,6)])
y <- X[-c(12,20,30,36),7]
lambdas <- 10^seq(-1, 5, by = 0.1)

lasso_reg <- cv.glmnet(x,y, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 10)

lambda_best <- lasso_reg$lambda.min 

lasso_model <- glmnet(x,y, alpha = 1, lambda = 5, standardize = TRUE)
c(lasso_model$a0,t(lasso_model$beta))

Y_pred = 25.2888 + 0.4579335*X_R[,2] + 1.0131185*X_R[,4]
plot(1:nrow(X_R),Y_R,type = "o",pch = 20)
lines(1:nrow(X_R),Y_pred,type = "o",pch = 20,col = "red",lty = 2)
abline(v = c(12,20,30,36),lty = 2,col = rgb(1,0,0,alpha = 0.3))

library("L1pack")
rmodel_l <- lad(formula = Y~A1+A2+A3+A4,data = X)
rmodel_l$coefficients

plot(1:nrow(X),X$Y,type = "o",pch = 20,ylab = "observed & predicted values",xlab = "Index")
lines(1:nrow(X),rmodel_l$fitted.values,type = "o",pch = 20,col = "red",lty = 2)
abline(v = 1:nrow(X),lty = 2,col = rgb(0,1,0,alpha = 0.3))
legend("topleft",legend = c("Obs","Pred"),fill = c("black","red"))

rmodel_l <- lqs(Y~A1+A2+A3+A4,data = X,method = "lms")
rmodel_l$coefficients

par(mfrow = c(1,1))
plot(1:nrow(X),X$Y,type = "o",pch = 20,ylab = "observed & predicted values",xlab = "Index")
lines(1:nrow(X),rmodel_l$fitted.values,type = "o",pch = 20,col = "red",lty = 2)
abline(v = 1:nrow(X),lty = 2,col = rgb(0,1,0,alpha = 0.3))
legend("topleft",legend = c("Obs","Pred"),fill = c("black","red"))

rmodel_l<-lqs(Y~A1+A2+A3+A4,data = X,method = "lts")
rmodel_l$coefficients

plot(1:nrow(X),X$Y,type = "o",pch = 20,ylab = "observed & predicted values",xlab = "Index")
lines(1:nrow(X),rmodel_l$fitted.values,type = "o",pch = 20,col = "red",lty = 2)
abline(v = 1:nrow(X),lty = 2,col = rgb(0,1,0,alpha = 0.3))
legend("topleft",legend = c("Obs","Pred"),fill = c("black","red"))

