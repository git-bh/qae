## OLS estimation using lm and martix (Rbasic4 파일의 부분을 보완)


#first, we generate x and y
n=1000; 
b0=0.1; b1=1 #true parameter

y=matrix(0,n,1)
x=matrix(rnorm(n),n,1)
eps = matrix(rnorm(n),n,1)

for (i in 1:n){
  y[i]=b0+b1*x[i]+eps[i]
}

plot(x,y)
plot(y~x)

# Estimate linear model with R default function
fit1 = lm(y~x)
summary(fit1)
fit1

# 95% confidence interval
confint(fit1, level=0.95)

# Plotting regression line
abline(fit1)

# fitted value of y
yhat = fitted(fit1)  
summary(yhat)
plot(yhat~x)

# residuals
uhat = resid(fit1)  
summary(uhat)
plot(uhat~x)


####################################
# OLS estimation using matrix
####################################

### Matrix facilities ###
# diag(k) : k x k identity matrix
# entry-wise multiplication : *
# matrix multiplication : %*%
# entry-wise division : /
# transpose : t(A)
# inverse : solve(A) if A is square matrix
# cross product(A'B) : crossprod(A,B)
# outer product(AB') : A%o%B


x=cbind(1,x) #the first column in x corresponds to intercept
head(x)

invx = solve(crossprod(x))    # (x'x)^(-1)
# 또는 invx = solve(t(x)%*%x)

olsb = invx%*%crossprod(x,y)  # (x'x)^(-1)x'y

olsb  #OLS estimate, compare this with fit1

#fit2 = lm(y~x-1) # regression without intercept
#summary(fit2)

u    = y-x%*%olsb   #residual
sig2 = crossprod(u)/(n-2)
varcov = matrix(sig2,2,2) * invx
varcov
se = rbind(sqrt(varcov[1,1]), sqrt(varcov[2,2]))  #standard error
# 또는 se = sqrt(diag(varcov)) 
tvalue = olsb/se # t-statistic
pvalue = 2*(1-pt(abs(tvalue),df=n-2)) # p-value
cbind(olsb, se, tvalue, pvalue)


########################
multiple regression case
########################

data = as.matrix(read.csv("DataHousingPrice.csv",header=T))

y = cbind(data[,1])   # or y = matrix(data[,1])
x = cbind(1,data[,2:4])
head(x)

n    = nrow(data)
k    = ncol(x) -1 

invx = solve(t(x)%*%x)
olsb = invx%*%t(x)%*%y
olsb

u    = y-x%*%olsb   #residual
sig2 = crossprod(u)/(n-k-1)
varcov = matrix(sig2,k+1,k+1) * invx
varcov

se = sqrt(diag(varcov)) #standard error
tvalue = olsb/se # t-statistic
pvalue = 2*(1-pt(abs(tvalue),df=n-k-1)) # p-value
cbind(olsb, se, tvalue, pvalue)


