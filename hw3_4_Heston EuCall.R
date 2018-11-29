myHestonReflection <- function(S0,K,r,T,sigma,beta,alpha,rho,v0,h,n){

Dt = T/h
S= matrix(nrow = n,ncol = h+1)
v= matrix(nrow = n,ncol = h+1)
B1 = matrix(data = rnorm(n*h),nrow = n,ncol = h)
B2 = matrix(data = rnorm(n*h),nrow = n,ncol = h)
W1= B1
W2= (rho*B1)+(sqrt(1-rho^2)*B2)
S[,1]=S0
v[,1]=v0

  for (j in 1:h){
    v [,j+1] = abs(v[,j])+alpha*(beta-abs(v[,j]))*Dt+sigma*sqrt(abs(v[,j]))*W2[,j]*sqrt(Dt)
    S [,j+1] = S[,j]+ r*S[,j]*Dt + (S[,j])*sqrt(abs(v[,j]))*W1[,j]*sqrt(Dt)
  }


myHestonReflection = mean(exp(-r*T)*(S[,h+1]>K)*(S[,h+1]-K))
myHestonReflection
}

myHestonPatialTrunc <- function(S0,K,r,T,sigma,beta,alpha,rho,v0,h,n){
  
  Dt = T/h
  S= matrix(nrow = n,ncol = h+1)
  v= matrix(nrow = n,ncol = h+1)
  B1 = matrix(data = rnorm(n*h),nrow = n,ncol = h)
  B2 = matrix(data = rnorm(n*h),nrow = n,ncol = h)
  W1= B1
  W2= (rho*B1)+(sqrt(1-rho^2)*B2)
  S[,1]=S0
  v[,1]=v0
  
    for (j in 1:h){
      v [,j+1] = v[,j] + alpha*(beta-v[,j])*Dt + sigma*sqrt(pmax(v[,j],0))*W2[,j]*sqrt(Dt)
      S [,j+1] = S[,j] + r*S[,j]*Dt + (S[,j])*sqrt(pmax(v[,j],0))*W1[,j]*sqrt(Dt)
    }
  
  
  myHestonPatialTrunc = mean(exp(-r*T)*(S[,h+1]>K)*(S[,h+1]-K))
  myHestonPatialTrunc
}

myHestonFullTrunc <- function(S0,K,r,T,sigma,beta,alpha,rho,v0,h,n){
  
  Dt = T/h
  S= matrix(nrow = n,ncol = h+1)
  v= matrix(nrow = n,ncol = h+1)
  B1 = matrix(data = rnorm(n*h),nrow = n,ncol = h)
  B2 = matrix(data = rnorm(n*h),nrow = n,ncol = h)
  W1= B1
  W2= (rho*B1)+(sqrt(1-rho^2)*B2)
  S[,1]=S0
  v[,1]=v0
  
    for (j in 1:h){
      v [,j+1] = v[,j]+alpha*(beta-pmax(v[,j],0))*Dt+sigma*sqrt(pmax(v[,j],0))*W2[,j]*sqrt(Dt)
      S [,j+1] = S[,j]+ r*S[,j]*Dt + (S[,j])*sqrt(pmax(v[,j],0))*W1[,j]*sqrt(Dt)
    }
  
  
  myHestonFullTrunc = mean(exp(-r*T)*(S[,h+1]>K)*(S[,h+1]-K))
  myHestonFullTrunc
}

myHestonReflection(S0 = 48,K = 50,r = 0.03,T = 0.5,sigma = 0.42,beta = 0.0625,alpha = 5.8,rho = -.6,v0 = 0.05,h = 1000,n = 1000)
myHestonPatialTrunc(S0 = 48,K = 50,r = 0.03,T = 0.5,sigma = 0.42,beta = 0.0625,alpha = 5.8,rho = -.6,v0 = 0.05,h = 1000,n = 1000)
myHestonFullTrunc(S0 = 48,K = 50,r = 0.03,T = 0.5,sigma = 0.42,beta = 0.0625,alpha = 5.8,rho = -.6,v0 = 0.05,h = 1000,n = 1000)
