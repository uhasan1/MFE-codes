# If you want a quicker but less precise result lower path number # 

############ functions ##########

zcbprice_CIR_explicit = function(r0,T,sigma,kappa,rbar){
  h1= sqrt(kappa^2+2*sigma^2)
  h2= (kappa+h1)/2
  h3= (2*kappa*rbar)/sigma^2
  
  A= ((h1*exp(h2*T))/(h2*(exp(h1*T)-1)+h1))^h3
  B= (exp(h1*T)-1)/(h2*(exp(h1*T)-1)+h1)
  
  P.0.T = A*exp(-B*r0)
  return(P.0.T)
}

price_numerix = function(T=30, WAC=0.08, notional=100000, r0=0.078, kappa=0.6, rbar=0.08, sigma=0.12){
  r = WAC/12
  N = T*12
  
  path= 500
  
  # Use one day as the time step
  n= round(T*365)
  dt =T/n
  
  rate= matrix(data = 0,nrow = path,ncol = n+1)
  rate[,1]= r0
  
  
  z_norm = matrix(data = rnorm(path*(n+1)),nrow = path,ncol = n+1)
  
  # Simulate short term rate using CIR model
  for (i in 1:n){
    rate[,i+1]= rate[,i]+ kappa*(rbar - rate[,i])*dt+ sigma * sqrt(dt) * sqrt(pmax(rate[,i],0))  * z_norm[,i]
  }
  
  # Compute the discount rate for each month
  #check
  DF = array(dim = N)
  for (i in 1:N){
    R =  dt*rowSums(rate[,2:(i*30)+1],dims = 1)
    DF[i] = 1/path * sum(exp(-R)) 
  }
  
  CPR= matrix(data = 0,nrow = N+1,ncol = 1)
  seasonality= c(0.94,0.76,0.74,0.95,0.98,0.92,0.98,1.10,1.18,1.22,1.23,0.98)
  
  PV=  matrix(data = 0,nrow = N+1,ncol = 1)
  PV[1]= notional
  
  rate_10y= matrix(data = 0,nrow = N+1,ncol = 1)
  rate_10y[1]= -1/10*log(zcbprice_CIR_explicit(r0 = mean(rate[,1]),T = 10,sigma = sigma,kappa = kappa,rbar = rbar))
  
  
  
  
  sum_CF= 0
  RI = array(data = 0,dim = N)
  BU = array(data = 0,dim = N)
  SG = array(data = 0,dim = N)
  SY = array(data = 0,dim = N)
  CPR = array(data = 0,dim = N)
  SP = array(data = 0,dim = N)
  IP = array(data = 0,dim = N)
  PP = array(data = 0,dim = N)
  c = array(data = 0,dim = N)
  for (t in 1:N){
    rate_10y[t+1]= -1/10*log(zcbprice_CIR_explicit(r0 = mean(rate[,t*30+1]),T = 10,sigma = sigma,kappa = kappa,rbar = rbar))
    RI[t+1]= 0.28+ 0.14*atan(-8.57+ 430*(12*r - rate_10y[t]))
    BU[t+1]= 0.3 + 0.7*(PV[t]/PV[1])
    SG[t+1]= min(1,(t-1)/30)
    SY[t+1]= seasonality[t%%12 +((t%%12)==0)*12]
    CPR[t+1]= RI[t+1]*BU[t+1]*SG[t+1]*SY[t+1]
    SP[t+1]= PV[t]*r*(1/(1-(1/((1+r)^(N-t+2))))-1)
    IP[t+1]= PV[t]*r
    PP[t+1]= (PV[t]-SP[t+1])*(1-(1-CPR[t+1])^(1/12))
    c [t+1]= SP[t+1]+PP[t+1]+IP[t+1]
    PV[t+1]= PV[t]-(SP[t+1]+PP[t+1])
    sum_CF= sum_CF + DF[t]*c[t+1]
  }
  
  return(sum_CF)
  
  
}


price_PSA = function(T, WAC, notional, r0, kappa, rbar, sigma){
  
  N= T*12
  
  path= 500
  
  # Use one day as the time step
  n= round(T*365)
  dt =T/n
  
  rate= matrix(data = 0,nrow = path,ncol = n+1)
  rate[,1]= r0
  
  
  z_norm = matrix(data = rnorm(path*(n+1)),nrow = path,ncol = n+1)
  
  # Simulate short term rate using CIR model
  for (i in 1:n){
    rate[,i+1]= rate[,i]+ kappa*(rbar - rate[,i])*dt+ sigma * sqrt(dt) * sqrt(pmax(rate[,i],0))  * z_norm[,i]
  }
  
  # Compute the discount rate for each month
  #check
  DF = array(dim = N)
  for (i in 1:N){
    R =  dt*rowSums(rate[,2:(i*30)+1],dims = 1)
    DF[i] = 1/path * sum(exp(-R)) 
  }
  
  CPR= matrix(data = 0,nrow = N+1,ncol = 1)
  CPR[1:31]= seq(from = 0, to = 0.06, by = 0.002)
  CPR[31:N+1]= CPR[31]
  
  r= WAC/12
  
  PV=  matrix(data = 0,nrow = N+1,ncol = 1)
  PV[1]= notional
  
  sum_CF= 0
  
  SP = array(data = 0,dim = N)
  IP = array(data = 0,dim = N)
  PP = array(data = 0,dim = N)
  c = array(data = 0,dim = N)
  
  for (t in 1:N){
    SP[t+1]= PV[t]*r*(1/(1-(1/((1+r)^(N-t+2))))-1)
    IP[t+1]= PV[t]*r
    PP[t+1]= (PV[t]-SP[t+1])*(1-(1-CPR[t+1])^(1/12))
    c [t+1]= SP[t+1]+PP[t+1]+IP[t+1]
    PV[t+1]= PV[t]-(SP[t+1]+PP[t+1])
    sum_CF= sum_CF + DF[t]*c[t+1]
  }
  
  return(sum_CF)
}

find_oas = function(oas){
  
  T= 30
  WAC= 0.08
  notional = 100000
  r0= 0.078
  kappa= 0.6
  rbar= 0.08
  sigma= 0.12
  
  r = WAC/12
  N = T*12
  
  path= 100
  
  # Use one day as the time step
  n= round(T*365)
  dt =T/n
  
  rate= matrix(data = 0,nrow = path,ncol = n+1)
  rate[,1]= r0
  
  
  z_norm = matrix(data = rnorm(path*(n+1)),nrow = path,ncol = n+1)
  
  # Simulate short term rate using CIR model
  for (i in 1:n){
    rate[,i+1]= rate[,i]+ kappa*(rbar - rate[,i])*dt+ sigma * sqrt(dt) * sqrt(pmax(rate[,i],0))  * z_norm[,i]
  }
  
  # Compute the discount rate for each month
  #check
  DF = array(dim = N)
  for (i in 1:N){
    R =  dt*rowSums(rate[,2:(i*30)+1],dims = 1)
    DF[i] = 1/path * sum(exp(-R)) 
  }
  
  CPR= matrix(data = 0,nrow = N+1,ncol = 1)
  seasonality= c(0.94,0.76,0.74,0.95,0.98,0.92,0.98,1.10,1.18,1.22,1.23,0.98)
  
  PV=  matrix(data = 0,nrow = N+1,ncol = 1)
  PV[1]= notional
  
  rate_10y= matrix(data = 0,nrow = N+1,ncol = 1)
  rate_10y[1]= -1/10*log(zcbprice_CIR_explicit(r0 = mean(rate[,1]),T = 10,sigma = sigma,kappa = kappa,rbar = rbar))
  
  disc_adj = array(dim = N)
  for (i in 1:N){
    R_adj =  dt*(rowSums(rate[,2:(i*30)+1],dims = 1) + oas*length(rate[2:(i*30)+1]))
    disc_adj[i] = 1/path * sum(exp(-R_adj)) 
  }
  
  
  sum_CF= 0
  RI = array(data = 0,dim = N)
  BU = array(data = 0,dim = N)
  SG = array(data = 0,dim = N)
  SY = array(data = 0,dim = N)
  CPR = array(data = 0,dim = N)
  SP = array(data = 0,dim = N)
  IP = array(data = 0,dim = N)
  PP = array(data = 0,dim = N)
  c = array(data = 0,dim = N)
  for (t in 1:N){
    rate_10y[t+1]= -1/10*log(zcbprice_CIR_explicit(r0 = mean(rate[,t*30+1]),T = 10,sigma = sigma,kappa = kappa,rbar = rbar))
    RI[t+1]= 0.28+ 0.14*atan(-8.57+ 430*(12*r - rate_10y[t]))
    BU[t+1]= 0.3 + 0.7*(PV[t]/PV[1])
    SG[t+1]= min(1,(t-1)/30)
    SY[t+1]= seasonality[t%%12 +((t%%12)==0)*12]
    CPR[t+1]= RI[t+1]*BU[t+1]*SG[t+1]*SY[t+1]
    SP[t+1]= PV[t]*r*(1/(1-(1/((1+r)^(N-t+2))))-1)
    IP[t+1]= PV[t]*r
    PP[t+1]= (PV[t]-SP[t+1])*(1-(1-CPR[t+1])^(1/12))
    c [t+1]= SP[t+1]+PP[t+1]+IP[t+1]
    PV[t+1]= PV[t]-(SP[t+1]+PP[t+1])
    sum_CF= sum_CF +  disc_adj[t]*c[t+1]
  }
  
  
  return(sum_CF)
  
  
}

price_numerix_with_oas = function(T=30, WAC=0.08, notional=100000, r0=0.078, kappa=0.6, rbar=0.08, sigma=0.12,oas){
  
  r = WAC/12
  N = T*12
  

  path= 500
  
  # Use one day as the time step
  n= round(T*365)
  dt =T/n
  
  rate= matrix(data = 0,nrow = path,ncol = n+1)
  rate[,1]= r0
  
  
  z_norm = matrix(data = rnorm(path*(n+1)),nrow = path,ncol = n+1)
  
  # Simulate short term rate using CIR model
  for (i in 1:n){
    rate[,i+1]= rate[,i]+ kappa*(rbar - rate[,i])*dt+ sigma * sqrt(dt) * sqrt(pmax(rate[,i],0))  * z_norm[,i]
  }
  
  CPR= matrix(data = 0,nrow = N+1,ncol = 1)
  seasonality= c(0.94,0.76,0.74,0.95,0.98,0.92,0.98,1.10,1.18,1.22,1.23,0.98)
  
  PV=  matrix(data = 0,nrow = N+1,ncol = 1)
  PV[1]= notional
  
  rate_10y= matrix(data = 0,nrow = N+1,ncol = 1)
  rate_10y[1]= -1/10*log(zcbprice_CIR_explicit(r0 = mean(rate[,1]),T = 10,sigma = sigma,kappa = kappa,rbar = rbar))
  
  disc_adj = array(dim = N)
  for (i in 1:N){
    R_adj =  dt*(rowSums(rate[,2:(i*30)+1],dims = 1) + oas*length(rate[2:(i*30)+1]))
    disc_adj[i] = 1/path * sum(exp(-R_adj)) 
  }
  
  
  sum_CF= 0
  RI = array(data = 0,dim = N)
  BU = array(data = 0,dim = N)
  SG = array(data = 0,dim = N)
  SY = array(data = 0,dim = N)
  CPR = array(data = 0,dim = N)
  SP = array(data = 0,dim = N)
  IP = array(data = 0,dim = N)
  PP = array(data = 0,dim = N)
  c = array(data = 0,dim = N)
  for (t in 1:N){
    rate_10y[t+1]= -1/10*log(zcbprice_CIR_explicit(r0 = mean(rate[,t*30+1]),T = 10,sigma = sigma,kappa = kappa,rbar = rbar))
    RI[t+1]= 0.28+ 0.14*atan(-8.57+ 430*(12*r - rate_10y[t]))
    BU[t+1]= 0.3 + 0.7*(PV[t]/PV[1])
    SG[t+1]= min(1,(t-1)/30)
    SY[t+1]= seasonality[t%%12 +((t%%12)==0)*12]
    CPR[t+1]= RI[t+1]*BU[t+1]*SG[t+1]*SY[t+1]
    SP[t+1]= PV[t]*r*(1/(1-(1/((1+r)^(N-t+2))))-1)
    IP[t+1]= PV[t]*r
    PP[t+1]= (PV[t]-SP[t+1])*(1-(1-CPR[t+1])^(1/12))
    c [t+1]= SP[t+1]+PP[t+1]+IP[t+1]
    PV[t+1]= PV[t]-(SP[t+1]+PP[t+1])
    sum_CF= sum_CF +  disc_adj[t]*c[t+1]
  }
  
  return(sum_CF)
  
  
}

price_numerix_IO_PO = function(T=30, WAC=0.08, notional=100000, r0=0.078, kappa=0.6, rbar=0.08, sigma=0.12){
  
  
  r = WAC/12
  N = T*12
  
  path= 500
  
  # Use one day as the time step
  n= round(T*365)
  dt =T/n
  
  rate= matrix(data = 0,nrow = path,ncol = n+1)
  rate[,1]= r0
  
  
  z_norm = matrix(data = rnorm(path*(n+1)),nrow = path,ncol = n+1)
  
  # Simulate short term rate using CIR model
  for (i in 1:n){
    rate[,i+1]= rate[,i]+ kappa*(rbar - rate[,i])*dt+ sigma * sqrt(dt) * sqrt(pmax(rate[,i],0))  * z_norm[,i]
  }
  
  # Compute the discount rate for each month
  DF = array(dim = N)
  for (i in 1:N){
    R =  dt*rowSums(rate[,2:(i*30)+1],dims = 1)
    DF[i] = 1/path * sum(exp(-R)) 
  }
  
  CPR= matrix(data = 0,nrow = N+1,ncol = 1)
  seasonality= c(0.94,0.76,0.74,0.95,0.98,0.92,0.98,1.10,1.18,1.22,1.23,0.98)
  
  PV=  matrix(data = 0,nrow = N+1,ncol = 1)
  PV[1]= notional
  
  rate_10y= matrix(data = 0,nrow = N+1,ncol = 1)
  rate_10y[1]= -1/10*log(zcbprice_CIR_explicit(r0 = mean(rate[,1]),T = 10,sigma = sigma,kappa = kappa,rbar = rbar))
  
  
  
  
  sum_CF_IO = 0
  sum_CF_PO = 0
  RI = array(data = 0,dim = N)
  BU = array(data = 0,dim = N)
  SG = array(data = 0,dim = N)
  SY = array(data = 0,dim = N)
  CPR = array(data = 0,dim = N)
  SP = array(data = 0,dim = N)
  IP = array(data = 0,dim = N)
  PP = array(data = 0,dim = N)
  c = array(data = 0,dim = N)
  for (t in 1:N){
    rate_10y[t+1]= -1/10*log(zcbprice_CIR_explicit(r0 = mean(rate[,t*30+1]),T = 10,sigma = sigma,kappa = kappa,rbar = rbar))
    RI[t+1]= 0.28+ 0.14*atan(-8.57+ 430*(12*r - rate_10y[t]))
    BU[t+1]= 0.3 + 0.7*(PV[t]/PV[1])
    SG[t+1]= min(1,(t-1)/30)
    SY[t+1]= seasonality[t%%12 +((t%%12)==0)*12]
    CPR[t+1]= RI[t+1]*BU[t+1]*SG[t+1]*SY[t+1]
    SP[t+1]= PV[t]*r*(1/(1-(1/((1+r)^(N-t+2))))-1)
    IP[t+1]= PV[t]*r
    PP[t+1]= (PV[t]-SP[t+1])*(1-(1-CPR[t+1])^(1/12))
    PV[t+1]= PV[t]-(SP[t+1]+PP[t+1])
    sum_CF_IO= sum_CF_IO + DF[t]*IP[t+1]
    sum_CF_PO= sum_CF_PO + DF[t]*(SP[t+1]+PP[t+1])
  }
  
  return(list(IO = sum_CF_IO, PO = sum_CF_PO))
  
}

############# writeup ###############

################# Q1 ################

price_1a = price_numerix(T=30, WAC=0.08, notional=100000, r0=0.078, kappa=0.6, rbar=0.08, sigma=0.12 )

kappa= seq(from = 0.3, to = 0.9, by = 0.1) 
price_1b = array(dim = length(kappa))
for (i in 1:length(kappa)){
  price_1b [i] = price_numerix(T=30, WAC=0.08, notional=100000, r0=0.078, kappa = kappa[i], rbar=0.08, sigma=0.12)
}

plot(kappa,price_1b,main = 'MBS Prices using Numerix Model',xlab = 'kappa', ylab ='Price',type = 'l',col = 2)

price_1c = array(dim = length(kappa))
kappa= 0.6
rbar= seq(from= 0.03, to = 0.09, by=0.01)
for (i in 1:length(rbar)){
  price_1c [i] = price_numerix(T=30, WAC=0.08, notional=100000, r0=0.078, kappa=0.6, rbar = rbar[i], sigma=0.12)
}

plot(x = rbar, y = price_1c, type = 'l',main = 'MBS Prices using Numerix Model' ,xlab = 'rbar', ylab ='Price', col = 3)

################# Q2 ################

price_2a = price_PSA(T=30, WAC=0.08, notional=100000, r0=0.078, kappa = 0.6, rbar=0.08, sigma=0.12)

kappa= seq(from = 0.3, to = 0.9, by = 0.1) 
price_2b = array(dim = length(kappa))
for (i in 1:length(kappa)){
  price_2b [i] = price_PSA(T=30, WAC=0.08, notional=100000, r0=0.078, kappa = kappa[i], rbar=0.08, sigma=0.12)
}

plot(kappa,price_2b,type = 'l',main = 'MBS Prices using PSA Model',xlab = 'kappa', ylab ='Price', col = 4)

################# Q3 ################


# I lowered the simulation steps to 100 for faster but not so concise result
# it still takes so much time
oasroot = function(X){
  Mortgage_value = 110000
  a = find_oas(X)
  return(a - Mortgage_value)
}
OAS = uniroot(oasroot,c(-0.1,0.1))
OAS = OAS$root

# OAS = -0.01268422

################# Q4 ################


y= 0.0005
P_pos = price_numerix_with_oas(T=30, WAC=0.08, notional=100000, r0=0.078, kappa = 0.6, rbar=0.08, sigma=0.12, oas = OAS-y)
P_neg = price_numerix_with_oas(T=30, WAC=0.08, notional=100000, r0=0.078, kappa = 0.6, rbar=0.08, sigma=0.12, oas = OAS+y)
P0 = price_numerix_with_oas(T=30, WAC=0.08, notional=100000, r0=0.078, kappa = 0.6, rbar=0.08, sigma=0.12, oas =0)

OAS_duration = (P_neg - P_pos)/(2*y*P0)
OAS_convexity = (P_pos + P_neg - 2*P0)/(2*P0*y^2)

################# Q5 ################

kappa= 0.6
rbar= seq(from= 0.03, to = 0.09, by=0.01)

IO = array(dim = length(rbar))
PO = array(dim = length(rbar))
for (i in 1:length(rbar)){
output = price_numerix_IO_PO(T=30, WAC=0.08, notional=100000, r0=0.078, kappa = 0.6, sigma=0.12, rbar = rbar[i])
IO[i] = output$IO
PO[i] = output$PO
}

plot(rbar,IO, type = 'l', main ='IO Tranch Prices using Numerix Model',xlab ='rbar',ylab = 'Price')
plot(rbar,PO, type = 'l', main ='PO Tranch Prices using Numerix Model',xlab ='rbar',ylab = 'Price')


