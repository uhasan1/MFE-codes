price_numerix = function(T=30, WAC=0.08, notional=100000, r0=0.078, k=0.6, r_bar=0.08, sigma=0.12){


# ============ generate path for short rate =========

path= 1000

# Use one day as the time step
n= round(T*360)
dt =T/n

short_rate= matrix(data = 0,nrow = path,ncol = n+1)
short_rate[,1]= r0


z = matrix(data = rnorm(path*(n+1)),nrow = path,ncol = n+1)

# Simulate the dynamics of the short term rate using CIR model
for (i in 1:n){
short_rate[,i+1]= short_rate[,i]+ k*(r_bar - short_rate[,i])*dt+ sigma * sqrt(dt) * sqrt(pmax(short_rate[,i],0))  * z[,i]
}

# Compute the discount rate for each month
#check
disc = array(dim = N)
for (i in 1:N){
R =  dt*rowSums(short_rate[,2:(i*30)+1],dims = 1)
disc[i] = 1/path * sum(exp(-R)) 
}

# sim_10= mean(short_rate(:,10*360+1))
# sim_20= mean(short_rate(:,20*360+1))
# sim_30= mean(short_rate(:,30*360+1))
# 
# expl_10= -1/10*log(zero_bond_CIR_explicit(k, r_bar, mean(short_rate(:,0*360+1)), sigma, 0, 10, 1))
# expl_20= -1/10*log(zero_bond_CIR_explicit(k, r_bar, mean(short_rate(:,10*360+1)), sigma, 0, 10, 1))
# expl_30= -1/10*log(zero_bond_CIR_explicit(k, r_bar, mean(short_rate(:,20*360+1)), sigma, 0, 10, 1))

# ================ Numerix model =====================
  
CPR= matrix(data = 0,nrow = N+1,ncol = 1)
seasonality= c(0.94,0.76,0.74,0.95,0.98,0.92,0.98,1.10,1.18,1.22,1.23,0.98)

PV=  matrix(data = 0,nrow = N+1,ncol = 1)
PV[1]= notional

rate_10y= matrix(data = 0,nrow = N+1,ncol = 1)
rate_10y[1]= -1/10*log(zcbprice_CIR_explicit(r0 = mean(short_rate[,1]),T = 10,sigma = sigma,kappa = k,rbar = r_bar))




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
rate_10y[t+1]= -1/10*log(zcbprice_CIR_explicit(r0 = mean(short_rate[,t*30+1]),T = 10,sigma = sigma,kappa = k,rbar = r_bar))
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
sum_CF= sum_CF + disc[t]*c[t+1]
}

return(sum_CF)


}