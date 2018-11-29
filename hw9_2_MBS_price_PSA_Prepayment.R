price_PSA = function(T, WAC, notional, r0, k, r_bar, sigma){

N= T*12

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
  sum_CF= sum_CF + disc[t]*c[t+1]
}

return(sum_CF)
}