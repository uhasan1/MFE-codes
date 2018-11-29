### I have assumes 252 for the day count in a year. (trading days) ###

# I have calculated all of the functions based on a face value of 1 dollar for future consistency, 
#in this case, K should be put as the percentage of Spot. 
#Then easily can multiply the final answer, either for bond or option, in the face value. 

r0 = 0.05
sigma = 0.1
rbar = 0.05
kappa = 0.82

n = 100    # number of simulations
T = 0.5

##### zcbprice_vasicek_MC #####
zcbprice_vasicek_MC = function(r0,T,sigma,kappa,rbar,n){

dt = 1/252 # time step is assumed to be a day
h = round(T/dt) # number of intervals

r = matrix(data = NA ,n,h+1)  # matrix to hold short rate paths
r[,1] = r0

W = matrix(data = rnorm(n*(h+1)),nrow = n,ncol = h+1) 
# the interest rate path
  for(j in 1:h){
    dr = kappa*(rbar-r[,j])*dt + sigma*sqrt(dt)*W[,j]
    r[,j+1] = r[,j] + dr
  }

P.0.T = mean(exp(-dt*apply(X = r,MARGIN = 1, FUN = sum)))

return((P.0.T))
}

##### zcbprice_vasicek_explicit #####
zcbprice_vasicek_explicit = function(r0,T,sigma,kappa,rbar){
    BT = (1/kappa)*(1-exp(-T*kappa)) 
    AT = exp((rbar-sigma^2/(2*kappa^2))*(BT-T)-(sigma^2)/(4*kappa)*BT^2)
    return(AT*exp(-BT*r0))
  }

zcbprice_vasicek_explicit(r0 = r0,T = T,sigma = sigma,kappa = kappa,rbar = rbar)
zcbprice_vasicek_MC(r0 = r0,T = T,sigma = sigma,kappa = kappa,rbar = rbar,n = 100)

##### coupon paying bond vasicek MC #####
# pricing a coupon paying bond with c/2 semi-annually and 1+c/2 at maturity
cbprice_vasicek_MC = function(maturity, coupon,r0,sigma,kappa,rbar,n,T_firstcoupon){
  cbprice_vasicek_MC = 0
  for (T in seq(T_firstcoupon,maturity,by = 0.5)){
    cbprice_vasicek_MC = cbprice_vasicek_MC + (coupon/2)*zcbprice_vasicek_MC(r0 = r0,sigma = sigma, kappa = kappa, rbar = rbar, n = n, T = T)
  }
  cbprice_vasicek_MC = cbprice_vasicek_MC + zcbprice_vasicek_MC(r0 = r0,sigma = sigma, kappa = kappa, rbar = rbar, n = n, T = T) #last payment
  return(cbprice_vasicek_MC)
}

##### coupon paying bond vasicek MC #####
cbprice_vasicek_explicit = function(maturity, coupon,r0,sigma,kappa,rbar,T_firstcoupon){
  cbprice_vasicek_explicit = 0
  for (T in seq(T_firstcoupon,maturity,by = 0.5)){
    cbprice_vasicek_explicit = cbprice_vasicek_explicit + (coupon/2)*zcbprice_vasicek_explicit(r0 = r0,sigma = sigma, kappa = kappa, rbar = rbar, T = T)
  }
  cbprice_vasicek_explicit = cbprice_vasicek_explicit + zcbprice_vasicek_explicit(r0 = r0,sigma = sigma, kappa = kappa, rbar = rbar, T = T) #last payment
  return(cbprice_vasicek_explicit)
}

cbprice_vasicek_MC(maturity = 4,coupon = 0.06,r0 = r0,sigma = sigma, kappa = kappa,rbar = rbar,T_firstcoupon = 0.5,n=100)
cbprice_vasicek_explicit(maturity = 4,coupon = 0.06,r0 = r0,sigma = sigma, kappa = kappa,rbar = rbar,T_firstcoupon = 0.6)


##### call on zero-coupon bond vasicek/option:MC/ bond:explicit #####
call_vasicek_MC_explicit = function(T, S, K, r0, sigma, kappa, rbar, n){
  if (S < T){
    cat ("bond maturity cannot be less than option maturity")
    return()
  }
  
  dt = 1/252 # time step is assumed to be a day
  h = round(T/dt) # number of intervals
  
  r1 = matrix(data = NA ,n,h+1)  # matrix to hold short rate paths
  r1[,1] = r0
  
  W = matrix(data = rnorm(n*(h+1)),nrow = n,ncol = h+1) 
  # the interest rate path
  for(j in 1:h){
    dr = kappa*(rbar-r1[,j])*dt + sigma*sqrt(dt)*W[,j]
    r1[,j+1] = r1[,j] + dr
  }
  
  zcbprice_vasicek_explicits = matrix(nrow=n,ncol = 1)
  
  for (i in 1:n){
    zcbprice_vasicek_explicits[i,1] = zcbprice_vasicek_explicit(r0 = r1[i,ncol(r1)],T = S-T,sigma = sigma,kappa = kappa,rbar = rbar)  
  }
  
  output = zcbprice_vasicek_explicits - K
  output_positive=(0<output)*output
  call.discounted = mean(exp(-dt*apply(X = r1,MARGIN = 1, FUN = sum))*output_positive)
  return(call.discounted)
}
  
call_vasicek_MC_explicit(T = 0.25,S = 0.5,r0 = r0,sigma = sigma,kappa = kappa,rbar = rbar,n = 100000,K = 0.98)
  
##### call_vasicek_MC_MC #####
call_vasicek_MC_MC = function(T, S, K, r0, sigma, kappa, rbar, n,m){
  if (S < T){
    cat ("bond maturity cannot be less than option maturity")
    return()
  }
  
  dt = 1/252 # time step is assumed to be a day
  h = round(T/dt) # number of intervals
  
  r1 = matrix(data = NA ,n,h+1)  # matrix to hold short rate paths
  r1[,1] = r0
  
  W = matrix(data = rnorm(n*(h+1)),nrow = n,ncol = h+1) 
  # the interest rate path
  for(j in 1:h){
    dr = kappa*(rbar-r1[,j])*dt + sigma*sqrt(dt)*W[,j]
    r1[,j+1] = r1[,j] + dr
  }
  
  cbprice_vasicek_MC = cbprice_vasicek_MC(maturity = S-T, T_firstcoupon = 0.25,coupon = 0.06, r0 = mean(r1[,ncol(r1)]), sigma = sigma,kappa = kappa,rbar = rbar,n=m)
  output = cbprice_vasicek_MC - K
  output_positive=(0<output)*output
  call.discounted = mean(exp(-dt*apply(X = r1,MARGIN = 1, FUN = sum)))*output_positive
  return(mean(call.discounted))
}
  
call_vasicek_MC_MC(T = 0.25,S = 4,K = 0.98,r0 = r0,sigma = sigma,kappa = kappa,rbar = rbar,n = 1000,m = 1000)

##### call_vasicek_MC_explicit_cb #####
call_vasicek_MC_explicit_cb = function(T, S, K, r0, sigma, kappa, rbar, n){
  if (S < T){
    cat ("bond maturity cannot be less than option maturity")
    return()
  }
  
  dt = 1/252 # time step is assumed to be a day
  h = round(T/dt) # number of intervals
  
  r1 = matrix(data = NA ,n,h+1)  # matrix to hold short rate paths
  r1[,1] = r0
  
  W = matrix(data = rnorm(n*(h+1)),nrow = n,ncol = h+1) 
  # the interest rate path
  for(j in 1:h){
    dr = kappa*(rbar-r1[,j])*dt + sigma*sqrt(dt)*W[,j]
    r1[,j+1] = r1[,j] + dr
  }
  
  cbprice_vasicek_explicit = cbprice_vasicek_explicit(maturity = S-T, T_firstcoupon = 0.25,coupon = 0.06, r0 = mean(r1[,ncol(r1)]), sigma = sigma,kappa = kappa,rbar = rbar)
  output = cbprice_vasicek_explicit - K
  output_positive=(0<output)*output
  call.discounted = mean(exp(-dt*apply(X = r1,MARGIN = 1, FUN = sum)))*output_positive
  return(mean(call.discounted))
}


zcbprice_vasicek_explicit(r0 = r0,T = T,sigma = sigma,kappa = kappa,rbar = rbar)
zcbprice_vasicek_MC(r0 = r0,T = T,sigma = sigma,kappa = kappa,rbar = rbar,n = 100)

cbprice_vasicek_MC(maturity = 4,coupon = 0.06,r0 = r0,sigma = sigma, kappa = kappa,rbar = rbar,T_firstcoupon = 0.5,n=100)
cbprice_vasicek_explicit(maturity = 4,coupon = 0.06,r0 = r0,sigma = sigma, kappa = kappa,rbar = rbar,T_firstcoupon = 0.6)

call_vasicek_MC_MC(T = 0.25,S = 4,K = 0.98,r0 = r0,sigma = sigma,kappa = kappa,rbar = rbar,n = 1000,m = 1000)
call_vasicek_MC_explicit(T = 0.25,S = 0.5,r0 = r0,sigma = sigma,kappa = kappa,rbar = rbar,n = 100000,K = 0.98)
call_vasicek_MC_explicit_cb(T = 0.25,S = 4,K = 0.98,r0 = r0,sigma = sigma,kappa = kappa,rbar = rbar,n = 100000)

# all the answers are pretty much the same. one of th emajor sources of difference can be due to trading days. 
# increasing the number of simulation will enhance the results but with the cost of more time spent.
