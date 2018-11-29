fdcn_X = function(S0=10, strike=10, r=0.02, T=2, sigma=0.2,m = 100 ,Dt = 0.002,type = c("call", "put"), style = c("european", "american")) {
  # m determines the convergance rate, don't set it too much!
  
  h = T / Dt
  m = 2*trunc((m+1)/2)                        
  X0 = log(S0)
  Xmin = X0 - 3*sigma*sqrt(T)
  Xmax = X0 + 3*sigma*sqrt(T)
  dX = (Xmax - Xmin)/m
  
  ## Set stock price limits to +/- 3 standard deviations.
  
  X.col = seq(from = Xmin, to = Xmax, length.out = m+1)
  f = matrix(NA, nrow=h+1,ncol =  m+1)
 
  if (type == "call"){
    f[h+1,] = pmax(exp(X.col)-strike,0)
  }
  else if (type == "put"){
    f[h+1,] = pmax(strike-exp(X.col),0)
  }
  else {
    cat("wrong input for option type! only 'call' or 'put' please! ")
  }
  
  
  pu = -Dt*0.25*(sigma^2/dX^2 + (r - 0.5*sigma^2)/dX)
  pm = 1 + Dt*(sigma^2/(2*dX^2) + r/2)
  pd = -Dt*0.25*(sigma^2/dX^2 - (r - 0.5*sigma^2)/dX)
  
  diagmat = rbind(0, cbind(diag(c(rep(pu, m-1), 1)), 0)) +
    diag(c(1, rep(pm, m-1), -1)) +
    rbind(cbind(0, diag(c(-1, rep(pd, m-1)))), 0)
  c = diagmat
  d = diag(2, m+1) - diagmat           
  for (i in (h:1)) {             # Iterate from last col to first col
    resmat = d %*% f[i+1,(m+1):1]
    if (type == 'call') {
      resmat[1] = exp(X.col[m+1]) - exp(X.col[m])   # boundary conditions
      resmat[m+1] = 0   # boundary conditions
    }
    else if (type == 'put') {
      resmat[1] = 0 # boundary conditions
      resmat[m+1] = exp(X.col[m+1]) - exp(X.col[m]) # bounary conditions
    }
    f[i,(m+1):1] = solve(c, resmat)
    
    if (type == 'put' && style == 'american')
      f[i,] = pmax(f[i,], strike - exp(X.col))
  }
  
  return(f[1, m/2+1])
}
