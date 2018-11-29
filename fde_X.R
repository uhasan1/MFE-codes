fde_X = function(S0=10, strike=10, r=0.04, T=0.5, sigma,m = 100 ,Dt=0.002,type = c("call", "put"), style = c("european", "american")) {
  # m determines the convergance rate, don't set it too much!
  
  #S0=10 ;strike=10 ;r=0.04; T=0.5; sigma= 0.2; m = 100 ;Dt = 0.002; type = "call"; style = "european";
  h = T/Dt

  m = 2*trunc((m+1)/2)                         
  X0 = log(S0)
  Xmin = X0 - 3*sigma*sqrt(T)
  Xmax = X0 + 3*sigma*sqrt(T)
  dX = (Xmax - Xmin)/m
  
  ## Set stock price limits to +/- 3 standard deviations.
  Dt = 0.001
  h = T/Dt
  X.col = seq(from = Xmin, to = Xmax, length.out = m+1)
  f = matrix(NA, nrow=h+1,ncol =  m+1)
  if (type == "call"){
    f[h+1,] <- pmax(exp(X.col)-strike,0)
  }
  else if(type == "put"){
    f[h+1,] <- pmax(strike-exp(X.col),0)
  }
  else {
    cat("wrong input for option type! only 'call' or 'put' please! ")
  }
  
  pu = (1 + r*Dt)^-1 * (-Dt/(2*dX)*(r - 0.5*sigma^2) + Dt/(2*dX^2)*sigma^2)
  pm = (1 + r*Dt)^-1 * (1 - Dt/dX^2*sigma^2)
  pd = (1 + r*Dt)^-1 * (Dt/(2*dX)*(r - 0.5*sigma^2) + Dt/(2*dX^2)*sigma^2)
  
  
  for (i in (h:1)) {             # Iterate from end to beginning.
    fillrow = rep(c(1:3), times=m-1) + rep(0:(m-2), each=3)
    f[i,c(2:m)] = matrix(f[i+1,fillrow], ncol=3, byrow=TRUE) %*% c(pu,pm,pd)
    
    if (type == 'call') {              
      f[i,m+1] = f[i,m] + exp(X.col[m+1]) - exp(X.col[m]) # when ITM
      f[i,1] = f[i,2]        # when OTM
    }
    else if (type == 'put') {            
      f[i,m+1] = f[i,m]      # when OTM
      f[i,1] = f[i,2] - (exp(X.col[m+1]) - exp(X.col[m]))   # when ITM
      if (style == 'american')
        f[i,] = pmax(f[i,], strike - exp(X.col))
    }
  }
   return(f[1, m/2+1])
}
