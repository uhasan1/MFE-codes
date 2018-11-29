fde_S = function(S0, strike, r, T, sigma, m = 100 ,h = 1000,type = c("call", "put"), style = c("european", "american")) {
  # m determines the convergance rate, don't set it too much!
  
  
  Dt = T / h
  m = 2*trunc((m+1)/2)                         
  
  S0.col = seq(from = 0, to = 2*S0, length.out = m+1)     
  
  f = matrix(data = NA, nrow=h+1,ncol = m+1)
  
  if (type == "call"){
    f[h+1,] = pmax(S0.col - strike, 0)
  }
  else if (type == "put"){
    f[h+1,] = pmax(strike - S0.col, 0)
  }
  else {
    cat("wrong input for option type! only 'call' or 'put' please! ")
  }
  
  
  for (i in h:1) {             
    for (j in (m-1):1) {
      pu = (1 + r*Dt)^-1 * (-0.5*r*j*Dt + 0.5*sigma^2*j^2*Dt)
      pm = (1 + r*Dt)^-1 * (1 - sigma^2*j^2*Dt)
      pd = (1 + r*Dt)^-1 * (0.5*r*j*Dt + 0.5*sigma^2*j^2*Dt)
      f[i,j+1] = t(c(pu,pm,pd)) %*% f[i+1,c(j,j+1,j+2)]
    }
    
    if (type == 'call') {              
      f[i,m+1] = f[i,m] + S0.col[m+1] - S0.col[m] # when deep in the money
      f[i,1] = f[i,2]        #when deep out of money
    }
    else if (type == 'put') {          
      f[i,m+1] = f[i,m]      # when deep out of the money
      f[i,1] = f[i,2] - (S0.col[2] - S0.col[1]) #when deep in money
      if (style == 'american')
        f[i,] = pmax(f[i,], strike - S0.col)
    }
  }
  
 return(f[1, m/2+1])
}

