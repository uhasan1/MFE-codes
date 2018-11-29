#comments are provided at the end of this program

# use the format below to get whatever information you want! You just need to change T,S,method
# my_Put_LSMC(S0 = 36, sigma = 0.2,strike = 40,r = 0.06,delta = 0,T = 1,n = 100000,h = 100,terms = 4,method = "monomials")
# my functions provides the price as well as upper and lower bounds of confidence interval. to use the price only add $price in the end 
my_Put_LSMC(S0 = 36, sigma = 0.2,strike = 40,r = 0.06,delta = 0,T = 1,n = 100000,h = 50,terms = 2,method = "laguerre")

# this one was an example


# constant factors in question
# strike = 40
# sigma = 0.2
# r = 0.06
# delta = 0
# h = 50
# n = 100000

monomials_results <- array(NA, dim=c(3,3,3),dimnames = c("S0","T","k"))

# generate data table inputs in the name format: "P_S0_T_k_first letter of method"
for (i in 1:3){
  for (j in 1:3){
    for (c in 1:3){
      S0 = ifelse(i==1,36,NA)
      S0 = ifelse(i==2,40,S0)
      S0 = ifelse(i==3,44,S0)
      T = ifelse(j==1, 0.5,NA)
      T = ifelse(j==2, 1,T)
      T = ifelse(j==3, 2,T)
      k = ifelse(c==1, 2,NA)
      k = ifelse(c==2, 3,k)
      k = ifelse(c==3, 4,k)
      monomials_results[i,j,c] = my_Put_LSMC(S0 = S0, sigma = 0.2,strike = 40,r = 0.06,delta = 0,T=T ,n = 10000,h = 50,terms = k,method = "monomials")$price
    }
  }
}

Hermite_results <- array(NA, dim=c(3,3,3),dimnames = c("S0","T","k"))

# generate data table inputs in the name format: "P_S0_T_k_first letter of method"
for (i in 1:3){
  for (j in 1:3){
    for (c in 1:3){
      S0 = ifelse(i==1,36,NA)
      S0 = ifelse(i==2,40,S0)
      S0 = ifelse(i==3,44,S0)
      T = ifelse(j==1, 0.5,NA)
      T = ifelse(j==2, 1,T)
      T = ifelse(j==3, 2,T)
      k = ifelse(c==1, 2,NA)
      k = ifelse(c==2, 3,k)
      k = ifelse(c==3, 4,k)
      Hermite_results[i,j,c] = my_Put_LSMC(S0 = S0, sigma = 0.2,strike = 40,r = 0.06,delta = 0,T=T ,n = 10000,h = 50,terms = k,method = "Hermite")$price
    }
  }
}

laguerre_results <- array(NA, dim=c(3,3,3),dimnames = c("S0","T","k"))

# generate data table inputs in the name format: "P_S0_T_k_first letter of method"
for (i in 1:3){
  for (j in 1:3){
    for (c in 1:3){
      S0 = ifelse(i==1,36,NA)
      S0 = ifelse(i==2,40,S0)
      S0 = ifelse(i==3,44,S0)
      T = ifelse(j==1, 0.5,NA)
      T = ifelse(j==2, 1,T)
      T = ifelse(j==3, 2,T)
      k = ifelse(c==1, 2,NA)
      k = ifelse(c==2, 3,k)
      k = ifelse(c==3, 4,k)
      laguerre_results[i,j,c] = my_Put_LSMC(S0 = S0, sigma = 0.2,strike = 40,r = 0.06,delta = 0,T=T ,n = 1000,h = 50,terms = k,method = "laguerre")$price
    }
  }
}


# the best result come in 2-3 terms, increasing terms or using different polynomials does not improve price recognition much
# increasing h does not improve the price recognition very much 
# increasing n improves the confidence interval
# increasing T increases the price of put option
# increasing K reduces price of put option
