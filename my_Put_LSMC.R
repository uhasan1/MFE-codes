my_Put_LSMC <- function (S0 , sigma, strike, r , delta = 0, T , n = 10000, h = 100, terms = 2, method = "monomial") 
{
  
Dt = T/h; 
method = "monomials"
ind_mat = matrix(data = 0, nrow = n, ncol = h+1)

onescol = matrix(data = 1,nrow = n, ncol = 1)
ranmat = matrix(data = rnorm(n*h/2),nrow = n/2,ncol = h)  
expranmat1= exp((r-0.5*sigma^2)*Dt+sigma*sqrt(Dt)*ranmat)
expranmat2= exp((r-0.5*sigma^2)*Dt-sigma*sqrt(Dt)*ranmat)

St = S0*cbind(onescol,rbind(t(apply(X = expranmat1,MARGIN = 1, FUN = cumprod)),t(apply(X = expranmat2,MARGIN = 1, FUN = cumprod))))

exercise_value <- pmax(strike - St,0) # take care that pmax preserves attributes of FIRST argument
# exercise_value = (strike>St)*(strike-St)

# For stability of faster convergence of the method, use in-the-money paths only in the least squares estimation step 
# as the goal is to estimate the expected continuation time
# So we are putting NA in OTM nodes.
X = matrix(ifelse(St < strike, St, 0),n,h+1) 

# we don't need the last column for regression so omit it
Xm <- X[,-ncol(X)]
k = terms
# define orthogonal functions
if (method == "monomials"){
  if (k == 2) {
    X <- Xm
  }
  else if (k == 3) {
    X <- Xm
    X2 <- Xm^2
  }
  else if (k == 4) {
    X <- Xm
    X2 <- Xm^2
    X3 <- Xm^3
  }
} else if (method == "Hermite"){
  if (k == 2) {
    X = 2*Xm
  }
  else if (k == 3){
    X = 2* Xm
    X2 = 4* (Xm^2) - 2
  }
  else if (k == 4){
    X = 2* Xm
    X2 = 4* (Xm^2) - 2
    X3 = 8* (Xm^3) - 12*Xm
  }
} else if (method == "laguerre"){
  if (k == 2){
    X = exp(-0.5* Xm)
  }
  else if (k ==3) {
    X = exp(-0.5 * Xm)
    X2 = exp(-0.5 * Xm) * (1 - Xm)
  }
  else if (k == 4){
    X = exp(-0.5 * Xm)
    X2 = exp(-0.5 * Xm) * (1 - Xm)
    X3 = exp(-0.5 * Xm) * (1 - (2 * Xm) + 0.5 * (Xm)^2)
  }
}else {
  cat("error! check your method again monomials/laguerre/Hermite")
}


Y <- matrix(NA, nrow = n, ncol = h+1)
cont_value <- matrix(NA, nrow = n, ncol = h) 

Y[,h+1] <-  exercise_value[,h+1]*exp(-r*Dt)

if (k==2){
  for (i in h:1) { # col1 contains the spot so cannot be used in regression: yields NA in the second coefficient
    regYX <- lm(Y[, i + 1] ~ X[, i])
    cont_value[, i] <- regYX$coefficients[1] + regYX$coefficients[2] * X[, i] 
    cont_value[, i] <- ifelse(is.na(cont_value[, i]), 0, cont_value[, i]) # 
    Y[, i] <- ifelse(exercise_value[, i] > cont_value[, i], exercise_value[,i] * exp (-r*Dt), Y[, i + 1] * exp(-r*Dt))  
  }
}

if (k==3){  
  for (i in h:1) { # col1 contains the spot so cannot be used in regression: yields NA in the second coefficient
    regYX <- lm(Y[, i + 1] ~ X[, i] + X2[, i])
    cont_value[, i] <- regYX$coefficients[1] + regYX$coefficients[2] * X[, i] + regYX$coefficients[3] * X2[, i]
    cont_value[, i] <- ifelse(is.na(cont_value[, i]), 0, cont_value[, i]) # 
    Y[, i] <- ifelse(exercise_value[, i] > cont_value[, i], exercise_value[,i] * exp (-r*Dt), Y[, i + 1] * exp(-r*Dt))  
  }
}

if (k==4){  
  for (i in h:1) { # col1 contains the spot so cannot be used in regression: yields NA in the second coefficient
    regYX <- lm(Y[, i + 1] ~ X[, i] + X2[, i] +X3[,i])
    cont_value[, i] <- regYX$coefficients[1] + regYX$coefficients[2] * X[, i] + regYX$coefficients[3] * X2[, i] + regYX$coefficients[3] * X3[, i]
    cont_value[, i] <- ifelse(is.na(cont_value[, i]), 0, cont_value[, i]) # 
    Y[, i] <- ifelse(exercise_value[, i] > cont_value[, i], exercise_value[,i] * exp (-r*Dt), Y[, i + 1] * exp(-r*Dt))  
  }
}

col0 = matrix(0,n,1)
cont_value = cbind(cont_value,col0) # adding zero col to the end of continuation value. there is no point in continuing

ind_mat = ifelse(exercise_value>cont_value, 1,0) # where you excercise

AM_PRICE =ind_mat*exercise_value
AM_PRICE = AM_PRICE[,-1]
ind_mat = ind_mat[,-1] 

final_price = matrix(data = NA, nrow = n, ncol = 1)

for (i in 1:n){
  j=1
  repeat{
    final_price[i] = AM_PRICE[i,j]*exp(-r*Dt*j)
    if ((AM_PRICE[i,j]!=0) | (j == h) ) break
    j = j+1
  }
}

price = mean(final_price)

s= sd(final_price)
Upper = price - 1.96*s/sqrt(n)
Lower = price + 1.96*s/sqrt(n)

AmPUT <- list("price" = price, "Upper.ConfidInt" = Upper, "Lower.confidInt" = Lower)

return(AmPUT)
}

