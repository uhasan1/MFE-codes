# values are stored in 3 matrices, the graphs are provided separately
# in order to reduce the time, reduce the number of simulations or steps in the main function


# constant factors in question
# V0 = 20000
# L0 = 22000
# mu = -0.1
# sigma = 0.2
# gamma = -0.4
 
# r0 = 0.02
# delta = 0.25
# alpha = 0.7
# epsilon = 0.95
# 
# h =100;# steps
# nsim = 10000; # number of simulations
library(plot3D)

default_option <- array(NA, dim=c(8,9,6),dimnames = c("landa1","landa2","T"))
default_probability <- array(NA, dim=c(8,9,6),dimnames = c("landa1","landa2","T"))
default_Exp_time <- array(NA, dim=c(8,9,6),dimnames = c("landa1","landa2","T"))


# generate data table inputs in the name format: "P_S0_T_k_first letter of method"
for (i in 1:8){
  for (j in 1:9){
    for (k in 1:6){
      landa1 = 0.05 + 0.05*(i-1)
      landa2 = 0 + 0.1*(j-1)
      T = 3 + 1*(k-1)
      my_default = mydefault(landa1,landa2,T)
      default_option[i,j,k] = my_default$D
      default_probability[i,j,k] = my_default$Prob
      default_Exp_time[i,j,k] = my_default$Et
    }
  }
}

par(mfrow=c(1,2))
landa2 <- (seq(0,0.8,length.out = 9))
T = (seq(3,8,length.out = 6))
persp3D(x = landa2, y = T, z = default_option[4,,],xlab = "landa2", ylab = "T", zlab = "option value")
persp3D(x = landa2, y = T, z = default_probability[4,,],xlab = "landa2", ylab = "T", zlab = "probability")

landa1 <- (seq(0.05,0.4,length.out = 8))
T = (seq(3,8,length.out = 6))
persp3D(x = landa1, y = T, z = default_option[,5,],xlab = "landa1", ylab = "T", zlab = "option value")
persp3D(x = landa1, y = T, z = default_probability[,5,],xlab = "landa1", ylab = "T", zlab = "probability")
