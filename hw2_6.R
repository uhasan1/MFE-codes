precision=1000
integral_start=0
integral_end=1
interval_length = (integral_end - integral_start)/precision

x=seq(integral_start,integral_end,by=interval_length)
y=4*sqrt(1-x^2)
estim1 = sum(interval_length*y[1:1000]) 
estim2 = sum(interval_length*y[2:1001])
# one is upper estimation, the other lower estimation; if the function is monotonic.
mean(c(estim1,estim2)) # use c! or else only returns the first number!

# m=2^31-1
# a=7^5
# b=0
# 
 n= 10000
# #first n uniform numbers
# U1= rep(NA,n)
# U1[1]=5
# for (i in 1:(length(U1)-1)){
#   U1[i+1]=(a*U1[i]+b)%%m
# }
# U1=U1/m

U1 =runif(n)
y=4*sqrt(1-U1^2)
mean(y)
