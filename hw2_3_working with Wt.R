#IBM algorithm parameters to generate uniform distribution
m=2^31-1
a=7^5
b=0

n= 100000
#first n uniform numbers
U1= rep(NA,n)
U1[1]=5
for (i in 1:(length(U1)-1)){
  U1[i+1]=(a*U1[i]+b)%%m
}
U1=U1/m
# second n uniform numbers
U2 = rep(NA,n)
U2[1]=6
for (i in 1:(length(U2)-1)){
  U2[i+1]=(a*U2[i]+b)%%m
}
U2=U2/m

Narray1_BM = rep(NA,n)
Narray2_BM = rep(NA,n)

Narray1_BM = sqrt(-2*log(U1))*cos(2*pi*U2)
Narray2_BM = sqrt(-2*log(U1))*sin(2*pi*U2)

W5=sqrt(5)*Narray1_BM
param_pos=W5^2+sin(W5)
param_neg=(-W5)^2+sin(-W5)
param_tot=c(param_pos,param_neg)

mean(param_pos)
var(param_pos)
avg+1.96*sqrt(s/length(param_pos))
avg-1.96*sqrt(s/length(param_pos))

mean_neg = mean(param_neg)
var_neg = var(param_neg)
mean_tot = mean(param_tot)
var_tot = var(param_tot)
mean_tot + 1.96*sqrt(var_tot/length(param_tot))
mean_tot - 1.96*sqrt(var_tot/length(param_tot))


#
for (t in c(0.5,3.2,6.5)){
  Wt=sqrt(t)*(Narray1_BM)
  print.default(mean(exp(t/2)*cos(Wt)))
  #print.default(sd(exp(t/2)*cos(Wt)))
}


