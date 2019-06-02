m=2^31-1
a=7^5
b=0

n= 100000
rho = 0.6
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

# check correlation for independence 
cor(U1,U2)

Narray1_BM = rep(NA,n)
Narray2_BM = rep(NA,n)

Narray1_BM = sqrt(-2*log(U1))*cos(2*pi*U2)
Narray2_BM = sqrt(-2*log(U1))*sin(2*pi*U2)

# ------ to test uncomment ------ #
# hist(Narray1_BM)
# mean(Narray1_BM)
# sd(Narray1_BM)
# hist(Narray2_BM)
# mean(Narray2_BM)
# sd(Narray2_BM)
# cor(Narray1_BM,Narray2_BM)

Narray1_BM_rho = rep(NA,n)
Narray2_BM_rho = rep(NA,n)

Narray1_BM_rho = Narray1_BM
Narray2_BM_rho = rho * Narray1_BM + sqrt(1-rho^2) * Narray2_BM
cor(Narray1_BM_rho,Narray2_BM_rho)

Narray1_BM_rho_mean= mean(Narray1_BM_rho)
Narray2_BM_rho_mean= mean(Narray2_BM_rho)

Narray12_BM_rho_cov = (1/(n-1))*sum((Narray1_BM_rho - Narray1_BM_rho_mean)*(Narray2_BM_rho - Narray2_BM_rho_mean))
Narray1_BM_rho_var =  (1/(n-1))*sum((Narray1_BM_rho - Narray1_BM_rho_mean)*(Narray1_BM_rho - Narray1_BM_rho_mean))
Narray2_BM_rho_var =  (1/(n-1))*sum((Narray2_BM_rho - Narray2_BM_rho_mean)*(Narray2_BM_rho - Narray2_BM_rho_mean))
pho_simul=Narray12_BM_rho_cov/(sqrt(Narray1_BM_rho_var*Narray2_BM_rho_var))
pho_simul

X=Narray1_BM_rho
Y=Narray2_BM_rho
output=(X^3+sin(Y)+(X^2*Y))
zerovec = rep(0,n)
output_positive=(zerovec<output)*output
mean(output_positive)