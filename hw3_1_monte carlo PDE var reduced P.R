
#X0 = 1;  
Y0 = 0.75; 

T = 2;
h =1000; # number of steps
Dt = T/h; 
n = 1000; # number of simulations
t=0;

W = matrix(rnorm(n*h),nrow = n,ncol = h)
#Z = matrix(rnorm(n*h),nrow = n,ncol = h)

Y_1 = matrix(nrow=n, ncol = h+1)
Y_2 = matrix(nrow=n, ncol = h+1)

Y_1[,1]=Y0
Y_2[,1]=Y0


  for (j in 1:h){
    Y_1[,j+1]=Y_1[,j]+((2/(1+t))*Y_1[,j]+(1+t^3)/3)*Dt+((1+t^3)/3)*sqrt(Dt)*W[,j]
    Y_2[,j+1]=Y_2[,j]+((2/(1+t))*Y_2[,j]+(1+t^3)/3)*Dt-((1+t^3)/3)*sqrt(Dt)*W[,j]
    t=t+Dt
  }

Y=rbind(Y_1,Y_2)
# expectedY = mean(Y[,h+1])
# sdY= sd(Y[,h+1])
# YLconfint = expectedY-1.96*sdY/sqrt(nrow(Y))
# YUconfint = expectedY+1.96*sdY/sqrt(nrow(Y))
# cat("E(Y3) is ", expectedY, " within the confidence interval [", YLconfint," , ",YUconfint, "]",sep = "")
sum(Y[,h+1]>5)/nrow(Y)
