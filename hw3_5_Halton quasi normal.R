m=2^31-1
a=7^5
b=0

n= 100
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

number2digits=function(n,base){ 
  #first digit in output is the least significant 
  digit=n%%base 
  if (n<base) 
    digit 
  else 
    c(digit,number2digits((n-digit)/base,base)) 
} 

digits2number=function(digits,base){ 
  #first digit in input should be the most significant 
  result=0 
  for (digit in digits) 
    result=(base*result)+digit 
  result 
} 

halton=function(n,base=2){ 
  #generate n first digits of the van der Corput sequence 
  output=NA*1:n 
  for(i in 1:n){ 
    digits=number2digits(i,base) 
    output[i]=digits2number(digits,base)/base^length(digits) 
  } 
  output 
} 

halton2 = halton(n = 100,base = 2)
halton7 = halton(n = 100,base = 7)
halton4 = halton(n = 100,base = 4)

plot(halton7)
plot(halton2,col="red")
points(halton4, col="blue")

halton2 = halton(n = 10000,base = 2)
halton7 = halton(n = 10000,base = 7)
halton4 = halton(n = 10000,base = 4)
halton5 = halton(n = 10000,base = 5)


result2_7= exp(-halton2*halton7)*(sin(6*pi*halton2)+(sign(cos(2*pi*halton7)) * abs(cos(2*pi*halton7))^(1/3)))
mean(result2_7)

result2_4= exp(-halton2*halton4)*(sin(6*pi*halton2)+(sign(cos(2*pi*halton4)) * abs(cos(2*pi*halton4))^(1/3)))
mean(result2_4)

result5_7= exp(-halton5*halton7)*(sin(6*pi*halton5)+(sign(cos(2*pi*halton7)) * abs(cos(2*pi*halton7))^(1/3)))
mean(result5_7)

cat("We only get correct answers when bases are coprime.")


