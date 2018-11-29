library(GUIDE)

############  credit var: uncorrelated loan default(var,expected loss,unexpected loss) P 361 - 365 #############

loan_number = 100 # N
loan_value = 15
portfolio = loan_number * loan_value # V
default_prob = 0.025   # p
LGD = 1

exp_loss = portfolio * default_prob * LGD

pbinom(4,loan_number,default_prob) # cumulative probability = probability of 4 defaults or less
dbinom(4,loan_number,default_prob) # probability of 4 defaults

var_level = 0.96

i = 0
while(pbinom(i,loan_number,default_prob)<var_level){
  i = i+1
}
# shows when cumulative prob exceeds specified var level

WCDR = i/loan_number # worst case default rate
unexp_loss = WCDR*portfolio - exp_loss

############  credit var P  364 #############

mu = portfolio * default_prob
sigma = mu * sqrt((1 - default_prob)/(loan_number*default_prob))

loss_var = qnorm(var_level,mu,sigma) # var
qnorm(var_level,mu,sigma) - exp_loss  # unexpected loss

# as N goes to infinity The VaR tends to the mean, and the Unexpected Loss VaR goes to zero

qnorm(0.2) # normal inverse cumulative -0.8416212 N^-1()
pnorm(0.2) # normal cumulative 0.5792597 N()
pnorm(qnorm(0.2)) # 0.2!

###### Simulating VaR for Independent Loans P 369 ######
Num=1000 #Number of loans
Size=15 #Dollar size of each loan
PD=0.025 #PD for one loan
alpha=0.95 #VaR alpha

N=10000 #Number of iterations
Iter_loss = array (0, dim=c(N)) #Distribution of losses per iteration 

for (iter in 1:N) {
  U = matrix(rnorm(Num,mean=0,sd=1), 1, Num) #Generate U_i
  Default = (U<qnorm(PD)) #Every loan, every iteration, did it default
  loan_loss = Default*Size #Total loss on each loan for this iteration (assuming LGD is 100%)
  Iter_loss[iter] = sum(loan_loss) #Total loss for this iteration
}

hist(Iter_loss,breaks = 1000) #Histogram of losses over iterations
EL = mean(Iter_loss) #Expected loss
# very important:
VaR = quantile(Iter_loss, alpha) #VaR 

###### Simulating VaR for correlated Loans P 373 ######
Num=1000 #Number of loans
Size=15 #Dollar size of each loan
PD=0.025 #PD for one loan
rho = 0.15 #Correlation between latent variables
alpha=0.95 #VaR alpha
N=10000 #Number of iterations
Iter_loss = array (0, dim=c(N)) #Distribution of losses per iteration 
for (iter in 1:N) {
  F = matrix(rnorm(1),1,Num) #One F Factor value per iteration
  Z = matrix(rnorm(Num,mean=0,sd=1), 1, Num) #idiosyncratic errors for all loans
  U = sqrt(rho)*F + sqrt(1-rho)*Z #Generate the U_i for all loans
  Default = (U<qnorm(PD)) #Every loan, every iteration, did it default
  loan_loss = Default*Size #Total loss on each loan for this iteration (assuming LGD is 100%)
  Iter_loss[iter] = sum(loan_loss) #Total loss for this iteration
}
hist(Iter_loss,breaks = 1000) #Histogram of losses over iterations
EL = mean(Iter_loss) #Expected loss
VaR = quantile(Iter_loss, alpha) #VaR

############  credit var: correlated loan default(var,expected loss,unexpected loss) P 376-379 #############

# P376
default_prob = 0.005   # p
rho = 0.3   # correlation
F = seq(-5,5,by=0.01)    #Factor
plot(F,pnorm((qnorm(default_prob)-sqrt(rho)*F)/(sqrt(1-rho))), ylab = "default conditional on F")    # default conditional on F 

# P378
default_prob = c(0.005,0.025,0.1)
rho = 0.3
X = 0.99 # X% worst case
WCDR_X = pnorm((qnorm(default_prob)-sqrt(rho)*qnorm(1-X))/(sqrt(1-rho)))  # qnorm(1-X) = -qnorm(X)

# P379
unexp_loss =  portfolio * LGD *(WCDR_X - default_prob) # unexpected loss credit var(X)

# Example P380
portfolio = 100
default_prob = 0.02
LGD = 0.4
rho = 0.1
X = 0.999
WCDR_X = pnorm((qnorm(default_prob)-sqrt(rho)*qnorm(1-X))/(sqrt(1-rho)))  # qnorm(1-X) = -qnorm(X)
unexp_loss =  portfolio * LGD *(WCDR_X - default_prob) # credit var(X)

# Example P383
portfolio = 325
default_prob = 0.01
LGD = 1.0
rho = 0.4
X = 0.962
max_loss = 19.5
WCDR_X = pnorm((qnorm(default_prob)-sqrt(rho)*qnorm(1-X))/(sqrt(1-rho)))  # qnorm(1-X) = -qnorm(X)
unexp_loss =  portfolio * LGD *(WCDR_X - default_prob) # credit var(X)


exp_loss_fun <- function (X) {
  WCDR_X = pnorm((qnorm(default_prob)-sqrt(rho)*qnorm(1-X))/(sqrt(1-rho)))
  return(portfolio * LGD *(WCDR_X)-max_loss)
}


# curve(exp_loss_fun, 0, 1)
## to drw a curve between 0,1

# IMPORTANT! TAKE CARE! IT IS EXPECTED LOSS SET TO ZERO! NOR UNEXPECTED LOSS!
uniroot(exp_loss_fun,c(0,1))
1 - uniroot(exp_loss_fun,c(0,1))$root


##### P385 #####
loan_number = 1000 # N
loan_value = 0.325
portfolio = loan_number * loan_value # V
default_prob = 0.01   # p

exp_loss = portfolio * default_prob
19.5/loan_value # haw many at least should default?
pbinom(59,loan_number,default_prob) # probability that less than that many default
1- pbinom(59,loan_number,default_prob) # how probable is for bank to pay?

mu = portfolio * default_prob
sigma = mu * sqrt((1 - default_prob)/(loan_number*default_prob))
pnorm(60,mu,sigma)
1- pnorm(60,mu,sigma)
################

###### MTM: mark to market P396 ######

######### finding yield with semiannual coupons & continuous compounding P137 ############
yield_find = function(y){
  Maturity = 3
  coupon_rate = 0.1
  bond_price = 0.94213
  
  bond = 0
  for (i in 1:(2*Maturity)){
    bond = bond + (coupon_rate/2)*exp(-i*y/2)
  }
  bond = bond + exp(-i*y/2)
  difference = bond - bond_price
  return(difference)
}

yield = uniroot(yield_find,c(-0.2,0.2))$root

################## Package guide for fixed income calculation #########################
library(GUIDE)

mybond_dur = function (faceval,couprate,discrate,maturity,frequency = "semi-annual",ratefreq = "continuous comp",durtype = "Modified") {
    if (frequency == "quarterly") {
      freq <- 4
      times <- seq(from = 0.25, by = 0.25, length.out = maturity * 
                     freq)
    }
    else if (frequency == "semi-annual") {
      freq <- 2
      times <- seq(from = 0.5, by = 0.5, length.out = maturity * 
                     freq)
    }
    else {
      freq <- 1
      times <- seq(from = 1, by = 1, length.out = maturity * 
                     freq)
    }
    if (ratefreq == "continuous comp") {
      pvfactors = exp(-discrate * times)
    }
    else if (ratefreq == "annual comp") {
      pvfactors = 1/(1 + discrate)^times
    }
    else {
      pvfactors = 1/(1 + discrate/freq)^(freq * times)
    }
    coupon <- couprate * faceval/(freq)
    cashflows <- rep(coupon, maturity * freq)
    cashflows[length(cashflows)] = cashflows[length(cashflows)] + 
      faceval
    price <- sum(cashflows * pvfactors)
    dur = sum(cashflows * pvfactors * times)/price
    if (durtype == "Modified") {
      dur <- dur/(1 + discrate/freq)
    }
    return(dur)
    
}

mybond_dur(faceval = 100,couprate = 0.10,discrate = 0.12,maturity = 3,frequency = "semi-annual",ratefreq = "continuous comp",durtype = "Modified")

#####################

mybond_price_panel = function () 
{
  my.draw <- function(panel) {
    faceval <- as.numeric(panel$facevalue)
    discrate = as.numeric(panel$discrate)/100
    maturity <- panel$maturity
    if (panel$frequency == "quarterly") {
      freq <- 4
      times <- seq(from = 0.25, by = 0.25, length.out = maturity * 
                     freq)
    }
    else if (panel$frequency == "semi-annual") {
      freq <- 2
      times <- seq(from = 0.5, by = 0.5, length.out = maturity * 
                     freq)
    }
    else {
      freq <- 1
      times <- seq(from = 1, by = 1, length.out = maturity * 
                     freq)
    }
    if (panel$ratefreq == "continuous comp") {
      pvfactors = exp(-discrate * times)
    }
    else if (panel$ratefreq == "annual comp") {
      pvfactors = 1/(1 + discrate)^times
    }
    else {
      pvfactors = 1/(1 + discrate/freq)^(freq * times)
    }
    coupon <- panel$couprate * faceval/(100 * freq)
    cashflows <- rep(coupon, maturity * freq)
    cashflows[length(cashflows)] = cashflows[length(cashflows)] + 
      faceval
    price <- sum(cashflows * pvfactors)
    price <- round(price, 2)
    plot(1:10, 1:10, type = "n", xlab = "", ylab = "", axes = FALSE, 
         frame = TRUE)
    text(5, 5, paste("Price: ", price), cex = 1.4)
    panel
  }
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, tkrp)
    panel
  }
  my.panel <- rp.control("Bond Price", frequency = "quarterly", 
                         couprate = 8, discrate = 10, maturity = 10)
  rp.textentry(panel = my.panel, variable = facevalue, labels = "Face Value:   ", 
               action = my.redraw, initval = "1000")
  rp.doublebutton(my.panel, variable = couprate, step = 0.25, 
                  title = "Coupon (%  p.a.)", initval = 10, range = c(0, 
                                                                      15), showvalue = TRUE, action = my.redraw)
  rp.doublebutton(my.panel, variable = discrate, step = 0.25, 
                  title = "Discount Rate(%  p.a.)", initval = 10, range = c(0, 
                                                                            15), showvalue = TRUE, action = my.redraw)
  rp.doublebutton(my.panel, variable = maturity, step = 0.5, 
                  title = "Maturity(Yrs)", initval = 10, range = c(0, 
                                                                   25), showvalue = TRUE, action = my.redraw)
  rp.radiogroup(panel = my.panel, variable = frequency, vals = c("quarterly", 
                                                                 "semi-annual", "annual"), action = my.redraw, title = "Coupon payments")
  rp.radiogroup(panel = my.panel, variable = ratefreq, vals = c("continuous comp", 
                                                                "same as coupon freq", "annual comp"), action = my.redraw, 
                title = "Frequency of discount rate")
  rp.tkrplot(panel = my.panel, name = tkrp, plotfun = my.draw)
}


mybond_dur_panel = function () 
{
  my.draw <- function(panel) {
    faceval <- as.numeric(panel$facevalue)
    discrate = as.numeric(panel$discrate)/100
    maturity <- panel$maturity
    if (panel$frequency == "quarterly") {
      freq <- 4
      times <- seq(from = 0.25, by = 0.25, length.out = maturity * 
                     freq)
    }
    else if (panel$frequency == "semi-annual") {
      freq <- 2
      times <- seq(from = 0.5, by = 0.5, length.out = maturity * 
                     freq)
    }
    else {
      freq <- 1
      times <- seq(from = 1, by = 1, length.out = maturity * 
                     freq)
    }
    if (panel$ratefreq == "continuous comp") {
      pvfactors = exp(-discrate * times)
    }
    else if (panel$ratefreq == "annual comp") {
      pvfactors = 1/(1 + discrate)^times
    }
    else {
      pvfactors = 1/(1 + discrate/freq)^(freq * times)
    }
    coupon <- panel$couprate * faceval/(100 * freq)
    cashflows <- rep(coupon, maturity * freq)
    cashflows[length(cashflows)] = cashflows[length(cashflows)] + 
      faceval
    price <- sum(cashflows * pvfactors)
    dur = sum(cashflows * pvfactors * times)/price
    if (panel$durtype == "Modified") {
      dur <- dur/(1 + discrate/freq)
    }
    dur <- round(dur, 2)
    plot(1:10, 1:10, type = "n", xlab = "", ylab = "", axes = FALSE, 
         frame = TRUE)
    text(5, 5, paste("Duration: ", dur), cex = 1.4)
    panel
  }
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, tkrp)
    panel
  }
  my.panel <- rp.control("Bond Duration", frequency = "quarterly", 
                         couprate = 8, discrate = 10, maturity = 10)
  rp.textentry(panel = my.panel, variable = facevalue, labels = "Face Value:   ", 
               action = my.redraw, initval = "1000")
  rp.doublebutton(my.panel, variable = couprate, step = 0.25, 
                  title = "Coupon (%  p.a.)", initval = 10, range = c(0, 
                                                                      15), showvalue = TRUE, action = my.redraw)
  rp.doublebutton(my.panel, variable = discrate, step = 0.25, 
                  title = "Discount Rate (%  p.a.)", initval = 10, range = c(1, 
                                                                             15), showvalue = TRUE, action = my.redraw)
  rp.doublebutton(my.panel, variable = maturity, step = 0.25, 
                  title = "Maturity (Yrs)", initval = 10, range = c(1, 
                                                                    25), showvalue = TRUE, action = my.redraw)
  rp.radiogroup(panel = my.panel, variable = frequency, vals = c("quarterly", 
                                                                 "semi-annual", "annual"), action = my.redraw, title = "Coupon payments")
  rp.radiogroup(panel = my.panel, variable = ratefreq, vals = c("continuous comp", 
                                                                "same as coupon freq", "annual comp"), action = my.redraw, 
                title = "Frequency of discount rate")
  rp.radiogroup(panel = my.panel, variable = durtype, vals = c("Macaulay", 
                                                               "Modified"), action = my.redraw, title = "Duration formula")
  rp.tkrplot(panel = my.panel, name = tkrp, plotfun = my.draw)
}


################ Quadratic/Delta-Gamma Model P174-179 ######################
S = 1500
sigma = 0.02
delta = 0.5
gamma = -0.07

mu_p = 0.5* S^2 * gamma * sigma^2
exp_p2 = (S*sigma*delta)^2 + 0.75*S^4*gamma^2*sigma^4
var_p = exp_p2 - mu_p^2
sigma_p = sqrt(var_p)
exp_p3 = 4.5* S^4 * delta^2*gamma*sigma^4 + 1.87*S^6*gamma^3*sigma^6
 # use moments to find skewness
skew_p = (exp_p3 -3*exp_p2 *mu_p +2*mu_p^3)/(sigma_p^3)
z_q = -qnorm(0.95) # confidence interval
w_q = z_q +(1/6)*(z_q^2 -1) *skew_p

abs(mu_p + w_q*sigma_p) ### VaR 95%

############### normal-linear approach: Duration Approach P 182  ###############3
portfolio = 820
duration = 5
dollar_dur = portfolio*duration

mu_y = 0
sigma_y = 0.001

var = qnorm(0.95,mean = mu_y,sd = sigma_y) * dollar_dur

######## Bond Portfolio Example P 185##############
corr_mat = matrix(data = c(1,0.9,0.6,0.9,1,0.7,0.6,0.7,1),nrow = 3,ncol = 3,byrow = T)
var_mat = diag(x = c(0.0006,0.001,0.002),nrow = 3)
cov_mat = var_mat %*% corr_mat %*% var_mat
size_mat = matrix(data = c(37397,331382,678074),nrow = 1)
weight_mat = size_mat/sum(size_mat)

var_p = size_mat %*% cov_mat %*% t(size_mat)
var_p = as.numeric(var_p)
sigma_p = sqrt(var_p)

T = 10 # 10 day
var_T = sqrt(T)*qnorm(0.99)*sigma_p


########### Linear Discriminant Analysis P 283 #############

p1 = 0.25
p0 = 0.75 # 1-p1
mu1 = 9
mu0 = 2
variance = 4

cutoff1 = 0.001
cutoff2 = 0.005

f = function(x){
  a = (mu1-mu0)*x/variance - 0.5* (mu1 + mu0)*(mu1-mu0)/variance + log(p1/p0) - log(cutoff1/(1-cutoff1))
  return(a)
}


g = function(x){
  a = (mu1-mu0)*x/variance - 0.5* (mu1 + mu0)*(mu1-mu0)/variance + log(p1/p0) - log(cutoff2/(1-cutoff2))
  return(a)
}

uniroot(f = f ,interval = c(2,9))
uniroot(f = g ,interval = c(2,9))

###### bond price explicit #######
mybond_price = function (faceval =100, couprate, discrate, maturity,freq = 2, compounding = c("same","continuous","annual")){
  
    times <- seq(from = (1/freq), by = (1/freq), length.out = maturity * freq)
    if (compounding == "continuous") {
      pvfactors = exp(-discrate * times)
    }
    else if (compounding == "annual") {
      pvfactors = 1/(1 + discrate)^times
    }
    else {
      pvfactors = 1/(1 + discrate/freq)^(freq * times)
    }
    coupon <- couprate * faceval/(freq)
    cashflows <- rep(coupon, maturity * freq)
    cashflows[length(cashflows)] = cashflows[length(cashflows)] + faceval
    price <- sum(cashflows * pvfactors)
    return(price)
}



################### Q1 excercise ##################


faceval = 100
couprate =0.07
freq = 2 # semi-annual coupons =2 ; quarterly =4, annually =1
rfr = 0.04
rfr_compounding = "semi-annual"
ytm = 0.05 
ytm_compounding = "semi-annual"
maturity = 3
recovery = 0.45 # considered as face value


mybond_price(faceval = faceval,couprate = couprate,discrate = rfr,compounding = rfr_compounding ,maturity = maturity,freq = freq)

T =seq(from = 1/freq, to = maturity, by=1/freq)
dat <- as.data.frame(matrix(ncol=10, nrow=length(T)))

colnames(dat)[1] = "T"
dat[,1] = T

for (i in 1:length(T)){
  dat$V2[i] = mybond_price(faceval = faceval, couprate = couprate, freq = freq , discrate = rfr, compounding = rfr_compounding, maturity = maturity-T[i]) + faceval*couprate/2
}
dat$V2[length(T)] = faceval + couprate* faceval/2 # last point
colnames(dat)[2] = "NO Default price"

colnames(dat)[3] = "recovery"
dat$recovery = faceval*recovery

colnames(dat)[4] = "loss"
dat$loss = dat$'NO Default price' - dat$recovery

DF = mybond_price(faceval = 1,couprate = 0, discrate = rfr, maturity = 1/freq,freq = freq, compounding = rfr_compounding)

colnames(dat)[5] = "PV loss"
dat$`PV loss` = dat$loss *(DF^(dat$T/dat$T[1]))

PV_risky = mybond_price(faceval = faceval, couprate = couprate, discrate = ytm, compounding = ytm_compounding, maturity = maturity, freq = freq)
PV_riskfree = mybond_price(faceval = faceval, couprate = couprate, discrate = rfr, compounding = rfr_compounding, maturity = maturity, freq = freq)
PV_exp_loss = PV_riskfree - PV_risky

Q = PV_exp_loss/sum(dat$`PV loss`)
cat("unconditionl probability of default in each frequency is", 100*Q, "percent") 

q = 0.002 # no matter what value you start with
colnames(dat)[6] = "PV loss conditional"
f = function(q){
  dat$`PV loss conditional` = dat$`PV loss`* q*((1-q)^((dat$T/dat$T[1])-1)) 
  diference = sum(dat$`PV loss conditional`) - PV_exp_loss
  return(diference)
}

uniroot(f, c(0,0.1))$root

################ Q2 excercise not finished ###################

faceval = 100
couprate = 0.08
freq = 1 # semi-annual coupons =2 ; quarterly =4, annually =1
rfr = 0.045
rfr_compounding = "continuous"
ytm = 0.06 
ytm_compounding = "continuous"
maturity = 1
recovery = 0.35 # considered as face value


mybond_price(faceval = faceval,couprate = couprate,discrate = rfr,compounding = rfr_compounding ,maturity = maturity,freq = freq)

T =seq(from = 1/2, to = maturity, by=1)
dat1 <- as.data.frame(matrix(ncol=10, nrow=length(T)))

colnames(dat1)[1] = "T"
dat1[,1] = T

for (i in 1:length(T)){
  dat1$V2[i] = mybond_price(faceval = faceval, couprate = couprate, freq = freq , discrate = rfr, compounding = rfr_compounding, maturity = maturity-T[i]) 
}
dat1$V2[length(T)] = faceval + couprate* faceval/2 # last point
colnames(dat1)[2] = "NO Default price"

colnames(dat1)[3] = "recovery"
dat1$recovery = faceval*recovery

colnames(dat1)[4] = "loss"
dat1$loss = dat1$`NO Default price` - dat1$recovery

DF = mybond_price(faceval = 1,couprate = 0, discrate = rfr, maturity = 1/freq,freq = freq, compounding = rfr_compounding)

colnames(dat1)[5] = "PV loss"
dat1$`PV loss` = dat1$loss *(DF^(dat1$T/dat1$T[1]))

PV_risky = mybond_price(faceval = faceval, couprate = couprate, discrate = ytm, compounding = ytm_compounding, maturity = maturity, freq = freq)
PV_riskfree = mybond_price(faceval = faceval, couprate = couprate, discrate = rfr, compounding = rfr_compounding, maturity = maturity, freq = freq)
PV_exp_loss = PV_riskfree - PV_risky

Q = PV_exp_loss/sum(dat1$`PV loss`)
cat("unconditionl probability of default in each frequency is", 100*Q, "percent") 


