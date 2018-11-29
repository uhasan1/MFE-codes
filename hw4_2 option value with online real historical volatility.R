# this is a great code that I have written! it fetches live stock data (from yahoo, though you can change the source) and option data (from google)
# as the risk free rate is set lower, the spread between the option and stock declines
# the remainder difference is due to the difference in implied volatility and historical volatility

library(quantmod)
library(fBasics)
library(RQuantLib)
r = 0.02
startdate = "2011-4-26"
enddate = "2016-4-26"
optionexpirydate = "2017-1-20" # this is expiration date available on google option

getSymbols(Symbols = "GOOG", from=startdate, to= enddate) #5 years data
GOOGLE=as.data.frame(GOOG)
GOOG.rtn=diff(log(GOOGLE$GOOG.Adjusted)) # Compute log returns
ann.vol= sqrt(252)*sd(GOOG.rtn) # annual volatility from daily volatility
S0 = GOOGLE$GOOG.Adjusted[nrow(GOOGLE)] #getting the last price
K = round(1.1*S0/20)*20  #strike price is the closest integer (divisible by 10) to 110% of the current price
T = as.numeric(as.Date("2017-1-20")-as.Date(enddate))/365 # life of option in years

GOOG170120C00820000_myprice = EuropeanOption(type = "call",underlying = S0,strike = K,dividendYield = 0,riskFreeRate = r,maturity = T, volatility = ann.vol)

# installs RCurl and jsonlite packages. will prompt you to select mirror for download
install.packages("RCurl")
install.packages("jsonlite")

library(RCurl)
library(jsonlite)

getOptionQuote <- function(symbol){
  output = list()
  url = paste('http://www.google.com/finance/option_chain?q=', symbol, '&output=json', sep = "")
  x = getURL(url)
  fix = fixJSON(x)
  json = fromJSON(fix)
  numExp = dim(json$expirations)[1]
  for(i in 1:numExp){
    # download each expirations data
    y = json$expirations[i,]$y
    m = json$expirations[i,]$m
    d = json$expirations[i,]$d
    expName = paste(y, m, d, sep = "_")
    if (i > 1){
      url = paste('http://www.google.com/finance/option_chain?q=', symbol, '&output=json&expy=', y, '&expm=', m, '&expd=', d, sep = "")
      json = fromJSON(fixJSON(getURL(url)))
    }
    output[[paste(expName, "calls", sep = "_")]] = json$calls
    output[[paste(expName, "puts", sep = "_")]] = json$puts
  }
  return(output)
}

fixJSON <- function(json_str){
  stuff = c('cid','cp','s','cs','vol','expiry','underlying_id','underlying_price',
            'p','c','oi','e','b','strike','a','name','puts','calls','expirations',
            'y','m','d')
  
  for(i in 1:length(stuff)){
    replacement1 = paste(',"', stuff[i], '":', sep = "")
    replacement2 = paste('\\{"', stuff[i], '":', sep = "")
    regex1 = paste(',', stuff[i], ':', sep = "")
    regex2 = paste('\\{', stuff[i], ':', sep = "")
    json_str = gsub(regex1, replacement1, json_str)
    json_str = gsub(regex2, replacement2, json_str)
  }
  return(json_str)
}
google_opt = getOptionQuote("GOOG")
jan2017call=google_opt$`2017_1_20_calls`
GOOG170120C00820000_price=jan2017call$p[which(jan2017call$strike=="780.00")]

GOOG170120C00820000_myprice$value
GOOG170120C00820000_price

EuropeanOptionImpliedVolatility(type = "call",value = as.numeric(GOOG170120C00820000_price),underlying = S0,strike = K,dividendYield = 0,riskFreeRate = r,maturity = T,volatility = 0.2)
