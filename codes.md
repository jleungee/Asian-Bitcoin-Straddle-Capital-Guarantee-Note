FINA-4354-final-project—source-codes
================
Group 4
12/7/2021

Bitcoin is one of the best incarnations of a “High Risk High Return”
asset. We understand that many retail investors would like to profit
from its high return yet are deterred by its tedious volatility. Hence,
we have designed a product, the Asian Bitcoin Straddle Capital
Guaranteed Note (ABC Note), which

1.  guarantees your initial capital;

2.  mitigates the risk greatly;

3.  has low entry barrier; and

4.  one can profit no matter price surges or plummets.

…………………………………………………………………….. ……………………………………………………………………..

**Annualized return of 15% in the investment period between 2019-2021**

…………………………………………………………………….. ……………………………………………………………………..

To briefly describe the ABC notes, it is compromised with **24 long
straddles** which will be executed **monthly**. For each straddle, its
**strike price will be set as the mean Bitcoin price in the previous
month** and all the monthly payoffs will be **accumulated** and given to
the investors after 24 months.

…………………………………………………………………….. ……………………………………………………………………..

This codebook have five sections:

1.  Set-up - Import necessary libraries and data

2.  Background research of Bitcoin - demonstrate Bitcoin’s
    characteristics quantitatively and graphically to support the
    stochastic model used in later sections.

3.  Monte Carlo Simulation pricing - optimize parameters in the
    simulation and use simulation to price ABC Note

4.  Backtesting of the ABC note - Backtest the return of the ABC note
    from 2017-2021

5.  Product variant - Explore the alternative structures of ABC Note and
    show their backtesting results

\*\*Please do not alter the file locations/working directories of the
data and the codebook

\*\*Please note that the five sections above are inter-related and
should be run one by one

\*\*The “current” date in this file is “2021-11-17” as it is when
Bitcoin option prices are retrieved from outer sources

\*\* It is not a must to run “Data from outer sources (csv and
quantmod)” as it is used to display how we get our raw data from
csv/quantmod. You can directly run “Data from outer sources (rds)”

# 1. Set-up

## Packages

``` r
# rm(list=ls()) # for debugging

# Installing packages (only run if needed)
# install.packages('quantmod')
# install.packages("tidyverse")
# So apparently reshape2 can only be run with an updated Rcpp
# install.packages("reshape2") 
# install.packages("Rcpp")
# update.packages("Rcpp")
# install.packages("doRNG")
# install.packages("doParallel")


# Libraries and libraries
library(magrittr) # the %>% library
library(tidyverse) # R packages designed for data science
library(quantmod) # the Quantitative Financial Modelling Framework
library(reshape2) # For melting dataframes
library(Rcpp) # For melting dataframes
library(doParallel) # For running parallel for loop 
library(doRNG) # For running parallel for loop with random numbers
registerDoParallel(cores=2) # standard

# For calculating how long this program needs
start.time <- Sys.time()

# Set the working directory for retrieving data
wd <- getwd()
wd <- paste(substr(wd,1,(nchar(wd)-4)),"data",sep="")
```

## Data from outer sources (csv and quantmod)

``` r
## The codes here are just a display of how the rds files are created
## Not a must to run it

# setwd(wd) # set working directory
## Getting 3-month Treasury-bill data from FRED
# dtb3 <- read.csv("DTB3.csv") %>% data.frame %>% na.omit
# dtb3$DTB3 <- dtb3$DTB3 %>% as.numeric # reformatting the price data as numbers
# dtb3$Date <- dtb3$Date %>% as.Date # reformatting dates to date data type 
## changing annualized price data to daily data
# dtb3 <- dtb3 %>% mutate(DTB3 = DTB3 / 36000) 

## Getting Bitcoin option price from Deribit
# Db.option <- read.csv("Bitcoin Option.csv") %>% data.frame %>% na.omit
## reformatting the date as date type format
# Db.option$Expiry.Date <- Db.option$Expiry.Date %>% as.Date 

## Getting bitcoin prices from quantmod
## you can add more crypto to the list
# tick <- tick.c <- c("BTC") # tickers that we want to get their price 
# year <- 10 # number of years of Bitcoin price data to be retrieved
# date <- as.Date("2021-11-17") # deemed as the "current date" 
# start.date <- date - 365*year # Bitcoin prices are retrieved from this date

# for (i in tick){
#     tryCatch({
#         # Getting price data
#         tick.usd <- paste(i,"-USD",sep="") # add "-USD" at the back of ticker
#         getSymbols(tick.usd, from = start.date) # get price data
#         # only get the adjusted price
#         ad.price <- tick.usd %>% get %>% Ad %>% data.frame %>% na.omit 
#         # Add a column of date which was originally the rownames
#         Date <- rownames(ad.price) %>% as.Date 
#         # assign the price variable name as the ticker without the "-USD" 
#         assign(i,cbind(Date,ad.price)) 
# 
#         # Calculating Return
#         tick.r <- paste(i,".r",sep="") # creating the variable name
#         r <- ad.price[,1] %>% log %>% diff %>% data.frame %>% na.omit # return 
#         # Add a column of date
#         r_df <- cbind(Date[2:length(Date)],r) %>% na.omit %>% data.frame 
#         assign(tick.r,r_df) # Assign the return dataframe with a name
# 
#     }, error=function(e){
#         # if there are problems accessing the data
#         tick.c <- tick.c[tick.c != i]
#         cat("Not Exist:",conditionMessage(e), "\n")}) # print error message
# }

## Saving data retrieved from csv files/ quantmod as rds
# saveRDS(dtb3,fil="dtb3.rds")
# saveRDS(Db.option,fil="Db.option.rds")
# saveRDS(BTC,fil="BTC.rds")
# saveRDS(BTC.r,fil="BTC.r.rds")
```

## Data from outer sources (rds)

``` r
setwd(wd)
# Getting data from rds files
dtb3 <- readRDS("dtb3.rds") # Risk Free Rate
Db.option <- readRDS("Db.option.rds") # Bitcoin Option
BTC <- readRDS("BTC.rds") # Bitcoin price
BTC.r <- readRDS("BTC.r.rds") # Bitcoin Return
```

# 2. Bitcoin background research

Bitcoin is a high risk high return cryptocurrency. As you can see in
this section, it can yield **&gt;1000% annual simple holding return**
yet it has an astonishing **standard error of 197%**. This marks part of
our economic motivations to create the ABC Notes - **capturing the high
yield of the asset while mitigating its risk**.

With this aim in mind, we tried to investigate its volatility and
distribution. We found that its volatility mean-reverting property, as
well as heavy tail property. These observations provide support for us
to use the Heston Model for simulating Bitcoin price in subsequent
sections.

## Get Bitcoin Prices

``` r
# This section of codes help retrieve the price data (USD) of a list of cryptos 
# You can add cryptos other than BTC in the tick vector to retrieve more cryptos

tick <- tick.c <- c("BTC") # A vector of tickers that we want to get their price 
year <- 10 # number of years of Bitcoin price data to be retrieved
date <- as.Date("2021-11-17") # deemed as the "current date" 
start.date <- date - 365*year # Bitcoin prices are retrieved from this date

for (i in tick){
    tryCatch({
        # Getting price data
        tick.usd <- paste(i,"-USD",sep="") # add "-USD" at the back of ticker
        getSymbols(tick.usd, from = start.date) # get price data
        # only get the adjusted price
        ad.price <- tick.usd %>% get %>% Ad %>% data.frame %>% na.omit 
        # Add a column of date which was originally the row names
        Date <- rownames(ad.price) %>% as.Date 
        # assign the price variable name as the ticker without the "-USD" 
        assign(i,cbind(Date,ad.price)) 

        # Calculating Return
        tick.r <- paste(i,".r",sep="") # creating the variable name
        r <- ad.price[,1] %>% log %>% diff %>% data.frame %>% na.omit # return 
        # Add a column of date
        r_df <- cbind(Date[2:length(Date)],r) %>% na.omit %>% data.frame 
        assign(tick.r,r_df) # Assign the return data frame with a variable name

    }, error=function(e){
        # if there are problems accessing the data
        tick.c <- tick.c[tick.c != i]
        cat("Not Exist:",conditionMessage(e), "\n")}) # print out error message
}
```

## Brief Review of Bitcoin

``` r
# Here provides a brief view of the price, return, and volatility of bitcoin
BTC.ad <- BTC$BTC.USD.Adjusted
BTC.sd <- sd(BTC.r$.)/sqrt(1/(252*year)) # annualized sd of bitcoin

cat("Max Price:",max(BTC.ad),"\t","Max Daily Return:",max(BTC.r$.),"\n")
```

    ## Max Price: 67566.83   Max Daily Return: 0.225119

``` r
cat("Min Price:",min(BTC.ad),"\t","Min Daily Return:",min(BTC.r$.),"\n")
```

    ## Min Price: 178.103    Min Daily Return: -0.4647302

``` r
cat("Mean Price:",mean(BTC.ad),"\t","Mean Daily Return:",mean(BTC.r$.),"\n\n")
```

    ## Mean Price: 10980.3   Mean Daily Return: 0.00168329

``` r
cat("Simple Holding Payoff:",
    as.numeric(BTC.ad[length(BTC.ad)])-as.numeric(BTC.ad[1]),"\n")
```

    ## Simple Holding Payoff: 41174.94

``` r
cat("Simple Holding Return: ",
    (as.numeric(BTC.ad[length(BTC.ad)])/as.numeric(BTC.ad[1])-1)*100,"%","\n\n",
    sep="")
```

    ## Simple Holding Return: 9003.253%

``` r
cat("Standard Error of Price:",sd(BTC.ad)/sqrt(1/(252*year)),"\n")
```

    ## Standard Error of Price: 797733.3

``` r
cat("Standard Error of Return:",sd(BTC.r$.)/sqrt(1/(252*year)),"\n")
```

    ## Standard Error of Return: 1.967353

``` r
# Super high return and fluctuations
```

``` r
# Visualizing the price (adjusted), return, and 30-day volatility movement

vol_30d <- NULL # initializing the vector for 30d volatility

# Calculating the 30d vol
for (i in (1:(length(BTC.r$.)-30))){
    vol_30d[i] <- (sd(BTC.r$.[i:(i+30)])/sqrt(1/252))^2 
}

# Merging ad.price, daily return, and 30d vol data
BTC.plot <- cbind(BTC[31:length(BTC.r$.),],
                  BTC.r$.[31:length(BTC.r$.)],vol_30d) %>% na.omit 
colnames(BTC.plot) <- c("Date","Ad","r","vol") # changing the column names

# Changing the dataframe formatting such that it is plottable
BTC.plot$Ad <- BTC.plot$Ad/10000 %>% unlist
BTC.plot.s <- subset(BTC.plot, Date > date-365) # only plot 1 year data
k <- melt(BTC.plot.s,id=c("Date"))

# plotting the data
ggplot(data=k,
       aes(x=Date, y=value, group=variable,colour=variable)) +
       geom_line()
```

![](codes_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# We can see that when Bitcoin price increases, volatility decreases, vice versa
```

``` r
# Visualizing only the movement of return and 30d volatility
BTC.plot2 <- subset(BTC.plot.s, select = -c(Ad))
k <- melt(BTC.plot2,id=c("Date"))

# plotting the data
ggplot(data=k,
       aes(x=Date, y=value, group=variable,colour=variable)) +
       geom_line()
```

![](codes_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
# It seems that return and volatility does not have apparent relationship
```

## Mean Reversion in Bitcoin volatility

``` r
# Plotting only volatility
BTC.plot3 <- subset(BTC.plot, select = -c(r,Ad),Date > date-365*3)
k <- melt(BTC.plot3,id=c("Date"))

# plotting the data
ggplot(data=k,
       aes(x=Date, y=value, group=variable,colour=variable)) +
       geom_line() +
       geom_smooth(method=lm,col=4)
```

    ## `geom_smooth()` using formula 'y ~ x'

![](codes_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# Despite great fluctuation in volatility, it reverts around a mean value
# which shows property of mean reversion
```

## Heavy tails in Bitcoin return

``` r
# QQplot of bitcoin market return, which shows significant heavy tails 
qqnorm(BTC.r$.)
qqline(BTC.r$., col = 2, lwd = 2)
```

![](codes_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

# Pricing ABC Notes with Monte Carlo Simulation and Heston Model

As ABC Note has an Asian Option feature and thus it has no closed form
solution. Hence, we need to use Monte Carlo Simulation to price the
option. As discussed in the previous section, Heston Model will be the
stochastic model used. This section will discuss how we should model the
Bitcoin price method.

Before delving into the actual pricing, we need to find the parameter
*κ* and *ξ* as they can’t be observed from historical data like other
parameters in the Heston model. Hence, we present two approaches to
estimate the two parameters:

1.  Minimizing the sum of squared difference between Monte Carlo
    Simulated European option prices and actual options prices
2.  Minimizing the sum of squared difference between simulated
    volatility and actual volatility

For both methods, Q-Q plots of the simulated price and correlation
between simulated and actual price are shown. The parameters calibrated
by method 1 has higher correlation between simulated an actual price and
create more significant heavy tails. Thus, **the *κ* and *ξ* calculated
by method 1 is chosen**.

After determining the parameters, we can start pricing the ABC Note by
Monte Carlo Simulation, which returns a value of **\~USD150,000**.
Undoubtedly, it is too pricy for an investor and thus reallocate the
capital and the final product worth USD1000 which:

1.  1% Management Fee
2.  79.2% Risk-free bond
3.  19.8% ABC Note

## Approximation of parameters: *κ*, *ξ*

Method 1: Minimizing the gap between actual option price and option
price predicted by Monte Carlo of Heston Model

We have taken around 200 Bitcoin options data from the crypto exchange
platform Deribit, including their prices, strikes, time to maturities,
and option types. With these information, as well as the Bitcoin price
data date till the options retrieval date (2021-11-17), we can run a
Monte Carlo simulation using the Heston model to price each option and
try to find the kappa and xi that minimize the sum of squared difference
between all simulated option prices and actual option prices.

``` r
# Cleaning the actual option data from Deribit 
Db <- Db.option[,c("Expiry.Date",
                   "Time.to.Maturity..Days.",
                   "Type",
                   "Strike",
                   "Last.Price..USD.",
                   "Volume")] %>% na.omit # subsetting the columns needed
colnames(Db) <- c("Date","T","Type","K","price","Volume") # renaming columns
Db <- subset(Db,Volume > 5) # remove options that aren't frequently traded
Db$T <- Db$T/365 # annualize the time to maturity
# risk-free return of each option
Db$r <- ((1 + dtb3$DTB3[length(dtb3$DTB3)])^Db$T-1)/T*365 
```

``` r
# The function "heston" can
# 1. return the sum of squared error between MC simulated option price 
#    calculated by the kappa and xi got with each method and actual price
# 2. Q-Q plot of the simulated Bitcoin price will also be returned
# 3. A vector of simulated price will be returned


heston <- function(para, # para[1] - kappa; para[2] - xi
                   data = Db, # data for actual option price
                   qq = FALSE,
                   plot = FALSE,
                   rf=dtb3,
                   voldf = BTC.plot,
                   retdf = BTC.r){
  
  n <- length(data$Date) # number of options
  kappa <- xi <- NULL # initializing the vector for kappa and xi   
  f.pre <- f.total <- NULL # initialize payoff vectors & squared error vector

  f.vec <- foreach (k = 1:n)%dorng%{ # parallelization
      tr <- data[k,] # for each option,
      K <- tr$K # we get its strike,
      T.day <- tr$T*365 # time to maturity (days)
      Type <- tr$Type # option type,
      price <- tr$price # and actual price.
      
      # then we set the Heston parameters
      d <- 50 # There is around 120 options so there are ~6000 trials
      r <- rf$DTB3[length(rf$DTB3)]*365 # current annualized risk free rate
      m <- 365 # number of days in a year
      dt <- 1/m # time step - day
      f <- rep(0,d) # initializing the payoff vector of each simulation
      s <- nu <- rep(0,T.day+1) # initialize price and volatility vector
      s[1] = voldf$Ad[length(voldf$Ad)]*10000 # start from current price
      theta <- (sd(retdf$.)/sqrt(1/(365)))^2 # annualized vol of return 
      nu[1] <- theta 

      # Next, we run the actual simulation
      for (j in 1:d){ # For each trials
        for (i in 2:(T.day+1)){ # For each time step
          # Initializing the stochastic terms           
          dWs <- rnorm(1)*sqrt(dt) # for stock
          dWt <- rnorm(1)*sqrt(dt) # for volatiltiy
          
          dnu <- para[1]*(theta-nu[i-1])*dt+para[2]*sqrt(nu[i-1])*dWt # dv 
          ds <- r*s[i-1]*dt + sqrt(nu[i-1])*s[i-1]*dWs # dS 

          # Recording the volatility and price of each time step
          nu[i] <- max(nu[i-1]+dnu,0)
          s[i] <- s[i-1]+ds
        }

        # Calculating the payoff for the options according to its type
        if (Type == "C"){ # call
          f[j] <- max(s[T.day+1]-K,0)*exp(-r*T.day/365)
        } else{ # put
          f[j] <- max(-(s[T.day+1]-K),0)*exp(-r*T.day/365)
        }
      }
      
      
      if(qq==FALSE & plot==FALSE){ 
          # squared difference between predicted and actual price
          (mean(f)-price)^2}else if(qq==TRUE & plot==FALSE){ 
              # simulated ABC Note payoffs
              mean(f)}else if(qq==TRUE & plot==TRUE){ 
                # simulated Bitcoin prices
                s}
  }

  if (qq == TRUE){
      return(f.vec %>% unlist) # return a vector of simulated prices
  }
  return(sum(f.vec %>% unlist)) # return the sum of squared error (SSE)
} 
```

``` r
# Then we find the kappa and xi which minimize the sum of squared error
opm <- optim(par = c(1.3, 9), # Initial Guess
             heston, # SSE between simulated and actual price
             data=Db) # data input
              
k.m1 <- opm$par[1] # optimized kappa by using method 1
x.m1 <- opm$par[2] # optimized xi by using method 1

cat("Kappa:",k.m1,"\n")
```

    ## Kappa: 1.3

``` r
cat("Vol of vol:",x.m1,"\n")
```

    ## Vol of vol: 9.9

QQPlot and Error plot of method 1

``` r
# Q-Q plot of a sample simulated Bitcoin price trials is shown. 

f.pre1 <- heston(para = c(k.m1, x.m1), # kappa and xi calculated by method 1
                data = Db,
                qq = TRUE)
s.pre1<- heston(para = c(k.m1, x.m1), # stock price sample by method 1
                data = Db,
                qq = TRUE,
                plot=TRUE) 
# Q-Q plot
# We plot QQ plot of the simulation 
s.r <- diff(log(s.pre1)) %>% na.omit
qqnorm(s.r)
qqline(s.r, col = 2, lwd = 2) 
```

![](codes_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
# The function plot actual option price against simulated (predicted) price,
# if actual and simulated prices are the same, all points should lie on
# the red line.
# Correlation between actual and simulated price is also returned 
act_pre_plot <- function(f.pre,
                         data,
                         p = 800000, # boundary of the simulated price 
                         a = 800000 # boundary of the simulated price
                         ){
  # Merging actual price and simulated price into a dataframe
  j <- cbind(f.pre, data$price) %>% data.frame 
  colnames(j) <- c("Predicted", "Actual") # renaming the dataframe
  j.s <- subset(j, Predicted < p & Actual < a) # Subsetting simulated price

  cat("Correlation:","\n")
  print(cor(j))
  # Plotting the actual price and predicted price
  ggplot(j.s,aes(x=Predicted,y=Actual)) + 
  geom_point(col = 1) +
  geom_abline(intercept = 0, slope = 1,col=2) +
  ggtitle("Actual Price against Predicted Price")+
  theme(plot.title = element_text(size = 20, hjust = .5)) 
  }
```

``` r
# Plotting actual price against predicted price, as well as their correlation
act_pre_plot(f.pre = f.pre1, data = Db) 
```

    ## Correlation: 
    ##           Predicted    Actual
    ## Predicted 1.0000000 0.4814905
    ## Actual    0.4814905 1.0000000

![](codes_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Method 2: Minimizing the actual 30 day volatility and predicted
volatility from Monte Carlo of Heston Model

As we can estimate the actual volatility by using 30-day volatility, we
tried to find the kappa and xi that minimize the sum of squared error
between simulated volatility by using Heston Model and the actual 30-day
volatility.

``` r
# Dataframes created here will be useful in
# approximating parameters, backtesting, and product variants

# Calculating the strike price in each month (mean of the previous month prices)
short_d <- substr(BTC$Date,1,7) # only get the year and month
btc.K <- cbind(BTC,short_d) 
colnames(btc.K) <- c("Date", "S","short_d")
group <- btc.K %>% group_by(short_d) %>% summarize(K = mean(S)) # monthly mean
group <- cbind(group[,"short_d"],lag(group[,"K"])) %>% na.omit 
# merge mean price with the price data, also omit NA
btc.K <- merge(x = btc.K, y = group, by = "short_d", all.x=TRUE) %>% na.omit 

# Merge strike, riskfree, price, return, and volatility datframe together
# note that BTC.plot is the dataframe for plotting price, ret, vol before
bt.btc <- list(btc.K,
               dtb3,
               BTC.plot[,c("Date","vol")]) %>% reduce(left_join, 
                                                      by = "Date") %>% na.omit
# Preparing dS and dv data in the Heston Model
sv_lag <- bt.btc[c("S","vol")] %>% lead 
colnames(sv_lag) <- c("S1","v1") # naming lagged data
bt.btc <- cbind(bt.btc,sv_lag) %>% na.omit # combine lag and non-lag data
bt.btc$dS <- bt.btc$S1 - bt.btc$S # their difference is dS
bt.btc$dv <- bt.btc$v1 - bt.btc$vol # and dv
n <- length(bt.btc$Date) # Checking the number of rows
bt.btc$dt <- rep(1/252,n) # creating dt column

theta <- (sd(BTC.r$.)/sqrt(1/(252)))^2 # long run volatility
bt.btc$theta <- rep(theta,n) # creating theta column
dt <- 1/252 # time step
bt.btc$dWs <- sqrt(dt)*rnorm(n)/sqrt(1/252) # modelling the stochastic term
bt.btc$dWv <- sqrt(dt)*rnorm(n)/sqrt(1/252) # modelling the stochastic term
colnames(bt.btc) = c("short_d","Date",
                     "S","K","r","v",
                     "S1","v1","dS","dv","dt",
                     "theta","dWs","dWv") # renaming the column names
# Only get data column needed
bt.btc <- bt.btc[c("Date","dS","r","S","K","v","dt","dv","theta","dWs","dWv")] 
bt.btc$r <- bt.btc$r*252 # annualizing the risk-free rate
```

``` r
# The sum of squared error between simulated and actual volatility
heston.v <- function(param,btc){
    sum((param[1]*(btc$theta-btc$v)*btc$dt +
         param[2]*sqrt(btc$v)*btc$dWv*sqrt(btc$dt)-btc$dv)^2)}
```

``` r
# Trying to Monte Carlo the parameters
d <- 10000 # trials
k <- xi <- rep(0,d) # initializing the list for kappa and xi
for (i in 1:d){
  bt.btc$dWs <- bt.btc$dWv <- rnorm(n)/sqrt(1/252) # stochastic term
  result <- optim(par = c(1.3, 9), 
                  heston.v, 
                  btc = bt.btc,
                  lower=c(0,0),
                  method = "L-BFGS-B") # optimization
  k[i] <- result$par[1] # result of the optimization: kappa
  xi[i] <- result$par[2] # result of the optimization: vol of vol
}

k.m2 <- mean(k) # optimized kappa by using method 2
x.m2 <- mean(xi) # optimized xi by using method 2
cat("Kappa:",round(k.m2,4),"\t","Standard Error:",sd(k)/sqrt(d),"\n")
```

    ## Kappa: 5.3106     Standard Error: 0.009679367

``` r
cat("Vol of vol:",round(x.m2,4),"\t","Standard Error:",sd(xi)/sqrt(d))
```

    ## Vol of vol: 0.0022    Standard Error: 3.158368e-05

QQPlot and Error plot of method 2

``` r
# Q-Q plot of a sample simulated Bitcoin price trials is shown. 

f.pre2 <- heston(para = c(k.m2, x.m2), # kappa and xi calculated by method 2
                data = Db,
                qq = TRUE)
s.pre2 <- heston(para = c(k.m2, x.m2), # stock price simulated by method 2
                data = Db,
                qq = TRUE,
                plot=TRUE)


# We plot the QQ plot of the simulation 
s.r2 <- diff(log(s.pre2)) %>% na.omit
qqnorm(s.r2)
qqline(s.r2, col = 2, lwd = 2) 
```

![](codes_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
# Plotting actual price against predicted price, as well as their correlation
act_pre_plot(f.pre = f.pre2, data = Db)
```

    ## Correlation: 
    ##           Predicted   Actual
    ## Predicted  1.000000 0.461307
    ## Actual     0.461307 1.000000

![](codes_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

## Pricing the ABC Note

For each month, we model the daily change of Bitcoin volatility and
prices by using the Heston Model. At the end of each month, we use the
mean price of the previous month as strike price (K), and the final
price of the month as stock price (S) to calculate the corresponding
straddle payoff in that month. We simulate the price path and payoff for
24 months and sum all the monthly straddle payoff together to get the
total payoff the ABC note. We repeat this process for 10000 trials and
take the mean of all trials to get the expected ABC Note payoff, which
is also the approximated option price.

``` r
# Function for pricing the ABC Note
pricing <- function(
    btc=bt.btc,
    ret = BTC.r,
    pricefn = "straddle", # monthly payoff structure: straddle or strangle
    d = 10000, # trials
    T = 2, # 2 years
    m = T*365,
    dt = T/m, # time step
    initial.S = btc$S[length(btc$S)], # most recent bitcoin price
    initial.K = btc$K[length(btc$K)], # most recent month mean price
    s1 = 1, # for strangle variant
    s2 = 1, # for strangle variant
    cnum = 1, # number of call
    pnum = 1, # number of put
    kappa = k.m1, # approximated through method 1
    x = x.m1, # approximated through method 1
    r = btc$r[length(btc$r)], # risk-free rate (3m t-bills)
    theta = (sd(ret$.)/sqrt(1/(365)))^2 # annualized volatility of return
){
  # Payoff Functions 
  call <- function(S,K) max(S-K,0) # European
  put <- function(S,K) max(K-S,0) # European
  # Long cnum call and pnum put at the same strike
  straddle <- function(S, K, cnum = 1, pnum = 1) {
      call(S,K)*cnum + put(S,K)*pnum} 
  # Long pnum put with strike K1 and cnum call with strike K2, whre K2 > K1
  strangle <- function(S, K1, K2, cnum = 1, pnum = 1){
      call(S,K2)*cnum + put(S,K1)*pnum} 

  # Initializing price, volatility, payoff path, and first strike price
  s <- nu <- NULL
  s[1] <- initial.S
  K <- initial.K 
  nu[1] <- theta # most recent 30d vol

  # The Monte Carlo Simulation of product price
  f <- foreach (i = 1:d)%dorng%{
      monthly.payoff <- rep(0,24)
      for (j in 1:24){ # 24 months
        for(h in 2:31){ # assuming 30 days in one month
            c <- h+30*(j-1)
          # stochastic terms
          dW1 <- sqrt(dt)*rnorm(1)
          dW2 <- sqrt(dt)*rnorm(1)

          # modelling the change of bitcoin price and volatility
          ds <- r*s[c-1]*dt+sqrt(nu[c-1])*s[c-1]*dW1
          dnu <- kappa *(theta-nu[c-1])*dt + x*sqrt(nu[c-1])*dW2
          s[c] <- s[c-1]+ds
          nu[c] <- max(nu[c-1]+dnu,0) # Ensure non-negative "nu"
        }
        # for each month. calculate the new K by the mean of last month price
        K <- mean(s[(1+30*(j-1)):(30+30*(j-1))])

        # Calculating the payoff of strangle and other product variants
        if (pricefn == "straddle"){ # Default
          monthly.payoff[j] <- straddle(s[30+30*(j-1)],
                                        K,
                                        cnum,
                                        pnum)*exp(-r*j/12)
                                        # Product Variant 1
                                        }else if(pricefn == "strangle"){
          monthly.payoff[j] <- strangle(s[30+30*(j-1)],
                                        K*s1,
                                        K*s2,
                                        cnum,
                                        pnum)*exp(-r*j/12)
                                        # Product Variant 3
                                        } else if(pricefn == "Asian"){
          monthly.payoff[j] <-straddle(mean((s[1+30*(j-1)])):(s[30+30*(j-1)]),
                                       K,
                                       cnum,
                                       pnum)*exp(-r*j/12)                         
                                        }
      }
      sum(monthly.payoff)
  }  
  return(f %>% unlist)
}
```

``` r
# Pricing ABC Note
f.vec <- pricing()
ABC.price <- mean(f.vec)
cat("Straddle Price:",round(ABC.price,4),"\n")
```

    ## Straddle Price: 275747.1

``` r
cat("Standard Error:",round(sd(f.vec)/sqrt(10000),4),"\n")
```

    ## Standard Error: 92156.08

``` r
# Raw straddle price is too high for retail investors, hence we restructure
# the product and run the backtest
```

## Restructuring the ABC Note

``` r
# Here with demonstrate for each 1000USD, where will it be invested
ca <-  function(target = 1000, # target price
                m.fee = 0.01, # management fee
                g.rate = 0.8, # capital guarantee rate
                rf = bt.btc$r[length(bt.btc$r)], # risk free rate
                rf.bond = 100*exp(-rf*2),
                price,
                printp = FALSE){ # risk free bond

# Participation rate: unit of straddle purchased
p.rate <- (target*(1-m.fee)-g.rate*target*exp(-2*rf))/price
# unit of bonds purchased
b.rate <- (target - p.rate*price)/rf.bond


if (printp == TRUE){
  cat("\n")
  cat("Participation rate:", p.rate,"\n")
  cat("unit in bond purchased:", b.rate,"\n")
  cat("unit in straddles:", p.rate,"\n\n")

  
}
return(c(p.rate,b.rate))
}
```

# Backtesting

In this section, we explore the return of ABC Note if it existed in
2017-2019. We find that the ABC Note (capital allocated) can return an
**annualized return of \~50%** in the investment period 2019-2021

## Backtesting return of ABC Note

``` r
bt.ret <- function(
  year = 4, # number of years before
  s1 = 1, # for strangle variation
  s2 = 1, # for strangle variation
  cnum = 1, # number of call
  pnum = 1, # number of put
  target = 1000, # target price of product
  pricefn = "straddle"){ # pricing function

  for(i in 2:year){ 
    start.d <- date-365*i # model start date
    end.d <- date - 365*(i-2)  # model end date

    # For pricing: simulate the data known at the start date
    bt.btc.i1 <- subset(bt.btc,Date <= start.d) 
    ret.i1 <- bt.btc.i1$S %>% log %>% diff %>% na.omit # calculate return

    # For payoff: only get the last price in each month
    btc.K.i <- subset(btc.K, Date <= end.d & Date>= start.d)
    bt.p <- btc.K.i %>% group_by(short_d) %>% summarise_all(last) 
    r <- bt.btc.i1$r[length(bt.btc.i1$r)] # intial risk free

    # Price the ABC Note at start date
    price <- pricing(
        d = 10000,
        pricefn = pricefn, # type of monthly payoff: straddle/strangle
        initial.S = bt.btc.i1$S[length(bt.btc.i1$S)], # initial Bitcoin price
        initial.K = bt.btc.i1$K[length(bt.btc.i1$K)], # initial strike
        r = r, # intial risk free
        cnum = cnum, # number of call
        pnum = pnum, # number of put
        theta = (sd(ret.i1)/sqrt(1/365))^2 # long run volatility
    ) %>% mean

    if(pricefn == "straddle"|pricefn=="strangle"){
      bt.p$putpayoff <- pmax(bt.p$K*s1-bt.p$S,0) # monthly put payoff
      bt.p$callpayoff <- pmax(bt.p$S-bt.p$K*s2,0) # monthly call payoff 

      # Total payoff = sum of call and put payoff, 
      # missing months in the data are approximated as the mean of other months
      payoff <- (sum(bt.p$putpayoff)+
                 sum(bt.p$callpayoff))*(24/length(bt.p$Date)) 
      ret <- (payoff)/price -1 # Calculating Investors' return
    }else if(pricefn == "Asian"){
        group2 <- btc.K.i %>% group_by(short_d) %>% summarize(S.asian = mean(S))
        bt.p2 <- merge(btc.K.i,group2) 
        bt.p2 <- bt.p2 %>% group_by(short_d) %>% summarise_all(last)
        bt.p2$putpayoff <- pmax(bt.p2$K-bt.p2$S.asian,0) # monthly put payoff
        bt.p2$callpayoff <- pmax(bt.p2$S.asian-bt.p2$K,0) # monthly call payoff 
        payoff <- (sum(bt.p2$putpayoff)+
                   sum(bt.p2$callpayoff))*(24/length(bt.p$Date)) 
        ret <- (payoff)/price -1 # Calculating Investors' return
    }

    print(paste("Investment Period:",
                as.Date(start.d,"%m/%d/%Y"),"-",as.Date(end.d,"%m/%d/%Y")))
    cat("straddle price:",price,"\n")
    cat("straddle payoff:",payoff,"\n")
    a.ret <- (1+ret)^0.5-1 # annualized raw return
    cat("Backtest annualized return of raw straddle: ",": ",
        a.ret*100,"%","\n",sep="")
    
    rates <- ca(rf = bt.btc.i1$r[length(bt.btc.i1$r)],
                price=price,
                printp =TRUE,
                target=target)

    # annualized capital allocated return
    a.ca.ret <- ((payoff*rates[1]+100*rates[2])/target)^0.5-1
    cat("Backtest annualized return of capital allocated ABC Note: ",": ",
        a.ca.ret*100,"%","\n",sep="")
    cat("\n\n")
  }

}
```

``` r
# Backtesting past returns
bt.ret()
```

    ## [1] "Investment Period: 2019-11-18 - 2021-11-17"
    ## straddle price: 18437.06 
    ## straddle payoff: 140210.8 
    ## Backtest annualized return of raw straddle: : 175.7687%
    ## 
    ## Participation rate: 0.01123082 
    ## unit in bond purchased: 8.102179 
    ## unit in straddles: 0.01123082 
    ## 
    ## Backtest annualized return of capital allocated ABC Note: : 54.43125%
    ## 
    ## 
    ## [1] "Investment Period: 2018-11-18 - 2020-11-17"
    ## straddle price: 14008.25 
    ## straddle payoff: 41470.88 
    ## Backtest annualized return of raw straddle: : 72.0599%
    ## 
    ## Participation rate: 0.0153808 
    ## unit in bond purchased: 8.103287 
    ## unit in straddles: 0.0153808 
    ## 
    ## Backtest annualized return of capital allocated ABC Note: : 20.34052%
    ## 
    ## 
    ## [1] "Investment Period: 2017-11-18 - 2019-11-18"
    ## straddle price: 337206.9 
    ## straddle payoff: 41966.14 
    ## Backtest annualized return of raw straddle: : -64.72222%
    ## 
    ## Participation rate: 0.000604935 
    ## unit in bond purchased: 8.10178 
    ## unit in straddles: 0.000604935 
    ## 
    ## Backtest annualized return of capital allocated ABC Note: : -8.590769%

# Product Variants

Here we explore alternatives of the ABC Note:

1.  **Changing all the monthly straddles into strangles:** With this
    change, the price of the product does decrease, however, the payoff
    of the ABC Note decreases even for a larger extent, hence resulting
    in a smaller return.

2.  **Changing ratio of call and put in the straddle:** With this
    change, price decreases while payoff remains in similar levels in
    all backtesting years. Hence, it may be profitable to change the
    ratio of call and put.

3.  **Using mean price as stock price:** This change increases price and
    decreases payoff, hence not recommended.

## Strangle version of ABC Note

``` r
# Function for backtesting return of strangle ver. ABC Note
ret.strng <- function(para,data){
    call <- function(S,K) max(S-K,0) # European
    put <- function(S,K) max(K-S,0) # European
    data$strangle.payoff <- call(data$S,
                                 data$K*(para[2]+para[1]))+
                                 put(data$S,data$K*(para[1]))

data <- subset(data,Date>date-365*2) # only get data after 2019
ret <- data$S %>% log %>% diff %>% na.omit

price <- pricing(
    ret = ret,
    d = 1000,
    pricefn = "strangle", # monthly payoff structure: straddle or strangle
    T = 2, # 2 years
    m = T*365,
    dt = 1/365, # time step
    initial.S = data$S[1], # most recent bitcoin price
    initial.K = data$K[1], # most recent month mean price
    s1 = para[1], # for strangle variant
    s2 = para[1]+para[2], # for strangle variant
    r = data$r[1], # risk-free rate (3m t-bills)
    theta = (sd(ret)/sqrt(1/(365)))^2 # annualized volatility of return
)

    return(sum(data$strangle.payoff)/mean(price)-1)
}

# Find the strike scalars that return the highest return
opm <- optim(par = c(1, 1), # Initial Guess
             ret.strng, # function for payoff for strangle ver. ABC Note
             data=bt.btc, 
             control=list(fnscale=-1), # maximize the function
             lower = c(0.5,0), # lower bound
             upper = c(1,0.5), # upper bound
             method = "L-BFGS-B") # for box constraint maximization


s1 <- opm$par[1]
s2 <- opm$par[1]+opm$par[2]
cat("Optimal Scalar for put strike price:",s1,"\n")
```

    ## Optimal Scalar for put strike price: 1

``` r
cat("Optimal Scalar for call strike price:",s2,"\n")
```

    ## Optimal Scalar for call strike price: 1.5

``` r
# Backtest variant result
bt.ret(year = 4, # backtest year
       s1 = s1, # put scalar obtained before
       s2 = s2, # call scalar obtained before
       pricefn = "strangle") # pricing function
```

    ## [1] "Investment Period: 2019-11-18 - 2021-11-17"
    ## straddle price: 21675.65 
    ## straddle payoff: 43555.16 
    ## Backtest annualized return of raw straddle: : 41.75349%
    ## 
    ## Participation rate: 0.009552811 
    ## unit in bond purchased: 8.102179 
    ## unit in straddles: 0.009552811 
    ## 
    ## Backtest annualized return of capital allocated ABC Note: : 10.73808%
    ## 
    ## 
    ## [1] "Investment Period: 2018-11-18 - 2020-11-17"
    ## straddle price: 18339.96 
    ## straddle payoff: 14789.61 
    ## Backtest annualized return of raw straddle: : -10.19943%
    ## 
    ## Participation rate: 0.01174801 
    ## unit in bond purchased: 8.103287 
    ## unit in straddles: 0.01174801 
    ## 
    ## Backtest annualized return of capital allocated ABC Note: : -0.7993347%
    ## 
    ## 
    ## [1] "Investment Period: 2017-11-18 - 2019-11-18"
    ## straddle price: 34826.2 
    ## straddle payoff: 26163.26 
    ## Backtest annualized return of raw straddle: : -13.3252%
    ## 
    ## Participation rate: 0.005857322 
    ## unit in bond purchased: 8.10178 
    ## unit in straddles: 0.005857322 
    ## 
    ## Backtest annualized return of capital allocated ABC Note: : -1.845806%

## Different ratio of call and put in the straddle

``` r
# Function for calculating strangle ver. ABC Note
ret.strng <- function(para,data){
    call <- function(S,K) max(S-K,0) # European
    put <- function(S,K) max(K-S,0) # European
    data$strangle.payoff <- call(data$S,data$K)*para+put(data$S,data$K)*(1-para)

data <- subset(data,Date>date-365*2) # only get data after 2019
ret <- data$S %>% log %>% diff %>% na.omit

price <- pricing(
    d = 1000,
    pricefn = "straddle", # monthly payoff structure: straddle or strangle
    T = 2, # 2 years
    m = T*365,
    dt = 1/365, # time step
    initial.S = data$S[1], # most recent bitcoin price
    initial.K = data$K[1], # most recent month mean price
    cnum = para,
    pnum = 1 - para,
    r = data$r[1], # risk-free rate (3m t-bills)
    theta = (sd(ret)/sqrt(1/(365)))^2 # annualized volatility of return
)

    return(sum(data$strangle.payoff)/mean(price))
}

# Find the strike scalars that return the highest return
opm <- optim(par = 0.5, # Initial Guess
             ret.strng, # function for payoff for strangle ver. ABC Note
             data=bt.btc, 
             control=list(fnscale=-1), # maximize the function
             lower = 0, # lower bound
             upper = 1, # upper bound
             method = "Brent") # one-dimensional constrained optimization


cnum <- opm$par[1] # ratio of call
pnum <- 1 - opm$par[1] # ratio of put
cat("Optimal ratio of call:",cnum,"\n")
```

    ## Optimal ratio of call: 0.7647753

``` r
cat("Optimal ratio of put:",pnum,"\n")
```

    ## Optimal ratio of put: 0.2352247

``` r
# Backtest variant result
bt.ret(year = 4, # backtest year
       cnum = cnum, # put scalar obtained before
       pnum = pnum, # call scalar obtained before
       pricefn = "straddle") # pricing function
```

    ## [1] "Investment Period: 2019-11-18 - 2021-11-17"
    ## straddle price: 11636.86 
    ## straddle payoff: 140210.8 
    ## Backtest annualized return of raw straddle: : 247.1146%
    ## 
    ## Participation rate: 0.01779375 
    ## unit in bond purchased: 8.102179 
    ## unit in straddles: 0.01779375 
    ## 
    ## Backtest annualized return of capital allocated ABC Note: : 81.79919%
    ## 
    ## 
    ## [1] "Investment Period: 2018-11-18 - 2020-11-17"
    ## straddle price: 6896.314 
    ## straddle payoff: 41470.88 
    ## Backtest annualized return of raw straddle: : 145.2241%
    ## 
    ## Participation rate: 0.03124251 
    ## unit in bond purchased: 8.103287 
    ## unit in straddles: 0.03124251 
    ## 
    ## Backtest annualized return of capital allocated ABC Note: : 45.12005%
    ## 
    ## 
    ## [1] "Investment Period: 2017-11-18 - 2019-11-18"
    ## straddle price: 8201.249 
    ## straddle payoff: 41966.14 
    ## Backtest annualized return of raw straddle: : 126.2088%
    ## 
    ## Participation rate: 0.02487283 
    ## unit in bond purchased: 8.10178 
    ## unit in straddles: 0.02487283 
    ## 
    ## Backtest annualized return of capital allocated ABC Note: : 36.16147%

## Mean stock price

Instead of calculating the monthly straddle payoff by using last price
in month n as stock price, we calculate straddle payoff by using the
mean price of month n as stock price.

``` r
# Backtest variant result
bt.ret(year = 4, # backtest year
       pricefn = "Asian") # pricing function
```

    ## [1] "Investment Period: 2019-11-18 - 2021-11-17"
    ## straddle price: 56085.9 
    ## straddle payoff: 106587.8 
    ## Backtest annualized return of raw straddle: : 37.85641%
    ## 
    ## Participation rate: 0.003691897 
    ## unit in bond purchased: 8.102179 
    ## unit in straddles: 0.003691897 
    ## 
    ## Backtest annualized return of capital allocated ABC Note: : 9.714599%
    ## 
    ## 
    ## [1] "Investment Period: 2018-11-18 - 2020-11-17"
    ## straddle price: 27305.25 
    ## straddle payoff: 28655.03 
    ## Backtest annualized return of raw straddle: : 2.441834%
    ## 
    ## Participation rate: 0.007890722 
    ## unit in bond purchased: 8.103287 
    ## unit in straddles: 0.007890722 
    ## 
    ## Backtest annualized return of capital allocated ABC Note: : 1.805576%
    ## 
    ## 
    ## [1] "Investment Period: 2017-11-18 - 2019-11-18"
    ## straddle price: 37492.54 
    ## straddle payoff: 32481.07 
    ## Backtest annualized return of raw straddle: : -6.922928%
    ## 
    ## Participation rate: 0.005440769 
    ## unit in bond purchased: 8.10178 
    ## unit in straddles: 0.005440769 
    ## 
    ## Backtest annualized return of capital allocated ABC Note: : -0.6571615%

``` r
# This variant, albeit not recommended, kinda exhibits scalability of
# this program as I just add a new "else if" in both the functions "pricing" 
# and "bt.ret" and this variant become a one-liner.
```

``` r
# Check how much time has been spent for running the program
cat("This program needs:","\n")
```

    ## This program needs:

``` r
Sys.time() - start.time 
```

    ## Time difference of 23.96163 mins
