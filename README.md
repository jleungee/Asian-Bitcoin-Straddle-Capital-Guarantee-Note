## Asian Bitcoin Straddles Capital Guaranteed Note (ABC Note)
Co-authors: Leung Sum Ming, Ho Tsz Chiu & Lin Jianan 

Bitcoin is one of the best incarnations of a "High Risk High Return" asset. 
We understand that many retail investors would like to profit from its high 
return yet are deterred by its tedious volatility. Hence, we have designed a 
product, the Asian Bitcoin Straddle Capital Guaranteed Note (ABC Note), which 

1. guarantees your initial capital;

2. mitigates the risk greatly;

3. has low entry barrier; and

4. one can profit no matter price surges or plummets.

5. Have an **Annualized return of 15% in the investment period between 2019-2021**

To briefly describe the ABC notes, it is compromised with **24 long straddles** 
which will be executed **monthly**. For each straddle, its 
**strike price will be set as the mean Bitcoin price in the previous month** 
and all the monthly payoffs will be **accumulated** and given to the investors 
after 24 months.

................................................................................
................................................................................

**Content of this project**
1. README.md
2. Summary Report that detailed the motivation, structure, hedging strategy of this project
3. Codebook that includes: Background Research, 
   Monte Carlo Simulation and Heston Model to price Bitcoin, backtesting, 
   and Product Variants
4. Graphs shown in the codebook (all the unnamed chunks)

................................................................................
................................................................................

**More details of the codebook:**

1. Set-up - Import necessary libraries and data

2. Background research of Bitcoin - demonstrate Bitcoin's characteristics 
quantitatively and graphically to support the use of Heston Model 
in later sections.

3. Monte Carlo Simulation pricing - optimize parameters in the simulation and 
use simulation to price ABC Note

4. Backtesting of the ABC note - Backtest the return of the ABC note from 
2017-2021

5. Product variant - Explore the alternative structures of ABC Note and show 
their backtesting results

**Here are a few things to be noted:**

1. The "current" date in this file is "2021-11-17" as it is when Bitcoin option 
prices are retrieved from outer sources (Deribit).

2. It is not a must to run "Data from outer sources (csv and quantmod)" as it 
is used to display how we get our raw data from csv/quantmod and you can 
directly run "Data from outer sources (rds)". In other words, the inclusion 
of csv files in the "data" file is completely extra and only for your reference.
It will not affect the codes even if you delete all the csv files.


