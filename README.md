## Asian Bitcoin Straddles Capital Guaranteed Note (ABC Note)

**This codebook have five sections:**

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


