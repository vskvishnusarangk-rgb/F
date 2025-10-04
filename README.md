This project is an interactive Shiny web application for portfolio optimization. 
It allows users to input stock symbols, define training and testing periods, and run different portfolio models such as Markowitz Optimization and Equal Weights Model. 
Users can also compare both models and analyze backtesting results.

🚀 Features

📈 Stock Data Fetching: Automatically retrieves stock data using tidyquant and quantmod.

⚖️ Portfolio Models:

Markowitz Mean-Variance Optimization,
 Equal Weights Model, 
Side-by-side comparison


🧮 Covariance Shrinkage Methods:

None (sample covariance), 
Eigenvalue shrinkage (nlshrink), 
Oracle Approximating Shrinkage (OAS)

📊 Outputs:

Optimized Portfolio Weights, 
Training Results (Annualized Return, Volatility, Reward-to-Risk), 
Backtesting Results (Annualized Return, Volatility, Reward-to-Risk)
