library(shiny)
library(shinyjs)
library(tidyverse)
library(tidyquant)
library(quadprog)
library(quantmod)

# User Interface
ui <- fluidPage(
  useShinyjs(),  
  titlePanel("Portfolio Optimization"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("stocks", "Enter Stock Symbols (comma-separated):", value = "RELIANCE.NS,TATAMOTORS.NS,SBIN.NS"),
      dateInput("train_from_date", "Training Period - From Date:", value = "2019-01-01"),
      dateInput("train_to_date", "Training Period - To Date:", value = "2021-01-01"),
      dateInput("test_from_date", "Backtesting Period - From Date:", value = "2021-01-01"),
      dateInput("test_to_date", "Backtesting Period - To Date:", value = "2022-01-01"),
      numericInput("rp", "Return Preference for Markowitz Model(monthly):", value = 0.02, min = 0, step = 0.001),
      selectInput("model", "Choose Portfolio Model:", 
                  choices = c("Markowitz Model", "Equal Weights Model", "Compare Both")),
      selectInput("shrinkage_method", "Covariance Shrinkage Method:",
                  choices = c("None", "eigenvalue shrinkage", "OAS")),
      actionButton("optimize", "Run Portfolio Optimization")
    ),
    
    mainPanel(
      tableOutput("weights"),
      tableOutput("results"),
      tableOutput("backtest_results")
    )
  )
)

# Backend
server <- function(input, output, session) {
  observe({
    # Disable shrinkage dropdown if Equal Weights Model is selected
    if (input$model == "Equal Weights Model") {
      disable("shrinkage_method")
    } else {
      enable("shrinkage_method")
    }
  })
  
  observeEvent(input$optimize, {
    
    stocks <- strsplit(input$stocks, ",")[[1]] %>% trimws()
    train_from_date <- as.character(input$train_from_date)
    train_to_date <- as.character(input$train_to_date)
    test_from_date <- as.character(input$test_from_date)
    test_to_date <- as.character(input$test_to_date)
    
    
    train_data <- tq_get(stocks, get = "stock.prices", complete_cases = TRUE, from = train_from_date, to = train_to_date)
    test_data <- tq_get(stocks, get = "stock.prices", complete_cases = TRUE, from = test_from_date, to = test_to_date)
    
    
    train_returns <- train_data %>%
      select(symbol, date, adjusted) %>%
      group_by(symbol) %>%
      tq_transmute(adjusted, mutate_fun = periodReturn, period = "monthly", col_rename = "Monthly_Return") %>%
      pivot_wider(names_from = 'symbol', values_from = 'Monthly_Return') %>%
      na.omit()
    
    test_returns <- test_data %>%
      select(symbol, date, adjusted) %>%
      group_by(symbol) %>%
      tq_transmute(adjusted, mutate_fun = periodReturn, period = "monthly", col_rename = "Monthly_Return") %>%
      pivot_wider(names_from = 'symbol', values_from = 'Monthly_Return') %>%
      na.omit()
    
    train_matrix <- data.matrix(train_returns[,-1])
    test_matrix <- data.matrix(test_returns[,-1])
    mean_returns <- colMeans(train_matrix, na.rm = TRUE)
    
    
    shrinkage_method <- input$shrinkage_method
    if (shrinkage_method == "None") {
      sig <- cov(train_matrix)
    } else if (shrinkage_method == "eigenvalue shrinkage") {
      library(nlshrink)
      sig <- nlshrink_cov(test_matrix)
    } else if (shrinkage_method == "OAS") {
      sample_cov <- cov(test_matrix)
      p <- ncol(train_matrix)
      n <- nrow(train_matrix)
      target_cov <- diag(mean(diag(sample_cov)), p)
      trace_sample <- sum(diag(sample_cov))
      norm_sample <- sum(sample_cov^2)
      lambda <- ((1 - (2 / p)) * trace_sample^2 + trace_sample^2) / 
        (n * (norm_sample - trace_sample^2 / p))
      lambda <- max(0, min(1, lambda))
      sig <- (1 - lambda) * sample_cov + lambda * target_cov
    }
    
    
    n <- ncol(train_matrix)
    results <- NULL
    weights_side_by_side <- data.frame(Stock = colnames(train_matrix))
    backtest_results <- NULL
    
    
    if (input$model == "Equal Weights Model" || input$model == "Compare Both") {
      weights_eq <- rep(1 / n, n)
      portfolio_returns_eq <- train_matrix %*% weights_eq
      cumulative_return_eq <- prod(1 + portfolio_returns_eq) - 1
      annualized_return_eq <- (1 + cumulative_return_eq)^(12 / nrow(train_matrix)) - 1
      stdev_eq <- sd(portfolio_returns_eq) * sqrt(12)  # Annualized standard deviation
      
      reward_risk_eq <- annualized_return_eq / stdev_eq  # Reward-to-Risk Ratio
      
      results <- rbind(results, data.frame(
        Model = "Equal Weights",
        Train_Annualized_Stdev = round(stdev_eq, 6),
        Train_Annualized_Return = paste0(round(annualized_return_eq * 100, 2), "%"),
        Train_Reward_Risk = round(reward_risk_eq, 4)
      ))
      
      
      backtest_eq <- test_matrix %*% weights_eq
      cumulative_return_test_eq <- prod(1 + backtest_eq) - 1
      annualized_return_test_eq <- (1 + cumulative_return_test_eq)^(12 / nrow(test_matrix)) - 1
      stdev_test_eq <- sd(backtest_eq) * sqrt(12)  # Annualized standard deviation
      reward_risk_test_eq <- annualized_return_test_eq / stdev_test_eq  # Reward-to-Risk Ratio for Backtest
      
      backtest_results <- rbind(backtest_results, data.frame(
        Model = "Equal Weights",
        Test_Annualized_Stdev = round(stdev_test_eq, 6),
        Test_Annualized_Return = paste0(round(annualized_return_test_eq * 100, 2), "%"),
        Test_Reward_Risk = round(reward_risk_test_eq, 4)
      ))
      
      weights_side_by_side$Equal_Weights <- round(weights_eq, 4)
    }
    
    
    if (input$model == "Markowitz Model" || input$model == "Compare Both") {
      rp <- input$rp
      dvec <- mean_returns
      Dmat <- sig
      Amat <- cbind(rep(1, n), mean_returns)
      bvec <- c(1, rp)
      solution <- solve.QP(Dmat, -dvec, Amat, bvec, meq = 1)
      weights_mk <- solution$solution
      portfolio_returns_mk <- train_matrix %*% weights_mk
      cumulative_return_mk <- prod(1 + portfolio_returns_mk) - 1
      annualized_return_mk <- (1 + cumulative_return_mk)^(12 / nrow(train_matrix)) - 1
      stdev_mk <- sd(portfolio_returns_mk) * sqrt(12)  # Annualized standard deviation
      
      reward_risk_mk <- annualized_return_mk / stdev_mk  # Reward-to-Risk Ratio
      
      results <- rbind(results, data.frame(
        Model = "Markowitz",
        Train_Annualized_Stdev = round(stdev_mk, 6),
        Train_Annualized_Return = paste0(round(annualized_return_mk * 100, 2), "%"),
        Train_Reward_Risk = round(reward_risk_mk, 4)
      ))
      
    
      backtest_mk <- test_matrix %*% weights_mk
      cumulative_return_test_mk <- prod(1 + backtest_mk) - 1
      annualized_return_test_mk <- (1 + cumulative_return_test_mk)^(12 / nrow(test_matrix)) - 1
      stdev_test_mk <- sd(backtest_mk) * sqrt(12)  # Annualized standard deviation
      reward_risk_test_mk <- annualized_return_test_mk / stdev_test_mk  # Reward-to-Risk Ratio for Backtest
      
      backtest_results <- rbind(backtest_results, data.frame(
        Model = "Markowitz",
        Test_Annualized_Stdev = round(stdev_test_mk, 6),
        Test_Annualized_Return = paste0(round(annualized_return_test_mk * 100, 2), "%"),
        Test_Reward_Risk = round(reward_risk_test_mk, 4)
      ))
      
      weights_side_by_side$Markowitz <- round(weights_mk, 4)
    }
    
    
    output$weights <- renderTable({ weights_side_by_side }, digits = 4)
    output$results <- renderTable({ results }, digits = 6)
    output$backtest_results <- renderTable({ backtest_results }, digits = 6)
  })
}

 
shinyApp(ui = ui, server = server)
