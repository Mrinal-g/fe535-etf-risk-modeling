library(quantmod)
library(ggplot2)

symbols <- c("IVV", "IJH", "IWF", "IJR", "IWM", "IWD", "IVW", "IWB", 
             "IVE", "IUSG", "IUSV","IWV", "IYW", "IWN", "IWO", 
             "IJK", "IJJ", "IJS", "IJT", "IYH", "IYF", "IYY", "IYK", "IYE", 
             "IYG", "IYJ", "IDU", "IYC", "IYM", "IYZ")

get_func <- function(x) get(getSymbols(x,from = "2010-01-01", to = "2023-09-30"))
list <- lapply(symbols,get_func) 
length(list)

adj_get <- function(x) x[,6]
adj_price <- lapply(list,adj_get)

head(adj_price[[1]])

adj_all <- Reduce(merge,adj_price)
adj <- adj_all["2010-01-01/2023-09-30"]

names(adj) <- symbols

log_ret <- na.omit(log(adj/lag(adj)))
log_ret <- log_ret["2010-01-01/2023-09-30"]

# Calculating the mean return, vol and sharpe ratio
Mu <- 252*apply(log_ret,2,mean)
sigma <- sqrt(252)*apply(log_ret,2,sd)
Sr <- Mu/sigma
results <- data.frame(cbind(Mu,sigma,Sr))
colnames(results) <- c("Mean_Return", "Volatility", "SR")
round(results,4)


summary_stats <- function(x) {
  return(c(mean = mean(x), Q1 = quantile(x, 0.25), median = median(x), Q3 = quantile(x, 0.75)))
}

mean_summary <- summary_stats(Mu)
volatility_summary <- summary_stats(sigma)
SR_summary <- summary_stats(Sr)

# Combine the summaries into a single data frame
summary_table <- data.frame(mean_summary, volatility_summary, SR_summary)
colnames(summary_table) <- c("Mean_Return", "Volatility", "Sharpe Ratio")
summary_table <- t(summary_table)

summary_table <- round(summary_table, 4)
print(summary_table)


plot <- ggplot(results, aes(x=Volatility, y=Mean_Return)) + 
  geom_point() +
  labs(title="Mean Returns vs Volatility", x="Volatility", y="Mean Return")  

print(plot)

##########################

Market_IVV <- log_ret[,1]
ETFs_all <- log_ret[,-1]

Beta <- sapply(ETFs_all,FUN=function(x) cov(x,Market_IVV)/var(Market_IVV))

jensen_alpha <- apply(ETFs_all,2,FUN=function(x) mean(x)*252)-Beta*mean(Market_IVV)*252

Treynor_ratio <- sapply(ETFs_all,FUN=function(x) mean(x)*252)/Beta

Track.Error <- sapply(ETFs_all,FUN=function(x) sqrt(var(x-Market_IVV)*sqrt(252)))

Info_ratio <- (sapply(ETFs_all,FUN=function(x) mean(x)*252)-mean(Market_IVV)*252) /Track.Error

results_R <- data.frame(cbind(Beta,jensen_alpha,Treynor_ratio,Track.Error,Info_ratio))
colnames(results_R) <- c("Beta", "Jensen_Alpha", "Treynor_Ratio","Tracking_Error","Information_Ratio")
round(results_R,4)

summary_stats_ratios <- function(x) {
  return(c(mean = mean(x), Q1 = quantile(x, 0.25), median = median(x), Q3 = quantile(x, 0.75)))
}

jensen_alpha_summary <- summary_stats_ratios(jensen_alpha)
Beta_summary <- summary_stats_ratios(Beta)
Treynor_ratio_summary <- summary_stats_ratios(Treynor_ratio)
Tracking_Error_summary <- summary_stats_ratios(Track.Error)
Info_ratio_summary <- summary_stats_ratios(Info_ratio)


summary_table_ratios <- data.frame(jensen_alpha_summary, Beta_summary,Treynor_ratio_summary, Tracking_Error_summary, Info_ratio_summary)


summary_table_ratios <- t(summary_table_ratios)


summary_table_ratios <- round(summary_table_ratios, 4)


print(summary_table_ratios)

# Relative Performance 
max(Treynor_ratio)
min(Treynor_ratio)
max(Track.Error)
min(Track.Error)
max(Info_ratio)
min(Info_ratio)


# Task 4 
data <- data.frame(Return = Mu[-1], Beta, Symbols = symbols[-1])

# Task 4 
data <- data.frame(Return = Mu[-1], Beta, Symbols = symbols[-1])

# Scatter plot 
plot_1 <- ggplot(data = data, aes(x=Beta, y=Return)) +
  geom_point(color = "blue") +
  labs(x = "Beta" , y = "Annualized Mean Return", title = "Mean Returns against Beta") +
  geom_text(aes(label = Symbols), hjust = 0.5, vjust = -0.5, size = 3) +
  geom_smooth(method = "lm", color = "yellow")

# Display the plot
print(plot_1)


################
# ETFs IVV, IYW and IYF

symbols <- c("IVV", "IYW", "IYF")
log_rets <- log_ret[,symbols]

Mu_2 <- Mu[symbols]
sigma_2 <- sigma[symbols]
Sr_2 <- Mu_2/sigma_2
result_2 <- data.frame(cbind(Mu_2,sigma_2,Sr_2))
colnames(result_2) <- c("Mean_ret", "Volatitliy", "SR")
round(result_2,4)

sig_mat <- data.matrix(var(log_rets)*252)


#############
w_function <- function(weights) {
  w_vec <- matrix(weights,3,1)
  mu_p <- t(w_vec)%*%Mu_2
  sig_p <- sqrt(t(w_vec)%*%sig_mat%*%w_vec)
  result <- c(mu_p,sig_p)
  return(result)
}

vec_ones <- rep(1,nrow(sig_mat))
w_0 <- solve(sig_mat)%*%vec_ones
w_0 <- w_0/sum(w_0)
B_mat <- solve(sig_mat)%*%( diag(vec_ones) - vec_ones%*%t(w_0))
w_1 <- B_mat%*%Mu_2
x<-as.numeric(t(vec_ones)%*%solve(sig_mat)%*%(Mu_2))
Sr_p<-(solve(sig_mat)%*%(Mu_2))/x
imp<-sapply(data.frame(w_0,w_1,Sr_p),w_function)
imp<-data.frame(imp)
imp[3,]<-imp[1,]/imp[2,]
rownames(imp) <- c("mu","sig","SR")

w_A_function <- function(A){
  w_vec <- w_0 + (1/A)*w_1
  mu_p <- t(w_vec)%*%Mu_2
  sig_p <- sqrt(t(w_vec)%*%sig_mat%*%w_vec)
  result <- c(mu_p,sig_p)
  return(result)
}

mu_0<-t(w_0)%*%Mu_2
mu_0<-as.numeric(mu_0)
m<-seq(mu_0,2*max(Mu_2),length.out = 100)
A_seq <-(m-t(w_0)%*%Mu_2)/(t(w_1)%*%Mu_2)
A_seq<-1/A_seq
sample(A_seq,5)

# compute optimal portfolios (MVEF) metrics
ds_A <- t(sapply(A_seq,w_A_function))
ds_A <- data.frame(ds_A)
names(ds_A) <- c("mu_p","sig_p")
ds_A$Sr <- (ds_A$mu_p)/ds_A$sig_p

plot(mu_p ~ sig_p,data = ds_A,
     type = "l", ylab = expression(mu[p]),
     xlab = expression(sigma[p]),
     xlim = range(ds_A$sig_p),
     ylim = range(ds_A$mu_p))
points(mu_p~sig_p,data = ds_A[which(ds_A$sig_p == imp['sig','w_0']),],
       col = 1,pch = 20,cex = 1.5)
points(mu_p~sig_p,data = ds_A[which.max(ds_A$Sr),],
       col = 1,pch = 20,cex = 1.5)
grid(10)

lambda<-seq(-1,1,by=0.001)
w_weight<-matrix(NA,2001,3)
for (i in lambda){
  w_weight[which(lambda==i),]<-i*w_0+(1-i)*Sr_p
}

mvef<-apply(w_weight, 1, w_function)
rownames(mvef) <- c("mu_p","sig_p")
mvef<-t(mvef)
mvef<-data.frame(mvef)

plot(mu_p ~ sig_p,data = ds_A,
     type = "l", ylab = expression(mu[p]),
     xlab = expression(sigma[p]),
     xlim = range(ds_A$sig_p),
     ylim = range(ds_A$mu_p))
points(mu_p~sig_p,data = ds_A[which(ds_A$sig_p == imp['sig','w_0']),],
       col = 1,pch = 20,cex = 1.5)
points(mu_p~sig_p,data = ds_A[which.max(ds_A$Sr),],
       col = 1,pch = 20,cex = 1.5)
grid(10)
lines(mu_p ~ sig_p,data = mvef,col = 'red',lty = 2,lwd = 2)

lambda<-seq(-1,1,by=0.001)
w_weight<-matrix(NA,2001,3)
for (i in lambda){
  w_weight[which(lambda==i),]<-(1-i)*Sr_p
}

mvef_0<-apply(w_weight, 1, w_function)
rownames(mvef_0) <- c("mu_p","sig_p")
mvef_0<-t(mvef_0)
mvef_0<-data.frame(mvef_0)

plot(mu_p ~ sig_p,data = ds_A,
     type = "l", ylab = expression(mu[p]),
     xlab = expression(sigma[p]),
     xlim = range(ds_A$sig_p),
     ylim = range(ds_A$mu_p))
points(mu_p~sig_p,data = ds_A[which(ds_A$sig_p == imp['sig','w_0']),],
       col = 1,pch = 20,cex = 1.5)
points(mu_p~sig_p,data = ds_A[which.max(ds_A$Sr),],
       col = 1,pch = 20,cex = 1.5)
grid(10)
lines(mu_p ~ sig_p,data = mvef_0,col = 'blue',lty = 2,lwd = 2)


model<-lm(mu_p ~ sig_p,mvef_0)
summary(model)



values <- seq(-1.5, 1.5, by = 0.01)
valid_combinations <- list()
for (i in 1:length(values)) {
  for (j in 1:length(values)) {
    for (k in 1:length(values)) {
      weight_vector <- c(values[i], values[j], values[k])
      if (sum(weight_vector) == 1) {
        valid_combinations <- c(valid_combinations, list(weight_vector))
      }
    }
  }
}
df <- as.matrix(do.call(rbind, valid_combinations))

# compute their mean returns and volatilities
results <- apply(df, 1, w_function)
rownames(results) <- c("mu_p","sig_p")
results<-as.matrix(results)
results<-t(results)
results<-data.frame(results)


plot(mu_p ~ sig_p, data = results,
     type = "p", ylab = expression(mu[p]),
     xlab = expression(sigma[p]),
     xlim = range(results$sig_p),
     ylim = range(results$mu_p),
     col = "yellow")

result_0 <- data.frame()
unique_values <- unique(round(results$mu_p,3))

for (value in unique_values) {
  subset_df <- results[round(results$mu_p,3) == value, ]
  min_value <- round(min(subset_df$sig_p),3)
  result_0 <- rbind(result_0, subset_df[round(subset_df$sig_p,3) == min_value, ])
}

plot(mu_p ~ sig_p, data = result_0,
     type = "b", ylab = expression(mu[p]),
     xlab = expression(sigma[p]),
     xlim = range(result_0$sig_p),
     ylim = range(result_0$mu_p),
     col = "red")
grid(10)

result_df2 <- data.frame()
unique_values <- unique(round(result_0$sig_p,2))

for (value in unique_values) {
  subset_df <- result_0[round(result_0$sig_p,2) == value, ]
  max_value <- round(max(subset_df$mu_p),2)
  result_df2 <- rbind(result_df2, subset_df[round(subset_df$mu_p,2) == max_value, ])
}

plot(mu_p ~ sig_p, data = result_df2,
     type = "p", ylab = expression(mu[p]),
     xlab = expression(sigma[p]),
     xlim = range(result_0$sig_p),
     ylim = range(result_0$mu_p),
     col = "blue") 
grid(10)


plot(mu_p ~ sig_p, data = ds_A,
     type = "l", ylab = expression(mu[p]),
     xlab = expression(sigma[p]),
     xlim = range(ds_A$sig_p),
     ylim = range(ds_A$mu_p),
     col = "blue")  
lines(mu_p ~ sig_p, data = result_df2, type = "p", col = "red")
grid(10)

#Q3 

#Part 1 
#### In this experiment we are trying to find the cost we bear per coin toss which is k. 

#### Firstly we declare a function coin_toss which runs a while loop and simulates a single run of the game,
#### it counts the number of tosses until three consecutive heads are seen.

coin_toss <- function(n) {
  heads<- 0
  tosses <- 0
  
  while (heads< 3) {
    toss <- sample(c(0, 1), size = 1)
    tosses <- tosses + 1
    
    if (toss == 1) {
      heads<- heads+ 1
    } else {
      heads<- 0
    }
  }
  
  return(tosses)
}

#### Simulations part here we use a function replicate which runs the coin_toss function 100000 times giving us large enough
#### sample set to estimate the value of X which is the average number of tosses to get three consecutive heads.

# simulations 
num_simulations <- 100000
tosses_data <- replicate(num_simulations, coin_toss(n = num_simulations))

#### In the follwoing part we compute the expected value of X and solve for k uisng 1 − E[X] × k = 0.

# expected value of X
E_X <- mean(tosses_data)
E_X
# Solve for k such that E[P&L] = 0
k <- 1 / E_X
print(k)

# Part 2 
#### We write a function which simulates a single run for the experiment. 
tur_exp <- function(num_turt) {
  turtles <- sample(num_turt)
  groups <- unique(cummin(turtles))
  return(length(groups))
}
# We chose the number of turtles to be 100 and to sample them we replicate it 10^5 times. Finally we take the mean to 
# obtain the average number of groups.
num_turt <- 100
num_groups <- replicate(10^5, tur_exp(num_turt))
mean_num_groups <- mean(num_groups)
mean_num_groups


# Part 3 
#### In the first line we create a vector of all the prices at P2 step. 
price_2 <- c(120,100,80)
#### We define the probabilities leading to P2.
prob <- c(0.55*0.55,0.45*0.55*2,0.45*0.45)
#### We calculate the expectation of P2 and Variance in the next steps.
exp <- sum(price_2*prob)
var <- sum((price_2-exp)^2*prob)
print(c(exp,var))


# Function to simulate the price change based on binomial distribution
simulate_price_change <- function(p_up, change_up, change_down) {
  if (rbinom(1, 1, p_up) == 1) {
    return(change_up)
  } else {
    return(change_down)
  }
}

# Function to simulate the price over 'n' steps
simulate_price_in_n_steps <- function(initial_price, n, p_up, change_up, change_down) {
  for (i in 1:n) {
    initial_price <- initial_price + simulate_price_change(p_up, change_up, change_down)
  }
  return(initial_price)
}

run_simulation <- function(initial_price, n, times, p_up, change_up, change_down) {
  final_prices <- replicate(times, simulate_price_in_n_steps(initial_price, n, p_up, change_up, change_down))
  mean_price <- mean(final_prices)
  var_price <- var(final_prices)
  return(list(mean = mean_price, variance = var_price))
}

# Parameters
initial_price <- 100
num_steps <- 2  # For M = 10 steps
num_simulations <- 10000
p_up <- 0.55  # Probability of price going up
change_up <- 10  # Price change when it goes up
change_down <- -10  # Price change when it goes down

# Run the simulation
simulation_results <- run_simulation(initial_price, num_steps, num_simulations, p_up, change_up, change_down)

# Output the results
cat("Expected price E[P10] =", simulation_results$mean, "\n")
cat("Variance of price V[P10] =", simulation_results$variance, "\n")

initial_price <- 100
num_steps <- 10  # For M = 10 steps
num_simulations <- 10000
p_up <- 0.55  # Probability of price going up
change_up <- 10  # Price change when it goes up
change_down <- -10  # Price change when it goes down

simulation_results <- run_simulation(initial_price, num_steps, num_simulations, p_up, change_up, change_down)

# Output the results
cat("Expected price E[P10] =", simulation_results$mean, "\n")
cat("Variance of price V[P10] =", simulation_results$variance, "\n")


# Q4 

# Part a 

library(quantmod)

getSymbols('IVV', from = "2010-01-01", to = "2023-09-30")

data <- to.monthly(IVV)
adj_p <- data$IVV.Adjusted
S0 <- data$IVV.Adjusted[[1]]
log_ret1 <-  na.omit(log(adj_p/lag(adj_p)))

log_ret <- mean(log_ret1) # log returns 

sigma <- sd(log_ret1)* sqrt(12)
mu <- mean(log_ret1) * 12 + sigma^2 / 2
table <- data.frame(mu,sigma)
print(table)

# Part b  

# Define the time increment and the GBM function
dt <- 1 / 12

gbm<- function(n) {
  drt <- rnorm(164, (mu - sigma^2/2) * dt, sigma * sqrt(dt))
  St <- S0 * exp(cumsum(drt))
  return(St)
}

# Simulate 1000 paths of GBM
s_matrix <- sapply(1:1000, gbm)

# Extract the simulated prices at month 164
simu <- s_matrix[164, ]

simu_expe <- mean(simu)


simu_expe_true <- S0 * exp(mu * 164 / 12)

simu_sigma <- sd(simu)

simu_sigma_true<- sqrt((exp((sigma^2 * 164 )/ 12) - 1) * S0^2 * exp((2 * mu * 164) / 12))


simu_vs_true <- data.frame(Mean = c(simu_expe, simu_expe_true), Sigma = c(simu_sigma, simu_sigma_true))
rownames(simu_vs_true) <- c('Simulation', 'True Value')
simu_vs_true


# Part c - 

norm<-function(p){
  norm<-c()
  for (i in 1:1000){
    norm<- c(norm, sum(abs(s_matrix[,i]-adj_p)^p)^(1/p))
  }
  norm<-matrix(norm,nrow=1000)  
  norm<-data.frame(norm)
  return(norm)
}

norm_2<-norm(p=2)
head(norm_2)


# Verifying the lowest distance from the true simulation
a<-sort(unlist(norm_2))[1]
i<-as.numeric(gsub('[A-z]',"",names(a)))
num_months <- length(adj_p)
# Define the GBM function with an argument for the number of periods
gbm <- function(S0, mu, sigma, dt, num_periods) {
  
  drt <- rnorm(num_periods - 1, mean = (mu - sigma^2 / 2) * dt, sd = sigma * sqrt(dt))
  
  St <- c(S0, S0 * exp(cumsum(drt)))
  return(St)
}


s_matrix <- sapply(1:1000, function(x) gbm(S0, mu, sigma, dt, num_months - 1))

s_matrix <- rbind(rep(S0, 1000), s_matrix)


# Creating the date sequence to match the number of simulation steps
date_seq <- seq.Date(from = as.Date("2010-02-01"), length.out = num_months, by = "months")

# Assuming you've correctly found the closest path index
closest_path_index <- which.min(unlist(norm_2))

# Now create the xts object
x <- xts(s_matrix[, closest_path_index], order.by = date_seq)

# Plot 
plot(x, type="l", main="Price Path: Simulated vs True", col = 'red', xlab = "Time", ylab = "Price")
lines(adj_p, col = 'green')

# Part d -

quant_99 <- quantile(simu*100,0.01)
quant_99
VaR <- mean(simu)*100 - quant_99
print(VaR)

# Part e - 


num_periods <- 165

# Define the GBM function 
gbm <- function(S0, mu, sigma, dt, num_periods) {
  drt <- rnorm(num_periods - 1, mean = (mu - sigma^2 / 2) * dt, sd = sigma * sqrt(dt))
  St <- S0 * exp(cumsum(c(0, drt)))  # c(0, drt) ensures that St starts with S0
  return(St)
}

# Define the range for sigma
sigma_range <- seq(0.10, 0.50, by = 0.01)

# Initialize a vector to store the VaR for each sigma
VaR_values <- numeric(length(sigma_range))

# Loop over the sigma range to calculate VaR for each value
for (j in 1:length(sigma_range)) {
  # Update sigma
  current_sigma <- sigma_range[j]
  # Simulate the GBM with the updated sigma
  s_matrix <- sapply(1:1000, function(x) gbm(S0, mu, current_sigma, dt, num_periods))
  
  # Compute the 99% VaR for each simulation
  simu_losses <- sort(S0 - s_matrix[num_periods, ])  # Assuming end of period is at num_periods
  VaR_values[j] <- -simu_losses[ceiling(0.01 * 1000)]  # 99% VaR
}

# Plot VaR against sigma
plot(sigma_range, VaR_values, type = 'l', col = 'blue', xlab = 'Sigma', ylab = '99% VaR', main = '99% VaR as a Function of Sigma')

# Part f - 
# Historical - 
quant <- quantile(log_ret1,0.01)
print(quant)

VaR <- log_ret - quant
print(VaR)

# Parametric - 

# Parametric approach

para<-matrix(NA,nrow=164,ncol=1000)
for (i in 1:1000){
  for (j in 2:165){
    para[j-1,i]<-log(s_matrix[j,i]/s_matrix[j-1,i])}}
quant <- quantile(para,0.01)
VaR_par <- mean(para) - quant
print(VaR_par)




simulated_log_returns <- matrix(NA, nrow = 164, ncol = 1000)
for (i in 1:1000) {
  simulated_log_returns[, i] <- diff(log(s_matrix[, i]))
}

simulated_log_returns_vector <- as.vector(simulated_log_returns)


historical_density <- density(log_ret1)

simulated_density <- density(simulated_log_returns_vector)


plot(historical_density, main = "Density of Historical vs Simulated Returns",
     xlab = "Log Returns", ylab = "Density", col = "blue", lwd = 2)
lines(simulated_density, col = "red", lwd = 2)
legend("topright", legend = c("Historical", "Simulated"), col = c("blue", "red"), lwd = 2)
