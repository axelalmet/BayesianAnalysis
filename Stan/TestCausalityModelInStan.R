library(rstan)
library(ggplot)
library(bayesplot)
library(rstanarm)
rstan_options(auto_write = TRUE)

### Bayesian Causality model proposed by Baldi and Shahbaba

# Import the stan model
model1 <- stan_model(file = "~/Documents/BayesianAnalysis/Stan/Models/causalitymodel.stan", verbose = FALSE)
model2 <- stan_model(file = "~/Documents/BayesianAnalysis/Stan/Models/contingencytables.stan", verbose = FALSE)

# Define the parameters
# table_data <- list(K1 = 12, N1 = 100, K2 = 30,  N = 200, N_rep = 100, u_prior = c(1, 1), v_prior = c(1, 1))
table_data1 <- list(K1 = 30, N1 = 100, K2 = 30,  N = 200, u_prior = c(1, 1), v_prior = c(1, 1))
table_data2 <- list(K1 = 30, N1 = 100, K2 = 30,  N = 200, theta1_prior = c(1, 1), theta2_prior = c(1, 1))

table1 <- matrix(c(70, 30, 70, 30), 2, 2, byrow=TRUE)
bf = contingencyTableBF(table1, sampleType = "indepMulti", fixedMargin = "rows")


# Fit the model now
fit1 <- sampling(model1,
            data = table_data1,    # named list of data
            chains = 4,             # number of Markov chains
            warmup = 1000,          # number of warmup iterations per chain
            iter = 2000,            # total number of iterations per chain
            cores = 1,              # number of cores (could use one per chain)
            refresh = 0)             # no progress shown

fit2 <- sampling(model2,
                 data = table_data2,    # named list of data
                 chains = 4,             # number of Markov chains
                 warmup = 1000,          # number of warmup iterations per chain
                 iter = 2000,            # total number of iterations per chain
                 cores = 1,              # number of cores (could use one per chain)
                 refresh = 0)             # no progress shown


output <- sbc(model, data = table_data, M = 500, cores = 1, refresh = 0)
posteriors <- extract(fit, pars = c("u", "v", "K1_pred", "K2_pred"))

### More classical treatment effects model
library(MASS)

# Simulate the data
# Basic setup
set.seed(123456)
N <- 500       # number of observations
alpha <- 1.0   # intercept in the Y model
tau <- 0.25    # treatment effect

# The assignment mechanism
N_t <- 200                                  # number of treated units
W <- sample(rep(c(0, 1), c(N - N_t, N_t)))  # binary treatment variable
ii_t <- which(W == 1); ii_c <- which(W == 0) # index arrays for treatment variable

# The science
mu_c <- alpha + 0*tau; sd_c <- 1   # mean and SD for the control
mu_t <- alpha + 1*tau; sd_t <- 1   # mean and SD for the treated

rho <- 0.0                                   # correlation between the errors
cov_mat <- rbind(c(sd_c^2, rho*sd_c*sd_t),     
                 c(rho*sd_c*sd_t, sd_t^2))   # variance-covariance matrix

science <- mvrnorm(n = N, mu = c(mu_c, mu_t), Sigma = cov_mat, empirical = TRUE)
Y0 <- science[, 1]        # potential outcome if W = 1
Y1 <- science[, 2]        # potential outcome if W = 0
tau_unit <- Y1 - Y0       # unit-level treatment effect

# The realization of potential outcomes by the assignment mechanism
Y_obs <- Y0 * (1 - W) + Y1 * W
Y_mis <- Y0 * W + Y1 * (1 - W)

# Actually simulate the model in Stan

# Collect data into a list format suitable for Stan
model3 <- stan_model(file = "~/Documents/BayesianAnalysis/Stan/Models/causaleffects.stan", verbose = FALSE)
stan_data <- list(N = N, y = Y_obs, w = W, rho = 0.0)

# Compile and run the stan model
fit_simdat <- sampling(model3,
                   data = stan_data,
                   chains = 4,             # number of Markov chains
                   warmup = 1000,          # number of warmup iterations per chain
                   iter = 2000,            # total number of iterations per chain
                   cores = 1,              # number of cores (could use one per chain)
                   refresh = 0)             # no progress shown
