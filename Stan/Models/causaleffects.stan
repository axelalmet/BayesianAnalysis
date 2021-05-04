data {

    int<lower=0> N; // Sample size
    vector[N] y; // Observed outcome
    vector[N] w; // Assigned treatments
    real<lower=-1, upper=1> rho; // Assumed correlation between potential outcomes (I think we can just say zero?)

}
parameters {

    real alpha; // Intercept
    real tau; // Super-population average treatment effect (what we really care about)
    real<lower=0> sigma_c; // residual SD for the control group
    real<lower=0> sigma_t; // residual SD for the treatment group

}
model {

    // Set the priors
    alpha ~ normal(0, 5);
    tau ~ normal(0, 5);
    sigma_c ~ normal(0, 5); 
    sigma_t ~ normal(0, 5);

    // Likelihood (regression)
    y ~ normal(alpha + tau * w, sigma_t * w + sigma_c * (1 - w));

}
generated quantities {

    real tau_fs; // Finite-sample average treatment effect
    real y0[N]; // Potential outcome if W = 0
    real y1[N]; // Potential outcome if W = 1
    real tau_unit[N]; // Unit-level treatment effect

    for (n in 1:N) {

        real mu_c = alpha;
        real mu_t = alpha + tau;

        if (w[n] == 1) {
            y0[n] = normal_rng(mu_c + rho * (sigma_c / sigma_t) * (y[n] - mu_t), sigma_c * sqrt(1 - rho^2));
            y1[n] = y[n];
        }
        else {
            y0[n] = y[n];
            y1[n] = normal_rng(mu_t + rho * (sigma_t / sigma_c) * (y[n] - mu_c), sigma_t * sqrt(1 - rho^2));
        }
        tau_unit[n] = y1[n] - y0[n];
    }

    tau_fs = mean(tau_unit);
}