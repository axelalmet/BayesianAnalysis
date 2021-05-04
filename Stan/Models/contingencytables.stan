data { 
    int<lower=0> K1; // # B -> C | A -> B
    int<lower=0> N1; // # B | A -> B
    int<lower=0> K2; // # B ->C | A -/-> B
    int<lower=0> N; // total number of B cells

    // Parameters for priors
    real<lower=0, upper=1> theta1_prior[2]; // Parameters for B -> C | A -> B
    real<lower=0, upper=1> theta2_prior[2]; // Parameters for B -> C | A -/-> B

}

parameters {

    // Model the rates as a function of the prior rate parameters (Beta distributions)
    real<lower=0, upper=1> theta1; 
    real<lower=0, upper=1> theta2;
    
}

model {
    // Specify the rates
    theta1 ~ beta(theta1_prior[1], theta1_prior[2]);
    theta2 ~ beta(theta2_prior[1], theta2_prior[2]);

    // Specify the model for K and L
    K1 ~ binomial(N1, theta1);
    K2 ~ binomial(N - N1, theta2);

}

generated quantities {
    real v = (theta2 - theta1) / (1 - theta1);
    real prob_v_affects_ab = v / ( theta1 + v - theta1 * v );
    real prob_only_v_affects_ab = (v - theta1 * v) / ( theta1 + v - theta1 * v );
    real odds_ratio = theta2 * (1 - theta1) / ((theta1 * (1 - theta2)));
}
