functions {

    real K2_loglikelihood(real u, real v, int K2, int N2) {

        real likelihood_sum = 0.0;
        real increment;
        real binom_coeffs;

        vector[( (K2 + 1) * (K2 + 2) )/ 2] elems;
        real temp = 0.0;

        int count = 1;

        for (m1 in 0:K2) {
            int m2_upper = K2 - m1;
            for (m2 in 0:m2_upper) {

                temp = lchoose(N2, m1) 
                        + lchoose(N2 - m1, m2)
                        + lchoose(N2 - m1 - m2, K2 - m1 - m2)
                        + (K2 - m2) * log(u)
                        + (N2 - K2 + m2) * log(1.0 - u)
                        + (K2 - m1) * log(v)
                        + (N2 - K2 + m1) * log(1.0 - v);

                elems[count] = temp;
                count += 1;
            }
        }
        
        return log_sum_exp(elems);
    }
}

data { 
    int<lower=0> K1; // # B -> C | A -> B
    int<lower=0> N1; // # B | A -> B
    int<lower=0> K2; // # B ->C | A -/-> B
    int<lower=0> N; // total number of B cells

    // Parameters for priors
    real<lower=0, upper=1> u_prior[2]; // Parameters for B -> C due to factors other than A -> B
    real<lower=0, upper=1> v_prior[2]; // Parameters for B -> C  due to A -> B

}

parameters {

    // Model the rates as a function of the prior rate parameters (Beta distributions)
    real<lower=0, upper=1> u; 
    real<lower=0, upper=1> v;
    
}

model {
    // Specify the priors
    target += beta_lpdf(u | u_prior[1] + K1, u_prior[2] + N1 - K1); // Beta prior for u
    target += beta_lpdf(v | v_prior[1], v_prior[2]); // Beta prior for v

    // Specify the likelihood for K2 (custom function for this)
    target += K2_loglikelihood(u, v, K2, N - N1);

}

generated quantities {

    // Some useful quantities to calculate
    real prob_ab_promotes_bc = v / (u + v - u * v); // The probability that A -> B promotes B -> C
    real prob_only_ab_promotes_bc = (v - u * v)/(u + v - u * v); // The probability that only    A -> B promotes B -> C
    real prob_ab_inhibits_bc = (1.0 - v) / (1.0 - u * v); // The probability that A -> B inhibits B -> C
    real prob_only_ab_inhibits_bc = (u - u * v) / (1.0 - u * v); // The probability that only A -> B inhibits B -> C

    // We'd like to sample from the likelihood to do some posterior predictive checks
    int K1_pred = binomial_rng(N1, u);
    int K2_pred = binomial_rng(N - N1, u + v - u * v);

}
