import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymc3 as pm
import scipy
import scipy.stats as st

#### Define the parameters and data
N = 200
N1 = 100
K1 = 12
K2 = 30

## Parameters for the Beta priors
a = 1
b = 1
c = 1
d = 1


with pm.Model() as model:
    u = pm.Beta('u', alpha = a + K1, beta = b + N1 - K1)
    v = pm.Beta('v', alpha = c, beta = d)


    # Let's define the causality model here to see if we can get it working
    def logp(u, v, K2, N2):

        num_elements = (K2 + 1) * (K2 + 2) / 2 # Number of values to consider
        elems = np.zeros(num_elements)

        count = 0

        # We define each element via its log for computational reasons
        for i in range(K2 + 1):

            j_upper = K2 + 1 - i

            for j in range(j_upper):
                
                temp = 

        return pm.logsumexp(elems)

    K2_likelihood = pm.DensityDist('K2_likelihood', logp, observed={''})