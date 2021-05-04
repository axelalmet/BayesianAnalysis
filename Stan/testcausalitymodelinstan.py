import pystan
import numpy as np
import matplotlib.pyplot as plt
import pickle
import gzip
import seaborn as sns

# model = pystan.StanModel('./Models/contingencytables.stan')
# with open('./Models/contingencytables.bin', 'wb') as f:
#     pickle.dump(model, f)

with open ('./Models/causalitymodel.bin', 'rb') as f:
    model1 = pickle.load(f)


with open ('./Models/contingencytables.bin', 'rb') as f:
    model2 = pickle.load(f)

# Define the data
data1 = {'N': 200, 'N1': 100, 'K1': 12, 'K2': 30, "u_prior": [1, 1], "v_prior": [1, 1]}
data2 = {'N': 200, 'N1': 100, 'K1': 12, 'K2': 30, "theta1_prior": [1, 1], "theta2_prior": [1, 1]}

fit1 = model1.sampling(data=data1, n_jobs=1)
fit2 = model2.sampling(data=data2, n_jobs=1)

# # Save the fit
# with open('./Output/baldishababamodel_stan.gz', 'wb') as output_f:
#     pickle.dump({'model': model, 'fit': fit}, output_f)

# Plot the fits
plt.figure(figsize=(5, 5))
sns.distplot(fit1['u'])
sns.distplot(fit1['v'])
plt.show()

plt.figure(figsize=(5, 5))
sns.distplot(fit2['theta1'])
sns.distplot(fit2['theta2'])
plt.show()
