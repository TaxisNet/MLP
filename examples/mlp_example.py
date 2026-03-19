import numpy as np

# add to PYTHONPATH
import sys, os
MPL_PATH = os.path.join(os.path.expanduser('~'), 'Documents', 'MLP', 'build')
sys.path.append(MPL_PATH)

import mlp_py


np.random.seed(100)  # For reproducibility
n = 100
dist = np.random.rand(n, n)*100  # Random distance matrix with values between 0 and 100
dist = (dist + dist.T) / 2.0
np.fill_diagonal(dist, 0.0)

weights = np.ones(n, dtype=np.float64)
weights = np.random.rand(n)  # Random weights for each node
weights /= np.sum(weights)  # Normalize weights to sum to 1
cost, solution = mlp_py.run(dist, weights, verbose=False)
print("Cost:", cost)
print("Solution:", solution)