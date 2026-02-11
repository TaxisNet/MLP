import numpy as np
import mlp_py


np.random.seed(0)  # For reproducibility
n = 10
dist = np.random.rand(n, n)
dist = (dist + dist.T) / 2.0
np.fill_diagonal(dist, 0.0)

weights = np.ones(n, dtype=np.float64)
weights[0] = 2.0
cost, solution = mlp_py.run(dist, weights, verbose=True)



# print(cost)
# print(solution)