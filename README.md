# Minimum Latency Problem (MLP) Solver

High‑performance C++ implementation of **GILS‑RVND** (Greedy Iterated Local Search with Randomized Variable Neighborhood Descent) for the Minimum Latency Problem (MLP), also known as the Traveling Repairman Problem (TRP).

## Overview

MLP seeks a tour that visits all nodes while minimizing the **sum of arrival times** (latencies) at each node. Unlike TSP (which minimizes total travel distance), MLP focuses on minimizing waiting time for all customers.

This solver combines:
- **Heuristic construction** (randomized cheapest insertion).
- **Local search** via **RVND** with multiple neighborhood moves.
- **Metaheuristic** **ILS** with perturbation (double‑bridge) for diversification.

## Solver Workings (Detailed)

### Objective
Given a permutation of nodes, the arrival time at node $v$ is the cumulative travel time from the start to $v$. The MLP objective is:

$$
\min \sum_{v=1}^{n} w_v \cdot T_v
$$

Where $w_v$ is the node weight (default 1.0) and $T_v$ is the arrival time at node $v$.

### 1) Construction Phase
The initial solution is built using **randomized cheapest insertion**:
- Start at the depot (the first vertex in the candidate list).
- Seed the solution with a few random vertices.
- Iteratively insert a vertex at the position that yields the lowest incremental cost, with randomization controlled by $\alpha$ (restricted candidate list size).

Implementation: [src/construction.cpp](src/construction.cpp)

### 2) RVND (Randomized Variable Neighborhood Descent)
RVND explores multiple neighborhoods in random order, applying the best improving move in each neighborhood. If an improvement is found, the neighborhood list is reset; otherwise the current neighborhood is removed.

Neighborhoods:
1. **Swap**: exchange two vertices.
2. **2‑opt**: reverse a segment.
3. **Reinsertion**: move a single vertex.
4. **Or‑opt 2**: move a length‑2 block.
5. **Or‑opt 3**: move a length‑3 block.

Implementation: [src/neighborhoods.cpp](src/neighborhoods.cpp)

### 3) ILS (Iterated Local Search)
After RVND converges to a local optimum, the solution is **perturbed** (double‑bridge) to escape local minima, then RVND is applied again. The best solution across ILS iterations is retained.

Implementation: [src/search.cpp](src/search.cpp), perturbation in [src/construction.cpp](src/construction.cpp)

### 4) Efficient Cost Updates
The solver maintains a **subsequence matrix** with:
- total travel time
- cumulative weighted cost
- sum of node weights
This enables **$O(1)$** evaluation of candidate moves during local search.

Implementation: [src/neighborhoods.cpp](src/neighborhoods.cpp)

## Indexing & Output Convention

Internally, the solver uses **1‑based node indices** for TSPLIB compatibility. If you need **0‑based output**, convert the final solution by subtracting 1 from each node index before returning/printing.

## Project Structure

- [src/](src): C++ sources.
- [instances/](instances): TSPLIB instances.
- [examples/](examples): Python usage examples.
- [benchmark/](benchmark): Benchmark scripts/results.
- [run.sh](run.sh): Benchmark runner.
- [CMakeLists.txt](CMakeLists.txt): Build configuration.

## Build & Run (C++)

### Prerequisites
- CMake 3.10+
- C++ compiler (C++11+)

### Build
```bash
cmake -S . -B build
cmake --build build -j
```

### Run
```bash
./build/mlp instances/berlin52.tsp
```

## Python Bindings (pybind11)

### Build
```bash
cmake -S . -B build
cmake --build build -j
```

### Use
```bash
PYTHONPATH=build python3 examples/mlp_example.py
```

Python API:
```python
import sys
sys.path.append("/path/to/build")
import mlp_py

cost, solution = mlp_py.run(distance_matrix, weights, iils=-1, seed=-1, verbose=False)
```

Notes:
- `distance_matrix` must be a **square** NumPy array.
- `weights` can be length $n$ (1‑based internally) or length $n+1$.
- Set `seed` for reproducible runs.

## Benchmarking
Use the provided script:

```bash
sh run.sh
```

## Reproducibility
For deterministic runs, fix the random seed. The C++ binary uses `srand(time(NULL))` by default. The Python binding accepts a `seed` parameter.

## License
See [LICENSE](LICENSE).