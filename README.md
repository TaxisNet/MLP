# Minimum Latency Problem (MLP) Solver

High-performance C++ implementation of **GILS-RVND** (Greedy Iterated Local Search with Randomized Variable Neighborhood Descent) for the **Minimum Latency Problem (MLP)**, also known as the **Traveling Repairman Problem (TRP)**. GILS-RVND Algorithm  introduced by [Silva et al. (2012)](https://doi.org/10.1016/j.ejor.2012.03.044).

## What this repository provides

- A fast C++ solver for TSPLIB-style instances.
- Optional Python bindings via pybind11 (`mlp_py`).
- Benchmark inputs and scripts.

MLP minimizes the **sum of arrival times** across customers (not tour length as in TSP).

---

## Algorithm summary

Given a permutation of nodes, arrival time at node $v$ is the cumulative travel time from the start to $v$. The objective is:

$$
\min \sum_{v=1}^{n} w_v \cdot T_v
$$

where $w_v$ is node weight (default $1.0$) and $T_v$ is arrival time.

The solver combines:

1. **Construction**: randomized cheapest insertion.  
	See [src/construction.cpp](src/construction.cpp).
2. **RVND local search**: Swap, 2-opt, Reinsertion, Or-opt(2), Or-opt(3).  
	See [src/neighborhoods.cpp](src/neighborhoods.cpp).
3. **ILS metaheuristic**: double-bridge perturbation + RVND restart.  
	See [src/search.cpp](src/search.cpp) and [src/construction.cpp](src/construction.cpp).
4. **Fast move evaluation** with subsequence information for efficient local search.  
	See [src/neighborhoods.cpp](src/neighborhoods.cpp).

---
## Algorithm parameters
Not all parameters are currently exposed in the Python API, but you can modify them in [src/search.cpp](src/search.cpp) (Notetation follows Silva et al.):
- `I_max`: max ILS iterations (default `10`)
- `I_ils`: max iterations without improvement before perturbation (default `min(100, diamention)`)
- `alpha`: randomness parameter for construction (default `0.3`)



## Requirements

### Core (C++)

- CMake >= 3.10
- C++ compiler with C++11 support (GCC/Clang/MSVC)

### Python bindings

- **Python 3.8+ recommended** (pybind11-based module)
- NumPy
- A Python environment matching the interpreter you want to import `mlp_py` from

> Note: pybind11 is resolved automatically from your system or fetched by CMake.

---

## Build instructions

### 1) Build C++ executable (and Python module by default)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This generates:

- `build/mlp` (C++ executable)
- `build/mlp_py.*` (Python extension module, when `BUILD_PYBIND=ON`)

### 2) Build only C++ executable (disable Python)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYBIND=OFF
cmake --build build -j
```

---

## Run instructions (C++)

### TSPLIB input

```bash
./build/mlp instances/berlin52.tsp
```

### JSON input

If the first argument ends with `.json`, the solver reads weighted data from JSON:

```bash
./build/mlp path/to/instance.json
```

JSON format is : 
```json
{
  "distance_matrix": [[0, d_12, ...], [d_21, 0, ...], ...],
  "weights": [w_1, w_2, ..., w_n] // Optional, defaults to 1.0 if missing
}
```
See [instances/burnma14.json](instances/burnma14.json) for an example.

---

## Run instructions (Python)

### Recommended setup (virtual environment)

```bash
python3 --version
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip numpy
```

Build the project (inside the same shell/environment):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Run example:

```bash
PYTHONPATH=build python examples/mlp_example.py
```

Alternative (interactive):

```python
import sys
sys.path.append("/absolute/path/to/MLP/build")
import mlp_py
```

If import fails, confirm:

1. You are using the same Python environment used during CMake configuration.
2. `PYTHONPATH` includes `build`.
3. The module file exists in `build` (`mlp_py*.so` on Linux).

---

## Python API

`mlp_py.run(distance_matrix, weights, iils=-1, seed=-1, verbose=False)`

### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `distance_matrix` | `numpy.ndarray` (2D, float64) | required | Square matrix of shape `(n, n)` |
| `weights` | `numpy.ndarray` (1D, float64) | required | Node weights of length `n` (or `n+1`) |
| `iils` | `int` | `-1` | If `-1`, auto-sets to `min(100, n)` |
| `seed` | `int` | `-1` | If `>= 0`, seeds RNG for reproducibility |
| `verbose` | `bool` | `False` | Enables verbose search logs |

### Returns

- `cost: float`
- `solution: list[int]` (0-based node indices)

---

## Indexing convention

- Internal solver indexing is **1-based**.
- Python API returns **0-based** indices.

---



## Benchmarking

Use:

```bash
sh run.sh
```

Benchmark assets are in [benchmark/](benchmark).

---

## Project layout

- [src/](src): core C++ implementation
- [instances/](instances): TSPLIB-style instances
- [examples/](examples): Python examples
- [benchmark/](benchmark): benchmark scripts/results
- [run.sh](run.sh): benchmark runner
- [CMakeLists.txt](CMakeLists.txt): build configuration

---

## Reproducibility

- C++ executable seeds with current time by default.
- Python binding supports deterministic runs with `seed >= 0`.

---

## License

See [LICENSE](LICENSE).