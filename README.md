# Minimum Latency Problem (MLP) Solver

A C++ implementation of the **GILS-RVND** (Greedy Iterated Local Search with Randomized Variable Neighborhood Descent) metaheuristic to solve the Minimum Latency Problem (MLP).

## Overview

The Minimum Latency Problem (MLP), also known as the Traveling Repairman Problem, searches for a tour that visits all nodes such that the sum of arrival times (latency) at each node is minimized. Unlike the standard TSP which minimizes total distance, MLP focuses on minimizing the waiting time for all customers.

This repository implements a high-performance solver that combines:
- **Heuristic Construction**: Randomized Cheapest Insertion.
- **Local Search**: RVND with multiple neighborhood structures.
- **Metaheuristic**: Iterated Local Search (ILS) with perturbation.

## Algorithms Implemented

The solver utilizes a robust local search strategy with the following neighborhood moves:
1. **Swap**: Exchange two vertices.
2. **2-opt**: Invert a subsequence of vertices.
3. **Reinsertion**: Move a vertex to a different position.
4. **Or-opt 2**: Move a sequence of 2 vertices.
5. **Or-opt 3**: Move a sequence of 3 vertices.

## Project Structure

- `src/`: Source code (`main.cpp`, `readData.cpp`).
- `instances/`: TSPLIB benchmark instances (e.g., `berlin52.tsp`, `kroA100.tsp`).
- `benchmark/`: Scripts and output files for results.
- `makefile`: Build configuration.
- `run.sh`: Automated benchmarking script.

## Getting Started

### Prerequisites
- G++ Compiler (supporting C++11 or later)
- Make

### Compilation

To compile the project, run the `make` command in the root directory:

```bash
make
```

This will generate the executable binary `mlp`.

### Running the Solver

To run the solver on a specific instance, provide the path to the instance file as an argument:

```bash
./mlp instances/<instance_filename>
```

**Example:**
```bash
./mlp instances/berlin52.tsp
```

### Benchmarking

A script is provided to automate running benchmarks across all available instances:

```bash
sh run.sh
```

*Note: Ensure the script targets the correct executable (`./mlp`). You may need to update `run.sh` if it references an older executable name.*


## Python Bindings (pybind11)

### Prerequisites
- CMake 3.10+
- A C++ compiler (C++11)
- Python 3 with `numpy`

### Build the Python Module

```bash
cmake -S . -B build
cmake --build build -j
```

### Using the Python Module
```python
# also include the build directory in your PYTHONPATH if necessary
# PYTHONPATH=build python3 python_example.py  
import mlp_py
```