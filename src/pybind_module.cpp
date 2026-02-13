#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include <stdexcept>
#include "search.h"

namespace py = pybind11;

double ** distanceMatrix = nullptr;
int dimension = 0;

static void free_distance_matrix() {
  if (!distanceMatrix) {
    return;
  }
  for (int i = 0; i <= dimension; i++) {
    delete[] distanceMatrix[i];
  }
  delete[] distanceMatrix;
  distanceMatrix = nullptr;
}

static void set_distance_matrix(const py::array_t<double, py::array::c_style | py::array::forcecast>& matrix) {
  auto buf = matrix.request();
  if (buf.ndim != 2) {
    throw std::invalid_argument("distance_matrix must be 2D");
  }
  if (buf.shape[0] != buf.shape[1]) {
    throw std::invalid_argument("distance_matrix must be square");
  }

  int n = static_cast<int>(buf.shape[0]);
  free_distance_matrix();
  dimension = n;

  distanceMatrix = new double*[n + 1];
  for (int i = 0; i <= n; i++) {
    distanceMatrix[i] = new double[n + 1];
  }

  for (int i = 0; i <= n; i++) {
    for (int j = 0; j <= n; j++) {
      distanceMatrix[i][j] = 0.0;
    }
  }

  double* data = static_cast<double*>(buf.ptr);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      distanceMatrix[i + 1][j + 1] = data[i * n + j];
    }
  }
}

static std::vector<double> build_weights(const py::array_t<double, py::array::c_style | py::array::forcecast>& weights, int n) {
  auto buf = weights.request();
  if (buf.ndim != 1) {
    throw std::invalid_argument("weights must be 1D");
  }

  int wlen = static_cast<int>(buf.shape[0]);
  std::vector<double> nodeWeights(n + 1, 1.0);
  double* data = static_cast<double*>(buf.ptr);

  if (wlen == n) {
    for (int i = 0; i < n; i++) {
      nodeWeights[i + 1] = data[i];
    }
    return nodeWeights;
  }

  if (wlen == n + 1) {
    for (int i = 0; i <= n; i++) {
      nodeWeights[i] = data[i];
    }
    return nodeWeights;
  }

  throw std::invalid_argument("weights length must be n or n+1");
}

static py::tuple run(const py::array_t<double, py::array::c_style | py::array::forcecast>& matrix,
                     const py::array_t<double, py::array::c_style | py::array::forcecast>& weights,
                     int iils,
                     int seed,
                     bool verbose) {
  set_distance_matrix(matrix);
  int n = dimension;
  if (n <= 0) {
    throw std::invalid_argument("distance_matrix size must be > 0");
  }

  if (iils < 0) {
    iils = (n >= 100) ? 100 : n;
  }

  if (seed >= 0) {
    srand(static_cast<unsigned int>(seed));
  }

  std::vector<double> nodeWeights = build_weights(weights, n);
  std::vector<int> solution;
  double cost = search(iils, n, nodeWeights, &solution, verbose);

  // Solution is 1-based, convert to 0-based for Python
  for (int& node : solution) {
    node -= 1;
  }
  return py::make_tuple(cost, solution);
}

PYBIND11_MODULE(mlp_py, m) {
  m.def("run", &run, py::arg("distance_matrix"), py::arg("weights"),
        py::arg("iils") = -1, py::arg("seed") = -1, py::arg("verbose") = false);
}