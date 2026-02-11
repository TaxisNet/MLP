#pragma once
#include <vector>

struct MlpResult {
  double cost;
  std::vector<int> solution;
  double elapsedSeconds;
};

MlpResult run_mlp(double** distanceMatrix,
                  int dimension,
                  const std::vector<double>& nodeWeights,
                  int iIls,
                  unsigned int seed,
                  int iterMax = 10);