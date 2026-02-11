#ifndef SEARCH_H
#define SEARCH_H

#include <vector>
using namespace std;

double search(int iIls, int dimension, vector<double> &nodeWeights,
              vector<int>* outSolution = nullptr, bool verbose = true);

#endif