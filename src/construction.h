#ifndef CONSTRUCTION_H
#define CONSTRUCTION_H


#include "structures.h"
#include <vector>
using namespace std;

vector <int> construction(vector <int> candidatesList, const vector<double>& nodeWeights,  double alpha);
vector <int> constructionSmith(vector <int> candidatesList, const vector<double>& nodeWeights, double alpha);
vector <int> pertub(vector <int> &solution);


#endif