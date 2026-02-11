#ifndef NEIGHBORHOODS_H
#define NEIGHBORHOODS_H

#include "structures.h"
#include <vector>
using namespace std;

neighborInfo swap(vector <int> &solution, vector <vector <subsequenceInfo>> &subsequenceMatrix , vector<double> &nodeWeights);
neighborInfo twoOpt(vector <int> &solution, vector <vector <subsequenceInfo>> &subsequenceMatrix, vector<double> &nodeWeights);
neighborInfo reinsertion(vector <int> &solution, vector <vector <subsequenceInfo>> &subsequenceMatrix, vector<double> &nodeWeights);
neighborInfo oropt2(vector <int> &solution, vector <vector <subsequenceInfo>> &subsequenceMatrix, vector<double> &nodeWeights);
neighborInfo oropt3(vector <int> &solution, vector <vector <subsequenceInfo>> &subsequenceMatrix, vector<double> &nodeWeights);
void RVND(vector <int> &solution, vector <vector <subsequenceInfo>> &subsequenceMatrix, vector<double> &nodeWeights);
void updatesMatrix(vector <vector <subsequenceInfo>> &subsequenceMatrix, vector <int> &solution, vector<double> &nodeWeights); 

#endif 