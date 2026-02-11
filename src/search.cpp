#include "search.h"
#include "structures.h"
#include "construction.h"
#include "neighborhoods.h"

#include <vector>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <ctime>

double search(int iIls, int dimension, vector<double> &nodeWeights,
              vector<int>* outSolution, bool verbose) {
  double bestCurrentCost = DBL_MAX, currentCost, finalCost = DBL_MAX;
  vector <int> vertices, bestCurrentSolution, currentSolution, finalSolution;
  vector <vector <subsequenceInfo>> subsequenceMatrix(dimension+1, vector <subsequenceInfo> (dimension+1));
  
  // time variables
  clock_t start = 0; clock_t end=0;

  if (verbose) {
    clock_t start = clock(); // Starts time counting
  }

  // Creates a list of vertices from 1 to dimension (so is also 1-indexed)
  for(int i = 0; i < dimension; i++){
    vertices.push_back(i+1);
  }

  for(int iterMax = 0; iterMax < 10; iterMax++){
    double alpha = (rand() % 90) / 100.0 + 0.1;

    currentSolution = construction(vertices, alpha); // Generates initial solution
    updatesMatrix(subsequenceMatrix, currentSolution, nodeWeights); // Updates subsequence matrix with new solution
    
    bestCurrentSolution = currentSolution;
    bestCurrentCost = subsequenceMatrix[0][currentSolution.size()-1].acumulateCost;

    int iterIls = 0;
    while(iterIls < iIls){
      RVND(currentSolution, subsequenceMatrix, nodeWeights);

      currentCost = subsequenceMatrix[0][currentSolution.size()-1].acumulateCost;

      if(currentCost < bestCurrentCost){
        bestCurrentSolution = currentSolution;
        bestCurrentCost = currentCost;
        iterIls = 0;
      }
      else{
        iterIls++;
      }

      currentSolution = pertub(bestCurrentSolution);
      updatesMatrix(subsequenceMatrix, currentSolution, nodeWeights);
    }

    if(bestCurrentCost < finalCost){
      finalSolution = bestCurrentSolution;
      finalCost = bestCurrentCost;
    }
  }

  if (outSolution) {
    *outSolution = finalSolution;
  }

  if (verbose) {
    clock_t end = clock(); // Ends time counting
    double time = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout << "Search time: " << setprecision( time<2.0? 3 : 2) <<time << " seconds." << endl;
    cout << "Final Cost: " << setprecision(finalCost<10e2 ? 0: 3) << finalCost << endl;
    setprecision(cout.precision());
    cout << "Solution: [";
    for(int i = 0; i < finalSolution.size(); i++){
      cout << finalSolution[i] << " ";
    }
    cout <<"]" << endl;

  }

  return finalCost;  
}