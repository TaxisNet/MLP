#include "search.h"
#include "structures.h"
#include "construction.h"
#include "neighborhoods.h"

#include <vector>
#include <cfloat>
#include <cstdlib>
#include <iostream>


double search(int iIls, int dimension){
  double bestCurrentCost = DBL_MAX, currentCost, finalCost = DBL_MAX;
  vector <int> vertices, bestCurrentSolution, currentSolution, finalSolution;
  vector <vector <subsequenceInfo>> subsequenceMatrix(dimension+1, vector <subsequenceInfo> (dimension+1));
  
  // nodeWeights[i] represents the weight of node i (so is also 1-indexed)
  vector <double> nodeWeights(dimension + 1, 1.0); // Initialize all node weights to 1.0
  // Creates vector with vertices
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

      // cout<<iterIls<<":"<<iIls<<": "<<bestCurrentCost<<" "<< currentCost<<endl;
      currentSolution = pertub(bestCurrentSolution);
      updatesMatrix(subsequenceMatrix, currentSolution, nodeWeights);

      // iterIls++;<<< is this a bug? added the else above
    }

    if(bestCurrentCost < finalCost){
      finalSolution = bestCurrentSolution;
      finalCost = bestCurrentCost;
    }
  }

  cout << endl << "Solution: ";
  for(int i = 0; i < finalSolution.size(); i++){
    cout << finalSolution[i] << " ";
  }

  return finalCost;  
}