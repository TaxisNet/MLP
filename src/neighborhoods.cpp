#include "neighborhoods.h"
#include "structures.h"
#include <vector>
#include <cfloat>
#include <cstdlib>

using namespace std;

extern double ** distanceMatrix; 
extern int dimension; 

neighborInfo swap(vector <int> &solution, vector <vector <subsequenceInfo>> &subsequenceMatrix , vector<double> &nodeWeights){
  double partialCost, cost, partialTime;
  int size = solution.size();

  neighborInfo bestNeighbor;
  bestNeighbor.bestCost = DBL_MAX;

  for(int i = 1; i < size-1; i++){
    for(int j = i+1; j < size; j++){
      if(j == i + 1){
      // Adjacent swap: positions i and i+1
        
        // Segment [0, i-1]: unchanged
        partialCost = subsequenceMatrix[0][i-1].acumulateCost;
        partialTime = subsequenceMatrix[0][i-1].totalTime;
        
        // Node j arrives at partialTime + d[i-1→j]
        double newEdge1 = distanceMatrix[solution[i-1]][solution[j]];
        partialTime += newEdge1;
        partialCost += nodeWeights[solution[j]] * partialTime;
        
        // Node i arrives after j
        double newEdge2 = distanceMatrix[solution[j]][solution[i]];
        partialTime += newEdge2;
        partialCost += nodeWeights[solution[i]] * partialTime;
        
        // Remaining segment [j+1, n]: shift by partialTime
        int n = solution.size() - 1;  // Last position in solution
        if (j< n){
          double edgeToNext = distanceMatrix[solution[i]][solution[j+1]];
          cost = partialCost + 
                 subsequenceMatrix[j+1][n].sumWeights * (partialTime + edgeToNext) +
                 subsequenceMatrix[j+1][n].acumulateCost;
        } else {
          cost = partialCost;
        }

      }else{
        // Non-adjacent swap
        
        // Segment [0, i-1]: unchanged
        partialCost = subsequenceMatrix[0][i-1].acumulateCost;
        partialTime = subsequenceMatrix[0][i-1].totalTime;
        
        // Node j arrives first
        double newEdge1 = distanceMatrix[solution[i-1]][solution[j]];
        partialTime += newEdge1;
        partialCost += nodeWeights[solution[j]] * partialTime;
        
        // Middle segment [i+1, j-1]: shifted by partialTime
        double edgeToMiddle = distanceMatrix[solution[j]][solution[i+1]];
        double middleShift = partialTime + edgeToMiddle;
        partialCost += subsequenceMatrix[i+1][j-1].sumWeights * middleShift + 
                       subsequenceMatrix[i+1][j-1].acumulateCost;
        partialTime = middleShift + subsequenceMatrix[i+1][j-1].totalTime;
        
        // Node i arrives after middle segment
        double edgeToI = distanceMatrix[solution[j-1]][solution[i]];
        partialTime += edgeToI;
        partialCost += nodeWeights[solution[i]] * partialTime;
        
        // Remaining segment [j+1, n]: shift by partialTime
        if (j+1 < size){
          double edgeToRest = distanceMatrix[solution[i]][solution[j+1]];
          double restShift = partialTime + edgeToRest;
          cost = partialCost + 
                 subsequenceMatrix[j+1][size-1].sumWeights * restShift +
                 subsequenceMatrix[j+1][size-1].acumulateCost;
        } else {
          cost = partialCost;
        }
      }

      

      if(cost < bestNeighbor.bestCost){    
        bestNeighbor.iBest = i;
        bestNeighbor.jBest = j;
        bestNeighbor.bestCost = cost;
	    } 
    }
  }

  return bestNeighbor;
}

neighborInfo twoOpt(vector <int> &solution, vector <vector <subsequenceInfo>> &subsequenceMatrix, vector<double> &nodeWeights){
  double partialCost, cost, partialTime;
  int size = solution.size();

  neighborInfo bestNeighbor;
  bestNeighbor.bestCost = DBL_MAX;
  
  for(int i = 1; i < size-1; i++){
    for(int j = i+1; j < size; j++){
      // 2-opt reverses segment [i, j]
      
      // Segment [0, i-1]: unchanged
      partialCost = subsequenceMatrix[0][i-1].acumulateCost;
      partialTime = subsequenceMatrix[0][i-1].totalTime;
      
      // Edge from i-1 to j (entering reversed segment)
      double edgeToReversed = distanceMatrix[solution[i-1]][solution[j]];
      double reversedShift = partialTime + edgeToReversed;
      
      // Reversed segment [j, i]: shifted by reversedShift
      // Use [j][i] which contains the reversed subsequence info
      partialCost += subsequenceMatrix[j][i].sumWeights * reversedShift + 
                     subsequenceMatrix[j][i].acumulateCost;
      partialTime = reversedShift + subsequenceMatrix[j][i].totalTime;
      
      // Remaining segment [j+1, n]: shifted by partialTime
      // Tail condition check
      if (j+1 < size){
      double edgeToRest = distanceMatrix[solution[i]][solution[j+1]];
      double restShift = partialTime + edgeToRest;
      cost = partialCost + 
             subsequenceMatrix[j+1][size-1].sumWeights * restShift +
             subsequenceMatrix[j+1][size-1].acumulateCost;
      }
       else {
        cost = partialCost;
      }


      if(cost < bestNeighbor.bestCost){    
        bestNeighbor.iBest = i;
        bestNeighbor.jBest = j;
        bestNeighbor.bestCost = cost;
      } 
    }
  }

  return bestNeighbor;
}

neighborInfo reinsertion(vector <int> &solution, vector <vector <subsequenceInfo>> &subsequenceMatrix, vector<double> &nodeWeights){
  double partialCost, cost, partialTime;
  int size = solution.size();

  neighborInfo bestNeighbor;
  bestNeighbor.bestCost = DBL_MAX;

  // Case 1 - Forward reinsertion: Move node i forward to position j (i < j)
  // Resulting path: [0, i-1] → [i+1, j] → i → [j+1, n]
  for(int i = 1; i < size-1; i++){
    for(int j = i + 1; j < size; j++){
      
      // Segment [0, i-1]: unchanged
      partialCost = subsequenceMatrix[0][i-1].acumulateCost;
      partialTime = subsequenceMatrix[0][i-1].totalTime;
      
      // Bridge: skip node i, go directly from i-1 to i+1
      double edgeSkip = distanceMatrix[solution[i-1]][solution[i+1]];
      double middleShift = partialTime + edgeSkip;
      
      // Segment [i+1, j]: shifted by middleShift
      // Formula: W[i+1,j] × T_shift + C[i+1,j]
      partialCost += subsequenceMatrix[i+1][j].sumWeights * middleShift + 
                     subsequenceMatrix[i+1][j].acumulateCost;
      partialTime = middleShift + subsequenceMatrix[i+1][j].totalTime;
      
      // Insert node i after position j
      double edgeToI = distanceMatrix[solution[j]][solution[i]];
      partialTime += edgeToI;
      partialCost += nodeWeights[solution[i]] * partialTime;
      
      // Remaining segment [j+1, n]: shifted by partialTime
      // Tail condition check
      if(j+1 < size){
        double edgeToRest = distanceMatrix[solution[i]][solution[j+1]];
        double restShift = partialTime + edgeToRest;
        cost = partialCost + 
              subsequenceMatrix[j+1][size-1].sumWeights * restShift +
              subsequenceMatrix[j+1][size-1].acumulateCost;
      } else {
        cost = partialCost;
      }


      if(cost < bestNeighbor.bestCost){    
        bestNeighbor.iBest = i;
        bestNeighbor.jBest = j;
        bestNeighbor.bestCost = cost;
      } 
    }
  }

  // Case 2 - Backwards reinsertion: Move node i backward to position j (j < i)
  // Resulting path: [0, j-1] → i → [j, i-1] → [i+1, n]
  for(int j = 1; j < size-1; j++){
    for(int i = j + 1; i < size; i++){
      
      // Segment [0, j-1]: unchanged
      partialCost = subsequenceMatrix[0][j-1].acumulateCost;
      partialTime = subsequenceMatrix[0][j-1].totalTime;
      
      // Insert node i at position j
      double edgeToI = distanceMatrix[solution[j-1]][solution[i]];
      partialTime += edgeToI;
      partialCost += nodeWeights[solution[i]] * partialTime;
      
      // Edge from i to j (where middle segment starts)
      double edgeToMiddle = distanceMatrix[solution[i]][solution[j]];
      double middleShift = partialTime + edgeToMiddle;
      
      // Segment [j, i-1]: shifted by middleShift
      // Formula: W[j,i-1] × T_shift + C[j,i-1]
      partialCost += subsequenceMatrix[j][i-1].sumWeights * middleShift + 
                     subsequenceMatrix[j][i-1].acumulateCost;
      partialTime = middleShift + subsequenceMatrix[j][i-1].totalTime;
      
      // Remaining segment [i+1, n]: shifted by partialTime
      // Bridge: skip node i, connect i-1 to i+1

      // Tail condition check
      if (i+1 < size){
      double edgeToRest = distanceMatrix[solution[i-1]][solution[i+1]];
      double restShift = partialTime + edgeToRest;
      cost = partialCost + 
             subsequenceMatrix[i+1][size-1].sumWeights * restShift +
             subsequenceMatrix[i+1][size-1].acumulateCost;
      } else {
        cost = partialCost;
      }


      if(cost < bestNeighbor.bestCost){    
        bestNeighbor.iBest = i;
        bestNeighbor.jBest = j;
        bestNeighbor.bestCost = cost;
      } 
    }
  }

  return bestNeighbor;
}

neighborInfo oropt2(vector <int> &solution, vector <vector <subsequenceInfo>> &subsequenceMatrix, vector<double> &nodeWeights){
  double partialCost, cost, partialTime;
  int size = solution.size();

  neighborInfo bestNeighbor;
  bestNeighbor.bestCost = DBL_MAX;

  for(int i = 1; i < size-2; i++){
    for(int j = 1; j <= size-3; j++){
      if(i != j){
        if(i < j){
          // Forward move: [0, i-1] → [i+2, j+1] → [i, i+1] → [j+2, n]
          
          // Segment [0, i-1]: unchanged
          partialCost = subsequenceMatrix[0][i-1].acumulateCost;
          partialTime = subsequenceMatrix[0][i-1].totalTime;
          
          // Bridge: skip nodes [i, i+1], connect i-1 to i+2
          double edgeSkip = distanceMatrix[solution[i-1]][solution[i+2]];
          double middleShift = partialTime + edgeSkip;
          
          // Segment [i+2, j+1]: shifted by middleShift
          partialCost += subsequenceMatrix[i+2][j+1].sumWeights * middleShift + 
                         subsequenceMatrix[i+2][j+1].acumulateCost;
          partialTime = middleShift + subsequenceMatrix[i+2][j+1].totalTime;
          
          // Insert first node (i) after segment [i+2, j+1]
          double edgeToI = distanceMatrix[solution[j+1]][solution[i]];
          partialTime += edgeToI;
          partialCost += nodeWeights[solution[i]] * partialTime;
          
          // Insert second node (i+1) after node i
          double edgeToI1 = distanceMatrix[solution[i]][solution[i+1]];
          partialTime += edgeToI1;
          partialCost += nodeWeights[solution[i+1]] * partialTime;
          
          // Remaining segment [j+2, n]: shifted by partialTime
          // Tail condition check
          if (j + 2 < size){
          double edgeToRest = distanceMatrix[solution[i+1]][solution[j+2]];
          double restShift = partialTime + edgeToRest;
          cost = partialCost + 
                 subsequenceMatrix[j+2][size-1].sumWeights * restShift +
                 subsequenceMatrix[j+2][size-1].acumulateCost;
          } else {
            cost = partialCost;
          }

                 
        }else{
          // Backward move: [0, j-1] → [i, i+1] → [j, i-1] → [i+2, n]
          
          // Segment [0, j-1]: unchanged
          partialCost = subsequenceMatrix[0][j-1].acumulateCost;
          partialTime = subsequenceMatrix[0][j-1].totalTime;
          
          // Insert first node (i) at position j
          double edgeToI = distanceMatrix[solution[j-1]][solution[i]];
          partialTime += edgeToI;
          partialCost += nodeWeights[solution[i]] * partialTime;
          
          // Insert second node (i+1) after node i
          double edgeToI1 = distanceMatrix[solution[i]][solution[i+1]];
          partialTime += edgeToI1;
          partialCost += nodeWeights[solution[i+1]] * partialTime;
          
          // Edge to middle segment [j, i-1]
          double edgeToMiddle = distanceMatrix[solution[i+1]][solution[j]];
          double middleShift = partialTime + edgeToMiddle;
          
          // Segment [j, i-1]: shifted by middleShift
          partialCost += subsequenceMatrix[j][i-1].sumWeights * middleShift + 
                         subsequenceMatrix[j][i-1].acumulateCost;
          partialTime = middleShift + subsequenceMatrix[j][i-1].totalTime;
          
          // Remaining segment [i+2, n]: shifted by partialTime
          // Bridge: skip nodes [i, i+1], connect i-1 to i+2
          // Tail condition check
          if (i + 2 < size){
          double edgeToRest = distanceMatrix[solution[i-1]][solution[i+2]];
          double restShift = partialTime + edgeToRest;
          cost = partialCost + 
                 subsequenceMatrix[i+2][size-1].sumWeights * restShift +
                 subsequenceMatrix[i+2][size-1].acumulateCost;
          } else {
            cost = partialCost;
          }
        }

        if(cost < bestNeighbor.bestCost){    
          bestNeighbor.iBest = i;
          bestNeighbor.jBest = j;
          bestNeighbor.bestCost = cost;
        } 
      }
    }
  }

  return bestNeighbor;
}
neighborInfo oropt3(vector <int> &solution, vector <vector <subsequenceInfo>> &subsequenceMatrix, vector<double> &nodeWeights){
  double partialCost, cost, partialTime;
  int size = solution.size();
  neighborInfo bestNeighbor;
  bestNeighbor.bestCost = DBL_MAX;

  for(int i = 1; i < size-3; i++){
    for(int j = 1; j <= size-4; j++){
      if(i != j){
        if(i < j){
          // Forward move: [0, i-1] → [i+3, j+2] → [i, i+1, i+2] → [j+3, n]
          
          // Segment [0, i-1]: unchanged
          partialCost = subsequenceMatrix[0][i-1].acumulateCost;
          partialTime = subsequenceMatrix[0][i-1].totalTime;
          
          // Bridge: skip nodes [i, i+1, i+2], connect i-1 to i+3
          double edgeSkip = distanceMatrix[solution[i-1]][solution[i+3]];
          double middleShift = partialTime + edgeSkip;
          
          // Segment [i+3, j+2]: shifted by middleShift
          partialCost += subsequenceMatrix[i+3][j+2].sumWeights * middleShift + 
                         subsequenceMatrix[i+3][j+2].acumulateCost;
          partialTime = middleShift + subsequenceMatrix[i+3][j+2].totalTime;
          
          // Insert the 3-node subsequence [i, i+1, i+2] after segment [i+3, j+2]
          double edgeToI = distanceMatrix[solution[j+2]][solution[i]];
          double threeNodeShift = partialTime + edgeToI;
          
          // The 3 nodes as a subsequence, shifted by threeNodeShift
          partialCost += subsequenceMatrix[i][i+2].sumWeights * threeNodeShift + 
                         subsequenceMatrix[i][i+2].acumulateCost;
          partialTime = threeNodeShift + subsequenceMatrix[i][i+2].totalTime;
          
          // Remaining segment [j+3, n]: shifted by partialTime
          // Tail condition check
          if (j + 3 < size){
          double edgeToRest = distanceMatrix[solution[i+2]][solution[j+3]];
          double restShift = partialTime + edgeToRest;
          cost = partialCost + 
                 subsequenceMatrix[j+3][size-1].sumWeights * restShift +
                 subsequenceMatrix[j+3][size-1].acumulateCost;
          } else {
            cost = partialCost;
          }
                 
        }else{
          // Backward move: [0, j-1] → [i, i+1, i+2] → [j, i-1] → [i+3, n]
          
          // Segment [0, j-1]: unchanged
          partialCost = subsequenceMatrix[0][j-1].acumulateCost;
          partialTime = subsequenceMatrix[0][j-1].totalTime;
          
          // Insert the 3-node subsequence [i, i+1, i+2] at position j
          double edgeToI = distanceMatrix[solution[j-1]][solution[i]];
          double threeNodeShift = partialTime + edgeToI;
          
          // The 3 nodes as a subsequence, shifted by threeNodeShift
          partialCost += subsequenceMatrix[i][i+2].sumWeights * threeNodeShift + 
                         subsequenceMatrix[i][i+2].acumulateCost;
          partialTime = threeNodeShift + subsequenceMatrix[i][i+2].totalTime;
          
          // Edge to middle segment [j, i-1]
          double edgeToMiddle = distanceMatrix[solution[i+2]][solution[j]];
          double middleShift = partialTime + edgeToMiddle;
          
          // Segment [j, i-1]: shifted by middleShift
          partialCost += subsequenceMatrix[j][i-1].sumWeights * middleShift + 
                         subsequenceMatrix[j][i-1].acumulateCost;
          partialTime = middleShift + subsequenceMatrix[j][i-1].totalTime;
          
          // Remaining segment [i+3, n]: shifted by partialTime
          // Bridge: skip nodes [i, i+1, i+2], connect i-1 to i+3
          // Tail condition check
          if (i + 3 < size){
          double edgeToRest = distanceMatrix[solution[i-1]][solution[i+3]];
          double restShift = partialTime + edgeToRest;
          cost = partialCost + 
                 subsequenceMatrix[i+3][size-1].sumWeights * restShift +
                 subsequenceMatrix[i+3][size-1].acumulateCost;
          }else {
            cost = partialCost;
          }
        }

        if(cost < bestNeighbor.bestCost){    
          bestNeighbor.iBest = i;
          bestNeighbor.jBest = j;
          bestNeighbor.bestCost = cost;
        }
      }
    }
  }

  return bestNeighbor;
}

void updatesMatrix(vector <vector <subsequenceInfo>> &subsequenceMatrix, vector <int> &solution, vector<double> &nodeWeights) { 
  int size = solution.size();

  // Updates total time
  for(int i = 0; i < size; i++){
    for(int j = i; j < size; j++){
      if(i == j){
        subsequenceMatrix[i][j].totalTime = 0;
        subsequenceMatrix[i][j].sumWeights = nodeWeights[solution[i]];
        subsequenceMatrix[i][j].vertices = 1;
      }else {
        subsequenceMatrix[i][j].totalTime = subsequenceMatrix[i][j-1].totalTime + distanceMatrix[solution[j-1]][solution[j]];
        subsequenceMatrix[j][i].totalTime = subsequenceMatrix[i][j].totalTime; 
        
        // Sum of weights from position i to j
        subsequenceMatrix[i][j].sumWeights = subsequenceMatrix[i][j-1].sumWeights + nodeWeights[solution[j]];
        subsequenceMatrix[j][i].sumWeights = subsequenceMatrix[i][j].sumWeights;
        
        // Number of vertices from position i to j
        subsequenceMatrix[i][j].vertices = j - i + 1;
        subsequenceMatrix[j][i].vertices = subsequenceMatrix[i][j].vertices;
      }
    }
  }

  // Updates accumulated cost - Forward direction
  for(int i = 0; i < size; i++){
    for(int j = i; j < size; j++){ 
      if(i == j){
        subsequenceMatrix[i][j].acumulateCost = 0;
      }else {
        // Weight of node at position j times its arrival time from position i
        //  nodeWeights[solution[j]] is the weight of the actual node
        subsequenceMatrix[i][j].acumulateCost = subsequenceMatrix[i][j-1].acumulateCost + nodeWeights[solution[j]] * subsequenceMatrix[i][j].totalTime;
      }      
    }
  }

  // Reverse direction
  for(int i = size-1; i >= 0; i--){
    for(int j = i; j >= 0; j--){
      if(i == j){
        subsequenceMatrix[i][j].acumulateCost = 0;
      }else{
        subsequenceMatrix[i][j].acumulateCost = subsequenceMatrix[i][j+1].acumulateCost + nodeWeights[solution[j]] * subsequenceMatrix[i][j].totalTime;
      }      
    }
  }
}

void RVND(vector <int> &solution, vector <vector <subsequenceInfo>> &subsequenceMatrix, vector<double> &nodeWeights){
  vector <int> neighborhoods = {0, 1, 2, 3, 4};
  neighborInfo neighbor;
  int size = solution.size();

  while(!neighborhoods.empty()){
    int choosen = rand() % neighborhoods.size();

    if(neighborhoods[choosen] == 0){
      neighbor = swap(solution, subsequenceMatrix, nodeWeights);

      if(neighbor.bestCost < subsequenceMatrix[0][size-1].acumulateCost){

        int aux = solution[neighbor.iBest];
	      solution[neighbor.iBest] = solution[neighbor.jBest];
	      solution[neighbor.jBest] = aux;

        // Updates subsequence matrix with new solution
        updatesMatrix(subsequenceMatrix, solution, nodeWeights);
        neighborhoods = {0, 1, 2, 3, 4};

      } else {
        neighborhoods.erase(neighborhoods.begin() + choosen);

      }

    }else if(neighborhoods[choosen] == 1){
      neighbor = reinsertion(solution, subsequenceMatrix, nodeWeights);

      if(neighbor.bestCost < subsequenceMatrix[0][size-1].acumulateCost){
        vector <int> solutionInicial = solution;

        solution.erase(solution.begin()+neighbor.iBest);
        solution.insert(solution.begin()+neighbor.jBest, solutionInicial[neighbor.iBest]);

        // Updates subsequence matrix with new solution
        updatesMatrix(subsequenceMatrix, solution, nodeWeights);
        neighborhoods = {0, 1, 2, 3, 4};

      } else {
        neighborhoods.erase(neighborhoods.begin() + choosen);

      }

    }else if(neighborhoods[choosen] == 2){
      neighbor = twoOpt(solution, subsequenceMatrix, nodeWeights);

      if(neighbor.bestCost < subsequenceMatrix[0][size-1].acumulateCost){

        int aux, k = neighbor.jBest - neighbor.iBest;

        if(k % 2 != 0){
          k = k + 1;
        }

        for(int q = 0; q < k/2; q++){
          aux = solution[neighbor.iBest+q];
          solution[neighbor.iBest+q] = solution[neighbor.jBest-q];
          solution[neighbor.jBest-q] = aux;
        }

        // Updates subsequence matrix with new solution
        updatesMatrix(subsequenceMatrix, solution, nodeWeights);
        neighborhoods = {0, 1, 2, 3, 4};

      } else {
        neighborhoods.erase(neighborhoods.begin() + choosen);

      }

    }else if(neighborhoods[choosen] == 3){

      neighbor = oropt2(solution, subsequenceMatrix, nodeWeights);

      if(neighbor.bestCost < subsequenceMatrix[0][size-1].acumulateCost){
        vector<int> block = {solution[neighbor.iBest], solution[neighbor.iBest + 1]};
  
        // If we remove elements BEFORE the target, the target index shifts left by 2 (for OR-Opt 2).
        int insertPos = neighbor.jBest;
        if (neighbor.iBest < neighbor.jBest) {
          insertPos -= 2; 
        }
        
        solution.erase(solution.begin() + neighbor.iBest, solution.begin() + neighbor.iBest + 2);
        
        solution.insert(solution.begin() + insertPos, block.begin(), block.end());

        // Updates subsequence matrix with new solution
        updatesMatrix(subsequenceMatrix, solution, nodeWeights);
        neighborhoods = {0, 1, 2, 3, 4};

      } else {
        neighborhoods.erase(neighborhoods.begin() + choosen);
      }   

    }else if(neighborhoods[choosen] == 4){
      neighbor = oropt3(solution, subsequenceMatrix, nodeWeights);

      if(neighbor.bestCost < subsequenceMatrix[0][size-1].acumulateCost){
        vector<int> block = {solution[neighbor.iBest], solution[neighbor.iBest + 1], solution[neighbor.iBest + 2]};

         // If we remove elements BEFORE the target, the target index shifts left by 3 (for OR-Opt 3).
        int insertPos = neighbor.jBest;
        if (neighbor.iBest < neighbor.jBest) {
          insertPos -= 3; 
        }

        solution.erase(solution.begin() + neighbor.iBest, solution.begin() + neighbor.iBest + 3);

        solution.insert(solution.begin() + insertPos, block.begin(), block.end());

        // Updates subsequence matrix with new solution
        updatesMatrix(subsequenceMatrix, solution, nodeWeights);
        neighborhoods = {0, 1, 2, 3, 4};

      } else {
        neighborhoods.erase(neighborhoods.begin() + choosen);

      }      
    }
  }
}