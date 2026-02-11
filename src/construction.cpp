#include "construction.h"
#include "structures.h"
#include <cstdlib>
#include <algorithm>
#include <cfloat>

extern double ** distanceMatrix; 
extern int dimension; 



bool compares(insertionInfo a, insertionInfo b){
  return a.cost < b.cost;
}



vector <int> construction(vector <int> candidatesList, double alpha){
  vector <int> initialSolution;

  // Insert depot
  initialSolution.push_back(candidatesList[0]);

  // Deletes it from candidates list
  candidatesList.erase(candidatesList.begin());

  // Chooses 3 random vertices
  for(int i = 0; i < 3; i++){
    int j = rand() % candidatesList.size();
    initialSolution.insert(initialSolution.begin()+1, candidatesList[j]); 
    candidatesList.erase(candidatesList.begin()+j);
  }

  // Calculates insertion cost of vertices on solution
  while(!candidatesList.empty()){
    vector <insertionInfo> insertionCost((initialSolution.size()) * candidatesList.size());

    for (int i = 0, j = 1, k = 0; i < initialSolution.size(); i++, j++) {
      for (auto l : candidatesList) {
        if (j == initialSolution.size()){
          insertionCost[k].cost = distanceMatrix[initialSolution[i]][l];
        }
        else{
          insertionCost[k].cost = 
          distanceMatrix[initialSolution[i]][l] +
          distanceMatrix[initialSolution[j]][l] - 
          distanceMatrix[initialSolution[i]][initialSolution[j]];
        }
        insertionCost[k].insertedNode = l;
        insertionCost[k].deletedEdge = i;
        k++;
      }
    }

    // Orders insertion costs
    sort(insertionCost.begin(), insertionCost.end(), compares);

    int elements = max(1, (int)(alpha * insertionCost.size())); // Bugfix: ensure at least one element is chosen
    int i = rand() % elements;

    // Inserts choosen node
    initialSolution.insert(initialSolution.begin() + (insertionCost[i].deletedEdge + 1), insertionCost[i].insertedNode);

    // Deletes inserted node from candidates list
    for(int j = 0; j < candidatesList.size(); j++){
      if(candidatesList[j] == insertionCost[i].insertedNode)
        candidatesList.erase(candidatesList.begin()+j);
    }
  }

  return initialSolution;
}

vector<int> constructionSmith(vector<int> candidatesList,
                              const vector<double>& nodeWeights) {
  vector<int> solution;

  // Start at depot (assumes candidatesList[0] is depot)
  int depot = candidatesList[0];
  solution.push_back(depot);
  candidatesList.erase(candidatesList.begin());

  int current = depot;

  while (!candidatesList.empty()) {
    int bestIndex = 0;
    double bestScore = DBL_MAX;

    vector<pair<double,int>> scored;
    scored.reserve(candidatesList.size());

    for (int i = 0; i < (int)candidatesList.size(); i++) {
      int v = candidatesList[i];
      double w = nodeWeights[v];
      double p = distanceMatrix[current][v];
      double score = (w > 0.0) ? (p / w) : DBL_MAX;
      scored.push_back({score, i});
    }

    sort(scored.begin(), scored.end(),
        [](pair<double,int>& a, pair<double,int>& b){ return a.first < b.first; });

    // pick random from top-K
    int K = max(1, (int)(0.2 * scored.size())); // 20% RCL
    int pick = rand() % K;
    int pickedIndex = scored[pick].second;

    current = candidatesList[pickedIndex];
    solution.push_back(current);
    candidatesList.erase(candidatesList.begin() + pickedIndex);
  }

  return solution;
}


vector <int> pertub(vector <int> &solution){
  vector <int> newSolution;
	vector <int> firstSubsequence;
	vector <int> secondSubsequence;
	int i, j, sizeFirst = 0, sizeSecond = 0, size = solution.size();
  // Double Bridge: Randomly selects two non-overlapping subsequences of size at least 2 and swaps them
  // Available positions exclude depot at 0
  int freeSlots = size - 1;
  if (freeSlots <= 3) {
    return solution; // not enough space for two non-overlapping length-2 blocks
  }

  if (size < 20) {
    sizeFirst = 2;
    sizeSecond = 2;
  } else {
    while(sizeFirst < 2 || sizeFirst > (freeSlots/2)) sizeFirst = rand() % freeSlots;
    while(sizeSecond < 2 || sizeSecond > (freeSlots/2)) sizeSecond = rand() % freeSlots;
  }

  // Ensure they fit
  while (sizeFirst + sizeSecond > freeSlots) {
    if (sizeFirst > 2) sizeFirst--;
    if (sizeSecond > 2) sizeSecond--;
  }

// i and j must be in [1, size-1] and non-overlapping
int maxI = freeSlots - sizeFirst + 1;
int maxJ = freeSlots - sizeSecond + 1;

i = 1 + rand() % maxI;

int leftMax = i - sizeSecond;
int rightMin = i + sizeFirst;

bool hasLeft = leftMax >= 1;
bool hasRight = rightMin <= maxJ;

if (!hasLeft && !hasRight) {
  // if sizeFirst and sizeSecond are too large to fit on either side of i, we need to choose a new i
  do { i = 1 + rand() % maxI;
       leftMax = i - sizeSecond;
       rightMin = i + sizeFirst;
       hasLeft = leftMax >= 1;
       hasRight = rightMin <= maxJ;
  } while (!hasLeft && !hasRight);
}

if (hasLeft && hasRight) {
  if (rand() % 2 == 0) j = 1 + rand() % leftMax;
  else j = rightMin + rand() % (maxJ - rightMin + 1);
} else if (hasLeft) {
  j = 1 + rand() % leftMax;
} else {
  j = rightMin + rand() % (maxJ - rightMin + 1);
}

	// Creates vector with first subsequence
	int iter = 0;
	for(int k = 0; k < size; k++){
		if(k >= i){
			firstSubsequence.push_back(solution[k]);
      iter++;
			if(iter == sizeFirst) break;
		}
	}
	
	// Creates vector with second subsequence
	iter = 0;
	for(int k = 0; k < size; k++){
		if(k >= j){
			secondSubsequence.push_back(solution[k]);
      iter++;
			if(iter == sizeSecond) break;
		}
	}

	if(j < i){	
    newSolution = solution;

    // Deletes first subsequence
    int deleted = 0;
    for(int k = 0; k < size; k++){
      if(k >= i){
        newSolution.erase(newSolution.begin() + k - deleted);
        deleted++;
        if(deleted == sizeFirst) break;
      }
    }

    // Deletes second subsequence
    deleted = 0;
    for(int k = 0; k < size; k++){
      if(k >= j){
        newSolution.erase(newSolution.begin() + k - deleted);
        deleted++;
        if(deleted == sizeSecond) break;
      }
    }

    // Inserts second subsequence
    int inserted = 0;
    for(int k = 0; k < secondSubsequence.size(); k++){
      newSolution.insert(newSolution.begin() + i + inserted - sizeSecond, secondSubsequence[k]);
      inserted++;
    }
        
    // Inserts first subsequence
    inserted = 0;
    for(int k = 0; k < firstSubsequence.size(); k++){
      newSolution.insert(newSolution.begin() + j + inserted, firstSubsequence[k]);
      inserted++;
    }
      
	} else {
    newSolution = solution;

    // Deletes second subsequence
    int deleted = 0;
    for(int k = 0; k < size; k++){
      if(k >= j){
        newSolution.erase(newSolution.begin() + k - deleted);
        deleted++;
        if(deleted == sizeSecond) break;
      }
    }

    // Deletes first subsequence
    deleted = 0;
    for(int k = 0; k < size; k++){
      if(k >= i){
        newSolution.erase(newSolution.begin() + k - deleted);
        deleted++;
        if(deleted == sizeFirst) break; 
      }
	  }

    // Inserts first subsequence
    int inserted = 0;
    for(int k = 0; k < firstSubsequence.size(); k++){
      newSolution.insert(newSolution.begin() + j + inserted - sizeFirst, firstSubsequence[k]);
      inserted++;
    }
        
    // Inserts second subsequence
    inserted = 0;
    for(int k = 0; k < secondSubsequence.size(); k++){
      newSolution.insert(newSolution.begin() + i + inserted, secondSubsequence[k]);
      inserted++;
    }
  }
		
	return newSolution;	
}