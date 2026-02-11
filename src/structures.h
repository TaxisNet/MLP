#ifndef STRUCTURES_H
#define STRUCTURES_H

struct insertionInfo{
  int insertedNode;
  int deletedEdge; 
  double cost;
};

struct neighborInfo{
  int iBest;
  int jBest;
  double bestCost;
};

struct subsequenceInfo{
  double totalTime;        // T(S): Sum of distances in segment [i,j]
  double acumulateCost;    // C(S): Already weighted cost within [i,j]
  double sumWeights;       // W(S): Sum of all node weights in segment [i,j]
  int vertices;            // Number of vertices in segment
};

#endif