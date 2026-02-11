#include "readData.h"
#include "search.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <cfloat>
#include <algorithm>

using namespace std;

double ** distanceMatrix; 
int dimension; 





int main(int argc, char** argv) {

  clock_t start = clock(); // Starts time counting
  

  // Make a fake command line argument list
  if (argc < 2) {
    char* fake_argv[] = {(char*)"main", (char*)"/home/taxis/Documents/MLP/instances/berlin52.tsp"};
    readData(2, fake_argv, &dimension, &distanceMatrix);
  }else {
    readData(argc, argv, &dimension, &distanceMatrix);
  }

  srand(time(NULL));

  int iIls;
    
  if(dimension >= 100){
    iIls = 100;
  } else {
    iIls = dimension;
  }
  // cout << "iIls: " << iIls << endl;
  double cost = search(iIls, dimension);

  // Ends time counting
  clock_t end = clock();
  double time = ((double) (end - start)) / CLOCKS_PER_SEC;

  cout << endl << "Cost: " << cost << endl;
  cout << "Time: " << time << endl << endl;
    
  return 0;

}