/*
 * This file stores some global constants and reusable temp arrays.
 *
 * Author: William Deng
 */

#ifndef HELPER_H_
#define HELPER_H_
#include "randomc.h"
#include "graph.h"
#include <algorithm>

class Helper {
public:
  Helper(Graph* graph, int32 s) {
    seed = s;
    randomNumberGenerator = new TRandomMersenne(seed);
    initPheromone = new double[graph->getNumEdges()];
    edgeIndex = new int[graph->getNumEdges()];
    numEdges = graph->getNumEdges();
    numVertexes = graph->getNumVertexes();
    pheromonePartialSum = new double*[graph->getNumVertexes()];
    for (int i = 0; i < graph->getNumVertexes(); ++i) {
      pheromonePartialSum[i] = new double[graph->getNumVertexes()];
    }
  }
  bool isEdgeIndexConsistent() {
    bool bitmap[numEdges];
    fill(bitmap, bitmap + numEdges, false);
    for (int i = 0; i < numEdges; ++i) {
      int k = edgeIndex[i];
      if (bitmap[k])
        return false;
      bitmap[k] = true;
    }
    for (int i = 0; i < numEdges; ++i) {
      if (!bitmap[i])
        return false;
    }
    return true;
  }

  ~Helper() {
    delete[] edgeIndex;
    delete randomNumberGenerator;
    delete[] initPheromone;
    for (int i = 0; i < numVertexes; ++i) {
      delete[] pheromonePartialSum[i];
    }
    delete[] pheromonePartialSum;
  }

  TRandomMersenne* randomNumberGenerator; //random number generator
  double minPheromone; 
  double maxPheromone;
  double *initPheromone;
  double** pheromonePartialSum; // store the pheromone partial sum for each vertex
  int* edgeIndex; //store the order of edges after sorting
  int32 seed; //random seed
private:
  int numEdges;
  int numVertexes;
};

#endif /* HELPER_H_ */
