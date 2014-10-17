/*
 * Implementation of the graph structure.
 *
 * Author: William Deng
 */
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <limits>
#include "graph.h"

/***   CREATES AN EDGE     ***/
Edge::Edge()
{
  cost = 0;
  pheromone = 0;
  v1 = -1;
  v2 = -1;
  index = -1;
  update = 0;
}

Graph* readEGraph(ifstream& inFile, double& maxCost, double& minCost) {
  int numVertexes;
  inFile >> numVertexes;
  Graph* graph = new Graph(numVertexes);
  double* xCord = new double[numVertexes];
  double* yCord = new double[numVertexes];
  for (int i = 0; i < numVertexes; i++) {
    inFile >> xCord[i] >> yCord[i];
  }
  maxCost = 0; // assuming that all costs are > 0
  minCost = std::numeric_limits<double>::max();
  for (int i = 0; i < numVertexes; i++) {
    for (int j = i + 1; j < numVertexes; j++) {
      double cost = sqrt((xCord[i] - xCord[j]) * (xCord[i] - xCord[j]) +
                         (yCord[i] - yCord[j]) * (yCord[i] - yCord[j]));
      if (maxCost < cost) {
        maxCost = cost;
      }
      if (minCost > cost) {
        minCost = cost;
      }
      Edge& edge = graph->getEdge(i,j);
      edge.v1 = i;
      edge.v2 = j;
      edge.cost = cost;
      edge.index = graph->edgeIndex(edge.v1, edge.v2);
    }
  }
  delete[] xCord;
  delete[] yCord;
  return graph;
}

Graph* readRGraph(ifstream& inFile, double& maxCost, double& minCost) {
  int numVertexes;
  inFile >> numVertexes;
  Graph* graph = new Graph(numVertexes);
  maxCost = 0; // assuming that all costs are > 0
  minCost = std::numeric_limits<double>::max();
  for (int i = 0; i < numVertexes; i++) {
    for (int j = 0; j < numVertexes; j++) {
      double cost;
      inFile >> cost;
      if (i == j) {
	assert(cost == 0);
      }
      if (i >= j) {
	continue;
      }
      assert(cost > 0);
      if (maxCost < cost) {
        maxCost = cost;
      }
      if (minCost > cost) {
        minCost = cost;
      }
      Edge& edge = graph->getEdge(i,j);
      edge.v1 = i;
      edge.v2 = j;
      edge.cost = cost;
      edge.index = graph->edgeIndex(edge.v1, edge.v2);
    }
  }
  return graph;
}
