/*
 * The header file for graph structure.
 *
 * Author: William Deng
 */
#ifndef GRAPH_H_
#define GRAPH_H_
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

class Edge {
public:
  Edge();

  double cost;
  double pheromone;
  int v1;
  int v2;
  int index;
  int update;
};

//Graph assume that complete graph
class Graph {
public:
  Graph(int nv) {
    edges = new Edge[nv * (nv - 1) / 2];
    numVertexes = nv;
    numEdges = nv * (nv - 1) / 2;
  }
  //get edge by two end vertexes
  //x and y are 0-based vertex index
  Edge& getEdge(int x, int y) {
    return edges[edgeIndex(x, y)];
  }
  Edge& getEdge(int index) {
    return edges[index];
  }
  ~Graph() {
    delete[] edges;
  }
  int getNumVertexes() const {
    return numVertexes;
  }
  int getNumEdges() const {
    return numEdges;
  }
  Edge* getEdges() {
    return edges;
  }
  // translates from two end points of an edge to the index
  int edgeIndex(int x, int y) const {
    assert (x >= 0 && x < numVertexes && y >= 0 && y < numVertexes);

    int idx;

    if (x <= y)
      idx = numVertexes * x - x * (x + 1) / 2 + (y - x - 1);
    else
      idx = numVertexes * y - y * (y + 1) / 2 + (x - y - 1);

    return idx;
  }

private:
  int numVertexes;
  int numEdges;
  Edge* edges;
};
// read euclidean graph
Graph* readEGraph(ifstream& inFile, double& max_cost, double& min_cost);
// read random graph
Graph* readRGraph(ifstream& inFile, double& max_cost, double& min_cost);
#endif /* GRAPH_H_ */
