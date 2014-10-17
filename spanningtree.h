/*
 * Data structure for spanning tree.
 *
 * Author: William Deng
 */

#ifndef SPANNINGTREE_H_
#define SPANNINGTREE_H_
#include <cassert>
#include "indexset.h"
#include "graph.h"
#include <algorithm>

class SpanningTree {
public:
  SpanningTree(Graph* graph);
  SpanningTree(const SpanningTree& source);
  void operator=(const SpanningTree& source);
  ~SpanningTree();
  void clear();
  void removeEdge(int index);
  void addEdge(int index);
  //access the index-th edge in the spanning tree
  int operator[](int index) const {
    return (*treeEdges)[index];
  }
  bool hasEdge(int edgeIndex) const {
    return treeEdges->contains(edgeIndex);
  }
  int size() const {
    return treeEdges->size();
  }
  bool isEmpty() const {
    return treeEdges->isEmpty();
  }
  bool isComplete() const {
    return treeEdges->size() + 1 == graph->getNumVertexes();
  }
  int getRoot() const {
    return root;
  }
  double getCost() const {
    return treeCost;
  }
  bool isNodeInTree(int node) {
    return node == root || depth[node] > 0;
  }
  int getDepthOfNode(int node) {
    return depth[node];
  }
  int getHeightOfNode(int node) {
    return height[node];
  }
  int getParentOfNode(int node) {
    return parent[node];
  }
  void setRoot(int r);
  bool repOK();
  IndexSet& getNeighbors(int node) const;

  friend std::ostream & operator <<(std::ostream& outs,
                                    const SpanningTree& tree);
private:
  void fixDepth(int node, int d);
  int computeHeight(int node);
  IndexSet* treeEdges;
  Graph* graph;
  double treeCost;
  int* depth;
  int* height;
  int* parent; // parent[i] is the parent of node i, -1 if i is a root
  IndexSet** children;
  int root;
  int numVertexes;
};

#endif /* SPANNINGTREE_H_ */
