/*
 * spanningtree.cpp
 *
 * Author: William Deng
 */

#include <cmath>
#include "spanningtree.h"

SpanningTree::SpanningTree(Graph* graph) {
  numVertexes = graph->getNumVertexes();
  treeEdges = new IndexSet(graph->getNumEdges());
  this->graph = graph;
  treeCost = 0;
  depth = new int[numVertexes];
  fill(depth, depth + numVertexes, 0);
  height = new int[numVertexes];
  fill(height, height + numVertexes, 0);
  root = -1;
  parent = new int[numVertexes];
  fill(parent, parent + numVertexes, -1);
  children = new IndexSet*[numVertexes];
  for (int i = 0; i < numVertexes; i++) {
    children[i] = new IndexSet(numVertexes);
  }
}

SpanningTree::SpanningTree(const SpanningTree& source) {
  graph = source.graph;
  numVertexes = source.numVertexes;
  treeEdges = new IndexSet(source.treeEdges->getCapacity());
  *treeEdges = *(source.treeEdges);
  treeCost = source.treeCost;

  depth = new int[numVertexes];
  copy(source.depth, source.depth + numVertexes, depth);
  height = new int[numVertexes];
  copy(source.height, source.height + numVertexes, height);
  root = source.root;
  parent = new int[numVertexes];
  copy(source.parent, source.parent + numVertexes, parent);
  children = new IndexSet*[numVertexes];
  for (int i = 0; i < numVertexes; i++) {
    children[i] = new IndexSet(numVertexes);
    *(children[i]) = *(source.children[i]);
  }
}

void SpanningTree::operator=(const SpanningTree& source) {
  assert(source.graph == graph);
  *treeEdges = *(source.treeEdges);
  treeCost = source.treeCost;
  copy(source.depth, source.depth + numVertexes, depth);
  copy(source.height, source.height + numVertexes, height);
  root = source.root;
  copy(source.parent, source.parent + numVertexes, parent);
  for (int i = 0; i < numVertexes; i++) {
    *(children[i]) = *(source.children[i]);
  }
}

SpanningTree::~SpanningTree() {
  delete treeEdges;
  delete[] depth;
  delete[] height;
  delete[] parent;
  for (int i = 0; i < numVertexes; i++) {
    delete children[i];
  }
  delete[] children;
}

void SpanningTree::clear() {
  treeEdges->clear();
  treeCost = 0.0;
  fill(depth, depth + numVertexes, 0);
  fill(height, height + numVertexes, 0);
  fill(parent, parent + numVertexes, -1);
  for (int i = 0; i < numVertexes; i++) {
    children[i]->clear();
  }
}

void SpanningTree::removeEdge(int index) {
  assert(0 <= index && index < treeEdges->getCapacity() && treeEdges->contains(index));
  treeEdges->remove(index);
  int v1 = graph->getEdge(index).v1;
  int v2 = graph->getEdge(index).v2;
  assert(graph->edgeIndex(v1, v2) == index);
  assert(parent[v1] == v2 || parent[v2] == v1);
  int p /*parent*/, c /*child*/;
  if (parent[v1] == v2) {
    p = v2;
    c = v1;
  } else {
    p = v1;
    c = v2;
  }
  parent[c] = -1;
  depth[c] = 0;
  children[p]->remove(c);
  treeCost -= graph->getEdge(index).cost;
  // fix heights of ancestors of p
  while (p != -1) {
    int maxChildHeight = -1;
    for (int i = 0; i < children[p]->size(); i++) {
      int child = children[p]->getIthElement(i);
      if (height[child] > maxChildHeight) {
        maxChildHeight = height[child];
      }
    }
    if (maxChildHeight + 1 == height[p]) { // no need to check parent of p anymore
      break;
    }
    assert(maxChildHeight + 1 < height[p]);
    height[p] = maxChildHeight + 1;
    p = parent[p];
  }
}

void SpanningTree::fixDepth(int node, int d) {
  depth[node] = d;
  for (int i = 0; i < children[node]->size(); i++) {
    fixDepth(children[node]->getIthElement(i), d + 1);
  }
}

void SpanningTree::addEdge(int index) {
  assert(0 <= index && index < treeEdges->getCapacity() && !treeEdges->contains(index));
  int v1 = graph->getEdge(index).v1;
  int v2 = graph->getEdge(index).v2;
  int inTreeNode, outTreeNode;
  if (parent[v1] == -1 && parent[v2] == -1) {
    assert(v1 == root || v2 == root);
    if (v1 == root) {
      inTreeNode = v1;
      outTreeNode = v2;
    } else {
      inTreeNode = v2;
      outTreeNode = v1;
    }
  } else {
    assert((parent[v1] == -1) ^ (parent[v2] == -1));
    if (parent[v1] == -1) {
      inTreeNode = v2;
      outTreeNode = v1;
    } else {
      inTreeNode = v1;
      outTreeNode = v2;
    }
  }
  parent[outTreeNode] = inTreeNode;
  depth[outTreeNode] = depth[inTreeNode] + 1;
  treeEdges->insert(index);
  assert(graph->edgeIndex(v1, v2) == index);
  children[inTreeNode]->insert(outTreeNode);
  fixDepth(outTreeNode, depth[inTreeNode] + 1);
  treeCost += graph->getEdge(index).cost;
  if (this->isComplete()) {
    if (height[root] == 0) {
      height[root] = computeHeight(root);
    } else {
      int p = inTreeNode;
      int c = outTreeNode;
      // invariant: c is a child of p and c's height increased
      while (p != -1) {
        if (height[p] >= height[c] + 1) {
          break;
        }
        height[p] = height[c] + 1;
        c = p;
        p = parent[p];
      }
    }
  }
}

// compute the height of node and compute and update
// the height array of all children of node.
int SpanningTree::computeHeight(int node) {
  int result = 0;
  for (int i = 0; i < children[node]->size(); i++) {
    int child = children[node]->getIthElement(i);
    height[child] = computeHeight(child);
    if (height[child] + 1 > result) {
      result = height[child] + 1;
    }
  }
  return result;
}

std::ostream & operator <<(std::ostream& outs,
                           const SpanningTree& tree) {
  outs << tree.numVertexes << " "
       << tree.graph->getNumEdges() << endl;
  for (int i = 0; i < tree.treeEdges->size(); ++i) {
    outs << "e " << tree.graph->getEdge(tree[i]).v1 + 1 << " "
         << tree.graph->getEdge(tree[i]).v2 + 1 << " "
         << tree.graph->getEdge(tree[i]).cost << endl;
  }
  return outs;
}

void SpanningTree::setRoot(int r) {
  root = r;
  parent[root] = -1;
}

bool SpanningTree::repOK() {
  if (root == -1) { // empty tree case
    if(treeEdges->size() != 0) {
      return false;
    }
    if (treeCost != 0) {
      return false;
    }
    for (int i = 0; i < numVertexes; i++) {
      if (depth[i] != 0) {
        return false;
      }
      if (height[i] != 0) {
        return false;
      }
      if (parent[i] != -1) {
        return false;
      }
      if (!children[i]->isEmpty()) {
        return false;
      }
    }
    return true;
  }
  // the tree is not empty case
  if (treeEdges->size() == 0) {
    return false;
  }
  double totalCost = 0;
  for (int i = 0; i < treeEdges->size(); i++) {
    Edge e = graph->getEdge((*treeEdges)[i]);
    totalCost += e.cost;
    if (parent[e.v1] != e.v2 && parent[e.v2] != e.v1) {
      return false;
    }
    if (parent[e.v1] == e.v2 && depth[e.v1] != depth[e.v2] + 1) {
      return false;
    }
    if (parent[e.v1] == e.v2 && !children[e.v2]->contains(e.v1)) {
      return false;
    }
    if (parent[e.v2] == e.v1 && depth[e.v2] != depth[e.v1] + 1) {
      return false;
    }
    if (parent[e.v2] == e.v1 && !children[e.v1]->contains(e.v2)) {
      return false;
    }
  }
  if (isComplete()) {
    for (int node = 0; node < graph->getNumVertexes(); node++) {
      int maxChildHeight = -1;
      for (int i = 0; i < children[node]->size(); i++) {
        int child = children[node]->getIthElement(i);
        if (height[child] > maxChildHeight) {
          maxChildHeight = height[child];
        }
      }
      if (maxChildHeight + 1 != height[node]) {
        return false;
      }
    }
  }
  double DELTA = 1e-6;
  if (fabs(totalCost - treeCost) > DELTA) {
    return false;
  }
  if (root != -1 && parent[root] != -1) {
    return false;
  }
  return true;
}

IndexSet& SpanningTree::getNeighbors(int node) const
{ // return a set containing the neighbors of node
  if (parent[node]==-1)
    return *children[node];

  IndexSet* neighbors;
  neighbors = new IndexSet(*children[node]);
  neighbors->insert(parent[node]);
  return *neighbors;
}//getNeighbors
