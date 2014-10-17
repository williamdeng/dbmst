#include "ant.h"
#include <algorithm>

Ant::Ant(int initLoc) {
  currentNode = initLoc;
  visitedNodes = new int[TABOO_LIST_SIZE];
  std::fill(visitedNodes, visitedNodes+TABOO_LIST_SIZE, -1);
  nextVisitedNode = 0;
}

void Ant::clearVisitedSet() {
  std::fill(visitedNodes, visitedNodes+TABOO_LIST_SIZE, -1);
  nextVisitedNode = 0;
}

bool Ant::visited(int node) const {
  for(int i=0; i<TABOO_LIST_SIZE; ++i) {
    if(node == visitedNodes[i]) {
      return true;
    }
  }
  return false;
}

void Ant::visit(int node) {
  visitedNodes[nextVisitedNode] = node;
  nextVisitedNode++;
  if(nextVisitedNode == TABOO_LIST_SIZE) {
    nextVisitedNode = 0;
  }
}
