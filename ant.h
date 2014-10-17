#ifndef _ant_h_
#define _ant_h_

class Ant
{
public:
  static const int TABOO_LIST_SIZE = 5;
  Ant(int initLoc);
  void clearVisitedSet();
  bool visited(int node) const;
  void visit(int node);
  ~Ant() { delete[] visitedNodes; }
  int location() { return currentNode; }
  void moveTo(int loc) { currentNode = loc; }
private:
  int currentNode;
  int *visitedNodes; // taboo list
  int nextVisitedNode;// taboo list front
};

#endif
