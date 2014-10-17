/*
 * The main function of bounded diameter minimal spanning tree using ant-based
 * algorithm.
 *
 * Author: William Deng
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <gflags/gflags.h>
#include "randomc.h"
#include "ant.h"
#include "graph.h"
#include "indexset.h"
#include "parameters.h"
#include "helper.h"
#include "spanningtree.h"
#include "queue.h"

DEFINE_string(file_name, "", "Input graph file name");
DEFINE_bool(euclidean, false, "Whether the input is euclid type");
DEFINE_bool(random, false, "Whether the input is random type");
DEFINE_int32(diameter, 10, "Max diameter of the spanning tree");
DEFINE_int64(random_seed, time(NULL), "Random seed");
DEFINE_bool(local_opt_selection_random, false, "Whether to pick a random edge or high cost edge (10% of the top) to be replaced");
DEFINE_int32(local_opt_replacement_strategy, 1, "0: random; 1: greedy with only positive improvement; 2: greedy with negative improvement allowed");

using namespace std;

//find the smallest
//invariant numbers[low] <= value < numbers[high]
inline
int binarySearch(double value, double numbers[], int high, int low) {
  assert(high >= low && value < numbers[high] && value >= numbers[low]);
  int mid;
  while (high != low + 1) {
    mid = (low + high) / 2;
    if (value < numbers[mid]) {
      high = mid;
    } else {
      low = mid;
    }
  }
  return high;
}

class compare_pheromone {
public:
  compare_pheromone(Graph* graph) {
    this->graph = graph;
  }

  bool operator()(int a, int b) {
    return (graph->getEdges())[a].pheromone > (graph->getEdges())[b].pheromone;
  }
private:
  Graph* graph;
};

class compare_cost {
public:

  compare_cost(Graph* graph) {
    this->graph = graph;
  }

  bool operator()(int a, int b) {
    return (graph->getEdges())[a].cost < (graph->getEdges())[b].cost;
  }
private:
  Graph* graph;
};

void setPheromoneInBound(Graph* graph, Edge& edge, const struct Helper& auxData) {
  if (edge.pheromone > auxData.maxPheromone) {
    edge.pheromone = auxData.maxPheromone
                     + auxData.initPheromone[graph->edgeIndex(edge.v1, edge.v2)];
  } else {
    if (edge.pheromone < auxData.minPheromone) {
      edge.pheromone = auxData.minPheromone
                       + auxData.initPheromone[graph->edgeIndex(edge.v1, edge.v2)];
    }
  }
}

void updatePheromone(Graph* graph, const struct Helper& auxData,
                     const Parameters& parameters) {
  for (int e = 0; e < graph->getNumEdges(); e++) {
    graph->getEdge(e).pheromone = (graph->getEdge(e).pheromone
                                   * parameters.decay) + graph->getEdge(e).update
                                  * (auxData.initPheromone)[e];

    setPheromoneInBound(graph, graph->getEdge(e), auxData);
    graph->getEdge(e).update = 0;
  }
}

void setupPheromonePartialSumArray(Graph * graph, double** pherPartialSum) {
  //for vertex i, partial sum of neighbor pheromone of vertex i: i, 0, 1, ..., i-1, i+1, ...., numVertexes-1
  for (int i = 0; i < graph->getNumVertexes(); ++i) {
    pherPartialSum[i][0] = 0.0;
    for (int j = 0; j < graph->getNumVertexes(); ++j) {
      if (j < i) {
        pherPartialSum[i][j + 1] = pherPartialSum[i][j] + graph->getEdge(i, j).pheromone;
      } else if (j > i) { // [i][i] represents i->i+1
        pherPartialSum[i][j] = pherPartialSum[i][j - 1] + graph->getEdge(i, j).pheromone;
      }
    }
  }
}

void antsMove(Ant** ants, const Parameters& parameters, struct Helper& auxData,
              Graph * graph) {
  setupPheromonePartialSumArray(graph, auxData.pheromonePartialSum);
  for (int j = 0; j < parameters.steps; j++) {
    if (j % parameters.pheromoneUpdatePeriod == 0 && j != 0) {
      updatePheromone(graph, auxData, parameters);
      setupPheromonePartialSumArray(graph, auxData.pheromonePartialSum);
    }
    for (int a = 0; a < graph->getNumVertexes(); a++) {
      int currentVertex = ants[a]->location();
      ants[a]->visit(currentVertex);
      int numTries = 0;
      int maxTries = 5;
      while (numTries < maxTries) {
        double rand_pher = auxData.randomNumberGenerator->Random()
                           * auxData.pheromonePartialSum[currentVertex][graph->getNumVertexes()
                               - 1];
        int nextVertex = binarySearch(rand_pher,
                                      auxData.pheromonePartialSum[currentVertex], graph->getNumVertexes()
                                      - 1, 0);
        if (nextVertex <= currentVertex)
          nextVertex--;
        assert(nextVertex >= 0 && nextVertex != currentVertex);
        if (!ants[a]->visited(nextVertex)) {
          /* New vertex, travel this edge				*/
          graph->getEdge(currentVertex, nextVertex).update++;
          ants[a]->moveTo(nextVertex);
          break;
        } else { // vertex is in the taboo list.
          numTries++;
        }
      }
    }
  }
}

void constructSpanningTree(Graph* graph, SpanningTree& currentSpanningTree,
                           int diameter, Helper & auxData) {
  for (int e = 0; e < graph->getNumEdges(); e++) {
    auxData.edgeIndex[e] = e;
  }
  compare_pheromone compare_obj_down(graph);

  sort(auxData.edgeIndex, auxData.edgeIndex + graph->getNumEdges(), compare_obj_down);
  currentSpanningTree.clear();
  int starting = 0;
  int ending = 0;
  int current = 0;
  compare_cost compareObj(graph);
  assert(auxData.isEdgeIndexConsistent());
  while (!currentSpanningTree.isComplete() && current < graph->getNumEdges()) {
    if (ending <= current) {
      starting = ending;
      ending += 5 * graph->getNumVertexes();

      assert(starting < graph->getNumEdges());

      if (ending > graph->getNumEdges())
        ending = graph->getNumEdges();
      sort(auxData.edgeIndex + starting, auxData.edgeIndex + ending, compareObj);
      assert(auxData.isEdgeIndexConsistent());
    }
    int v1 = graph->getEdge(auxData.edgeIndex[current]).v1;
    int v2 = graph->getEdge(auxData.edgeIndex[current]).v2;
    if (currentSpanningTree.isEmpty()) {
      currentSpanningTree.setRoot(v1);
      currentSpanningTree.addEdge(auxData.edgeIndex[current]);
    } else {
      bool suitable = false;
      if (currentSpanningTree.isNodeInTree(v1) ^ currentSpanningTree.isNodeInTree(v2)) {
        suitable = (currentSpanningTree.isNodeInTree(v1) && (currentSpanningTree.getDepthOfNode(v1) + 1) <= diameter / 2) ||
	  (currentSpanningTree.isNodeInTree(v2) && (currentSpanningTree.getDepthOfNode(v2) + 1) <= diameter / 2);
      }
      if (suitable) {
        currentSpanningTree.addEdge(auxData.edgeIndex[current]);
      }
    }
    current++;
  }
}

void updateEdgePheromone(Graph* graph, Edge& edge,
                         const Helper& auxData, double factor) {
  edge.pheromone *= factor;
  setPheromoneInBound(graph, edge, auxData);
}

void restart(Graph* graph, SpanningTree& tree, Helper& auxData, Parameters& parameters) {
  for (int e = 0; e < tree.size(); e++) {
    double factor = auxData.randomNumberGenerator->IRandom(10, 30)
                    / 100.0; /* Evaporate all but 10-30% of pheromone  */
    updateEdgePheromone(graph, graph->getEdge(tree[e]), auxData, factor);
  }
}

int main(int argc, char *argv[]) {
  google::ParseCommandLineFlags(&argc, &argv, true);
  if (FLAGS_local_opt_replacement_strategy < 0 || 
      FLAGS_local_opt_replacement_strategy > 2) {
    std::cout << "local_opt_replacement_strategy can only be 0, 1, or 2\n";
    return 1;
  }

  void localOneOpt(SpanningTree& tree, Graph* graph, Helper& auxData, int diameter);

  const int ITERATIONS = 1;
  const int MAX_NO_IMPROVEMENT_ITERATIONS = 2500;
  const int RESTART_THRESHOLD = 100;

  time_t startTime = time(NULL), endTime;
  ifstream inFile;
  inFile.open(FLAGS_file_name.c_str());
  Graph* graph;
  double maxCost = 0, minCost = 0;
  if (FLAGS_euclidean) {
    graph = readEGraph(inFile, maxCost, minCost);
  } else if (FLAGS_random) {
    graph = readRGraph(inFile, maxCost, minCost);
  } else {
    cout << "Unknown file type" << endl;
    return 1;
  }

  Parameters parameters(graph, FLAGS_diameter);
  Helper auxData(graph, FLAGS_random_seed);
  auxData.minPheromone = (maxCost - minCost) / 3.0;
  auxData.maxPheromone = (maxCost - minCost) / 3.0 * 4000;
  for (int e = 0; e < graph->getNumEdges(); e++) {
    auxData.initPheromone[e] = (maxCost - graph->getEdge(e).cost)
                               + (maxCost - minCost) / 3.0;
    graph->getEdge(e).pheromone = auxData.initPheromone[e];
  }

  // each vertex has one ant in the beginning.
  Ant** ants = new Ant*[graph->getNumVertexes()];
  for (int i = 0; i < graph->getNumVertexes(); i++) {
    ants[i] = new Ant(i);
  }

  int iterationBest = 0;
  int iterationRestart = 0;
  SpanningTree bestTree(graph);
  SpanningTree currentTree(graph);
  int i = 0;
  for (i = 0; i < ITERATIONS; i++) {
    if (((i + 1) % 1000) == 0) {
      parameters.decay = parameters.decay * 1.05;
      if (parameters.decay > 0.98) {
        parameters.decay = 0.98;
      }
      parameters.pheromoneEnhance = parameters.pheromoneEnhance * 1.05;
    }

    for (int a = 0; a < graph->getNumVertexes(); a++) {
      if (i > 0) {
        //half of the ants remain the same places and other half move to random vertexes
        if (auxData.randomNumberGenerator->Random() > 0.5) {
          ants[a]->moveTo(auxData.randomNumberGenerator->IRandom(0, graph->getNumVertexes() - 1));
        }
      }
      ants[a]->clearVisitedSet();
    }

    antsMove(ants, parameters, auxData, graph);
    updatePheromone(graph, auxData, parameters);

    currentTree.clear();
    constructSpanningTree(graph, currentTree, parameters.diameter, auxData);
    localOneOpt(currentTree, graph, auxData, parameters.diameter);
    assert(currentTree.isComplete() && currentTree.repOK());
    if (bestTree.size() == 0 || bestTree.getCost() > currentTree.getCost()) {
      bestTree = currentTree;
      iterationBest = i;
    }

    for (int e = 0; e < bestTree.size(); e++) {
      updateEdgePheromone(graph, graph->getEdge(bestTree[e]),
                          auxData, parameters.pheromoneEnhance);
    }
    if (i - iterationBest > MAX_NO_IMPROVEMENT_ITERATIONS)
      break;
    if ((i - iterationBest) > RESTART_THRESHOLD && (i - iterationRestart)
        > RESTART_THRESHOLD) {//restart

      restart(graph, bestTree, auxData, parameters);

      iterationRestart = i;
    }
  }
  endTime = time(NULL);
  cout << bestTree.getCost() << "\t" << i << "\t" << difftime(endTime, startTime)
       << "\t" << auxData.seed << endl;
  const char* fn = "tree.txt";
  if(argc == 5)
    fn = argv[4];
  ofstream outTree(fn);
  outTree << bestTree;
  outTree.close();
  for (int a = 0; a < graph->getNumVertexes(); a++) {
    delete ants[a];
  }
  delete[] ants;
  delete graph;
  return 0;
}

// check whether edge (sroot, cand) can reconnect the subtree rooted at sroot back
// to the big tree.
bool isFeasible(SpanningTree& tree, Graph* graph, int diameter, int sroot, int cand) {
  Edge& e = graph->getEdge(sroot, cand);
  if (tree.getDepthOfNode(cand) < diameter / 2 - tree.getHeightOfNode(sroot)  
      && !tree.hasEdge(e.index)) {
    // checks whether cand is in the subtree of sroot.
    if (tree.getDepthOfNode(cand) <= tree.getDepthOfNode(sroot)) {
      return true;
    } else {
       int node = cand;
       while (tree.getDepthOfNode(node) > tree.getDepthOfNode(sroot)) {
         node = tree.getParentOfNode(node);
       }
       return (node != sroot);
     }
  } 
  return false;
}

void localOneOpt(SpanningTree& tree, Graph* graph, Helper& auxData, int diameter) {
  assert(tree.getHeightOfNode(tree.getRoot()) <= diameter / 2);
  int MAX_NO_IMPROVEMENTS = 6;
  int numNoImprov = 0;
  // for non-random method, pick an edge with cost >= COST_PERCENTILE of all tree edges.
  double COST_PERCENTILE = 90; 
  int greedyStartIndex = (COST_PERCENTILE) / 100 * tree.size();
  compare_cost costComparator(graph);
  bool randomSelection = FLAGS_local_opt_selection_random;
  if (!randomSelection) {
    for (int i = 0; i < tree.size(); i++) {
      auxData.edgeIndex[i] = tree[i];
    }
  }
  SpanningTree bestTree(tree);
  while (numNoImprov < MAX_NO_IMPROVEMENTS) {
    int remEdgeIndex;
    int toRemoveIndexInAuxData = -1;
    if (randomSelection) {
      remEdgeIndex= tree[auxData.randomNumberGenerator->IRandom(0,
	  graph->getNumVertexes() - 2)];
    } else {
      toRemoveIndexInAuxData = auxData.randomNumberGenerator->IRandom(
	  greedyStartIndex, tree.size() - 1);
      nth_element(auxData.edgeIndex, auxData.edgeIndex +
	  toRemoveIndexInAuxData, auxData.edgeIndex + tree.size(), 
	  costComparator);
      remEdgeIndex = auxData.edgeIndex[toRemoveIndexInAuxData];
    }
    Edge remEdge = graph->getEdge(remEdgeIndex);
    if (tree.getDepthOfNode(remEdge.v1) == 0 || tree.getDepthOfNode(remEdge.v2) == 0) {
      continue;
    }
    int vert; 
    if (tree.getDepthOfNode(remEdge.v1) > tree.getDepthOfNode(remEdge.v2)) {
      vert = remEdge.v1;
    } else {
      vert = remEdge.v2;
    }
    int candidates[graph->getNumVertexes()];
    int numCandidates = 0;
    for (int i = 0; i < graph->getNumVertexes(); i++) {
      if (i == vert) { 
	continue;
      }
      Edge& e = graph->getEdge(vert, i);
      if (isFeasible(tree, graph, diameter, vert, i) && (e.cost < remEdge.cost
          || FLAGS_local_opt_replacement_strategy == 2)) {
        candidates[numCandidates++] = e.index;
      }
    }
    if (numCandidates != 0) {
      tree.removeEdge(remEdgeIndex);
      int replacementIndex;
      if (FLAGS_local_opt_replacement_strategy == 0) { // random 
	replacementIndex = auxData.randomNumberGenerator->IRandom(0, numCandidates - 1);
      } else { // greedy
	replacementIndex = 0;
	for (int k = 1; k < numCandidates; k++) {
	  if (graph->getEdge(candidates[k]).cost < graph->getEdge(
		candidates[replacementIndex]).cost) {
	    replacementIndex = k;
	  }
	}
      }
      tree.addEdge(candidates[replacementIndex]);
      if (!randomSelection) {
        auxData.edgeIndex[toRemoveIndexInAuxData] = candidates[replacementIndex];
      }
      assert(tree.repOK());
      assert(tree.getHeightOfNode(tree.getRoot()) <= diameter / 2);
    }
    if (tree.getCost() < bestTree.getCost()) {
      numNoImprov = 0;
      bestTree = tree;
    } else {
      numNoImprov++;
    }
    if (numNoImprov == MAX_NO_IMPROVEMENTS && !randomSelection) {
      randomSelection = true;
      numNoImprov = 0;
    }
  }
  if (tree.getCost() > bestTree.getCost()) {
    tree = bestTree;
  }
}

int findTreeDiameter(SpanningTree& tree, Graph& graph)
{
  int numVertices = tree.size()+1;  //tree.size() returns number of edges in tree
  Queue queue(numVertices+1);
  int levelMarker = -1;  //this is used to mark the end of each level as we do the BFS of the tree
  int numVerticesVisited = 0;
  int newRoot;
  bool done = false;

  bool* visited = new bool[numVertices];
  for (int i = 0; i < numVertices ; i++)
    visited[i] = false;

  queue.enqueue(tree.getRoot());
  visited[tree.getRoot()] = true;
  numVerticesVisited++;
  queue.enqueue(levelMarker);

  IndexSet neighbors(numVertices);
  //First BFS to find a vertex furthest from the root of tree.
  //This vertex will be called new root
  while (!queue.isEmpty() && !done){
    int u = queue.dequeue();
    if (u == levelMarker){
      queue.enqueue(levelMarker);  //mark the start of a new level
      continue;
    }
    if (numVerticesVisited == numVertices){// all vertices have been visited
      // let the newRoot be the last visited vertex
      newRoot = u;
      done = true;   //terminate BFS right away
      continue;
    }
    neighbors = tree.getNeighbors(u);
    int v;
    // add all unvisited neighbors of u into queue
    for (int i=0; i < neighbors.size(); i++){
      v = neighbors.getIthElement(i);//neighbors[i]
      if (!visited[v]){
	queue.enqueue(v);
	visited[v] = true;
	numVerticesVisited++;
      }
    }//for
  }//while

  //clean up in preparation for the 2nd BFS
  for (int i = 0; i < numVertices; i++){
    visited[i] = false;
  }
  queue.clear();
  done = false;
  
  // Start 2nd BFS from newRoot
  queue.enqueue(newRoot);
  visited[newRoot] = true;
  numVerticesVisited = 1;

  queue.enqueue(levelMarker);
  neighbors.clear();

  int diameter = 0;

  while (!queue.isEmpty() && !done){
    int u = queue.dequeue();
    if (u == levelMarker){
      queue.enqueue(levelMarker);
      diameter++;
      continue;
    }
    visited[u] = true;
    numVerticesVisited++;
    if (numVerticesVisited == numVertices){
      done = true;  //terminate 
      continue;
    }
    neighbors = tree.getNeighbors(u);
    int v;
    // add all unvisited neighbors of u into queue
    for (int i=0; i < neighbors.size(); i++){
      v = neighbors.getIthElement(i);
      if (!visited[v]){
	queue.enqueue(v);
	visited[v] = true;
	}
    }//for
  }//while
 
  return diameter;
}//findTreeDiameter
