/*
 * This file stores global parameters.
 *
 * Author: William Deng
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_
struct Parameters {
  int diameter;
  int steps; /* # edges ants travel per iteration	*/
  int pheromoneUpdatePeriod;
  double decay;
  double pheromoneEnhance;

  Parameters(Graph* graph, int d) {
    diameter = d;
    steps = (int) ((2.0 / 3.0) * (graph->getNumVertexes()));
    pheromoneUpdatePeriod = steps / 3;
    decay = 0.7;
    pheromoneEnhance = 1.5;
  }
};

#endif /* PARAMETERS_H_ */
