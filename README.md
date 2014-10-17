This package contains a C++ implementation of ant-based algorithm for solving
the bounded diameter minimum spanning tree (DBMST) problem. DBMST is a variation
of standard minimum spanning tree (http://http://en.wikipedia.org/wiki/Minimum_spanning_tree) with an additional constraint that the diameter of the spanning
tree is less than or equal to some input constant d. DBMST is shown to be
NP-hard by Carey and Johnson in their "Computers and intractability: a guide
to the theory of NP-completeness" book.

This package tries to provide a ant-based algorithm for an approximate DBMST.
For more information about ant-based algorithm, please refer to wiki page:
http://en.wikipedia.org/wiki/Ant_colony_optimization_algorithms.
This code is a part of an ongoing research project by William Deng at
Google and Thang Bui at Penn State at Harrisburg.

The project includes following files:
Makefile: the build file (using Unix make)
ant.h/cpp: the ants that walk around the graph
bdmst.cpp: the main file
euclid100-1.txt: sample input Euclidean graph
graph.h/cpp: the graph data structure
helper.h: stores some global constants and reusable temp arrays
indexmap.h/cpp: a map with constant get ith position operation
indexset.h/cpp: a set with constant get ith position operation
mersenne.cpp and randomc.h: random number generator (NOT implemented by us)
parameters.h: some parameters used in the code
queue.h/cpp: circular array queue (by Thang Bui)
spanningtree.h/cpp: spanning tree  
There are two tests files indexmap_test.cpp and spanningtree_test.cpp which
requires gtest (https://code.google.com/p/googletest/).

To compile the project, in additional to Make, shell, and C++ compiler,
gflags library (https://code.google.com/p/gflags/) is also required.

