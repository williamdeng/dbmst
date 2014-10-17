#include "spanningtree.h"
#include "gtest/gtest.h"

class SpanningTreeTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    ifstream inFile("egraph.test");
    double maxCost, minCost;
    graph = readEGraph(inFile, maxCost, minCost);
  }

  virtual void TearDown() {
    delete graph;
  }

  Graph* graph;
};
TEST_F(SpanningTreeTest, Empty) {
  SpanningTree tree(graph);
  EXPECT_EQ(0, tree.size());
  EXPECT_TRUE(tree.repOK());
}

TEST_F(SpanningTreeTest, NonEmpty) {
  SpanningTree tree(graph);
  tree.setRoot(0);
  tree.addEdge(graph->getEdge(0, 1).index);
  EXPECT_TRUE(tree.repOK());
  tree.addEdge(graph->getEdge(1, 2).index);
  EXPECT_TRUE(tree.repOK());
  tree.addEdge(graph->getEdge(1, 3).index);
  EXPECT_EQ(3, tree.size());
  EXPECT_EQ(0, tree.getDepthOfNode(0));
  EXPECT_EQ(1, tree.getDepthOfNode(1));
  EXPECT_EQ(2, tree.getDepthOfNode(2));
  EXPECT_EQ(2, tree.getDepthOfNode(3));

  EXPECT_EQ(2, tree.getHeightOfNode(0));
  EXPECT_EQ(1, tree.getHeightOfNode(1));
  EXPECT_EQ(0, tree.getHeightOfNode(2));
  EXPECT_EQ(0, tree.getHeightOfNode(3));
  ASSERT_TRUE(tree.repOK());
  tree.clear();
  EXPECT_TRUE(tree.isEmpty());
  EXPECT_EQ(0, tree.size());
}

TEST_F(SpanningTreeTest, RemoveEdge) {
  SpanningTree tree(graph);
  tree.setRoot(0);
  tree.addEdge(graph->getEdge(0, 1).index);
  tree.addEdge(graph->getEdge(1, 2).index);
  tree.addEdge(graph->getEdge(1, 3).index);
  tree.removeEdge(graph->getEdge(0, 1).index);
  EXPECT_EQ(2, tree.size());
  EXPECT_EQ(0, tree.getHeightOfNode(0));
  EXPECT_EQ(1, tree.getHeightOfNode(1));
  EXPECT_EQ(0, tree.getHeightOfNode(2));
  EXPECT_EQ(0, tree.getHeightOfNode(3));
} 

TEST_F(SpanningTreeTest, RemoveEdgeNoChange) {
  SpanningTree tree(graph);
  tree.setRoot(0);
  tree.addEdge(graph->getEdge(0, 1).index);
  tree.addEdge(graph->getEdge(1, 2).index);
  tree.addEdge(graph->getEdge(1, 3).index);
  EXPECT_TRUE(tree.isComplete());
  tree.removeEdge(graph->getEdge(1, 2).index);
  EXPECT_EQ(2, tree.size());
  EXPECT_EQ(2, tree.getHeightOfNode(0));
  EXPECT_EQ(1, tree.getHeightOfNode(1));
  EXPECT_EQ(0, tree.getHeightOfNode(2));
  EXPECT_EQ(0, tree.getHeightOfNode(3));
} 

TEST_F(SpanningTreeTest, RemoveEdgeAddEdge) {
  SpanningTree tree(graph);
  tree.setRoot(0);
  tree.addEdge(graph->getEdge(0, 1).index);
  tree.addEdge(graph->getEdge(1, 2).index);
  tree.addEdge(graph->getEdge(1, 3).index);
  tree.removeEdge(graph->getEdge(1, 2).index);
  tree.addEdge(graph->getEdge(2, 3).index);
  EXPECT_EQ(3, tree.size());
  
  EXPECT_EQ(0, tree.getDepthOfNode(0));
  EXPECT_EQ(1, tree.getDepthOfNode(1));
  EXPECT_EQ(3, tree.getDepthOfNode(2));
  EXPECT_EQ(2, tree.getDepthOfNode(3));

  EXPECT_EQ(3, tree.getHeightOfNode(0));
  EXPECT_EQ(2, tree.getHeightOfNode(1));
  EXPECT_EQ(0, tree.getHeightOfNode(2));
  EXPECT_EQ(1, tree.getHeightOfNode(3));
} 
