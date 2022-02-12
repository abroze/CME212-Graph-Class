#include <gtest/gtest.h>
#include <fstream>
#include <iostream>

#include "CME212/Util.hpp"
#include "Graph.hpp"
#include "shortest_path.hpp"
#include "subgraph.hpp"


class GraphPointFixture : public ::testing::Test {
 protected:
   //Define types
  using GraphType = Graph<int, int>;
  using NodeType  = typename GraphType::node_type;
  using NodeIter = typename GraphType::node_iterator;
  using EdgeType  = typename GraphType::edge_type;

  //Set up Graph and Points
  GraphType graph;
  std::vector<Point> points;
  virtual void SetUp() {
    for(int i = 0; i < 10; i++)
      points.push_back(Point(i));
  }
  
};

// Test adding node with default value
TEST_F(GraphPointFixture, DefaultNodeVal){
  graph.add_node(points[0]);
  EXPECT_EQ( graph.node(0).value(), 0 ) << "add_node does not intalize node vale with a default 0 value";
}

// Test degree function
TEST_F(GraphPointFixture, Degree){
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  graph.add_edge(n0, n1);

  EXPECT_EQ(n2.degree(),0)  << "n2 degree is 0";
  EXPECT_EQ(n1.degree(), 1) << "n1 degree is 1";
}

// Test node iterator
TEST_F(GraphPointFixture, NodeIter){
  graph.add_node(points[0]);
  graph.add_node(points[1]);
  graph.add_node(points[2]);
  
  int iter = 0;
  for(auto ni = graph.node_begin(); ni != graph.node_end(); ++ni){
    ++iter;
  }
  EXPECT_EQ(iter, 3) << " error in node iteration " ;
}

// Test incident iterator
TEST_F(GraphPointFixture, IncidentIter){
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);
  NodeType n4 = graph.add_node(points[4]);
  NodeType n5 = graph.add_node(points[5]);
  graph.add_edge(n0, n1);
  graph.add_edge(n0, n2);
  graph.add_edge(n0, n3);
  graph.add_edge(n0, n4);
  graph.add_edge(n1, n2);
  graph.add_edge(n1, n3);
  graph.add_edge(n2, n3);
  graph.add_edge(n3, n4);

  for(auto ni = graph.node_begin(); ni != graph.node_end(); ++ni){
    NodeType n = *ni;
    for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
      EdgeType e = *ii;
      std::cout << "node1 index : " << e.node1().index() << std::endl;
      EXPECT_EQ(e.node1(), n) << " error in incident iteration " << e.node1().position() << "  " << n.position();
    }
  }
}


//Test nearest node
TEST_F(GraphPointFixture, ShortestPath){
  graph.add_node(points[0]);
  graph.add_node(points[1]);
  graph.add_node(points[2]);
  
  NodeIter nearest = nearest_node(graph, Point(0));
  EXPECT_EQ( *nearest, graph.node(0)) << " error finding nearest node " ;
}
