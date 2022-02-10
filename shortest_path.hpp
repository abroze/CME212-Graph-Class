/**
 * @file shortest_path.cpp
 * Implimentation file for using our templated Graph to determine shortest paths.
 */


#include <vector>
#include <fstream>
#include <queue>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"

// Define our types
using GraphType = Graph<int>;
using NodeType  = typename GraphType::node_type;
using NodeIter  = typename GraphType::node_iterator;

/** Find the node with the minimum euclidean distance to a point.
 * @param g  The graph of nodes to search.
 * @param point  The point to use as the query.
 * @return An iterator to the node of @a g with the minimun Eucliean
 *           distance to @a point.
 *           graph.node_end() if graph.num_nodes() == 0.
 *
 * @post For all i, 0 <= i < graph.num_nodes(),
 *          norm(point - *result) <= norm(point - g.node(i).position())
 */

struct CompareNodes {
private:
  const Point& point_;

public:
  CompareNodes(const Point& point) :
  point_(point) {}

  bool operator()(const NodeType& n1, const NodeType& n2) {
    Point pt_diff_1 = n1.position() - point_;
    Point pt_diff_2 = n2.position() - point_;

    double dist1 = norm_2(pt_diff_1);
    double dist2 = norm_2(pt_diff_2);

    return (dist1 < dist2);
  }
};

NodeIter nearest_node(const GraphType& g, const Point& point) {
  NodeIter node_it = std::min_element(g.node_begin(), g.node_end(), CompareNodes(point));
  return node_it;
}


/** Update a graph with the shortest path lengths from a root node.
 * @param[in,out] g     Input graph
 * @param[in,out] root  Root node to start the search.
 * @return The maximum path length found.
 *
 * @post root.value() == 0
 * @post Graph has modified node values indicating the minimum path length
 *           to the root.
 * @post Graph nodes that are unreachable from the root have value() == -1.
 *
 * This sets all nodes' value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(GraphType& g, NodeType& root) {

  NodeIter node_it = g.node_begin();
  while(node_it != g.node_end()){
    (*node_it).value() = -1;
    ++node_it;
  }

  std::queue<NodeType> Q;

  root.value() = 0;
  Q.push(root);

  while(!(Q.empty())) {
    NodeType current = Q.front();
    Q.pop();

    auto edge_it = current.edge_begin();
    while (edge_it != current.edge_end()) {
      if ((*edge_it).node2().value() == -1) {
        (*edge_it).node2().value() = (*edge_it).node1().value() + 1;
        Q.push((*edge_it).node2());
      }
      ++edge_it;
    }
  }

  int longest_path = 0;

  auto ni = g.node_begin();
  while (ni != g.node_end()) {
    NodeType node = *ni;
    if (node.value() > longest_path)
      longest_path = node.value();
    ++ni;
  }
  return longest_path;
}

