/**
 * @file shortest_path.cpp
 * Implimentation file for using our templated Graph to determine shortest paths.
 */


#include <vector>
#include <fstream>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"

// Define our types
using GraphType = Graph<double>;
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
  Point p_;
  CompareNodes(const Point& p) : p_(p) {
  };

  template <typename NODE>
  bool operator()(const NODE& node1, const NODE& node2) const {
    Point diff1 = node1.position() - p_;
    Point diff2 = node2.position() - p_;
    if (norm(diff1) < norm(diff2)) return true;
    return false;
    }
};


Graph<int>::Node nearest_node(const GraphType& g, const Point& point)
{
  // HW1 #3: YOUR CODE HERE
  CompareNodes cn = CompareNodes(point);

  // Find the closest node to the given Point as root.
  Graph<int>::node_iterator nifirst = g.node_begin();
  Graph<int>::node_iterator nilast = g.node_end();
  Graph<int>::node_iterator niroot = std::min_element(nifirst, nilast, cn);
  Graph<int>::node_type root = *niroot;

  return root;
  
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
int shortest_path_lengths(GraphType& g, NodeType& root)
{
  /*
  // HW1 #3: YOUR CODE HERE
  // Set all the nodes' default values to -1 and the root's value to 0.
  for(; nifirst != nilast; ++nifirst){
    (*nifirst).value() = -1; 
  }
  root.value() = 0;

  // Set the current longest distance from root.
  int max = 0;

  /** The queue is to store the nodes needed to be evaluated.
   *  For each node n in the queue, we iterate the adjent nodes and set their values,
   *  push them into the queue and then pop n. In this process, update max.
   
  std::queue<Graph<int, int>::node_type> waiting;
  waiting.push(root);
  while(!waiting.empty()){
    Graph<int, int>::node_type r = waiting.front();
    int cur = r.value();
    Graph<int, int>::incident_iterator rbegin = r.edge_begin();
    Graph<int, int>::incident_iterator rend = r.edge_end();
    for(; rbegin != rend; ++rbegin){
      if((*rbegin).node2().value() == -1){
        (*rbegin).node2().value() = cur + 1;
        if((*rbegin).node2().value() > max) max = (*rbegin).node2().value();
        waiting.push((*rbegin).node2());
      }
      if((*rbegin).node2().value() > cur + 1){
        (*rbegin).node2().value() = cur + 1;
        if((*rbegin).node2().value() > max) max = (*rbegin).node2().value();
      }
    }
    waiting.pop();
  }
  return max;
  
*/
  return 0;
}

