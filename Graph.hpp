#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <iostream>
#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <chrono>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph: private totally_ordered <Graph<V>> {
 private:

  // predeclare the internal struct for node and edge
  struct internal_node;
  struct internal_edge; 
  

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of a node value. */
  using node_value_type = V;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Constructs an empty graph. */
  Graph() : nodes_(), edges_(), internal_nodes_(), internal_edges_(), 
  nodes_to_edge_(), nodes_adj_edges_() {};

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node: private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */

    /** Constructs an invalid node. */
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_.at(id_).position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id_;
    }

    node_value_type& value() {
      return (graph_->nodes_.at(id_).value_);
    }

    const node_value_type& value() const {
      return (graph_->nodes_.at(id_).value_);
    }

    size_type degree() const {
      return (graph_->nodes_.at(id_).degree_);
    }

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, id_, 0);
    }

    incident_iterator edge_end() const {
      return IncidentIterator(graph_, id_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (graph_ == n.graph_ && id_ == n.id_)
        return true;
      return false;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      if (graph_ == n.graph_ && id_ < n.id_)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    Graph* graph_;
    size_type id_;

    /** private Node constructor */
    Node(const graph_type* graph, size_type id): 
      graph_(const_cast<graph_type*>(graph)),
      id_(id)
      {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& a = node_value_type()) {
    internal_node new_internal_node;

    size_type new_internal_node_id = num_nodes();

    new_internal_node.position_ = position;
    new_internal_node.value_ = a;
    new_internal_node.degree_ = 0;

    nodes_.push_back(new_internal_node);
    std::pair<size_type, internal_node> pair(new_internal_node_id, new_internal_node);
    internal_nodes_.insert(pair);

    nodes_adj_edges_[new_internal_node_id];

    return Node(this, new_internal_node_id);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (this == n.graph_ && n.id_ < num_nodes())
      return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i >= 0);
    assert(i < num_nodes());
    return Node(this, i);
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(graph_->edges_[id_].node1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(graph_->edges_[id_].node2_id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ == e.graph_ && id_ == e.id_)
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_ && id_ < e.id_)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* graph_;
    size_type id_;

    /** private Edge constructor */
    Edge(const graph_type* graph, size_type id):
    graph_(const_cast<graph_type*>(graph)),
    id_(id)
    {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i >= 0);
    assert(i < num_edges());
    return Edge(this, i);
  }

  std::string create_edge_key(const Node& a, const Node& b) const {
    // get the nodes ids
    std::string a_id = std::to_string(a.id_);
    std::string b_id = std::to_string(b.id_);

    // create key
    std::string edge_key = a_id + "," + b_id;

    return edge_key;
  }

  Edge find_edge(const Node& a, const Node& b){
    std::string key1 = create_edge_key(a, b);
    std::string key2 = create_edge_key(b, a);

    if (nodes_to_edge_.find(key1) == nodes_to_edge_.end()) {
      size_type edge_id = nodes_to_edge_.find(key2)->second;
      return edge(edge_id);
    }
    else if (nodes_to_edge_.find(key2) == nodes_to_edge_.end()) {
      size_type edge_id = nodes_to_edge_.find(key1)->second;
      return edge(edge_id);
    }
    assert(false);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(!(a==b));
    std::string key1 = create_edge_key(a, b);
    std::string key2 = create_edge_key(b, a);
    if (nodes_to_edge_.find(key1) == nodes_to_edge_.end() &&
        nodes_to_edge_.find(key2) == nodes_to_edge_.end()) {
      return false;
    }
    return true;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    assert(!(a==b));

    if (has_edge(a, b))
      return find_edge(a, b);

    size_type new_internal_edge_id = num_edges();

    internal_edge new_internal_edge;

    new_internal_edge.node1_id_ = a.index();
    new_internal_edge.node2_id_ = b.index();
    new_internal_edge.edge_id_ = new_internal_edge_id;

    edges_.push_back(new_internal_edge);

    // update degree of both nodes
    nodes_[a.index()].degree_++;
    nodes_[b.index()].degree_++;

    // create pair for unordered map
    std::string key = create_edge_key(a, b);
    std::pair<std::string, size_type> pair(key, new_internal_edge_id);
    nodes_to_edge_.insert(pair);

    nodes_adj_edges_[a.index()].push_back(new_internal_edge_id);
    nodes_adj_edges_[b.index()].push_back(new_internal_edge_id);

    return Edge(this, new_internal_edge_id);
    }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    internal_nodes_.clear();
    internal_edges_.clear();
    nodes_to_edge_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    Node operator*() const {
      return Node(graph_, node_iter_id_);
    }

    NodeIterator& operator++() {
      node_iter_id_++;
      return (*this);
    }

    bool operator==(const NodeIterator& n) const {
      if (graph_ == n.graph_ && node_iter_id_ == n.node_iter_id_)
        return true;
      return false;
    }

   private:
    friend class Graph;

    Graph* graph_;
    size_type node_iter_id_;

    NodeIterator(const Graph* graph, size_type node_iter_id):
    graph_(const_cast<Graph*>(graph)),
    node_iter_id_(node_iter_id)
    {};

  };

  node_iterator node_begin() const {
    return NodeIterator(const_cast<Graph*>(this), 0);
  }

  node_iterator node_end() const {
    return NodeIterator(const_cast<Graph*>(this), num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    Edge operator*() const {
      return edge_val_;
    }

    IncidentIterator& operator++() {
      inc_edge_id_++;
      if (inc_edge_id_ != graph_->nodes_.at(node_id_).degree_) {
        edge_id_ = graph_->nodes_adj_edges_.at(node_id_).at(inc_edge_id_);
        edge_val_ = graph_->edge(edge_id_);
      }
      return (*this);
    }

    bool operator==(const IncidentIterator& i) const {
      if (graph_ == i.graph_ && node_id_ == i.node_id_ &&
          inc_edge_id_ == i.inc_edge_id_){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    friend class Node;

    Graph* graph_;
    size_type node_id_;
    size_type inc_edge_id_;
    size_type edge_id_;
    value_type edge_val_;

    IncidentIterator(const Graph* graph, size_type node_id, size_type inc_edge_id):
    
    graph_(const_cast<Graph*>(graph)),
    node_id_(node_id),
    inc_edge_id_(inc_edge_id)
    {
      if (inc_edge_id_ != graph_->nodes_.at(node_id_).degree_) {

        edge_id_ = graph_->nodes_adj_edges_.at(node_id).at(inc_edge_id_);
        edge_val_ = graph_->edge(edge_id_);

        if (graph_->edges_.at(edge_id_).node1_id_ != node_id_) {
          
          // swap node1_id and node2_id  
          graph_->edges_.at(edge_id_).node2_id_ = graph_->edges_.at(edge_id_).node1_id_;
          graph_->edges_.at(edge_id_).node1_id_ = node_id_;
        }
      }
    };

  };


  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    Edge operator*() const {
      return Edge(graph_, edge_iter_id_);
    }

    EdgeIterator& operator++() {
      edge_iter_id_++;
      return (*this);
    }

    bool operator==(const EdgeIterator& e) const {
      if (graph_ == e.graph_ && edge_iter_id_ == e.edge_iter_id_)
        return true;
      return false;
    }

   private:
    friend class Graph;
    
    Graph* graph_;
    size_type edge_iter_id_;

    EdgeIterator(Graph* graph, size_type edge_iter_id) :
    graph_(const_cast<graph_type*>(graph)),
    edge_iter_id_(edge_iter_id) {};
  };

  edge_iterator edge_begin() const {
    return EdgeIterator(const_cast<Graph*>(this), 0);
  }

  edge_iterator edge_end() const {
    return EdgeIterator(const_cast<Graph*>(this), num_edges());
  }

 private:

  struct internal_node {
    Point position_;
    node_value_type value_;
    size_type degree_;
  };

  struct internal_edge {
    size_type node1_id_;
    size_type node2_id_;
    size_type edge_id_;
  };

  // vectors to store nodes and edges in Graph
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;

  // maps to store internal nodes and internal edges in Graph
  std::map<size_type, internal_node> internal_nodes_;
  std::map<size_type, internal_edge> internal_edges_;

  // unordered map to store internal nodes and internal edges in Graph
  std::unordered_map<std::string, size_type> nodes_to_edge_;

  // unordered map to store nodes and vector with adjacent edges
  std::unordered_map<size_type, std::vector<size_type>> nodes_adj_edges_;

};

#endif // CME212_GRAPH_HPP