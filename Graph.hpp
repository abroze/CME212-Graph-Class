#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template <typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct NodeElement;
  struct EdgeElement;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  using node_value_type = V;

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

  /** Construct an empty graph. */
  Graph() {
    n_nodes = 0;
    n_edges = 0;
  }



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

  class Node : private totally_ordered<Node> {
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
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return (graph_->node_position(uid_));
    }

    Point& position() {
      return (graph_->node_position(uid_));
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return (graph_->nodes_.at(uid_).uid_);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;


    node_value_type& value() {
      return (graph_->node(uid_));
    }

    const node_value_type& value() const {
      return (graph_->node(uid_));
    }

    // return the number of incident edges.
    size_type degree() const {
      return (graph_->adj_).at(uid_).size(); 
    }

    // Start of the incident iterator.
    incident_iterator edge_begin() const {
      return IncidentIterator((*this), 0);
    }

    // End of incident iterator.
    incident_iterator edge_end() const {
      return IncidentIterator((*this), degree());
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return(graph_ == n.graph_ && uid_ == n.uid_);
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
      if (graph_ == n.graph_ && uid_ < n.uid_){
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // Pointer back to the Graph container
    Graph* graph_;
    // This element's unique identification number
    size_type uid_;


    /** Private Constructor */
    Node(const Graph* graph, size_type uid) :
    graph_(const_cast<Graph*>(graph)),
    uid_(uid) {}

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return n_nodes;
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
  Node add_node(const Point& position, 
                const node_value_type& value = node_value_type()) {

    Node n = Node(this, next_uid());
    nodes_.push_back(NodeElement(position, value, n_nodes));
    idx2uid_.push_back(n.uid_);

    n_nodes++;

    adj_.push_back(std::vector<EdgeElement>());

    return n;

  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_ == this && n.index() < n_nodes)
      return (n.uid_ == idx2uid_[n.index()]);
    else
      return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    size_type uid = idx2uid_.at(i);
    return Node(this, uid);
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_id_);  
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_id_);     
    }

    /** Return the distance between the two nodes */
    double length() const {
      return norm(node1().position() - node2().position());
    }



    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if(graph_ != e.graph_) {
        return false;
      }

      return (( node1_id_ == e.node1_id_ && node2_id_ == e.node2_id_)
          ||  ( node1_id_ == e.node2_id_ && node2_id_ == e.node1_id_));
    }


    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      size_type m1 = std::min(node1_id_, node2_id_);
      size_type m2 = std::min(e.node1_id_, e.node2_id_);

      if(graph_ != e.graph_) {
        return (graph_ < e.graph_);
      }
      else if(m1 != m2) {
        return (m1 < m2);
      }
      else
        return (std::max(node1_id_, node2_id_) < std::max(e.node1_id_, e.node2_id_));
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to the Graph container
    Graph* graph_;
    // This nodes' unique identification number
    size_type node1_id_;
    size_type node2_id_;

    /** Private constructor through Nodes, for valid Edge objects. */
    Edge(const Graph* graph, Node n1, Node n2) :
      graph_(const_cast<Graph*>(graph)),
      node1_id_(n1.uid_),
      node2_id_(n2.uid_) {}

    /** Private constructor through Node Ids, for valid Edge objects. */
    Edge(const Graph* graph, size_type id1, size_type id2) :
    graph_(const_cast<Graph*>(graph)),
    node1_id_(id1),
    node2_id_(id2) {}

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  /** Return the number of edges in the graph.
   *
   * Complexity: O(1).
   */
  size_type edges_size() const {
    return n_edges;
  }

  /** Synonym for size(). */
  size_type num_edges() const {
    return edges_size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return *std::next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
 
 bool has_edge(const Node& a, const Node& b) const {
    //search on node with smallest degree
    if(a.degree() < b.degree()) {
      return (find_adj_node(b.uid_, a.uid_) != a.degree());
    }

    else
      return (find_adj_node(a.uid_, b.uid_) != b.degree());
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

    Edge e(this, a, b);
    if(has_edge(a,b)){
      return e;
    }

    else{
      // add each node in the adjacency list of the other one
      // twin_idx is the degree before insertion
      adj_.at(a.uid_).push_back(EdgeElement(b.uid_, b.degree()));
      adj_.at(b.uid_).push_back(EdgeElement(a.uid_, a.degree() - 1));

      ++n_edges;

      return e;
    }
  }


  /** Remove the element at position @a i, moving later elements down. */
  // remove_edge to be implemented


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    n_nodes = 0;
    n_edges = 0;
    nodes_.clear();
    idx2uid_.clear();
    adj_.clear();
  }

  



  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy


    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    NodeIterator& operator++() {
      node_idx_ += 1;
      return(*this);
    }

    bool operator==(const NodeIterator node_iter) const {
      return( node_iter.node_idx_ == this->node_idx_ &&
              node_iter.graph_ == this->graph_);
    }

    Node operator*() const {
      size_type node_uid = (graph_->nodes_).at(node_idx_).uid_;
      return Node(graph_, node_uid);
    }



   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type node_idx_;

    //Private constructor
    NodeIterator(const Graph* graph, size_type node_idx) :
      graph_(const_cast<Graph*>(graph)),
      node_idx_(node_idx) {}

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  NodeIterator node_begin() const {
      return NodeIterator(this, 0);
    }

  NodeIterator node_end() const {
      return NodeIterator(this, num_nodes());
  }

  


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const


    Edge operator*() const {
      size_type idx2 = (graph_->adj_).at(node_id_).at(idx_).adj_uid_;
      return Edge(graph_, node_id_, idx2);
    }

    IncidentIterator& operator++() {
      ++idx_;
      return *this;
    }

    bool operator==(const IncidentIterator& ii) const {
      return((graph_ == ii.graph_) && (node_id_ == ii.node_id_) && (idx_ == ii.idx_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    /* A pointer back to the graph. */
    Graph* graph_;
    /* The uid of the node whose edges this iterates. */
    size_type node_id_;
    /* Index of adjacent node. */
    size_type idx_;

    IncidentIterator(const Node& n, size_type i) :
      graph_(const_cast<Graph*>(n.graph_)),
      node_id_(n.uid_),
      idx_(i) {}

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    Edge operator*() const {
      size_type idx2 = (graph_->adj_).at(edge_idx1_).at(edge_idx2_).adj_uid_;
      return Edge(graph_, edge_idx1_, idx2);
    }

    EdgeIterator& operator++() {
      ++ edge_idx2_;
      fix_indices();
      return (*this);
    }

    bool operator==(const EdgeIterator& ei) const {
      return ((edge_idx1_ == ei.edge_idx1_) && (edge_idx2_== ei.edge_idx1_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    /* A pointer to the graph whose edges it iterates through. */
    Graph* graph_;
    /* First index */
    size_type edge_idx1_;
    /* Second index */
    size_type edge_idx2_;

    EdgeIterator(const Graph* graph, size_type edge_idx1, size_type edge_idx2) :
    graph_(const_cast<Graph*>(graph)),
    edge_idx1_(edge_idx1),
    edge_idx2_(edge_idx2)
    {
      fix_indices();
    }

    void fix_indices() {
      if(edge_idx1_ >= graph_->next_uid()){
        // then it should be edge_end()
        edge_idx1_ = graph_->next_uid();
        edge_idx2_ = 0;
      }
      else if(edge_idx2_ >= (graph_->adj_).at(edge_idx1_).size()){
        ++ edge_idx1_;
        edge_idx2_ = 0;
        fix_indices();
      }
      else if(edge_idx1_ >= (graph_->adj_).at(edge_idx1_).at(edge_idx2_).adj_uid_){
        ++ edge_idx2_;
        fix_indices();
      }
      // else everything is fine
    }

  };


  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  EdgeIterator edge_begin() const {
    EdgeIterator ei(this, 0, 0);
    return ei;
  }

  EdgeIterator edge_end() const {
    return EdgeIterator(this, next_uid(), 0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.


  // Internal type for set elements
  struct NodeElement {
    Point position_;          // The position held by a node
    node_value_type value_;
    size_type uid_;           // The unique identifcation for a node

    NodeElement(Point position, node_value_type value, size_type uid) :
      position_(position),
      value_(value),
      uid_(uid) {}
  };


  // Internal type for set elements
  struct EdgeElement {
    size_type adj_uid_;       // The first node held by an edge
    size_type twin_uid_;      // The second node held by an edge

    EdgeElement(size_type adj_uid, size_type twin_uid) :
      adj_uid_(adj_uid),
      twin_uid_(twin_uid) {}
  };

  std::vector<NodeElement> nodes_;
  std::vector<std::vector<EdgeElement>> adj_;
  std::vector<size_type> idx2uid_;

  size_type n_nodes;
  size_type n_edges;

  size_type next_uid() const {
    return nodes_.size();
  }

  Point& node_position(size_type id) {
    return nodes_.at(id).position_;
  }

  node_value_type& node_value(size_type id) {
      return nodes_.at(id).value_;
  }


  size_type find_adj_node(size_type uid, size_type from_uid) const {
    size_type idx = 0;
    size_type d = adj_.at(from_uid).size();
    while(idx < d && adj_.at(from_uid).at(idx).adj_uid_ != uid) {
      ++idx;
    }
    return idx;
  }
  


  // Disable copy and assignment of a Graph
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

};



#endif // CME212_GRAPH_HPP