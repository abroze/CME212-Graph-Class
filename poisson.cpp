/**
 * @file poisson.cpp
 * Test script for treating for using the GraphSymmetricMatrix class
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SFML_Viewer to visualize the solution.
 */

#include "CME212/SFML_Viewer.hpp"
#include "GraphSymmetricMatrix.hpp"

#include "Graph.hpp"


/* * Traits that MTL uses to determine properties of our GraphSymmetricMatrix . */
namespace mtl {
  namespace ashape {

    /* * Define GraphSymmetricMatrix to be a non - scalar type . */
    template <>
    struct ashape_aux <GraphSymmetricMatrix> {
        typedef nonscal type;
    };
  } // end namespace ashape

  /* * GraphSymmetricMatrix implements the Collection concept
      * with value_type and size_type */
  template < >
  struct Collection <GraphSymmetricMatrix> {
    typedef double value_type;
    typedef unsigned size_type;
  };
} // end namespace mtl


namespace itl {
  template <class Real>

  // Class to visualize the iterations of the solution
  class visual_iteration : public cyclic_iteration<Real> {

    typedef cyclic_iteration<Real> super;
    typedef visual_iteration self;

    void update_viewer() {
      viewer_->add_nodes(graph_->node_begin(), graph_->node_end(), nc_, NodePosition(*x_), node_map_);
      viewer_->add_edges(graph_->edge_begin(), graph_->edge_end(), node_map_);
      viewer_->set_label(this->i);
    }

  public:

    visual_iteration( mtl::vec::dense_vector<double>& r0,
                      mtl::vec::dense_vector<double>* x,
                      std::map<NodeType, unsigned> node_map,
                      int max_iterations,
                      Real tolerance,
                      CME212::SFML_Viewer* viewer,
                      GraphType* graph,
                      NodeColor nc,
                      int cycle = 10,
                      Real atolerance = Real(0)) :
    
    super(r0, max_iterations, tolerance, atolerance, cycle),
    viewer_(viewer),
    graph_(graph),
    nc_(nc),
    x_(x),
    node_map_(node_map)
    {update_viewer();}

    bool finished() {
      update_viewer();
      return super::finished();
    }

    template <typename T>
    bool finished(const T& r) {
       bool ret = super::finished(r);
       update_viewer();

       // Add delay for better visualization
       std::this_thread::sleep_for(std::chrono::milliseconds(2000));
       return ret;
    }

  private:
    CME212::SFML_Viewer* viewer_;
    GraphType* graph_;
    NodeColor nc_;
    mtl::vec::dense_vector<double>* x_;
    std::map<NodeType, unsigned> node_map_;
  };
}


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

  {
    // Create a nodes_file from the first input argument
    std::ifstream nodes_file(argv[1]);
    // Interpret each line of the nodes_file as a 3D Point and add to the Graph
    std::vector<NodeType> node_vec;
    Point p;
    while (CME212::getline_parsed(nodes_file, p))
      node_vec.push_back(graph.add_node(2*p - Point(1,1,0)));

    // Create a tets_file from the second input argument
    std::ifstream tets_file(argv[2]);
    // Interpret each line of the tets_file as four ints which refer to nodes
    std::array<int,4> t;
    while (CME212::getline_parsed(tets_file, t)) {
      graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
      graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
      graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
      graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
    }
  }

  // Get the edge length, should be the same for each edge
  auto it = graph.edge_begin();
  assert(it != graph.edge_end());
  double h = norm((*it).node1().position() - (*it).node2().position());

  // Make holes in our Graph
  remove_box(graph, Box3D(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));


  size_t num_nodes = graph.size();
  mtl::vec::dense_vector<double> b(num_nodes, 0.0);
  
  for (size_t i = 0; i < num_nodes; ++i) {
    auto n = graph.node(i);
    auto x = n.position();
    if (boundary(n)) {
      b[i] = g(x);
    }
    else {
      double gx = h * h * f(x);
      for (auto ni = n.edge_begin(); ni != n.edge_end(); ++ni) {
        auto n2 = (*ni).node2();
        if (boundary(n2))
          gx -= g(n2.position());
      }
      b[i] = gx;
    }
  }

  GraphSymmetricMatrix A(&graph);
  A.make_sparse();

  mtl::vec::dense_vector<double> x(num_nodes, 0.0);

  // Launch Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();

  // Simulation thread to visualize  solution
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&]() {
    itl::visual_iteration<double> iter(b, &x, node_map, 1000, 1.e-11, &viewer, &graph, NodeColor());
    itl::cg(A, x, b, iter);
  });

  viewer.event_loop();

  // Interrupting window if returning from event loop
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
