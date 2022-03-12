/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include "CME212/SFML_Viewer.hpp"

#include <fstream>
#include <chrono>
#include <thread>

#include "thrust/for_each.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/system/omp/execution_policy.h"

#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"

#include "Graph.hpp"
#include "SpaceSearcher.hpp"

#include "mass_spring.hpp"


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct an empty graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  std::vector<typename GraphType::node_type> nodes;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(nodes[t[0]], nodes[t[1]]);
    graph.add_edge(nodes[t[0]], nodes[t[2]]);

    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);

    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  /*
   // Set initial conditions for your nodes, if necessary.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    (*it).value().mass = 1.0/double(graph.size());
  }
  */

  std::vector<Node> fixed_corners;
  std::vector<Point> fixed_positions;

  for (auto node_it = graph.node_begin(); node_it != graph.node_end(); ++node_it) {
    auto n = *node_it;
    n.value().mass = (double) 1/graph.size();
    n.value().vel = Point(0, 0, 0);

    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      fixed_corners.push_back(n);
      fixed_positions.push_back(n.position());
    }
  }


  // Parallel version to set initial conditions for all edges 
  thrust::for_each(thrust::omp::par,
                   //thrust::seq,
                   graph.node_begin(),
                   graph.node_end(),
                   ThrustNodeInitialization(graph.num_nodes(),
                   fixed_corners,
                   fixed_positions));

  
  for (auto eit = graph.edge_begin(); eit != graph.edge_end(); ++eit) {
    (*eit).value().K = 100.0;
    (*eit).value().L = (*eit).length();
  }
  

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // We want viewer interaction and the simulation at the same time
  // Viewer is thread-safe, so launch the simulation in a child thread
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&](){

      CME212::Clock clock;

      // Begin the mass-spring simulation
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;

      clock.start();

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        GravityForce F1 = GravityForce();
        MassSpringForce F2 = MassSpringForce();
        DampingForce F3 = DampingForce((double)1/graph.size());
        CombinedForces combined_force_fn = make_combined_force(F1, F2, F3);

        Point center = Point(0.5, 0.5, -0.5);
        double radius = 0.15;
        double plane_thresh = -0.75;

        PinConstraint C1 = PinConstraint(fixed_corners, fixed_positions);
        PlaneConstraint C2 = PlaneConstraint(plane_thresh);
        SphereConstraint C3 = SphereConstraint(center, radius);
        RemoveSphereConstraint C4 = RemoveSphereConstraint(center, radius);
        SelfCollisionConstraint C5 = SelfCollisionConstraint();
        CombinedConstraints combined_constraint_fn = make_combined_constraint(C1, C2, C5);

        // Apply all forces and constraints to to the graph
        symp_euler_step(graph, t, dt, combined_force_fn, combined_constraint_fn);

        // Clear the viewer ’s nodes and edges
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        if (graph.size() < 100)
          std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }

      double total_time = clock.seconds();
      std::cout << "Time: " << total_time << " s" << std::endl;

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}








