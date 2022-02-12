/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include <fstream>
#include <chrono>
#include <thread>

#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
    Point vel;                       // Node velocity
    double mass;                   // Node mass
    NodeData() : vel(0), mass(1) {}
};


/** Custom structure of data to store with Edges */
struct EdgeData {
    double K;     // Edge spring constant
    double L;     // Edge spring L
    EdgeData() : K(100), L(0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


struct Force {
    Force() {}

    virtual Point operator()(Node n, double t) = 0;
    virtual ~Force(){}
};


struct GravityForce: public Force {

    virtual Point operator()(Node n, double t) {
      (void) t;
      return n.value().mass * Point(0, 0, -grav);
    }
};


struct MassSpringForce: public Force {
    MassSpringForce(){}

    virtual Point operator()(Node n, double t) {

      (void) t;
      Point F_spring = Point(0, 0, 0);
      
      for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
        Edge e = *ii;
        Point xi = e.node1().position();
        Point xj = e.node2().position();
        
        double edge_length = e.length();
        double K = e.value().K;
        double L = e.value().L;

        F_spring += -K * (xi - xj) * (edge_length - L)/edge_length;
      }
      return F_spring;
    }
    virtual ~MassSpringForce(){}

};


struct DampingForce: public Force {
    DampingForce(): damp_const_(0) {}
    DampingForce(double c): damp_const_(c) {}

    virtual Point operator()(Node n, double t) {
      (void) t;
      Point F_damp = n.value().vel * -damp_const_;
      return F_damp;
    }
    virtual ~DampingForce(){}

    private:
      double damp_const_;
};


/** Force function object for HW2 #1. */

/** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force.
   */
struct Problem1Force {

    template <typename NODE>
    Point operator()(NODE n, double t) {
      (void) t; 

      if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
        return Point(0, 0, 0);
      }

      MassSpringForce msf = MassSpringForce();
      GravityForce gf = GravityForce();

      Point F_spring = msf(n, t);
      Point F_grav = gf(n, t);

      return (F_spring + F_grav);
    }
};

struct CombinedForces {

    CombinedForces(std::vector<Force*> force_vec) :
    force_vec_(force_vec) {}

    Point operator()(Node n, double t) {

        Point p = Point(0, 0, 0);

        for (auto it = force_vec_.begin(); it != force_vec_.end(); ++it){
            p += (*(*it))(n,t);
        }
        return p;
    }

    private:
    std::vector<Force*> force_vec_;

};


template <typename F1, typename F2>
CombinedForces make_combined_force(F1& f1, F2& f2){
    std::vector<Force*> v = {&f1, &f2};
    return CombinedForces(v);
}

template <typename F1, typename F2, typename F3>
CombinedForces make_combined_force(F1& f1, F2& f2, F3& f3){
    std::vector<Force*> v = {&f1, &f2, &f3};
    return CombinedForces(v);
}


struct Constraint {
    Constraint() {}
    virtual void operator()(GraphType& g, double t) = 0;
    virtual ~Constraint(){}
};


struct PinConstraint: public Constraint {

    PinConstraint(std::vector<Node>& fixed_nodes, std::vector<Point>& fixed_positions) :
    fixed_nodes_(fixed_nodes),
    fixed_positions_(fixed_positions) {}

    virtual void operator()(GraphType& g, double t){
      (void) g;
      (void) t;
      for(Point::size_type k = 0; k < fixed_positions_.size(); k++){
        Node fix_node = fixed_nodes_[k];
        fix_node.position() = fixed_positions_[k];
        fix_node.value().vel = Point(0, 0, 0);
        }
    };

    virtual ~PinConstraint(){}

    private:
    std::vector<Node> fixed_nodes_;
    std::vector<Point> fixed_positions_;
};


struct PlaneConstraint : public Constraint {
    PlaneConstraint(double thresh) :
    threshold_(thresh){}

    virtual void operator()(GraphType& g, double t){
      (void) t;
      for(auto i = g.node_begin(); i != g.node_end(); ++i){
          Node n = *i;
          if (dot(n.position(), Point(0, 0, 1)) < threshold_){
              n.position().z = threshold_;
              n.value().vel.z = 0;
          }
      }

    }
    private:
    double threshold_;
};


struct SphereConstraint : public Constraint {
    SphereConstraint(Point c, double r): 
    center_(c), radius_(r){}

    virtual void operator()(GraphType& g, double t){
      (void) t;
      for (auto i = g.node_begin(); i != g.node_end(); ++i){
        Node n = *i;

        Point vec_r = n.position() - center_;

        if (norm(vec_r) < radius_) {
            Point vec_unit_radius = vec_r/norm(vec_r);
            n.position() = vec_unit_radius * radius_ + center_;
            n.value().vel -= dot(n.value().vel, vec_unit_radius) * vec_unit_radius;
        }
      }
    }

    private:
    Point center_;
    double radius_;
};

struct RemoveSphereConstraint : public Constraint {
    RemoveSphereConstraint(const Point center, double r) :
    center_(center), radius_(r) {}

    virtual void operator()(GraphType& g, double t){
        (void) t;

        for (auto ni = g.node_begin(); ni!= g.node_end(); ++ni){
            auto n = *ni;
            double dist = norm(n.position() - center_);
            if (dist < radius_){
                ni = g.remove_node(ni);
            }
        }
    }
    virtual ~RemoveSphereConstraint(){};

    private:
    const Point center_;
    const double radius_;
};


struct CombinedConstraints {
    CombinedConstraints(std::vector<Constraint*> constr_vec):
    constr_vec_(constr_vec){}

    void operator()(GraphType& g, double t){
      for(unsigned int i = 0; i < constr_vec_.size(); ++i){
          (*constr_vec_[i])(g, t);
      }
    }
    private:
      std::vector<Constraint*> constr_vec_;

};

template <typename C1, typename C2>
CombinedConstraints make_combined_constraint(C1& c1, C2& c2){
    std::vector<Constraint*> constr_vec = {&c1, &c2};
    return CombinedConstraints(constr_vec);
}


template <typename C1, typename C2, typename C3>
CombinedConstraints make_combined_constraint(C1& c1, C2& c2, C3& c3){
    std::vector<Constraint*> constr_vec = {&c1, &c2, &c3};
    return CombinedConstraints(constr_vec);
}



/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */


template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {

  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}














