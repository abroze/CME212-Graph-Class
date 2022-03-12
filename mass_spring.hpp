/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include "Graph.hpp"
#include <fstream>
#include <chrono>
#include <thread>

#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include <thrust/for_each.h>
#include <thrust/execution_policy.h>
#include "thrust/iterator/transform_iterator.h"
#include "thrust/system/omp/execution_policy.h"


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




struct CutDownVelocity {

    CutDownVelocity(Node& n, double ra) :
    n1(n), radius(ra) {}

    void operator()(Node n2) {

        Point r = n1.position() - n2.position();
        double l2 = normSq(r);

        if (n1 != n2 && l2 < radius) {

            // Prevent hit by removing velocity component in r
            n1.value().vel -= (dot(r, n1.value().vel) / l2) * r;
        }
    }

    Node n1;
    double radius;
};



Box3D FindBoxIntersection(Box3D& smallbox, Box3D& bigbox) {

    // Obtain maximum of each box
    Point max_smallbox = smallbox.max();
    Point max_bigbox = bigbox.max();

    // Obtain minimum of each box
    Point min_smallbox = smallbox.min();
    Point min_bigbox = bigbox.min();

    // Declare points of resulting box and update them
    Point resulting_max = max_smallbox;
    Point resulting_min = min_smallbox;

    for (Point::size_type i = 0; i < min_smallbox.size(); ++i) {
        
        if (max_smallbox[i] > max_bigbox[i]) 
            resulting_max[i] = max_bigbox[i];

        if (min_smallbox[i] < min_bigbox[i])
            resulting_min[i] = min_bigbox[i];
    }

    // Comparison
    for (Point::size_type i = 0; i < min_smallbox.size(); ++i) {
        
        if (resulting_max[i] < min_bigbox[i])
            resulting_max[i] = min_bigbox[i];

        if (resulting_min[i] > max_bigbox[i]) 
            resulting_min[i] = max_bigbox[i];

    }
    return Box3D(resulting_min,resulting_max);
}


struct DetermineInfluence {

    DetermineInfluence(SpaceSearcher<Node>& ss): searcher_(ss){}

    void operator()(Node n) {

        const Point& center = n.position();
        double radius2 = std::numeric_limits<double>::max();

        // Determine squared radius
        for (auto eit = n.edge_begin(); eit != n.edge_end(); ++eit) {
            radius2 = std::min(radius2, normSq((*eit).node2().position() - center));
        }
        radius2 *= 0.9;

        // Construct small bounding box around node using scaled radius
        Point upper = center - sqrt(radius2);
        Point lower = center + sqrt(radius2);
        Box3D SmallBB(lower, upper);

        // Get the big bounding box from searcher
        Box3D BigBB = searcher_.bounding_box();

        // Get an intersection box contained by the big bounding box
        Box3D result_bb = FindBoxIntersection(SmallBB, BigBB);

        assert(searcher_.bounding_box().contains(result_bb));

        // Constraint is enforced within the bounding box in parallel
        thrust::for_each(searcher_.begin(result_bb), searcher_.end(result_bb), CutDownVelocity(n, radius2));
    }

    private:
    SpaceSearcher<Node>& searcher_;
};



struct SelfCollisionConstraint : public Constraint {

    void operator()(GraphType& g, double t) {

        auto n2p = [](const Node& n) {return n.position();};
        Box3D bigbb = Box3D(thrust::make_transform_iterator(g.node_begin(), n2p),
                            thrust::make_transform_iterator(g.node_end(), n2p));

        Point extended_lower = bigbb.min();
        Point extended_upper = bigbb.max();

        for (unsigned i = 0; i < extended_lower.size(); ++i) {
            extended_lower[i] = -abs(extended_lower[i])*1.5;
            extended_upper[i] = abs(extended_upper[i])*1.5;
        }

        bigbb = Box3D(extended_lower, extended_upper);

        SpaceSearcher<Node> searcher(bigbb, g.node_begin(), g.node_end(), n2p);

        thrust::for_each(g.node_begin(), g.node_end(), DetermineInfluence(searcher));
    }
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



struct ThrustUpdateNodePos  {

    ThrustUpdateNodePos(double difftime) :
    dt_(difftime) {}

    void operator()(Node n) {
        n.position() += n.value().vel * dt_;
    }

    double dt_;
};


template <typename F>
struct ThrustUpdateNodeVel  {

    ThrustUpdateNodeVel(double time, double diff_time, F& force) :
    f_(force), t_(time), dt_(diff_time) {}

    void operator()(Node n) {
        n.value().vel += f_(n, t_) * (dt_ / n.value().mass);
    }

    F f_; 
    double t_; 
    double dt_; 
};



template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {

//  // Compute the t+dt position
//  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
//    auto n = *it;
//
//    // Update the position of the node according to its velocity
//    // x^{n+1} = x^{n} + v^{n} * dt
//    n.position() += n.value().vel * dt;
//  }

  thrust::for_each(thrust::omp::par, g.node_begin(), g.node_end(), ThrustUpdateNodePos(dt));

  constraint(g, t);

  // Compute the t+dt velocity
  //for (auto it = g.node_begin(); it != g.node_end(); ++it) {
  //  auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
  //  n.value().vel += force(n, t) * (dt / n.value().mass);
  //}

  thrust::for_each(thrust::omp::par, g.node_begin(), g.node_end(), ThrustUpdateNodeVel<F>(t,dt,force));

  return t + dt;
}


struct ThrustNodeInitialization {

    ThrustNodeInitialization(Point::size_type num_nodes,
                             std::vector<Node>& fixed_nodes,
                             std::vector<Point>& fixed_positions) :
        num_nodes_(num_nodes),
        fixed_nodes_(fixed_nodes),
        fixed_positions_(fixed_positions){}

    void operator()(Node n) {

        // Setting initial conditions for nodes
        n.value().mass = 1.0/num_nodes_; 
        n.value().vel = Point(0, 0, 0);

        // Fixing points on corners
        if (n.position() == Point(0, 0, 0) or n.position() == Point(1, 0, 0)) {
            fixed_nodes_.push_back(n);
            fixed_positions_.push_back(n.position());
        }
    }

    Point::size_type num_nodes_; 
    std::vector<Node>& fixed_nodes_; 
    std::vector<Point>& fixed_positions_;
};












