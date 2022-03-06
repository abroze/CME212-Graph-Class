/**
 * @file GraphSymmetricMatrix.hpp
 * Implimentation file for treating the Graph as a MTL Matrix
 */

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

#include "CME212/Color.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"

#include <fstream>

#include "Graph.hpp"


using GraphType = Graph<char,char>;  //<  DUMMY Placeholder
using NodeType  = typename GraphType::node_type;


/** Remove all the nodes in graph @a g whose posiiton is within Box3D @a bb.
 * @param[in,out] g  The Graph to remove nodes from
 * @param[in]    bb  The BoundingBox, all nodes inside this box will be removed
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const Box3D& bb) {
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    NodeType n = *it;
    if (bb.contains(n.position())) {
      g.remove_node(n);
    }
  }
  return;
}

// check whether node is in boundary
bool boundary(const NodeType& n) {
  const Point p = n.position();
  
  if (norm_inf(p) == 1)
    return true;

  if (norm_inf(p - Point(0.6, 0.6, 0)) < 0.2 or
      norm_inf(p - Point(-0.6, 0.6, 0)) < 0.2 or
      norm_inf(p - Point(0.6, -0.6, 0)) < 0.2 or
      norm_inf(p - Point(-0.6, -0.6, 0)) < 0.2) {

    return true;
  } 

  Box3D box(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1));
  if(box.contains(p))
    return true;
  
  return false;
}


// Color functor
struct NodeColor {
  NodeColor(){}

  CME212::Color operator()(GraphType::node_type n) {
    double x = n.position().z;
    double c = std::abs(x) / norm(n.position());
    return CME212::Color::make_heat(c);
  }
};


// Position functor
struct NodePosition
{
  NodePosition(mtl::vec::dense_vector<double> x) : x_(x){}

  Point& operator()(NodeType& n)
  {
    Point& p = n.position();
    p.z = x_[n.index()];
    return p;
  }
  mtl::vec::dense_vector<double> x_;
};


// function f(x)
double f(Point& p) {return 5 * cos(norm_1(p));}

// function g(x) 
double g(Point& p){
  if (norm_inf(p) == 1)
    return 0.0;

  if (norm_inf(p - Point(0.6, 0.6, 0)) < 0.2 or
      norm_inf(p - Point(-0.6, -0.6, 0)) < 0.2 or
      norm_inf(p - Point(0.6, -0.6, 0)) < 0.2 or
      norm_inf(p - Point(-0.6, 0.6, 0)) < 0.2 ) {
      
      return -0.2;
  }

  Box3D box(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1));

  if (box.contains(p))
    return 1.0;

  return 1.0;

}


/* GraphSymmetricMatrix which is implemented by graph. */
class GraphSymmetricMatrix {

public:
  GraphSymmetricMatrix(GraphType* graph) : graph_(graph) {}

  // Elements of the symmetric matrix
  double element(size_t i, size_t j) {
    if (i == j && boundary(graph_->node(i)))
      return 1.0;

    if (i != j && (boundary(graph_->node(i)) or boundary(graph_->node(j))))
      return 0;

    else {
      if (i == j)
        return -double(graph_->node(i).degree());

      if (graph_->has_edge(graph_->node(i), graph_->node(j)))
        return 1.0;

      else
        return 0.0;
    }
  }

  /** Helper function to perform multiplication . Allows for delayed
    * evaluation of results .
    * Assign :: apply (a , b ) resolves to an assignment operation such as
    * a += b , a -= b , or a = b .
    * @pre @a size ( v ) == size ( w ) */

  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult(const VectorIn& v, VectorOut& w, Assign) const{
    size_t size = graph_->size();

    for (size_t i = 0; i < size; ++i) {
      double count = 0;

      for (size_t j = pt_idx_vec_[i]; j < pt_idx_vec_[i+1]; ++j) {
        count +=  elements_[j] * v[idx_vec_[j]];
      }

      Assign::apply(w[i], count);
    }
  }

  /** Matvec forwards to MTL’s lazy mat_cvec_multiplier operator */
    template <typename Vector>
    mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>
    operator *(const Vector& v) const{
        return {*this, v};
    }

  // Construct a sparse matrix using the graph's feature
  void make_sparse() {
    pt_idx_vec_.push_back(0);
    for (size_t i = 0; i < num_rows(); ++i) {
      for (size_t j = 0; j < num_cols(); ++j) {
        if (this->element(i, j) != 0) {
          elements_.push_back(this->element(i, j));
          idx_vec_.push_back(j);
        }
      }
      pt_idx_vec_.push_back(idx_vec_.size());
    }
  }


  size_t num_rows() const{
    return graph_->size();
  }

  size_t num_cols() const{
    return graph_->size();
  }

  size_t m_size() const{
    return graph_->size() * graph_->size();
  }

private:
  GraphType* graph_;
  std::vector<double> elements_;
  std::vector<size_t> pt_idx_vec_;
  std::vector<size_t> idx_vec_;
};


/* * The number of elements in the matrix . */
inline std::size_t size(const GraphSymmetricMatrix& A){
  return A.m_size();
}

/* * The number of rows in the matrix . */
inline std::size_t num_rows(const GraphSymmetricMatrix& A){
  return A.num_rows();
}

/* * The number of columns in the matrix . */
inline std::size_t num_cols(const GraphSymmetricMatrix& A){
  return A.num_cols();
}








