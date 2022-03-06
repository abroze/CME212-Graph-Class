/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

#include "IdentityMatrix.hpp"

/* * The number of elements in the matrix . */
inline std::size_t size(const IdentityMatrix& A){return A.m_size();}

/* * The number of rows in the matrix . */
inline std::size_t num_rows(const IdentityMatrix& A){return A.num_rows();}

/* * The number of columns in the matrix . */
inline std::size_t num_cols(const IdentityMatrix& A){return A.num_cols();}

/* * Traits that MTL uses to determine properties of our IdentityMatrix . */
namespace mtl {
namespace ashape {

/* * Define IdentityMatrix to be a non - scalar type . */
template <>
struct ashape_aux <IdentityMatrix> {
    typedef nonscal type ;
};
} // end namespace ashape

/* * IdentityMatrix implements the Collection concept
    * with value_type and size_type */
template < >
struct Collection <IdentityMatrix> {
    typedef double value_type ;
    typedef unsigned size_type ;
};
} // end namespace mtl




int main()
{
  typedef mtl::vec::dense_vector<double> Vector;
  typedef IdentityMatrix IdentityMatrix;

  const size_t size = 40;
  const size_t N = size * size;

  IdentityMatrix IM(N);

  //Create and Identity Preconditioner (there is no conditioning)
  itl::pc::identity<IdentityMatrix, double> P(IM);

  //Set b such that x == 1 is the solution; start with x == 0
  Vector x(N, 1.0), b(N);
  b = IM * x;
  x = 0.0;

  //Termination criterion: r < 1e-6 * b or N iterations
  itl::noisy_iteration<double> iter(b, 500, 1.e-6);

  //Solve Ax == b with left preconditioner P
  itl::bicgstab(IM, x, b, P, iter);

  return 0;
}

