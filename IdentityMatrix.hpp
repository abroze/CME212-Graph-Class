/**
 * @file mtl_test.cpp
 * Implementation file for interfacing with MTL4 and it's linear solvers.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/itl/itl.hpp>
#include <boost/numeric/mtl/mtl.hpp>


// HW3: YOUR CODE HERE
// Define a IdentityMatrix that interfaces with MTL

class IdentityMatrix {

public:
    IdentityMatrix(size_t size) :
    size_(size) {}

    /** Helper function to perform multiplication . Allows for delayed
        * evaluation of results .
        * Assign :: apply (a , b ) resolves to an assignment operation such as
        * a += b , a -= b , or a = b .
        * @pre @a size ( v ) == size ( w ) */

    template <typename VectorIn, typename VectorOut, typename Assign>
    void mult(const VectorIn& v, VectorOut& w, Assign) const {
        for (size_t i = 0; i < size_; ++i) {Assign::apply(w[i], v[i]);}
    }

    /* * Matvec forwards to MTL’s lazy mat_cvec_multiplier operator */
    template <typename Vector>
    mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector>
    operator *(const Vector& v) const{
        return {*this, v};
    }

    size_t num_rows() const {return size_;}

    size_t num_cols() const {return size_;}

    size_t m_size() const {return size_ * size_;}

private:
    size_t size_;

};