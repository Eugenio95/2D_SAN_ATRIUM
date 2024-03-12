#ifndef MAT2D_HPP
#define MAT2D_HPP

#include <algorithm>
#include <cassert>

// A simple implementation of a 2D matrix represented as a C-style
// array of rows.
template <typename T>
class mat2D {
private:
    T *mat;
    const int r, c;
public:
    // create a matrix with r rows and c cols
    mat2D(size_t rows, size_t cols) : r(rows), c(cols)
    {
        assert(r > 0);
        assert(c > 0);
        mat = new T[r * c];
        assert( mat != 0 );
    }

    // create a matrix with r rows and c cols, and fill the matrix
    // with the initial value val
    mat2D(size_t rows, size_t cols, T val) : mat2D(rows, cols)
    {
        std::fill(mat, mat + r*c, val);
    }

    virtual ~mat2D( )
    {
        delete mat;
    }

    // get a pointer to the internal representation of the matrix
    T *data( void )
    {
        return mat;
    }

    // get a pointer to the beginning of row r
    T *operator[]( size_t row )
    {
        assert(row >= 0 && row < rows);
        return mat + row*c;
    }

    // get the number of rows
    size_t rows( void ) const { return r; }

    // get the number of columns
    size_t cols( void ) const { return c; }
};

#endif
