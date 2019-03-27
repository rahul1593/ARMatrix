/*
 * ARMATRIX.H
 * This is Multi-dimensional matrix operations library.It contains the required functions 
 * for operating on multi-dimensional matrices. It also contains some other mathematical 
 * functions like convolutions, correlations, etc.
 * 
 * Author       : Rahul Bhartari
 * Date Created : Dec 1, 2018
 * License      : LGPL2.1
 */

#ifndef _MATRICS_H_
#define _MATRICS_H_

#include <stdarg.h>

#ifdef MODE_C
    #include <stdlib.h>
    #include <stdio.h>
    #include <stdint.h>
    #include <math.h>

    #define mtx_type    double
#else
    #include <iostream>
    #include <vector>
    #include <cmath>
    #include <cstdlib>
    #include <cstdint>

    //using namespace std;
#endif

// operation definitions
#define OP_ADD      1
#define OP_SUB      2
#define OP_MUL      3
#define OP_DIV      4
#define OP_EXP      5
#define OP_POW      6
#define OP_SIN      7
#define OP_COS      8
#define OP_SINH     9
#define OP_COSH     10
#define OP_LOG      11
#define OP_SQR      12
#define OP_SQRT     13
#define OP_FACT     14
#define OP_MOD      15
#define OP_REM      16
#define OP_RELU     17

// matrix data type definitions
#define AMTX_UI8    0   // unsigned and signed 8/16/32/64-bit types
#define AMTX_UI32   1
#define AMTX_UI64   2
#define AMTX_I8     3
#define AMTX_I32    4
#define AMTX_I64    5
#define AMTX_F32    6   // single and double precision floating types
#define AMTX_F64    7
#define AMTX_ANY    8   // user defined type, input as class object

typedef struct{
    int x, y;
} Point2D;

typedef struct{
    int x, y, z;
} Point3D;

#ifndef MODE_C


template <typename mtx_type> class NDMat{
private:
    int                     type;           // data type of the matrix
    std::vector<int>        shape;          // matrix dimensions, outermost to innermost
    std::vector<int>        increments;     // increment in data array for each element in respective dimension
    int                     dlen;           // length of data array
    std::vector<mtx_type>   data;           // vector of any specified data type, data is stored in row major vector, type depends on usage
/*
    bool                    is_submat;          // is a child submatrix or not
    bool                    has_submat;         // does it has subtrices or not
    NDMat<mtx_type>   *parent_mat=NULL;         // pointer to parent matrix of this matrix
    std::vector<NDMat<mtx_type>*> submat_list;  // list of pointers to submatrices
    std::vector<Point2D>    join_indx;          // set of indexes to join disjoint data, used in referencing sliced matrices
*/
    bool                    is_sparse;          // is a sparse matrix with a lot of zeros, eg, identity matrix
    std::vector<std::vector<int>> sp_data;      // index data for sparse matrix
    // Sparse optimization support will be added later on
protected:
    // check feasibility for his matrix to be converted and stored as a sparse matrix
    //bool __checkSparseFeasibility();
    // convert and store this matrix as a sparse matrix
    //bool __converToSparse();

    // check dimension compatibility for broadcasting and matrix dot product
    std::vector<int> __checkCompatibleBroadcastDim(NDMat<mtx_type> A, NDMat<mtx_type> B);    // if yes return normalized dimension for matrix with less number of dimensions
    bool __checkCompatibleDotProductDim(NDMat<mtx_type> A, NDMat<mtx_type> B);
    static void __recursivePaddedVectorCopy(mtx_type* src_data, int* src_shape, int* src_inc, mtx_type* tgt_data, int* tgt_shape, int* tgt_inc, int dim_cnt);
public:
    /*  Constructor with shape as a vector, default type is AMTX_F32, i.e, float32.
        Arguments:
            shape       : Vector defining the dimensionality of the n-D matrix, from higher to lower dimension
            padding     : (optional) The size of extra padding (of zeros) on each side for each of the dimensions
            type        : (optional) Data type as defined in this header
            init_value  : (optional) Pointer to value of given (or default) type, which is used to initialise the matrix
            type_object : (optional) Custom type of object to be stored in matrix
        
        Return: Object of class NDMat having given (or default) attributes.
    */
    NDMat(const NDMat<mtx_type> &src_mat);
    NDMat(std::vector<int> shape, mtx_type init_value = 0);
    
    // destructor for NDMat object
    ~NDMat();

    // generate an identity matrix of given shape and data type
    static NDMat<mtx_type> identity(int side);

    // accessor and mutator methods
    // get shape of this matrix
    std::vector<int> getShape();
    
    // check if this matrix is sparse
    bool isSparse();
    
    // get total number of elements in this matrix
    int len();

    // Reshape the matrix
    void reshape(std::vector<int> shape);

    // change the data type of the matrix
    NDMat<mtx_type> asType(int type);

    // print this matrix on console
    void print();

    // create a duplicate matrix by copying the current matrix, with and without padding
    NDMat<mtx_type> copy();
    NDMat<mtx_type> copy(std::vector<int> padding);

    // create zero initialized matrix similar to current matrix
    static NDMat<mtx_type> zeros_like(NDMat<mtx_type> src_mat);

    /*  Create a submatrix by slicing the current matrix
        Arguments:
            axis_min_max : Vector containing ranges for each dimension to be sliced
        Return: Object of class NDMat representing the submatrix
        Variations:
            slice   : This function creates a completely new object and copies data from source matrix
            vslice  : This function creates an object in which data points to source matrix and uses 'joinIndex' variable
                      for stiching the data together.
    */
    NDMat<mtx_type> slice(std::vector<Point2D> axis_min_max);  // duplicate data in new matrix
    NDMat<mtx_type> vslice(std::vector<Point2D> axis_min_max); // slice but refer to original data

    /*  Assign values for a slice of this matrix
        Arguments:
            axis_min_max : Vector containing ranges for each dimension to be sliced
    */
   void assignSlice(std::vector<Point2D> axis_min_max, NDMat<mtx_type> set_val);

    /*  Mask the matrix and return the masked copy
        Arguments:
            condition   : masking condition
            axis        : axis on which masking needs to be done, default to '-1' for all values
            val         : pointer to the value to be set where 'condition' is true
    */
    NDMat<mtx_type> maskDup(NDMat<bool> condition, mtx_type val, int axis = -1);
    void mask(NDMat<bool> condition, mtx_type val, int axis = -1);

    /*  Perform the operation defined by user on each element, axis or dimension of given matrix
        Arguments:
            src_mat     : source matrix to be operated on
            op_function : function to operate on data
            axis        : axis on which operation needs to be done, default to '-1' for all values
    */
    NDMat<mtx_type> operateDup(mtx_type *(func)(mtx_type *), int axis = -1);
    void operate(mtx_type *(func)(mtx_type *), int axis = -1);

    // reverse a matrix around the given axis (default to '-1' for all values)
    NDMat<mtx_type> reverseDup(int axis = -1);    // returns a duplicate matrix
    void reverse(int axis = -1);        // in place reverse

    // Transpose of this matrix with copied data or in place
    NDMat<mtx_type> transposeDup();
    void transpose();

    // Inverse of this matrix
    NDMat<mtx_type> inverseDup();
    void inverse();

    // dot product of two matrices
    NDMat<mtx_type> dot(NDMat<mtx_type> B);

    // add two matrices
    NDMat<mtx_type> add(NDMat<mtx_type> B);

    // subtract matrix B from A
    NDMat<mtx_type> sub(NDMat<mtx_type> B);

    // elementwise product
    NDMat<mtx_type> mul(NDMat<mtx_type> B);

    // elementwise division
    NDMat<mtx_type> div(NDMat<mtx_type> B);

    // scale a matrix, i.e, multiply by a constant pointed by 'c'
    NDMat<mtx_type> scaleDup(mtx_type c);
    void scale(mtx_type c);
    
    // Overloaded arithmetic operators for easy arithmetic
    // assign a matrix to another
    NDMat<mtx_type> &operator=(NDMat<mtx_type> &m2);
    NDMat<mtx_type> &operator=(mtx_type val);

    // add two matrices
    NDMat<mtx_type> operator+(NDMat<mtx_type> &m2);
    NDMat<mtx_type> operator+(mtx_type val);
    NDMat<mtx_type> &operator+=(NDMat<mtx_type> &m2);
    NDMat<mtx_type> &operator+=(mtx_type val);

    // subtract one matrix from another matrix
    NDMat<mtx_type> operator-(NDMat<mtx_type> &m2);
    NDMat<mtx_type> operator-(mtx_type val);
    NDMat<mtx_type> &operator-=(NDMat<mtx_type> &m2);
    NDMat<mtx_type> &operator-=(mtx_type val);

    // elementwise multiply two matrices
    NDMat<mtx_type> operator*(NDMat<mtx_type> &m2);
    NDMat<mtx_type> operator*(mtx_type val);
    NDMat<mtx_type> &operator*=(NDMat<mtx_type> &m2);
    NDMat<mtx_type> &operator*=(mtx_type val);
    
    // elementwise divide one matrix by another
    NDMat<mtx_type> operator/(NDMat<mtx_type> &m2);
    NDMat<mtx_type> operator/(mtx_type val);
    NDMat<mtx_type> &operator/=(NDMat<mtx_type> &m2);
    NDMat<mtx_type> &operator/=(mtx_type val);

    // elementwise mod one matrix by another
    NDMat<mtx_type> operator%(NDMat<mtx_type> &m2);
    NDMat<mtx_type> operator%(mtx_type val);
    NDMat<mtx_type> &operator%=(NDMat<mtx_type> &m2);
    NDMat<mtx_type> &operator%=(mtx_type val);

    // Overloaded comparison operators
    friend NDMat<bool> operator<(NDMat<mtx_type> &m1, NDMat<mtx_type> &m2);
    friend NDMat<bool> operator>(NDMat<mtx_type> &m1, NDMat<mtx_type> &m2);
    friend NDMat<bool> operator>=(NDMat<mtx_type> &m1, NDMat<mtx_type> &m2);
    friend NDMat<bool> operator<=(NDMat<mtx_type> &m1, NDMat<mtx_type> &m2);
    friend NDMat<bool> operator!=(NDMat<mtx_type> &m1, NDMat<mtx_type> &m2);
    friend NDMat<bool> operator==(NDMat<mtx_type> &m1, NDMat<mtx_type> &m2);

    friend NDMat<bool> operator<(NDMat<mtx_type> &m1, mtx_type val);
    friend NDMat<bool> operator>(NDMat<mtx_type> &m1, mtx_type val);
    friend NDMat<bool> operator>=(NDMat<mtx_type> &m1, mtx_type val);
    friend NDMat<bool> operator<=(NDMat<mtx_type> &m1, mtx_type val);
    friend NDMat<bool> operator!=(NDMat<mtx_type> &m1, mtx_type val);
    friend NDMat<bool> operator==(NDMat<mtx_type> &m1, mtx_type val);


    // Overloaded stream operator
    friend NDMat<mtx_type> &operator << (std::ostream &out, NDMat<mtx_type> &src_mat);

    // Overloaded function operator to get any value
    // Takes dimension index as arguments, returns the value or submatrix address
    NDMat<mtx_type> &operator()(...);

    // operator overloading for masking
    NDMat<mtx_type> operator[](NDMat<bool> condition);
};

#else
/*
 * Data structure to store the matrix
 */
typedef struct{
    int* shape;         // matrix dimensions, outermost to innermost
    int* increments;    // store length of each dimension
    int dim_cnt;        // number of dimensions
    int dlen;           // length of data array
    mtx_type* data;     // data can be of stated type, integer by default, data is stored in row major array
    int is_sparse;		// is a sparse matrix with a lot of zeroes, eg., identity matrix
    int** coordinates;  // store the coordinates for values in sparse matrix
} Matrix;


/*
 * Data structure to store tree whMtxVectorich has data on leaf nodes for efficient multidimensional add/sub/mul/dot operation
 */
typedef union {
    struct MtxNode* children;   // all the child nodes
    mtx_type* data;             // contains single data value only in leaf node
    uint32_t dlen, vlen;        // vector(child nodes) length
} MtxNode;

/*
 * -------- Matrix Operation functions --------
 * User must take care of allocating or freeing space for all the matrices.
 */

//* pretty-print the matrix
void mtx_print(Matrix* A);

//* create a matrix initialized with given value
Matrix* mtx_create(int* shape, int dim_cnt, mtx_type init_val);

//* create a matrix similar to the given matrix
Matrix* mtx_create_as(Matrix M);

//* destroy the input matrix and free the space allocated
void mtx_destroy(Matrix** A);

//* reshape the matrix in place
int mtx_reshape(Matrix*A, int* shape, int dim_cnt);

//* create copy of a matrix
void mtx_copy(Matrix* src, Matrix* tgt);

// matrix addition, C = A+B. The addition is done according to the dimensionality of the matrices
void mtx_add(Matrix* A, Matrix* B, Matrix* C);

// matrix subtraction, C = A-B
void mtx_sub(Matrix* A, Matrix* B, Matrix* C);

// matrix multiplication, C = A*B
void mtx_mul(Matrix* A, Matrix* B, Matrix* C);

// matrix division, C = A/B
void mtx_div(Matrix* A, Matrix* B, Matrix* C);

//* scale a matrix, i.e, multiply by a constant 'c'
void mtx_scale(Matrix* A, mtx_type c);

// matrix dot product, C = A.B
void mtx_dot(Matrix* A, Matrix* B, Matrix* C);

// matrix transpose
void mtx_transpose(Matrix* A, Matrix *T);

// matrix invert
void mtx_invert(Matrix* A, Matrix* Ai);

// identity matrix of dimensions as specified by 'shape' (array of length 'dim')
void mtx_identity(int* shape, int dim, Matrix* I);

/* 
 * -- Planned to be implemented later ......

// convolution 2d, padding='s'(same) or 'v'(valid)
void mtx_conv2d(Matrix* A, Matrrix* F, int* strides, char padding, Matrix* C);

// correlation 2d
void mtx_corr2d(Matrix* A, Matrrix* F, Matrix* C);

// 2d max pooling
void mtx_maxpool2d(Matrix* A, int* kernel_dim, int* strides, char padding, Matrix* C);

// 2d avg pooling
void mtx_avgpool2d(Matrix* A, int* kernel_dim, int* strides, char padding, Matrix* C);

// activations
void mtx_activation(Matrix* A, int activation_type, Matrix* C);

// fast fourier transform
void mtx_fft(Matrix* A, Matrix* C);

*/
#endif

#endif
