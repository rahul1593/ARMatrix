/*
 * ARMATRIX
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

#ifdef MODE_C
    #include <stdlib.h>
    #include <stdio.h>
    #include <stdint.h>
    #include <math.h>
#else
    #include <iostream>
    #include <vector>
    #include <cmath>
    #include <cstdlib>
    #include <cstdint>

    using namespace std;
#endif

#define mtx_type    double

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
#define AMTX_UI16   1
#define AMTX_UI32   2
#define AMTX_UI64   3
#define AMTX_I8     4
#define AMTX_I16    5
#define AMTX_I32    6
#define AMTX_I64    7
#define AMTX_F32    8   // single and double precision floating types
#define AMTX_F64    9
#define AMTX_UCR    10  // character type
#define AMTX_ANY    11  // user defined type, input as class object

typedef struct{
    int x,y;
} Index2D;

typedef struct{
    int x,y,z;
} Index3D;

class NDMat{
private:
    int                 type;           // data type of the matrix
    vector<int>         shape;          // matrix dimensions, outermost to innermost
    int                 dlen;           // length of data array
    void                *data;          // pointer to the arrray of any specified data type
    vector<Index2D>     joinIndx;       // set of indexes to join disjoint data, used in referencing sliced matrices
                                        // data is stored in row major array, type depends on user
    bool                is_sparse;      // is a sparse matrix with a lot of zeros, eg, identity matrix
    vector<vector<int>> sp_data;    // index data for sparse matrix
public:
    // constructor with shape as a vector, default type is AMTX_F32, i.e, float32
    NDMat(vector<int> shape, int type=AMTX_F32, void* type_object=(void*)NULL);
    NDMat(vector<int> shape, int type=AMTX_F32, void* init_value=(void*)NULL, void* type_object=(void*)NULL);
    NDMat(vector<int> shape, vector<int> padding, int type=AMTX_F32, void* type_object=(void*)NULL);
    NDMat(vector<int> shape, vector<int> padding, int type=AMTX_F32, void* init_value=(void*)NULL, void* type_object=(void*)NULL);
    // destructor for NDMat object
    ~NDMat();

    // create using current object, with given padding in each dimension
    NDMat NDMatFrom(NDMat src_mat);
    NDMat NDMatFrom(NDMat src_mat, vector<int> padding);

    // create a duplicate matrix by copying the current matrix
    NDMat copy();

    // change the data type of the matrix
    void toType(int type, void* type_object=(void*)NULL);

    // create zero initialized matrix similar to current matrix
    NDMat zeros();

    // create a submatrix by slicing the current matrix
    NDMat slice(vector<Index2D> axis_min_max);  // duplicate data in new matrix
    NDMat vslice(vector<Index2D> axis_min_max); // slice but refer to original data

    // reshape the matrix
    void reshape(vector<int> shape);

    // create a copy of cu
};


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
