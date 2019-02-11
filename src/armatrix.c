/*
 * ARMATRIX
 * C Source file for functions in "armatrix.h".
 * 
 * Author       : Rahul Bhartari
 * Date Created : Dec 1, 2018
 * License      : LGPL2.1
 */

#include "armatrix.h"

/*
 * Internal function declarations
 */
//void _int_vector_op(mtx_type* data_A, mtx_type *data_B, int dlen, int start_pos, int op);

int _int_print_mtx(mtx_type* data, int *shape, int dim_cnt, int start_pos, int lvl);
int* _int_check_compatible_broadcast_dim(Matrix *A, Matrix *B);
//int _int_check_compatible_mul_dim(Matrix *A, Matrix *B);
void __mtx_op(Matrix* A, Matrix* B, Matrix* C, int op);
void __internal_mat_op(Matrix* A, Matrix* B, Matrix* C, int dim_diff, int op);
void __vector_op(mtx_type* vec1, mtx_type* vec2, int ln1, int ln2, int op, mtx_type* ovec);

/*
 *  Definitions for internal functions of library
 */

int _int_print_mtx(mtx_type* data, int *shape, int dim_cnt, int start_pos, int lvl){
    int i, d, new_pos;
    int *new_shape;
    
    if(dim_cnt > 1){
        new_shape = (int*)malloc(sizeof(int)*(dim_cnt-1));
        for(i = 1; i < dim_cnt; ++i){
            new_shape[i-1] = shape[i];
        }
        // update the start position
        new_pos = start_pos;
        for(i = 0; i < shape[0]; ++i){
            if(i > 0){// print the space denoting the level
                printf("\n");
                for(d = 0; d < lvl; ++d) printf(" ");
            }
            // print the inner dimensions
            printf("[");
            new_pos = _int_print_mtx(data, new_shape, dim_cnt-1, new_pos, lvl+1);
            printf("]");
        };
        free(new_shape);
    } else {// print the innermost dimension
        for(i = start_pos; i < start_pos+shape[0]; ++i){
            printf("%0.4f ", data[i]);
        }
        printf("\b");
        // update the start position for upcoming prints
        new_pos = start_pos + shape[0];
    }
    return new_pos;
}

int* _int_check_compatible_broadcast_dim(Matrix *A, Matrix *B){
    int i, dmc, df, f = 0;
    Matrix *M1, *M2;
    int *cshape;
    if(A->dim_cnt >= B->dim_cnt){
        M1 = A;
        M2 = B;
    } else {
        M1 = B;
        M2 = A;
    }
    
    // check if dimensions match after equalising the dimensionality by prepending 1's to Matrix with less dimensions
    dmc = M1->dim_cnt;			// maximum dimension
    cshape = (int*)malloc(sizeof(int)*dmc);
	df = M1->dim_cnt - M2->dim_cnt;	// dimensionality difference
	for(i=0; i < df; ++i) cshape[i] = M1->shape[i];
	
	for(i=df; i < dmc; ++i){
		if(M1->shape[i] != M2->shape[i - df]){		// if dimension does NOT match and none of the dimensions is 1
			if(M1->shape[i] != 1 && M2->shape[i - df] != 1){
                free(cshape);
                return NULL;
            }
            cshape[i] = (M1->shape[i] == 1) ? M2->shape[i] : M1->shape[i];
		} else {
            cshape[i] = M1->shape[i];
        }
	}
    return cshape;
}


void __vector_op(mtx_type* vec1, mtx_type* vec2, int ln1, int ln2, int op, mtx_type* ovec){
    mtx_type *v1, *v2;
    int s1, s2;
    int i;
    
    if(ln1 == 1){
        v2 = vec1;
        s2 = ln1;
        v1 = vec2;
        s1 = ln2;
    } else {
        v1 = vec1;
        s1 = ln1;
        v2 = vec2;
        s2 = ln2;
    }
    
    switch(op){
        case OP_ADD:
            if(s2 == 1){
                for(i = 0; i < s1; ++i){
                    ovec[i] = v1[i] + v2[0];
                }
            } else {
                for(i = 0; i < s1; ++i){
                    ovec[i] = v1[i] + v2[i];
                }
            }
            break;
        case OP_SUB:
            if(s2 == 1){
                for(i = 0; i < s1; ++i){
                    ovec[i] = v1[i] - v2[0];
                }
            } else {
                for(i = 0; i < s1; ++i){
                    ovec[i] = v1[i] - v2[i];
                }
            }
            break;
        case OP_MUL:
            if(s2 == 1){
                for(i = 0; i < s1; ++i){
                    ovec[i] = v1[i] * v2[0];
                }
            } else {
                for(i = 0; i < s1; ++i){
                    ovec[i] = v1[i] * v2[i];
                }
            }
            break;
        case OP_DIV:
            if(s2 == 1){
                for(i = 0; i < s1; ++i){
                    ovec[i] = v1[i] / v2[0];
                }
            } else {
                for(i = 0; i < s1; ++i){
                    ovec[i] = v1[i] / v2[i];
                }
            }
            break;
        default:
            printf("Error: Invalid Operation for %s\n", __func__);
            break;
    }
}

/*
typedef struct Matrix{
    int* shape;         // matrix dimensions, outermost to innermost
    int* increments;    // store length of each dimension
    int dim_cnt;        // number of dimensions
    int dlen;           // length of data array
    mtx_type* data;     // data can be of stated type, integer by default, data is stored in row major array
	int is_sparse;		// is a sparse matrix with a lot of zeroes, eg., identity matrix
	int** coordinates;  // store the coordinates for values in sparse matrix
};
 */

void __internal_mat_op(Matrix* A, Matrix* B, Matrix* C, int dim_diff, int op){
    Matrix M1, M2, Mr;
    Matrix *M2_p;
    int i;
    if(A->dim_cnt > 1){
        // update submatrix dimensions
        M1.shape = A->shape + 1;
        M1.increments = A->increments + 1;
        M1.dim_cnt = A->dim_cnt - 1;
        M1.dlen = A->increments[0];
        if(dim_diff == 0){
            M2.shape = B->shape + 1;
            M2.increments = B->increments + 1;
            M2.dim_cnt = B->dim_cnt-1;
            M2.dlen = B->increments[0];
            M2_p = &M2;
        } else {
            M2_p = B;
            --dim_diff;
        }
        Mr.shape = C->shape + 1;
        Mr.increments = C->increments + 1;
        Mr.dim_cnt = C->dim_cnt - 1;
        Mr.dlen = C->increments[0];
        for(i=0; i < C->shape[0]; ++i){
            if(A->shape[0] == 1){
                M1.data = A->data;
            } else {
                M1.data = A->data + i*A->increments[0];
            }
            if(B->shape[0] == 1){
                M2.data = B->data;
            } else {
                M2.data = B->data + i*B->increments[0];
            }
            Mr.data = C->data + i*C->increments[0];
            __internal_mat_op(&M1, M2_p, &Mr, dim_diff, op);
        }
    } else {
        // do the vector addition operation
        __vector_op(A->data, B->data, A->dlen, B->dlen, op, C->data);
    }
}

void __mtx_op(Matrix* A, Matrix* B, Matrix* C, int op){
    int i, n, lvl, oln, cA, cB, cC;
    Matrix *M1, *M2;
    int *cshape;
    
    // check the dimensions of the matrices
    cshape = _int_check_compatible_broadcast_dim(A, B);
    if(cshape == NULL){
        printf("Error: Matrix dimension mismatch.\n");
        return;
    }
    // set the matrices as M1 with maximum dimension count
    if(A->dim_cnt >= B->dim_cnt){
        M1 = A;
        M2 = B;
    } else {
        M1 = B;
        M2 = A;
    }
    /* M2 is the matrix with less number of dimensions
     * Copy matrix with less dimensions according to dimensionality
     */
    
    // align the output matrix with the output
    for(i=0, oln=1; i < M1->dim_cnt; ++i) oln *= cshape[i];
    free(C->shape);
    C->shape = cshape;
    if(C->dim_cnt != M1->dim_cnt || C->dlen != oln){
        C->data = (mtx_type*)realloc(C->data, sizeof(mtx_type)*oln);
        C->dim_cnt = M1->dim_cnt;
        C->dlen = oln;
    }
    // pass submatrices in recursive fashion
    // submatrices must have the pointer for data from the original matrices
    __internal_mat_op(M1, M2, C, M1->dim_cnt - M2->dim_cnt, op);
}

/*
 * -------- Library function definitions --------
 */

/*
 * pretty-print the matrix
 */
void mtx_print(Matrix* A){
    // get the total number of items in the array
    printf("[");
    _int_print_mtx(A->data, A->shape, A->dim_cnt, 0, 1);
    printf("]\n");
}

/*
 * create a matrix initialized with given value
 */
Matrix* mtx_create(int* shape, int dim_cnt, mtx_type init_val){
    int i, dl = 1;
    Matrix *A;
    A = (Matrix*)malloc(sizeof(Matrix));
    A->shape = (int*)malloc(sizeof(int)*dim_cnt);
    A->increments = (int*)malloc(sizeof(int)*dim_cnt);
    for(i = 0; i < dim_cnt; ++i){
        A->shape[i] = shape[i];
        A->increments[dim_cnt-i-1] = dl;
        dl *= shape[dim_cnt-i-1];
    }
    
    A->dim_cnt = dim_cnt;
    A->dlen = dl;
    A->data = (mtx_type*)malloc(sizeof(mtx_type)*dl);
    for(i = 0; i < dl; ++i){
        A->data[i] = init_val;
    }
    return A;
}

/*
 * Create a matrix similar to given matrix
 */
Matrix* mtx_create_as(Matrix M){
    return mtx_create(M.shape, M.dim_cnt, 1);
}

/*
 * destroy the input matrix and free the space allocated
 */
void mtx_destroy(Matrix** A){
    free((*A)->shape);
    free((*A)->increments);
    free((*A)->data);
    
    (*A)->shape = (int*)NULL;
    (*A)->increments = (int*)NULL;
    (*A)->data = (mtx_type*)NULL;
    
    free((*A));
    *A = (Matrix*)NULL;
}

/*
 * Reshape the matrix in place
 * Returns -1 if invalid shape
 */
int mtx_reshape(Matrix* A, int* shape, int dim_cnt){
    int i, dl_new = 1;
    // check total number of entries as per new shape
    for(i = 0; i < dim_cnt; ++i){
        dl_new *= shape[i];
    }
    // check if array lengths are equal
    if(dl_new != A->dlen)
        return -1;
    // reallocate old space for shape
    A->shape = (int*)realloc(A->shape, sizeof(int)*dim_cnt);
    // store new shape
    for(i = 0; i < dim_cnt; ++i){
        A->shape[i] = shape[i];
    }
    // update the dimensionality
    A->dim_cnt = dim_cnt;
    A->dlen = dl_new;
    return 0;
}

/*
 * Create copy of a matrix
 * Arguments:
 *      src : pointer to source matrix
 *      tgt : pointer to target matrix
 */
void mtx_copy(Matrix* src, Matrix* tgt){
    int i;
    for(i = 0; i < src->dlen; ++i){
        tgt->data[i] = src->data[i];
    }
}

/* 
 * Matrix addition, C = A + B. The addition is done according to the dimensionality of the matrices
 */
void mtx_add(Matrix* A, Matrix* B, Matrix* C){
    __mtx_op(A, B, C, OP_ADD);
}

/*
 * Matrix subtraction, C = A - B
 */
void mtx_sub(Matrix* A, Matrix* B, Matrix* C){
    __mtx_op(A, B, C, OP_SUB);
}

/*
 * Matrix scalar multiplication, C = A * B
 */
void mtx_mul(Matrix* A, Matrix* B, Matrix* C){
    __mtx_op(A, B, C, OP_MUL);
}

/*
 * Matrix scalar multiplication, C = A / B
 */
void mtx_div(Matrix* A, Matrix* B, Matrix* C){
    __mtx_op(A, B, C, OP_DIV);
}

/* 
 * Scale a matrix, i.e, multiply by a constant 'c' of type same as of matrix A
 */
void mtx_scale(Matrix *A, mtx_type c){
    int i;
    for(i = 0; i < A->dlen; ++i){
        A->data[i] *= c;
    }
}

/*
 * Matrix dot product, C = A.B
 * Arguments:
 * 	C = A.B
 */
void mtx_dot(Matrix* A, Matrix* B, Matrix* C){
	
}

/*
 * Matrix transpose
 */
void mtx_transpose(Matrix *A, Matrix *T){
    
}

/* 
 * Matrix invert
 */
void mtx_invert(Matrix* A, Matrix* Ai){
    
}

/*
 * Identity matrix of dimensions as specified by 'shape' (array of length 'dim')
 */
void mtx_identity(int* shape, int dim, Matrix* I){
    
}
