/*
 * ARMATRIX.CPP
 * C++ Source file for functions in "matrics.h".
 * 
 * Author       : Rahul Bhartari
 * Date Created : Feb 10, 2019
 * License      : LGPL2.1
 */

#include "armatrix.h"

template<class mtx_type>
std::vector<int> NDMat<mtx_type>::__checkCompatibleBroadcastDim(NDMat<mtx_type> A, NDMat<mtx_type> B){

}



/*
    Constructor definitions for template class NDMat
*/
template<class mtx_type>
NDMat<mtx_type>::NDMat(const NDMat &src_mat){
    // copy the shape
    this->type = src_mat.type;
    this->shape = src_mat.shape;
    this->increments = src_mat.increments;
    this->dlen = src_mat.dlen;
    // copy data
    this->data = src_mat.data;
    // copy other metadata
    this->is_sparse = src_mat.is_sparse;
    this->sp_data = src_mat.sp_data;
}

template<class mtx_type>
NDMat<mtx_type>::NDMat(std::vector<int> mshape, mtx_type init_value){
    int dl = 1;
    int dim_cnt = mshape.size();
    this->shape = mshape;
    this->dlen = mshape[0];
    for(int i = 0; i < dim_cnt; ++i){
        this->increments[dim_cnt-i-1] = dl;
        dl *= mshape[dim_cnt-i-1];
    }
    this->dlen = dl;
    // copy data
    std::vector<mtx_type> vdata(this->dlen, init_value);
    this->data = vdata;
    // copy other metadata
    this->is_sparse = false;
}

/*
    Destructor definition for template class NDMat<mtx_type>
*/
template<class mtx_type>
NDMat<mtx_type>::~NDMat(){
    
}


// generate an identity matrix of given shape and data type
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::identity(int side){
    std::vector<int> __shape[2];
    __shape[0] = __shape[1] = side;
    
    NDMat<mtx_type> m1(__shape);
    int stride = side + 1;

    for(int i = 0; i < side*side; i += stride){
        m1.data[i] = 1;
    }
    return m1;
}

// get shape of this matrix
template<class mtx_type>
std::vector<int> NDMat<mtx_type>::getShape(){
    return this->shape;
}

// check if this matrix is sparse
template<class mtx_type>
bool NDMat<mtx_type>::isSparse(){
    return this->is_sparse;
}

// get total number of elements in this matrix
template<class mtx_type>
int NDMat<mtx_type>::len(){
    return this->dlen;
}

// Reshape the matrix
template<class mtx_type>
void NDMat<mtx_type>::reshape(std::vector<int> nshape){
    int i, dl_new = 1;
    // check total number of entries as per new shape
    for(i = 0; i < dim_cnt; ++i){
        dl_new *= nshape[i];
    }
    // check if array lengths are equal
    if(dl_new != this->dlen) return;
    // else update the shape
    this->shape = nshape;
    this->increments = nshape;  // assign array similar to nshape
    for(i = 0, dl_new = 1; i < dim_cnt; ++i){
        this->increments[dim_cnt-i-1] = dl_new;
        dl_new *= mshape[dim_cnt-i-1];
    }
}

// change the data type of the matrix
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::asType(int type){
    switch(type){
        case    AMTX_UI8:
            NDMat<uint8_t> new_mat(this->shape);
            std::vector<uint8_t> vdata(this->data.begin(), this->data.end());
            new_mat.data = vdata;
            new_mat.is_sparse = this->is_sparse;
            new_mat.sp_data = this->sp_data;
            return new_mat;
        case    AMTX_UI32:
            NDMat<uint32_t> new_mat(this->shape);
            std::vector<uint32_t> vdata(this->data.begin(), this->data.end());
            new_mat.data = vdata;
            new_mat.is_sparse = this->is_sparse;
            new_mat.sp_data = this->sp_data;
            return new_mat;
        case    AMTX_UI64:
            NDMat<uint64_t> new_mat(this->shape);
            std::vector<uint64_t> vdata(this->data.begin(), this->data.end());
            new_mat.data = vdata;
            new_mat.is_sparse = this->is_sparse;
            new_mat.sp_data = this->sp_data;
            return new_mat;
        case    AMTX_I8 :
            NDMat<int8_t> new_mat(this->shape);
            std::vector<int8_t> vdata(this->data.begin(), this->data.end());
            new_mat.data = vdata;
            new_mat.is_sparse = this->is_sparse;
            new_mat.sp_data = this->sp_data;
            return new_mat;
        case    AMTX_I32:
            NDMat<int32_t> new_mat(this->shape);
            std::vector<int32_t> vdata(this->data.begin(), this->data.end());
            new_mat.data = vdata;
            new_mat.is_sparse = this->is_sparse;
            newnew_matMat.sp_data = this->sp_data;
            return new_mat;
        case    AMTX_I64:
            NDMat<int64_t> new_mat(this->shape);
            std::vector<int64_t> vdata(this->data.begin(), this->data.end());
            new_mat.data = vdata;
            new_mat.is_sparse = this->is_sparse;
            new_mat.sp_data = this->sp_data;
            return new_mat;
        case    AMTX_F32:
            NDMat<float> new_mat(this->shape);
            std::vector<float> vdata(this->data.begin(), this->data.end());
            new_mat.data = vdata;
            new_mat.is_sparse = this->is_sparse;
            new_mat.sp_data = this->sp_data;
            return new_mat;
        case    AMTX_F64:
            NDMat<double> new_mat(this->shape);
            std::vector<double> vdata(this->data.begin(), this->data.end());
            new_mat.data = vdata;
            new_mat.is_sparse = this->is_sparse;
            new_mat.sp_data = this->sp_data;
            return new_mat;
        default         :
            break;
    }
    return *this;
}

// print this matrix on console ****
template<class mtx_type>
void NDMat<mtx_type>::print(){

}

// create a duplicate matrix by copying the current matrix, with and without padding
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::copy(){
    NDMat<mtx_type> new_mat(*this);
    return new_mat;
}

template<class mtx_type>
void NDMat<mtx_type>::__recursivePaddedVectorCopy(mtx_type* src_data, int* src_shape, int* src_inc, 
                                                mtx_type* tgt_data, int* tgt_shape, int* tgt_inc, int dim_cnt){
    // copy data from source with given shape data pointer to target with given shape and data pointer
    int padding=(tgt_shape[0]-src_shape[0])/2;
    if(dim_cnt > 1){
        // recursively traverse inwards
        for(int i = padding; i < src_shape[0]+padding; ++i){
            NDMat<mtx_type>::__recursivePaddedVectorCopy(src_data + ((i - padding) * src_inc[0]),
                                                        src_shape + 1,
                                                        src_inc + 1,
                                                        tgt_data + (i * tgt_inc[0]),
                                                        tgt_shape + 1,
                                                        tgt_inc + 1,
                                                        dim_cnt - 1);
        }
    } else {
        // copy the innermost dimension
        for(int i = 0; i < src_shape[0]; ++i){
            tgt_data[padding + i] = src_data[i];
        }
    }
}

template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::copy(std::vector<int> padding){
    // create new shape vector
    int dc = this->shape.size();    // dimension count
    std::vector<int> nshape(this->shape.begin(), this->shape.end());        // new shape vector
    std::vector<int> v_off(dc, 1);     // index offset for different dimensions
    
    // add padding for each dimensions
    for(int i = 0; i < dc; ++i){
        nshape[i] += 2*padding[i];
    }
    // new padded matrix
    NDMat<mtx_type> new_mat(nshape);

    for(int i = padding[0]; i < this->shape[0]+padding[0]; ++i){
        NDMat<mtx_type>::__recursivePaddedVectorCopy(this->data.begin() + ((i - padding[0]) * this->increments[0]),
                                                    this->shape.begin()+1,
                                                    this->increments.begin()+1,
                                                    new_mat.data.begin() + (i * new_mat.increments[0]),
                                                    nshape.begin()+1,
                                                    new_mat.increments.begin()+1,
                                                    dc-1);
    }
    return new_mat;
}

// create zero initialized matrix similar to current matrix
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::zeros_like(NDMat<mtx_type> src_mat){
    return NDMat<mtx_type>(src_mat.getShape());
}

template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::slice(std::vector<Point2D> axis_min_max){  // duplicate data in new matrix
    // new shape for the matrix
    std::vector<mtx_type> nshape(this->shape.size());
    for(int i=0; i < this->shape.size(); i++){
        nshape[i] = axis_min_max[i].y - axis_min_max[i].x;
    }

    NDMat<mtx_type> new_mat(nshape);
    // prepare offsets
    for(int i = axis_min_max[0].x; i < axis_min_max[0].y; ++i){
        NDMat<mtx_type>::__recursivePaddedVectorCopy(this->data.begin() + ((i - padding[0]) * this->increments[0]),
                                                    this->shape.begin()+1,
                                                    this->increments.begin()+1,
                                                    new_mat.data.begin() + (i * new_mat.increments[0]),
                                                    nshape.begin()+1,
                                                    new_mat.increments.begin()+1,
                                                    dc-1);
    }
    return new_mat;
}

template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::vslice(std::vector<Point2D> axis_min_max){ // slice but refer to original data

}

template<class mtx_type>
void NDMat<mtx_type>::assignSlice(std::vector<Point2D> axis_min_max, NDMat<mtx_type> set_val){

}

template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::maskDup(NDMat<bool> condition, mtx_type val, int axis){

}

template<class mtx_type>
void NDMat<mtx_type>::mask(NDMat<bool> condition, mtx_type val, int axis){

}

template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::operateDup(mtx_type *(func)(mtx_type *), int axis){

}

template<class mtx_type>
void NDMat<mtx_type>::operate(mtx_type *(func)(mtx_type *), int axis){

}

// reverse a matrix around the given axis (default to '-1' for all values)
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::reverseDup(int axis){    // returns a duplicate matrix

}

template<class mtx_type>
void NDMat<mtx_type>::reverse(int axis){        // in place reverse

}

// Transpose of this matrix with copied data or in place
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::transposeDup(){

}

template<class mtx_type>
void NDMat<mtx_type>::transpose(){

}

// Inverse of this matrix
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::inverseDup(){

}

template<class mtx_type>
void NDMat<mtx_type>::inverse(){

}

// dot product of two matrices
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::dot(NDMat<mtx_type> m2){

}

// add two matrices
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::add(NDMat<mtx_type> m2){

}

// subtract matrix B from A
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::sub(NDMat<mtx_type> m2){

}

// elementwise product
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::mul(NDMat<mtx_type> m2){

}

// elementwise division
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::div(NDMat<mtx_type> m2){

}

// scale a matrix, i.e, multiply by a constant pointed by 'c'
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::scaleDup(mtx_type c){

}

template<class mtx_type>
void NDMat<mtx_type>::scale(mtx_type c){

}

// Overloaded arithmetic operators for easy arithmetic
// assign a matrix to another
template<class mtx_type>
NDMat<mtx_type> & NDMat<mtx_type>::operator=(NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<mtx_type> & NDMat<mtx_type>::operator=(mtx_type val){

}

// add two matrices
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::operator+(NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::operator+(mtx_type val){

}

template<class mtx_type>
NDMat<mtx_type> & NDMat<mtx_type>::operator+=(NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<mtx_type> & NDMat<mtx_type>::operator+=(mtx_type val){

}

// subtract one matrix from another matrix
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::operator-(NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::operator-(mtx_type val){

}

template<class mtx_type>
NDMat<mtx_type> & NDMat<mtx_type>::operator-=(NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<mtx_type> & NDMat<mtx_type>::operator-=(mtx_type val){

}

// elementwise multiply two matrices
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::operator*(NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::operator*(mtx_type val){

}

template<class mtx_type>
NDMat<mtx_type> & NDMat<mtx_type>::operator*=(NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<mtx_type> & NDMat<mtx_type>::operator*=(mtx_type val){

}

// elementwise divide one matrix by another
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::operator/(NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::operator/(mtx_type val){

}

template<class mtx_type>
NDMat<mtx_type> & NDMat<mtx_type>::operator/=(NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<mtx_type> & NDMat<mtx_type>::operator/=(mtx_type val){

}

// elementwise mod one matrix by another
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::operator%(NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::operator%(mtx_type val){

}

template<class mtx_type>
NDMat<mtx_type> & NDMat<mtx_type>::operator%=(NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<mtx_type> & NDMat<mtx_type>::operator%=(mtx_type val){

}

// Overloaded comparison operators
template<class mtx_type>
NDMat<bool> operator<(NDMat<mtx_type> &m1, NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<bool> operator>(NDMat<mtx_type> &m1, NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<bool> operator>=(NDMat<mtx_type> &m1, NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<bool> operator<=(NDMat<mtx_type> &m1, NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<bool> operator!=(NDMat<mtx_type> &m1, NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<bool> operator==(NDMat<mtx_type> &m1, NDMat<mtx_type> &m2){

}

template<class mtx_type>
NDMat<bool> operator<(NDMat<mtx_type> &m1, mtx_type val){

}

template<class mtx_type>
NDMat<bool> operator>(NDMat<mtx_type> &m1, mtx_type val){

}

template<class mtx_type>
NDMat<bool> operator>=(NDMat<mtx_type> &m1, mtx_type val){

}

template<class mtx_type>
NDMat<bool> operator<=(NDMat<mtx_type> &m1, mtx_type val){

}

template<class mtx_type>
NDMat<bool> operator!=(NDMat<mtx_type> &m1, mtx_type val){

}

template<class mtx_type>
NDMat<bool> operator==(NDMat<mtx_type> &m1, mtx_type val){

}


// Overloaded stream operator
template<class mtx_type>
NDMat<mtx_type> & operator << (std::ostream &out, NDMat<mtx_type> &src_mat){

}

// Overloaded function operator to get any value
// Takes dimension index as arguments, returns the value or submatrix address
template<class mtx_type>
NDMat<mtx_type> & NDMat<mtx_type>::operator()(...){

}

// operator overloading for masking
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::operator[](NDMat<bool> condition){

}
