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
    this->dlen = src_mat.dlen;
    // copy data
    this->data = src_mat.data;
    // copy other metadata
    this->is_sparse = src_mat.is_sparse;
    this->sp_data = src_mat.sp_data;
}

template<class mtx_type>
NDMat<mtx_type>::NDMat(std::vector<int> mshape, mtx_type init_value){
    // check for data type
    // to be done
    this->shape = mshape;
    this->dlen = mshape[0];
    for(int i = 1; i < mshape.size(); i++){
        this->dlen *= mshape[i];
    }
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
NDMat<mtx_type> NDMat<mtx_type>::copy(std::vector<int> padding){
    // create new shape vector
    int dc = this->shape.size();
    std::vector<int> nshape(this->shape.begin(), this->shape.end());
    std::vector<int> vcount_t(nshape.begin(), nshape.end());
    std::vector<int> vcount_s(nshape.begin(), nshape.end());
    std::vector<int> vi(dc, 1);
    // add padding for each dimensions
    nshape[dc-1] += padding[i];
    vcount_t[dc-1] = 1;
    vcount_s[dc-1] = 1;
    for(int i = dc-2; i > -1; --i){
        nshape[i] += padding[i];
        vcount_t[i] = vcount_t[i+1] * nshape[i+1];
        vcount_s[i] = vcount_s[i+1] * this->shape[i];
    }
    NDMat<mtx_type> new_mat(nshape);
    
}

// create zero initialized matrix similar to current matrix
template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::zeros_like(NDMat<mtx_type> src_mat){
    return NDMat<mtx_type>(src_mat.getShape());
}

template<class mtx_type>
NDMat<mtx_type> NDMat<mtx_type>::slice(std::vector<Point2D> axis_min_max){  // duplicate data in new matrix

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
