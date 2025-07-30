#ifndef MATRIX_H
#define MATRIX_H

#include "utils.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <seal/util/numth.h>
#include <iomanip>


// using namespace seal;
using namespace seal::util;
// using namespace std;




// begin class Matrix

template <typename T>
class Matrix
{
protected:
    std::size_t n, d;
    std::vector<std::vector<T>> M;

public:
    // matrix();

    // empty matrix
    Matrix() : n(0), d(0) {}

    //~matrix();

    //  rows * cols, all elements are initialized with the default constructor of T
    Matrix(std::size_t rows, std::size_t cols) : n(0), d(0)
    {
        resize(rows, cols);
    }

    void clear()
    {
        n = d = 0;
        M.clear();
    }

    void resize(std::size_t rows, std::size_t cols)
    {
        std::size_t j;
        n = rows;
        d = cols;
        M.resize(d);
        for (j = 0; j < d; j++)
        {
            M[j].resize(n);
        }
    } // resize

    // return the number of rows
    size_t get_rows() const
    {
        return n;
    }

    // return the number of columns
    size_t get_cols() const
    {
        return d;
    }

    // a reference to the element (i, j)
    T &operator()(const std::size_t i, const std::size_t j) { return M[j][i]; }

    const T &operator()(const std::size_t i, const std::size_t j) const { return M[j][i]; }

    inline T &get(const std::size_t i, const std::size_t j) { return M[j][i]; }

    inline void set(const std::size_t i, const std::size_t j, const T a) { M[j][i] = a; }

    // the transpose of the matrix
    Matrix<T> transpose() const
    {
        Matrix<T> B(d, n);
        for (std::size_t i = 0; i < n; i++)
            for (std::size_t j = 0; j < d; j++)
                B(j, i) = M[j][i];
        return B;
    }

    // return the i-th row of the matrix; indices start from 0.
    std::vector<T> get_row(std::size_t i)
    {
        std::vector<T> v;
        v.resize(d);
        for (std::size_t j = 0; j < d; j++)
            v[j] = M[j][i];

        return v;
    }

    std::vector<T> get_row(std::size_t i) const
    {
        std::vector<T> v;
        v.resize(d);
        for (std::size_t j = 0; j < d; j++)
            v[j] = M[j][i];

        return v;
    }

    // return the matirx by rows
    std::vector<std::vector<T>> get_matrix_by_rows()
    {
        std::vector<std::vector<T>> destination;
        std::vector<T> row;
        for(std::size_t i=0;i<n;i++){
            row=get_row(i);
            destination.push_back(row);
        }
        return destination;
    }

    //generate matrix by rows vector
    void gen_matrix_by_rows(std::vector<std::vector<T>> &data){
        n=data.size();
        d=data[0].size();
        resize(n,d);
        for(size_t i=0;i<n;i++){
            for(size_t j=0;j<d;j++){
                M[j][i]=data[i][j];
            }
        }
    }

    // flatten matrix to a vector by rows 
    std::vector<T> flatten_matrix_to_rows_vector() {
        std::vector<T> flat_vector(n*d);

        for (std::size_t i = 0; i < n; i++) {
            for (std::size_t j = 0; j < d; j++) {
                flat_vector[i*d+j] = M[j][i] ;
            }
        }
        return flat_vector;
    }

    // flatten matrix to a vector by cols 
    std::vector<T> flatten_matrix_to_cols_vector() {
        std::vector<T> flat_vector(n*d);

        for (std::size_t i = 0; i < d; i++) {
            for (std::size_t j = 0; j < n; j++) {
                flat_vector[i*n+j] = M[i][j] ;
            }
        }
        return flat_vector;
    }

    // return the last row of the matrix
    std::vector<T> get_last_row()
    {
        return get_row(n - 1);
    }

    std::vector<T> get_last_row() const
    {
        return get_row(n - 1);
    }

    // return the i-th column of the matrix; indices start from 0.
    std::vector<T> get_col(std::size_t j)
    {
        std::vector<T> v;
        v.resize(n);
        for (std::size_t i = 0; i < n; i++)
            v[i] = M[j][i];

        return v;
    }

    std::vector<T> get_col(std::size_t j) const
    {
        std::vector<T> v;
        v.resize(n);
        for (std::size_t i = 0; i < n; i++)
            v[i] = M[j][i];

        return v;
    }

    // return the last row of the matrix
    std::vector<T> get_last_col()
    {
        return get_col(d - 1);
    }

    std::vector<T> get_last_col() const
    {
        return get_col(d - 1);
    }

    std::vector<T> diagonal_vector(std::size_t ell){
        std::vector<T> v;
        v.clear();
        if (d % n == 0){
            if ((ell >=0)&&(ell < n)){
                for(std::size_t i = 0; i < d; i++){
                    v.push_back(M[(ell + i) % d][i % n]);
                }
            }
            else{
                std::cerr << "the input parameter must be in [0, #rows)"  << std::endl;
            }
        }
        else{
            std::cerr << "#rows must divide #cols"  << std::endl;
        }
        return v;
    }


    std::vector<T> diagonal_encoding(){
        std::vector<T> v;
        if (std::gcd(n, d) == 1){
            
            for (std::size_t i = 0; i < n*d; i++){
                v.push_back(M[i % d][i % n]);
            }
            
        }
        else{
            std::cerr << "diagonal_encoding: the input dimension must be coprime. " << std::endl;
        }
        return v;

    }

    //Jiang function generate matrix
    //generate u_sigma matrix
    void generate_u_sigma(std::size_t rows, std::size_t cols) {
        resize(rows * cols, rows * cols);
        for (std::size_t i = 0; i < cols; i++) {
            for (std::size_t j= 0; j < rows; j++) {
                M[i * rows + (j + i) % rows][i * rows + j] = 1;
            }
        }
    }

    //generate u_tau matrix
    void generate_u_tau(std::size_t rows, std::size_t cols) {
        resize(rows * cols, rows * cols);
        for (std::size_t i = 0; i < cols; i++) {
            for (std::size_t j= 0; j < rows; j++) {
                M[((rows + 1) * j + cols * i) % (rows * cols)][i * rows + j] = 1;
                // std::cout<<"rows:"<<i * rows + j<<"  cols:"<<((rows + 1) * j + cols * i) % (rows * cols)<<std::endl;
            }
        }
    }

    // get i-th diag vector
    std::vector<T> diag_vector(std::size_t i,size_t slot_conunt) {
        std::vector<T> diag_vector;
        if (n != d) {
            std::cout << "matrix is not a squre" << std::endl;
            return diag_vector;
        }
        for (std::size_t j = 0; j < n; j++) {
            diag_vector.push_back(M[(i + j) % d][j]);
        }
        diag_vector.resize(slot_conunt);
        return diag_vector;
    }
    /*------------------------------------------------------------*/

    void scaling_inplace(T c){
        for (std::size_t i=0; i < n; i++){
            for (std::size_t j = 0; j<d; j++){
                M[j][i] *= c;
            }
        }
    }

    Matrix<T> scaling(T c){
        Matrix<T> R(n, d);
        for (std::size_t i=0; i < n; i++){
            for (std::size_t j = 0; j<d; j++){
                R.set(i, j, M[j][i] * c);
            }
        }
        return R;
    }

    // return if the matrix is the zero matrix in the sense 
    // that the abs of each elements is less than pow(10, a)
    bool is_zero(int64_t a = -5){
        double large_error=0;
        for (std::size_t i=0; i < n; i++){
            for (std::size_t j = 0; j<d; j++){
                if (abs(double((M[j][i]))) > (double)pow(10, a)){
                    std::cout << "Largest absolute error > "<<pow(10,a) << std::endl;
                    return false;
                }
                if(abs(double((M[j][i])))>large_error){
                    large_error=abs(double((M[j][i])));
                }
            }
        }
        std::cout<<"Largest absolute error: |"<<std::setprecision(10)<<large_error<<"|"<<std::endl;
        return true;
    }


     Matrix<T> add(Matrix<T> & N){
        std::size_t m = N.get_rows();
        std::size_t p = N.get_cols();
        Matrix<T> C;
        C.resize(n, d);
        if ((m != n) || (d != p)){
            std::cerr << "matrix::add: the dimensions do not match. " << std::endl;
        }
        else{
            for (std::size_t i=0; i < n; i++){
                for (std::size_t j = 0; j<d; j++){
                    C.set(i, j, M[j][i] + N(i, j));
                }
            }
        }
        return C;
    }

    Matrix<T> multiply(Matrix<T> & N){
        std::size_t m = N.get_rows();
        std::size_t p = N.get_cols();
        Matrix<T> C;
        C.resize(n, p);
        if (m != d){
            std::cerr << "matrix::multiply: the dimensions do not match. " << std::endl;
        }
        else{
            for (std::size_t i=0; i < n; i++){
                for (std::size_t j = 0; j<p; j++){
                    T t = 0;
                    for(std::size_t k = 0; k < m; k++){
                        t += M[k][i] * N(k, j);
                    }
                    C.set(i, j, t);
                }
            }
        }
        return C;
    }

    Matrix<T> submatrix(std::size_t p, std::size_t q)
    {
        Matrix<T> A(p, q);
        for(std::size_t i = 0; i < p; i++){
            for (std::size_t j = 0; j < q; j++){
                A.set(i, j, M[j][i]);
            }
        }
        return A;
    }

    // augumenting the original matrix to a matrix with larger dimensions 
    // by fill with zeros. 
    void augumented_matrix(std::size_t nn, std::size_t mm, Matrix<T> &B)
    {
        if ((nn < n) || (mm < d)){
            std::cerr << "augumented_matrix: the input dimensions must be at least as the original" << std::endl;
            return;
        }
        B.resize(nn, mm);
        for (std::size_t i = 0; i < n; i++){
            for(std::size_t j = 0; j < d; j++){
                B.set(i, j, M[j][i]);
            }
        }
        for (std::size_t i = n; i < nn; i++){
            for(std::size_t j = d; j < mm; j++){
                B.set(i, j, 0);
            }
        }
    }

    void Merage_matrix(size_t start_row,size_t end_row,size_t start_col,size_t end_col,Matrix<T> B){
        for(size_t i=start_row;i<end_row;i++){
            for(size_t j=start_col;j<end_col;j++){
                M[j][i]=B.get(i-start_row,j-start_col);
            }
        }
    }

    void print(std::size_t rows = 6, std::size_t cols = 6)
    {


        if ((n > 2*rows) && (d > 2*cols)){
            size_t r = rows / 2;
            size_t c = cols / 2;

            for (size_t i = 0; i < r; i++){
                std::cout << "    [";
                for (size_t j = 0; j < c; j++){
                    std::cout << std::setw(6) << std::right << std::right << M[j][i] << ",";
                }
                std::cout << std::setw(6) << std::right << std::right << "..."
                      << ",";
                for (size_t j = d - cols + c; j < d - 1; j++){
                    std::cout << std::setw(6) << std::right << std::right << M[j][i] << ",";
                }
                std::cout << std::setw(6) << std::right << std::right << M[d - 1][i] << "]" << std::endl;
            }

            std::cout << "    [";
            for (size_t j = 0; j < c; j++){
                std::cout << std::right << std::setw(6) << std::right << "..."
                      << ",";
            }
            std::cout << std::setw(6) << std::right << std::right << "..."
                  << ",";
            for (size_t j = d - cols + c; j < d - 1; j++){
                std::cout << std::setw(6) << std::right << std::right << "..."
                      << ",";
            }
            std::cout << std::setw(6) << std::right << std::right << "..."
                  << "]" << std::endl;

            for (size_t i = n - rows + r; i < n; i++){
                std::cout << "    [";
                for (size_t j = 0; j < c; j++){
                    std::cout << std::setw(6) << std::right << std::right << M[j][i] << ",";
                }
                std::cout << std::setw(6) << std::right << std::right << "..."
                      << ",";
                for (size_t j = d - cols + c; j < d - 1; j++){
                    std::cout << std::setw(6) << std::right << std::right << M[j][i] << ",";
                }
                std::cout << std::setw(6) << std::right << std::right << M[d - 1][i] << "]" << std::endl;
            }
        }
        else{
            for (size_t i = 0; i < n; i++){
                std::cout << "    [";
                for (size_t j = 0; j < d-1; j++){
                    std::cout << std::setw(6) << std::right << std::right << M[j][i] << ",";
                }
                std::cout << std::setw(6) << std::right << std::right << M[d - 1][i] << "]" << std::endl;
            }
        }
    }

}; // end of class matrix



void random_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A);
void order_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A);









// Some useful functions related to the matrix class

template <typename T>
Matrix<T> diagonal_decoding(std::vector<T> &diag_ecd, std::size_t n, std::size_t m){
    Matrix<T> A(n, m);
    if (std::gcd(n, m) != 1){
        std::cerr << "diagonal_decoding: the input dimension must be coprime. " << std::endl;
    }
    else{
        for (std::size_t k = 0; k < n*m; k++){
            A.set(k%n, k%m, diag_ecd[k]);
        }
    }
    return A;
}






template <typename T>
Matrix<T> diagonal_decoding_alternative(std::vector<T> &diag_ecd, std::size_t n, std::size_t m){
        
    // tuple<uint64_t, int64_t, int64_t> gcd;
    auto gcd = seal::util::xgcd(n, m);
    uint64_t g = get<0>(gcd);
    int64_t s = get<1>(gcd);
    int64_t t = get<2>(gcd);
    std::cout << "s: " << s << std::endl;
    std::cout << "t: " << t << std::endl;
    std::cout << "(-3) % 12 = " << mod(-3, 12) <<  std::endl;

    Matrix<T> A(n, m);
    if (g != 1){
        std::cerr << "diagonal_decoding: the input dimension must be coprime. " << std::endl;
    }
    else{
        for (std::size_t i = 0; i < n; i++){
            for (std::size_t j = 0; j < m; j++){
                std::size_t k = (i*t*m + j*s*n);
                // std::cout << "(i, j): (" << i << ", " << j << "); k: " << k << std::endl;
                k = mod(k, n*m);
                // std::cout << "(i, j): (" << i << ", " << j << "); k: " << k << std::endl;
                A.set(i, j, diag_ecd[k]);
            }
        }
    }
    return A;
}

// return the tranformation matrix for the (k, l) submatrix of a (n, m) matrix
void submatrix_transformation(std::size_t n, std::size_t m, std::size_t k, std::size_t l, Matrix<double> &T);


// return a vector consisting of the rotated nonzero diagonal vectors of the transformation matrix
//! note that the returned vector is of length slot_count, so when it is not suitable for packing 
//! different matrices in one ciphertext!
void rotated_diagonal_vectors_for_submatrix(std::size_t n, std::size_t m, std::size_t p, std::size_t q, 
                                            std::size_t slot_count,
                                            Matrix<double> &T, //the transformation matrix
                                            std::vector<std::size_t> &ind, // the indices of nonzero diagonal vectors of T
                                            std::vector<std::vector<double>> &rotated_diagonal_vectors);


// split matrix of dimension (n, m) into blocks, each of which has dimension (p, q)
// if p (q) does not divides n (m), the original matrix will be automatedly augumented by filling with zeros.
template <typename T>
void split_matrix(Matrix<T> &AA, std::size_t p, std::size_t q, std::vector<std::vector<Matrix<T>>> &S)
{
    std::size_t n = AA.get_rows();
    std::size_t m = AA.get_cols();
    uint64_t nn = n;
    uint64_t mm = m;
    if (nn % p != 0){
        nn += (p - nn % p); 
    }
    if (mm % q != 0){
        mm += (q - mm % q); 
    }
    // std::cout << "augumented (n, m): (" << nn <<", " << mm << ")" << std::endl;
    Matrix<T> A;
    AA.augumented_matrix(nn, mm, A);

    // A.print(10);
    
    S.resize(nn/p);
    for(std::size_t i = 0; i < nn/p; i++){
        // cout << " 1st level loop " << endl;
        S[i].resize(mm/q);
        for (std::size_t j = 0; j < mm/q; j++){
            // cout << " 2nd level loop " << endl;
            Matrix<T> B(p, q);
            for(std::size_t k = 0; k < p; k++){
                // cout << " 3rd level loop " << endl;
                for (std::size_t l = 0; l < q; l ++){
                    // cout << " 4th level loop " << endl;
                    // cout << "(i*p + k, j*q+l)): " << i*p + k << ", " << j*q+l << endl;
                    B.set(k, l, A(i*p + k, j*q+l));
                }
            }
            S[i][j] = B;
        }
    }
}



/*
    Given a vector<vector<vector<T>>>, which comes from decrypt_decode a blocked matrix with size (n, p),
    each blocked submatrix is of size (b_n, b_p), i.e., X[i][j] is the blocked submatrix, return the original 
    (n, p) matrix.
*/

template <typename T>
void join_matrix(std::vector<std::vector<std::vector<T>>> &X, std::size_t n, std::size_t p, std::size_t b_n, std::size_t b_p,  Matrix<T> & A)
{
    A.resize(n, p);
    std::size_t n_r = X.size();
    std::size_t n_c = X[0].size();
    // std::cout << "(n_r, n_c): (" << n_r << ", " << n_c << ")" << std::endl;
    // std::cout << "(n, p): (" << n << ", " << p << ")" << std::endl;
    // std::cout << "(b_n, b_p): (" << b_n << ", " << b_p << ")" << std::endl;
    for (std::size_t i = 0; i < n_r; i++){
        for (std::size_t j = 0; j < n_c; j++){
            Matrix<T> a;
            a = diagonal_decoding(X[i][j], b_n, b_p);
            for(std::size_t k = 0; k < b_n; k++){
                for (std::size_t l = 0; l < b_p; l++){
                    A(b_n*i + k, b_p*j + l) = a(k, l);
                }
            }
        }  
    }
}


#endif // matrix.h