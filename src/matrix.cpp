#include "utils.h"
#include "matrix.h"


#include <vector>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <seal/util/numth.h>


// using namespace seal;
using namespace seal::util;
// using namespace std;



void random_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A){
    A.resize(n, m);
    for (size_t i = 0; i < m; i++){
        for(size_t j = 0; j < n; j++){
            A.set(j, i, pow(-1, i+j)*rand()/(pow(2, 30)));
        }
    }
}
void order_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A){
    A.resize(n, m);
    for (size_t i = 0; i < m; i++){
        for(size_t j = 0; j < n; j++){
            A.set(j, i, double(i*n+j));
        }
    }
}






// Some useful functions related to the matrix class






void submatrix_transformation(std::size_t n, std::size_t m, std::size_t k, std::size_t l, Matrix<double> &T)
{
    if ((n < k) || (m < l)){
        std::cerr << "submatrix_transformation: " << n << " must be >= " << k 
                  << " and " << m << " must be >= " << l << std::endl;
        return;
    }
    
    // tuple<uint64_t, int64_t, int64_t> gcd;
    auto gcd = seal::util::xgcd(n, m);
    uint64_t g = get<0>(gcd);
    if ((g != 1) || (std::gcd(k, l) != 1)){
        std::cerr << "submatrix_transformation: the dimension (" << n << ", " << m 
             << ") or (" << k << ", " << l << ") are not coprime " << std::endl;
        return;
    }


    int64_t s = get<1>(gcd);
    int64_t t = get<2>(gcd);
    T.resize(n*m, n*m);
    for (std::size_t r = 0; r < k * l; r++){
        int64_t rk = mod(r, k);
        int64_t rl = mod(r, l);
        std::size_t j = mod(rk * t * m + rl * s * n, n * m);
        T.set(r, j, 1.);
    }
}  




void rotated_diagonal_vectors_for_submatrix(std::size_t n, std::size_t m, std::size_t p, std::size_t q, std::size_t slot_count,
                                            Matrix<double> &T, //the transformation matrix
                                            std::vector<std::size_t> &ind, // the indices of nonzero diagonal vectors of T
                                            std::vector<std::vector<double>> &rotated_diagonal_vectors)
{
    // Matrix<double> T; //the transformation matrix
    // vector<size_t> ind; // the indices of nonzero diagonal vectors of T
    std::size_t n_r; // the number of different values of r
    // submatrix_transformation(n, m, p, q, T);
    // nonzero_index_set(n, m, p, q, ind);
    n_r = ind.size();
    std::size_t r = ceil(sqrt(double(n_r)));
    std::size_t l, k;
    k = r;
    if (is_square(n_r)){
        l = r;
    }
    else{
        l = n_r/r + 1;
    }
    
    std::cout << "(k, l)= (" << k << ", " << l <<")" << std::endl;

    for(std::size_t i = 0; i < l; i++){
        std::cout << " i = " << i << std::endl; 
        for(std::size_t j = 0; j < k; j++){
            std::cout << " j = " << j << std::endl;
            if ((i < l-1) || ((i == l -1) && ((l-1)*k + j < n_r))){
                std::vector<double> d, t;
                d = T.diagonal_vector(ind[k*i+j]);
                std::cout << "before rotate: ";
                print_vector(d, 15);
                rotate_vector(d, -1*(int64_t)(ind[k*i]), slot_count, t);
                std::cout << "after rotate " << -1*(int64_t)(ind[k*i]) << ": ";
                print_vector(t, 15);
                rotated_diagonal_vectors.push_back(t);
            }  
            else{
                std::vector<double> t(slot_count, 10e-10); // CKKS::Plaintext does not allow an all-zero vector.
                rotated_diagonal_vectors.push_back(t);
            }
        }
    }
}



