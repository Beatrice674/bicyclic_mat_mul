#ifndef ENCRYPTED_MATRIX_H
#define ENCRYPTED_MATRIX_H


#include "matrix.h"
#include <iostream>
#include <seal/seal.h>


// using namespace std;
// using namespace seal;


//! note that this is not the 2*sqrt() algorithm, 
//! since it  may not hold that 
//! ind[k*i + j] = ind[k*i] + ind[j]
void submatrix_extracting(seal::Decryptor &decryptor, 
                          seal::CKKSEncoder &encoder, seal::Encryptor &encryptor, seal::Evaluator &evaluator, 
                          seal::RelinKeys &relin_keys, seal::GaloisKeys &galois_keys, 
                          seal::Ciphertext &ct_mat, std::size_t n, std::size_t m, std::size_t p, std::size_t q, seal::Ciphertext &ct_r);





// Input two matrices A and B, and the number of slot
template <typename T>
void packing_for_matrix_multiplication(Matrix<T> &A, Matrix<T> &B, 
                                       std::vector<T> &dA, std::vector<T> &dB, 
                                       std::size_t slot_count)
{
    std::size_t n, m, p;
    n = A.get_rows();
    m = A.get_cols();
    if (m != B.get_rows()){
        std::cerr << "packing: #cols of A does not match with #rows of B" << std::endl;
    }
    else{
        p = B.get_cols();
        if ((slot_count < n*m*(std::ceil(double(p)/m) + 1)) || (slot_count < m*p*(std::ceil(double(n)/m)+1))) {
            std::cerr << "packing: the number of slot is not enough for the  (" << n << ", " 
                 << m <<", " <<  p << ") matrix multiplication" << std::endl;    
        }
        else{
            dA.resize(slot_count);
            dB.resize(slot_count);
            std::vector<T> da, db;
            da = A.diagonal_encoding();
            db = B.diagonal_encoding();
            for (std::size_t i = 0; i < n*m; i++){
            /*********************************************************************
            * in fact, -i*n mod n*m = (m-i)*n, which is not greater than (m-1)*n. 
            * after that, the first n*p elements will be taken. 
            * so, j*n*m >= (m-1)*n + n*p.
            * so, j is at least 1 + (p-1)/m.
            * i.e., da will be repeated 1 + (p-1)/m times at least.
            * simply repeat ceil(p/m) + 1 times, always enough.
            **********************************************************************/
                for (std::size_t j = 0; j < std::ceil(double(p)/m)+1; j++){
                    dA[i + n*m*j] = da[i];
                }
            }

            for (std::size_t i = 0; i < p*m; i++){
            // Similarily, we need to repeatedly store db 1 + (n-1)/m times at least.
                for (std::size_t j = 0; j < std::ceil(double(n)/m)+1; j++){
                    dB[i + p*m*j] = db[i];
                }
            }
        } 
    }
}

// Input two matrices A and B, and the number of slot
template <typename T>
void packing_for_matrix_multiplication_v2(Matrix<T> &A, Matrix<T> &B, 
                                       std::vector<T> &dA, std::vector<T> &dB, 
                                       std::size_t slot_count)
{
    std::size_t n, m, p;
    n = A.get_rows();
    m = A.get_cols();
    if (m != B.get_rows()){
        std::cerr << "packing: #cols of A does not match with #rows of B" << std::endl;
    }
    else{
        p = B.get_cols();
        if ((slot_count < n*m*(std::ceil(double(p)/m) )) || (slot_count < m*p*(std::ceil(double(n)/m)))) {
            std::cerr << "packing: the number of slot is not enough for the  (" << n << ", " 
                 << m <<", " <<  p << ") matrix multiplication" << std::endl;    
        }
        else{
            dA.resize(slot_count);
            dB.resize(slot_count);
            std::vector<T> da, db;
            da = A.diagonal_encoding();
            db = B.diagonal_encoding();
            dA=da;
            dB=db;
        } 
    }
}



// Packing the right matrix of size (m, p) to be multiplied by another (n, m) matrix from left,
// returns a vector of size slot_count. 
template <typename T>
void packing_right_matrix(Matrix<T> &B, std::vector<T> &dB, std::size_t n, std::size_t slot_count)
{
    std::size_t m, p;
    m = B.get_rows();
    p = B.get_cols();
        if ((slot_count < n*m*(std::ceil(double(p)/m) + 1)) || (slot_count < m*p*(std::ceil(double(n)/m)+1))) {
            std::cerr << "packing: the number of slot is not enough for the  (" << n << ", " 
                 << m <<", " <<  p << ") matrix multiplication" << std::endl;    
        }
        else{
            // dA.resize(slot_count);
            dB.resize(slot_count);
            std::vector<T> db;
            // da = A.diagonal_encoding();
            db = B.diagonal_encoding();
            // for (std::size_t i = 0; i < n*m; i++){
            // /*********************************************************************
            // * in fact, -i*n mod n*m = (m-i)*n, which is not greater than (m-1)*n. 
            // * after that, the first n*p elements will be taken. 
            // * so, j*n*m >= (m-1)*n + n*p.
            // * so, j is at least 1 + (p-1)/m.
            // * i.e., da will be repeated 1 + (p-1)/m times at least.
            // * simply repeat ceil(p/m) + 1 times, always enough.
            // **********************************************************************/
            //     for (std::size_t j = 0; j < std::ceil(double(p)/m)+1; j++){
            //         dA[i + n*m*j] = da[i];
            //     }
            // }

            for (std::size_t i = 0; i < p*m; i++){
            // Similarily, we need to repeatedly store db 1 + (n-1)/m times at least.
                for (std::size_t j = 0; j < std::ceil(double(n)/m)+1; j++){
                    dB[i + p*m*j] = db[i];
                }
            }
        } 
}

template <typename T>
void packing_right_matrix_v2(Matrix<T> &B, std::vector<T> &dB, std::size_t n, std::size_t slot_count)
{
    std::size_t m, p;
    m = B.get_rows();
    p = B.get_cols();
        if ((slot_count < n*m*(std::ceil(double(p)/m) )) || (slot_count < m*p*(std::ceil(double(n)/m)))) {
            std::cerr << "packing: the number of slot is not enough for the  (" << n << ", " 
                 << m <<", " <<  p << ") matrix multiplication" << std::endl;    
        }
        else{
            // dA.resize(slot_count);
            dB.resize(slot_count);
            std::vector<T> db;
            // da = A.diagonal_encoding();
            db = B.diagonal_encoding();
            for (std::size_t i = 0; i < p*m; i++){
            // Similarily, we need to repeatedly store db 1 + (n-1)/m times at least.
                for (std::size_t j = 0; j < std::ceil(double(n)/m); j++){
                    dB[i + p*m*j] = db[i];
                }
            }
        } 
}

// Packing the left matrix of size (n, m) to be multiplied by another (m, p) matrix from right,
// returns a vector of size slot_count. 
template <typename T>
void packing_left_matrix(Matrix<T> &A, std::vector<T> &dA, std::size_t p, std::size_t slot_count)
{
    std::size_t n, m;
    n = A.get_rows();
    m = A.get_cols();

        if ((slot_count < n*m*(std::ceil(double(p)/m) + 1)) || (slot_count < m*p*(std::ceil(double(n)/m)+1))) {
            std::cerr << "packing: the number of slot is not enough for the  (" << n << ", " 
                 << m <<", " <<  p << ") matrix multiplication" << std::endl;    
        }
        else{
            dA.resize(slot_count);
            // dB.resize(slot_count);
            std::vector<T> da;
            da = A.diagonal_encoding();
            // db = B.diagonal_encoding();
            for (std::size_t i = 0; i < n*m; i++){
            /*********************************************************************
            * in fact, -i*n mod n*m = (m-i)*n, which is not greater than (m-1)*n. 
            * after that, the first n*p elements will be taken. 
            * so, j*n*m >= (m-1)*n + n*p.
            * so, j is at least 1 + (p-1)/m.
            * i.e., da will be repeated 1 + (p-1)/m times at least.
            * simply repeat ceil(p/m) + 1 times, always enough.
            **********************************************************************/
                for (std::size_t j = 0; j < std::ceil(double(p)/m)+1; j++){
                    dA[i + n*m*j] = da[i];
                }
            }

            // for (std::size_t i = 0; i < p*m; i++){
            // // Similarily, we need to repeatedly store db 1 + (n-1)/m times at least.
            //     for (std::size_t j = 0; j < std::ceil(double(n)/m)+1; j++){
            //         dB[i + p*m*j] = db[i];
            //     }
            // }
        } 
}

template <typename T>
void packing_left_matrix_v2(Matrix<T> &A, std::vector<T> &dA, std::size_t p, std::size_t slot_count)
{
    std::size_t n, m;
    n = A.get_rows();
    m = A.get_cols();

        if ((slot_count < n*m*(std::ceil(double(p)/m))) || (slot_count < m*p*(std::ceil(double(n)/m)))) {
            std::cerr << "packing: the number of slot is not enough for the  (" << n << ", " 
                 << m <<", " <<  p << ") matrix multiplication" << std::endl;    
        }
        else{
            dA.resize(slot_count);
            // dB.resize(slot_count);
            std::vector<T> da;
            da = A.diagonal_encoding();
            // dA=da;
            for (std::size_t i = 0; i < n*m; i++){
            /*********************************************************************
            * in fact, -i*n mod n*m = (m-i)*n, which is not greater than (m-1)*n. 
            * after that, the first n*p elements will be taken. 
            * so, j*n*m >= (m-1)*n + n*p.
            * so, j is at least 1 + (p-1)/m.
            * i.e., da will be repeated 1 + (p-1)/m times at least.
            * simply repeat ceil(p/m) + 1 times, always enough.
            **********************************************************************/
                for (std::size_t j = 0; j < std::ceil(double(p)/m); j++){
                    dA[i + n*m*j] = da[i];
                }
            }
        } 
}







// The following only works for m = max(n, m, p) for the moment
void encrypted_matrix_multiplication(//seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, 
                                     seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, 
                                     seal::GaloisKeys &galois_keys, uint64_t n, uint64_t m, uint64_t p,
                                     seal::Ciphertext &encrypted_A, seal::Ciphertext &encrypted_B, 
                                     seal::Ciphertext &encrypted_C);

void rotate_cipher_vector(seal::CKKSEncoder &encoder, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys,seal::RelinKeys &relin_keys,
                          int64_t move_step, int64_t slot_all,
                          seal::Ciphertext &cipher, seal::Ciphertext &destination);

void encrypted_matrix_multiplication_v2(seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, 
                                     seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, 
                                     seal::GaloisKeys &galois_keys, uint64_t n, uint64_t m, uint64_t p,
                                     seal::Ciphertext &encrypted_A, seal::Ciphertext &encrypted_B, 
                                     seal::Ciphertext &encrypted_C);


// Input A and B (two matrices to be multiplied), it first splits the two matrices into blocks, then encodes
// and encrypts each block into a ciphertext. The structure of the output is vector<vector<Ciphertext>>.
// NOTE: the resulting ciphertext is for matrix multiplication, so the plainvector is repeated serveral times.
void encode_encrypt_block_matrix(seal::CKKSEncoder &encoder, seal::Encryptor &encryptor,  double scale, Matrix<double> &A, 
                                 Matrix<double> &B, std::size_t b_n, std::size_t b_m, std::size_t b_p, 
                                 std::vector<std::vector<seal::Ciphertext>> &A_encrypted, 
                                 std::vector<std::vector<seal::Ciphertext>> &B_encrypted);
void encode_encrypt_block_matrix_v2(seal::CKKSEncoder &encoder, seal::Encryptor &encryptor,  double scale, 
                                Matrix<double> &A, Matrix<double> &B, std::size_t b_n, std::size_t b_m, std::size_t b_p, 
                                 std::vector<std::vector<seal::Ciphertext>> &A_encrypted, 
                                 std::vector<std::vector<seal::Ciphertext>> &B_encrypted);


void encrypted_block_matrix_multiplication(//seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, 
                                     seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, 
                                     seal::GaloisKeys &galois_keys, uint64_t n, uint64_t m, uint64_t p,
                                     std::vector<std::vector<seal::Ciphertext>> &encrypted_A, 
                                     std::vector<std::vector<seal::Ciphertext>> &encrypted_B, 
                                     std::vector<std::vector<seal::Ciphertext>> &encrypted_C);
void encrypted_block_matrix_multiplication_v2(seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, 
                                     seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, 
                                     seal::GaloisKeys &galois_keys, uint64_t n, uint64_t m, uint64_t p,
                                     std::vector<std::vector<seal::Ciphertext>> &encrypted_A, 
                                     std::vector<std::vector<seal::Ciphertext>> &encrypted_B, 
                                     std::vector<std::vector<seal::Ciphertext>> &encrypted_C);







#endif // encrypted_matrix.h