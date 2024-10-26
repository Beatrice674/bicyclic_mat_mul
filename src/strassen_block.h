#include <iostream>
#include <seal/seal.h>
#include "matrix.h"
#include "utils.h"
#include "coprime_enc_cipher.h"

using namespace std;
using namespace seal;

class matrix_block{
private:
    size_t n;//rows
    size_t m;//cols
    size_t b_rows;//block rows
    size_t b_cols;//block cols
    vector<vector<seal::Ciphertext>> block_matrix_cipher;//block matrix
public:
    matrix_block();
    matrix_block(vector<vector<seal::Ciphertext>> data,size_t n,size_t m,size_t b_rows,size_t b_cols);
    matrix_block(Matrix<double> data,double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    matrix_block(Matrix<double> data,size_t n,size_t m,size_t repeate,double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    matrix_block extract_matrix(size_t start_row,size_t end_row,size_t start_col,size_t end_col);
    void dec_block_matrix_jiang(Matrix<double> &destination,seal::CKKSEncoder &encoder,seal::Decryptor &decryptor);
    void dec_block_matrix_coprime(Matrix<double> &destination,size_t s_n,size_t s_m,seal::CKKSEncoder &encoder,seal::Decryptor &decryptor);
    ~matrix_block();
    size_t get_rows();
    size_t get_cols();
    size_t get_block_rows();
    size_t get_block_cols();
    vector<vector<seal::Ciphertext>> get_cipher_matrix();
    vector<vector<seal::Ciphertext>> get_sub_cipher_matrix(seal::Evaluator &evaluator,seal::CKKSEncoder &encoder);
};

inline matrix_block::matrix_block()
{
    n=0;
    m=0;
}

inline matrix_block::matrix_block(vector<vector<seal::Ciphertext>> data, size_t n, size_t m, size_t b_rows, size_t b_cols)
{
    block_matrix_cipher=data;
    this->n=n;
    this->m=m;
    this->b_rows=b_rows;
    this->b_cols=b_cols;
}

// Jiang block
inline matrix_block::matrix_block(Matrix<double> data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    std::size_t slot_count = encoder.slot_count();
    std::size_t d=sqrt(slot_count);
    
    std::vector<std::vector<Matrix<double>>> Sdata;
    split_matrix(data, d, d, Sdata);

    n=data.get_rows();
    m=data.get_cols();
    b_rows = Sdata.size();
    b_cols = Sdata[0].size();

    block_matrix_cipher.resize(b_rows);

    for(std::size_t i = 0; i < b_rows; i++){
        // A_encrypted[i].resize(b_cols_a);
        for(std::size_t j = 0; j < b_cols; j++){
            //! to be written!
            std::vector<double> dA;
            dA=Sdata[i][j].flatten_matrix_to_rows_vector();
            seal::Ciphertext tmp;
            encode_encrypt(encoder, encryptor, dA, scale, tmp);
            block_matrix_cipher[i].push_back(tmp);
        }
    }
}

inline matrix_block::matrix_block(Matrix<double> data, size_t s_n, size_t s_m,size_t repeate, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    std::vector<std::vector<Matrix<double>>> Sdata;
    split_matrix(data, s_n, s_m, Sdata);

    
    n=data.get_rows();
    m=data.get_cols();
    b_rows = Sdata.size();
    b_cols = Sdata[0].size();

    // cout<<"n:"<<n<<"  m:"<<m<<"  b_rows"<<b_rows<<"  b_cols"<<b_cols<<endl;

    block_matrix_cipher.resize(b_rows);

    for(std::size_t i = 0; i < b_rows; i++){
        // A_encrypted[i].resize(b_cols_a);
        for(std::size_t j = 0; j < b_cols; j++){
            //! to be written!
            std::vector<double> dA,da;
            da=Sdata[i][j].diagonal_encoding();
            dA.resize(encoder.slot_count());
            //repeate
            for (size_t s_i = 0; s_i < s_n * s_m; s_i++)
            {
                for (size_t s_j = 0; s_j < repeate; s_j++)
                {
                    dA[s_i + s_n * s_m * s_j] = da[s_i];
                }
            }
            seal::Ciphertext tmp;
            encode_encrypt(encoder, encryptor, dA, scale, tmp);
            block_matrix_cipher[i].push_back(tmp);
        }
    }
}

inline matrix_block matrix_block::extract_matrix(size_t start_row, size_t end_row, size_t start_col, size_t end_col)
{
    vector<vector<seal::Ciphertext>> submatrix;

    for(size_t i=start_row;i<end_row;i++){
        vector<seal::Ciphertext> cipher_tmp;
        for(size_t j=start_col;j<end_col;j++){
            cipher_tmp.push_back(block_matrix_cipher[i][j]);
        }
        submatrix.push_back(cipher_tmp);
    }
    return matrix_block(submatrix,n,m,end_row-start_row,end_col-start_col);
}

inline void matrix_block::dec_block_matrix_jiang(Matrix<double> &destination, seal::CKKSEncoder &encoder, seal::Decryptor &decryptor)
{
    
    size_t d=sqrt(encoder.slot_count());
    destination.resize(b_rows*d,b_cols*d);
    Matrix<double> tmp_matrix;
    tmp_matrix.resize(d,d);
    for(size_t i=0;i<b_rows;i++){
        for(size_t j=0;j<b_cols;j++){
            cipher_matrix_jiang dec_tmp(block_matrix_cipher[i][j],d,d,d);
            dec_tmp.dec_matrix_cipher(tmp_matrix,encoder,decryptor);
            destination.Merage_matrix(i*d,(i+1)*d,j*d,(j+1)*d,tmp_matrix);
        }
    }
    destination.resize(n,m);
    // cout<<"n:"<<n<<"  m:"<<m<<endl;
}

inline void matrix_block::dec_block_matrix_coprime(Matrix<double> &destination,size_t s_n,size_t s_m, seal::CKKSEncoder &encoder, seal::Decryptor &decryptor)
{
    // destination.resize(n,m);
    // size_t d=sqrt(encoder.slot_count());
    // cout<<n<<"  "<<m<<endl;
    destination.resize(b_rows*s_n,b_cols*s_m);
    Matrix<double> tmp_matrix;
    tmp_matrix.resize(s_n,s_m);
    for(size_t i=0;i<b_rows;i++){
        for(size_t j=0;j<b_cols;j++){
            cipher_matrix_coprime dec_tmp(block_matrix_cipher[i][j],s_n,s_m);
            dec_tmp.dec_matrix_cipher(tmp_matrix,encoder,decryptor);
            // tmp_matrix.print(6,6);
            destination.Merage_matrix(i*s_n,(i+1)*s_n,j*s_m,(j+1)*s_m,tmp_matrix);
        }
    }
}

matrix_block::~matrix_block()
{
    ;
}

inline size_t matrix_block::get_rows()
{
    return n;
}

inline size_t matrix_block::get_cols()
{
    return m;
}

inline size_t matrix_block::get_block_rows()
{
    return b_rows;
}

inline size_t matrix_block::get_block_cols()
{
    return b_cols;
}

inline vector<vector<seal::Ciphertext>> matrix_block::get_cipher_matrix()
{
    return block_matrix_cipher;
}

inline vector<vector<seal::Ciphertext>> matrix_block::get_sub_cipher_matrix(seal::Evaluator &evaluator, seal::CKKSEncoder &encoder)
{
    vector<vector<seal::Ciphertext>> sub_block_matrix_cipher=block_matrix_cipher;

    vector<double> vec_tmp(encoder.slot_count(),-1);
    Plaintext plain_tmp;
    encoder.encode(vec_tmp,1,plain_tmp);
    evaluator.mod_switch_to_inplace(plain_tmp,sub_block_matrix_cipher[0][0].parms_id());

    for(size_t i=0;i<b_rows;i++){
        for(size_t j=0;j<b_cols;j++){
            evaluator.multiply_plain_inplace(sub_block_matrix_cipher[i][j],plain_tmp);
        }
    }
    return sub_block_matrix_cipher;
}

// Add block matrix
matrix_block ADD(matrix_block matrix_1,matrix_block matrix_2,seal::Evaluator &evaluator){
    if(matrix_1.get_block_rows()!=matrix_2.get_block_rows() ||matrix_1.get_block_cols()!=matrix_2.get_block_cols())
    {
        cerr<<"Matrix blocks not equal"<<endl;
    }
    size_t block_rows=matrix_1.get_block_rows();
    size_t block_cols=matrix_1.get_block_cols();
    vector<vector<seal::Ciphertext>> block_matrix_cipher,block_matrix_cipher_1,block_matrix_cipher_2;
    block_matrix_cipher.resize(block_rows,std::vector<seal::Ciphertext>(block_cols));
    block_matrix_cipher_1=matrix_1.get_cipher_matrix();
    block_matrix_cipher_2=matrix_2.get_cipher_matrix();

    for(size_t i=0;i<block_rows;i++){
        for(size_t j=0;j<block_cols;j++){
            if(block_matrix_cipher_1[i][j].scale()!=block_matrix_cipher_2[i][j].scale()){
                block_matrix_cipher_1[i][j].scale()=block_matrix_cipher_2[i][j].scale();
            }
            evaluator.add(block_matrix_cipher_1[i][j],block_matrix_cipher_2[i][j],block_matrix_cipher[i][j]);
        }
    }

    return matrix_block(block_matrix_cipher,matrix_1.get_rows(),matrix_1.get_cols(),block_rows,block_cols);
}

// Subtract block matrix
matrix_block SUB(matrix_block matrix_1,matrix_block matrix_2,seal::Evaluator &evaluator,seal::CKKSEncoder &encoder){
    if(matrix_1.get_block_rows()!=matrix_2.get_block_rows() ||matrix_1.get_block_cols()!=matrix_2.get_block_cols())
    {
        cerr<<"Matrix blocks not equal"<<endl;
    }

    size_t block_rows=matrix_1.get_block_rows();
    size_t block_cols=matrix_1.get_block_cols();
    vector<vector<seal::Ciphertext>> block_matrix_cipher,block_matrix_cipher_1,block_matrix_cipher_2;
    block_matrix_cipher.resize(block_rows,std::vector<seal::Ciphertext>(block_cols));
    block_matrix_cipher_1=matrix_1.get_cipher_matrix();
    // cout<<"get matrix1 success"<<endl;
    block_matrix_cipher_2=matrix_2.get_sub_cipher_matrix(evaluator,encoder);
    // cout<<"get matrix2 success"<<endl;

    // cout<<block_matrix_cipher_1.size()<<"  "<<block_matrix_cipher_1[0].size()<<endl;
    // cout<<block_matrix_cipher_2.size()<<"  "<<block_matrix_cipher_2[0].size()<<endl;

    for(size_t i=0;i<block_rows;i++){
        for(size_t j=0;j<block_cols;j++){
            if(block_matrix_cipher_1[i][j].scale()!=block_matrix_cipher_2[i][j].scale()){
                block_matrix_cipher_1[i][j].scale()=block_matrix_cipher_2[i][j].scale();
            }
            evaluator.add(block_matrix_cipher_1[i][j],block_matrix_cipher_2[i][j],block_matrix_cipher[i][j]);
        }
    }

    return matrix_block(block_matrix_cipher,matrix_1.get_rows(),matrix_1.get_cols(),block_rows,block_cols);
}

// Merage block matrix
matrix_block Merage_matrix(matrix_block c11,matrix_block c12,matrix_block c21,matrix_block c22)
{
    vector<vector<seal::Ciphertext>> block_matrix_cipher,block_matrix_cipher_11,block_matrix_cipher_12,block_matrix_cipher_21,block_matrix_cipher_22;
    block_matrix_cipher_11=c11.get_cipher_matrix();
    block_matrix_cipher_12=c12.get_cipher_matrix();
    block_matrix_cipher_21=c21.get_cipher_matrix();
    block_matrix_cipher_22=c22.get_cipher_matrix();

    for(size_t i=0;i<block_matrix_cipher_11.size();i++){
        for(size_t j=0;j<block_matrix_cipher_12[i].size();j++){
            block_matrix_cipher_11[i].push_back(block_matrix_cipher_12[i][j]);
        }
    }

    for(size_t i=0;i<block_matrix_cipher_21.size();i++){
        for(size_t j=0;j<block_matrix_cipher_22[i].size();j++){
            block_matrix_cipher_21[i].push_back(block_matrix_cipher_22[i][j]);
        }
    }

    for(size_t i=0; i<block_matrix_cipher_21.size();i++){
        block_matrix_cipher_11.push_back(block_matrix_cipher_21[i]);
    }

    return matrix_block(block_matrix_cipher_11,c11.get_rows(),c11.get_cols(),block_matrix_cipher_11.size(),block_matrix_cipher_11[0].size());
}

//Strassen block jinag multiply
matrix_block matrix_block_jiang_mul(matrix_block encrypted_A,matrix_block encrypted_B,double scale, seal::CKKSEncoder& encoder,
                                    seal::Evaluator &evaluator,seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys)
{
    seal::Ciphertext cipher_A,cipher_B;
    cipher_A=encrypted_A.get_cipher_matrix()[0][0];
    cipher_B=encrypted_B.get_cipher_matrix()[0][0];

    size_t d= sqrt(encoder.slot_count());
    cipher_matrix_jiang matrix_A(cipher_A,d,d,d);
    cipher_matrix_jiang matrix_B(cipher_B,d,d,d);
    cipher_matrix_jiang matrix_C;

    matrix_A.change_cipher_to_matrix_A(scale,encoder,evaluator,galois_keys,relin_keys);
    matrix_B.change_cipher_to_matrix_B(scale,encoder,evaluator,galois_keys,relin_keys);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,galois_keys,matrix_A,matrix_B,matrix_C);

    vector<vector<seal::Ciphertext>> cipher_matrix(1,vector<seal::Ciphertext>(1));
    cipher_matrix[0][0]=matrix_C.get_cipher_matrix();
    return matrix_block(cipher_matrix,encrypted_A.get_rows(),encrypted_B.get_cols(),encrypted_A.get_block_rows(),encrypted_B.get_block_cols());
}

// Strassen block by jiang
matrix_block Strassen_block_jiang(matrix_block encrypted_A,matrix_block encrypted_B,double scale, seal::CKKSEncoder& encoder,
                                  seal::Evaluator &evaluator,seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys)
{
    size_t block_rows_A=encrypted_A.get_block_rows();
    size_t block_cols_A=encrypted_A.get_block_cols();
    size_t block_rows_B=encrypted_B.get_block_rows();
    size_t block_cols_B=encrypted_B.get_block_cols();

    // cout<<encrypted_A.get_cols()<<"  "<<encrypted_A.get_rows()<<endl;
    matrix_block A11,A12,A21,A22;
    matrix_block B11,B12,B21,B22;
    matrix_block M1,M2,M3,M4,M5,M6,M7;
    matrix_block C11,C12,C21,C22;

    A11=encrypted_A.extract_matrix(0,block_rows_A/2,0,block_cols_A/2);
    A12=encrypted_A.extract_matrix(0,block_rows_A/2,block_cols_A/2,block_cols_A);
    A21=encrypted_A.extract_matrix(block_rows_A/2,block_rows_A,0,block_cols_A/2);
    A22=encrypted_A.extract_matrix(block_rows_A/2,block_rows_A,block_cols_A/2,block_cols_A);

    B11=encrypted_B.extract_matrix(0,block_rows_B/2,0,block_cols_B/2);
    B12=encrypted_B.extract_matrix(0,block_rows_B/2,block_cols_B/2,block_cols_B);
    B21=encrypted_B.extract_matrix(block_rows_B/2,block_rows_B,0,block_cols_B/2);
    B22=encrypted_B.extract_matrix(block_rows_B/2,block_rows_B,block_cols_B/2,block_cols_B);

    // stop
    if(block_rows_A/2<=1 && block_cols_A/2<=1){
        M1=matrix_block_jiang_mul(ADD(A11,A22,evaluator),ADD(B11,B22,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
        M2=matrix_block_jiang_mul(ADD(A21,A22,evaluator),B11,scale,encoder,evaluator,galois_keys,relin_keys);
        M3=matrix_block_jiang_mul(A11,SUB(B12,B22,evaluator,encoder),scale,encoder,evaluator,galois_keys,relin_keys);
        M4=matrix_block_jiang_mul(A22,SUB(B21,B11,evaluator,encoder),scale,encoder,evaluator,galois_keys,relin_keys);
        M5=matrix_block_jiang_mul(ADD(A11,A12,evaluator),B22,scale,encoder,evaluator,galois_keys,relin_keys);
        M6=matrix_block_jiang_mul(SUB(A21,A11,evaluator,encoder),ADD(B11,B12,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
        M7=matrix_block_jiang_mul(SUB(A12,A22,evaluator,encoder),ADD(B21,B22,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
    }
    else{
        M1=Strassen_block_jiang(ADD(A11,A22,evaluator),ADD(B11,B22,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
        M2=Strassen_block_jiang(ADD(A21,A22,evaluator),B11,scale,encoder,evaluator,galois_keys,relin_keys);
        M3=Strassen_block_jiang(A11,SUB(B12,B22,evaluator,encoder),scale,encoder,evaluator,galois_keys,relin_keys);
        M4=Strassen_block_jiang(A22,SUB(B21,B11,evaluator,encoder),scale,encoder,evaluator,galois_keys,relin_keys);
        M5=Strassen_block_jiang(ADD(A11,A12,evaluator),B22,scale,encoder,evaluator,galois_keys,relin_keys);
        M6=Strassen_block_jiang(SUB(A21,A11,evaluator,encoder),ADD(B11,B12,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
        M7=Strassen_block_jiang(SUB(A12,A22,evaluator,encoder),ADD(B21,B22,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
    }

    C11=ADD(ADD(M1,M4,evaluator),SUB(M7,M5,evaluator,encoder),evaluator);
    C12=ADD(M3,M5,evaluator);
    C21=ADD(M2,M4,evaluator);
    C22=ADD(SUB(M1,M2,evaluator,encoder),ADD(M3,M6,evaluator),evaluator);

    return  Merage_matrix(C11,C12,C21,C22);
}

// Naive block by jiang
matrix_block Naive_block_jiang(matrix_block encrypted_A,matrix_block encrypted_B,double scale, seal::CKKSEncoder& encoder,
                               seal::Evaluator &evaluator,seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys)
{
    size_t block_rows_A=encrypted_A.get_block_rows();
    size_t block_cols_A=encrypted_A.get_block_cols();
    size_t block_rows_B=encrypted_B.get_block_rows();
    size_t block_cols_B=encrypted_B.get_block_cols();

    size_t d= sqrt(encoder.slot_count());

    if(block_cols_A!=block_rows_B){
        cerr<<"Matrix blocks not equal"<<endl;
    }

    vector<vector<seal::Ciphertext>> cipher_matrix(block_rows_A,vector<seal::Ciphertext>(block_cols_B));
    seal::Ciphertext cipher_tmp;
    for(size_t i=0;i<block_rows_A;i++){
        for(size_t j=0;j<block_cols_B;j++){
            vector<seal::Ciphertext> cipher_vector;
            for(size_t k=0;k<block_cols_A;k++){
                cipher_matrix_jiang matrix_A(encrypted_A.get_cipher_matrix()[i][k],d,d,d);
                cipher_matrix_jiang matrix_B(encrypted_B.get_cipher_matrix()[k][j],d,d,d);
                cipher_matrix_jiang matrix_C;

                matrix_A.change_cipher_to_matrix_A(scale,encoder,evaluator,galois_keys,relin_keys);
                matrix_B.change_cipher_to_matrix_B(scale,encoder,evaluator,galois_keys,relin_keys);
                encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,galois_keys,matrix_A,matrix_B,matrix_C);
                cipher_vector.push_back(matrix_C.get_cipher_matrix());
            }
            evaluator.add_many(cipher_vector,cipher_tmp);
            cipher_matrix[i][j]=cipher_tmp;
        }
    }
    return matrix_block(cipher_matrix,encrypted_A.get_rows(),encrypted_B.get_cols(),block_rows_A,block_cols_B);
}


//Strassen block coprime multiply
matrix_block matrix_block_coprime_mul(size_t s_n,size_t s_m,size_t s_p,matrix_block encrypted_A,matrix_block encrypted_B,double scale, seal::CKKSEncoder& encoder,
                                      seal::Evaluator &evaluator,seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys)
{
    seal::Ciphertext cipher_A,cipher_B;
    cipher_A=encrypted_A.get_cipher_matrix()[0][0];
    cipher_B=encrypted_B.get_cipher_matrix()[0][0];

    cipher_matrix_coprime matrix_A(cipher_A,s_n,s_m);
    cipher_matrix_coprime matrix_B(cipher_B,s_m,s_p);
    cipher_matrix_coprime matrix_C;

    encrypted_coprime_matrix_multiplication(encoder,scale,evaluator,relin_keys,galois_keys,matrix_A,matrix_B,matrix_C);

    vector<vector<seal::Ciphertext>> cipher_matrix(1,vector<seal::Ciphertext>(1));
    cipher_matrix[0][0]=matrix_C.get_cipher_matrix();
    return matrix_block(cipher_matrix,encrypted_A.get_rows(),encrypted_B.get_cols(),encrypted_A.get_block_rows(),encrypted_B.get_block_cols());
}

//Strassen block coprime multiply imporve log
matrix_block matrix_block_coprime_mul_imporve_log(size_t s_n,size_t s_m,size_t s_p,matrix_block encrypted_A,matrix_block encrypted_B,double scale, seal::CKKSEncoder& encoder,
                                      seal::Evaluator &evaluator,seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys)
{
    seal::Ciphertext cipher_A,cipher_B;
    cipher_A=encrypted_A.get_cipher_matrix()[0][0];
    cipher_B=encrypted_B.get_cipher_matrix()[0][0];

    cipher_matrix_coprime matrix_A(cipher_A,s_n,s_m);
    cipher_matrix_coprime matrix_B(cipher_B,s_m,s_p);
    cipher_matrix_coprime matrix_C;

    encrypted_coprime_matrix_multiplication_imporve_log(encoder,scale,evaluator,relin_keys,galois_keys,matrix_A,matrix_B,matrix_C);

    vector<vector<seal::Ciphertext>> cipher_matrix(1,vector<seal::Ciphertext>(1));
    cipher_matrix[0][0]=matrix_C.get_cipher_matrix();
    return matrix_block(cipher_matrix,encrypted_A.get_rows(),encrypted_B.get_cols(),encrypted_A.get_block_rows(),encrypted_B.get_block_cols());
}

// Strassen block by coprime
matrix_block Strassen_block_coprime(size_t s_n,size_t s_m,size_t s_p,matrix_block encrypted_A,matrix_block encrypted_B,double scale, seal::CKKSEncoder& encoder,
                                    seal::Evaluator &evaluator,seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys)
{
    size_t block_rows_A=encrypted_A.get_block_rows();
    size_t block_cols_A=encrypted_A.get_block_cols();
    size_t block_rows_B=encrypted_B.get_block_rows();
    size_t block_cols_B=encrypted_B.get_block_cols();

    matrix_block A11,A12,A21,A22;
    matrix_block B11,B12,B21,B22;
    matrix_block M1,M2,M3,M4,M5,M6,M7;
    matrix_block C11,C12,C21,C22;

    A11=encrypted_A.extract_matrix(0,block_rows_A/2,0,block_cols_A/2);
    A12=encrypted_A.extract_matrix(0,block_rows_A/2,block_cols_A/2,block_cols_A);
    A21=encrypted_A.extract_matrix(block_rows_A/2,block_rows_A,0,block_cols_A/2);
    A22=encrypted_A.extract_matrix(block_rows_A/2,block_rows_A,block_cols_A/2,block_cols_A);

    B11=encrypted_B.extract_matrix(0,block_rows_B/2,0,block_cols_B/2);
    B12=encrypted_B.extract_matrix(0,block_rows_B/2,block_cols_B/2,block_cols_B);
    B21=encrypted_B.extract_matrix(block_rows_B/2,block_rows_B,0,block_cols_B/2);
    B22=encrypted_B.extract_matrix(block_rows_B/2,block_rows_B,block_cols_B/2,block_cols_B);

    // stop
    if(block_rows_A/2<=1 && block_cols_A/2<=1){
        M1=matrix_block_coprime_mul(s_n,s_m,s_p,ADD(A11,A22,evaluator),ADD(B11,B22,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
        M2=matrix_block_coprime_mul(s_n,s_m,s_p,ADD(A21,A22,evaluator),B11,scale,encoder,evaluator,galois_keys,relin_keys);
        M3=matrix_block_coprime_mul(s_n,s_m,s_p,A11,SUB(B12,B22,evaluator,encoder),scale,encoder,evaluator,galois_keys,relin_keys);
        M4=matrix_block_coprime_mul(s_n,s_m,s_p,A22,SUB(B21,B11,evaluator,encoder),scale,encoder,evaluator,galois_keys,relin_keys);
        M5=matrix_block_coprime_mul(s_n,s_m,s_p,ADD(A11,A12,evaluator),B22,scale,encoder,evaluator,galois_keys,relin_keys);
        M6=matrix_block_coprime_mul(s_n,s_m,s_p,SUB(A21,A11,evaluator,encoder),ADD(B11,B12,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
        M7=matrix_block_coprime_mul(s_n,s_m,s_p,SUB(A12,A22,evaluator,encoder),ADD(B21,B22,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
    }
    else{
        M1=Strassen_block_coprime(s_n,s_m,s_p,ADD(A11,A22,evaluator),ADD(B11,B22,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
        M2=Strassen_block_coprime(s_n,s_m,s_p,ADD(A21,A22,evaluator),B11,scale,encoder,evaluator,galois_keys,relin_keys);
        M3=Strassen_block_coprime(s_n,s_m,s_p,A11,SUB(B12,B22,evaluator,encoder),scale,encoder,evaluator,galois_keys,relin_keys);
        M4=Strassen_block_coprime(s_n,s_m,s_p,A22,SUB(B21,B11,evaluator,encoder),scale,encoder,evaluator,galois_keys,relin_keys);
        M5=Strassen_block_coprime(s_n,s_m,s_p,ADD(A11,A12,evaluator),B22,scale,encoder,evaluator,galois_keys,relin_keys);
        M6=Strassen_block_coprime(s_n,s_m,s_p,SUB(A21,A11,evaluator,encoder),ADD(B11,B12,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
        M7=Strassen_block_coprime(s_n,s_m,s_p,SUB(A12,A22,evaluator,encoder),ADD(B21,B22,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
    }

    C11=ADD(ADD(M1,M4,evaluator),SUB(M7,M5,evaluator,encoder),evaluator);
    C12=ADD(M3,M5,evaluator);
    C21=ADD(M2,M4,evaluator);
    C22=ADD(SUB(M1,M2,evaluator,encoder),ADD(M3,M6,evaluator),evaluator);

    return  Merage_matrix(C11,C12,C21,C22);
}

// Strassen block imporve log by coprime
matrix_block Strassen_block_coprime_imporve_log(size_t s_n,size_t s_m,size_t s_p,matrix_block encrypted_A,matrix_block encrypted_B,double scale, seal::CKKSEncoder& encoder,
                                    seal::Evaluator &evaluator,seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys)
{
    size_t block_rows_A=encrypted_A.get_block_rows();
    size_t block_cols_A=encrypted_A.get_block_cols();
    size_t block_rows_B=encrypted_B.get_block_rows();
    size_t block_cols_B=encrypted_B.get_block_cols();

    matrix_block A11,A12,A21,A22;
    matrix_block B11,B12,B21,B22;
    matrix_block M1,M2,M3,M4,M5,M6,M7;
    matrix_block C11,C12,C21,C22;

    A11=encrypted_A.extract_matrix(0,block_rows_A/2,0,block_cols_A/2);
    A12=encrypted_A.extract_matrix(0,block_rows_A/2,block_cols_A/2,block_cols_A);
    A21=encrypted_A.extract_matrix(block_rows_A/2,block_rows_A,0,block_cols_A/2);
    A22=encrypted_A.extract_matrix(block_rows_A/2,block_rows_A,block_cols_A/2,block_cols_A);

    B11=encrypted_B.extract_matrix(0,block_rows_B/2,0,block_cols_B/2);
    B12=encrypted_B.extract_matrix(0,block_rows_B/2,block_cols_B/2,block_cols_B);
    B21=encrypted_B.extract_matrix(block_rows_B/2,block_rows_B,0,block_cols_B/2);
    B22=encrypted_B.extract_matrix(block_rows_B/2,block_rows_B,block_cols_B/2,block_cols_B);

    // stop
    if(block_rows_A/2<=1 && block_cols_A/2<=1){
        M1=matrix_block_coprime_mul_imporve_log(s_n,s_m,s_p,ADD(A11,A22,evaluator),ADD(B11,B22,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
        M2=matrix_block_coprime_mul_imporve_log(s_n,s_m,s_p,ADD(A21,A22,evaluator),B11,scale,encoder,evaluator,galois_keys,relin_keys);
        M3=matrix_block_coprime_mul_imporve_log(s_n,s_m,s_p,A11,SUB(B12,B22,evaluator,encoder),scale,encoder,evaluator,galois_keys,relin_keys);
        M4=matrix_block_coprime_mul_imporve_log(s_n,s_m,s_p,A22,SUB(B21,B11,evaluator,encoder),scale,encoder,evaluator,galois_keys,relin_keys);
        M5=matrix_block_coprime_mul_imporve_log(s_n,s_m,s_p,ADD(A11,A12,evaluator),B22,scale,encoder,evaluator,galois_keys,relin_keys);
        M6=matrix_block_coprime_mul_imporve_log(s_n,s_m,s_p,SUB(A21,A11,evaluator,encoder),ADD(B11,B12,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
        M7=matrix_block_coprime_mul_imporve_log(s_n,s_m,s_p,SUB(A12,A22,evaluator,encoder),ADD(B21,B22,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
    }
    else{
        M1=Strassen_block_coprime(s_n,s_m,s_p,ADD(A11,A22,evaluator),ADD(B11,B22,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
        M2=Strassen_block_coprime(s_n,s_m,s_p,ADD(A21,A22,evaluator),B11,scale,encoder,evaluator,galois_keys,relin_keys);
        M3=Strassen_block_coprime(s_n,s_m,s_p,A11,SUB(B12,B22,evaluator,encoder),scale,encoder,evaluator,galois_keys,relin_keys);
        M4=Strassen_block_coprime(s_n,s_m,s_p,A22,SUB(B21,B11,evaluator,encoder),scale,encoder,evaluator,galois_keys,relin_keys);
        M5=Strassen_block_coprime(s_n,s_m,s_p,ADD(A11,A12,evaluator),B22,scale,encoder,evaluator,galois_keys,relin_keys);
        M6=Strassen_block_coprime(s_n,s_m,s_p,SUB(A21,A11,evaluator,encoder),ADD(B11,B12,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
        M7=Strassen_block_coprime(s_n,s_m,s_p,SUB(A12,A22,evaluator,encoder),ADD(B21,B22,evaluator),scale,encoder,evaluator,galois_keys,relin_keys);
    }

    C11=ADD(ADD(M1,M4,evaluator),SUB(M7,M5,evaluator,encoder),evaluator);
    C12=ADD(M3,M5,evaluator);
    C21=ADD(M2,M4,evaluator);
    C22=ADD(SUB(M1,M2,evaluator,encoder),ADD(M3,M6,evaluator),evaluator);

    return  Merage_matrix(C11,C12,C21,C22);
}



// Naive block
matrix_block Naive_block_coprime(size_t s_n,size_t s_m,size_t s_p, matrix_block encrypted_A,matrix_block encrypted_B,double scale, seal::CKKSEncoder& encoder,
                                 seal::Evaluator &evaluator,seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys)
{
    size_t block_rows_A=encrypted_A.get_block_rows();
    size_t block_cols_A=encrypted_A.get_block_cols();
    size_t block_rows_B=encrypted_B.get_block_rows();
    size_t block_cols_B=encrypted_B.get_block_cols();


    if(block_cols_A!=block_rows_B){
        cerr<<"Matrix blocks not equal"<<endl;
    }

    vector<vector<seal::Ciphertext>> cipher_matrix(block_rows_A,vector<seal::Ciphertext>(block_cols_B));
    seal::Ciphertext cipher_tmp;
    for(size_t i=0;i<block_rows_A;i++){
        for(size_t j=0;j<block_cols_B;j++){
            vector<seal::Ciphertext> cipher_vector;
            for(size_t k=0;k<block_cols_A;k++){
                // cout<<"i:"<<i<<"  j:"<<j<<"  k:"<<k<<endl;
                cipher_matrix_coprime matrix_A(encrypted_A.get_cipher_matrix()[i][k],s_n,s_m);
                cipher_matrix_coprime matrix_B(encrypted_B.get_cipher_matrix()[k][j],s_m,s_p);
                cipher_matrix_coprime matrix_C;

                encrypted_coprime_matrix_multiplication(encoder,scale,evaluator,relin_keys,galois_keys,matrix_A,matrix_B,matrix_C);
                cipher_vector.push_back(matrix_C.get_cipher_matrix());
            }
            evaluator.add_many(cipher_vector,cipher_tmp);
            cipher_matrix[i][j]=cipher_tmp;
        }
    }
    return matrix_block(cipher_matrix,encrypted_A.get_rows(),encrypted_B.get_cols(),block_rows_A,block_cols_B);
}


// Naive block imporve log
matrix_block Naive_block_coprime_imporve_log(size_t s_n,size_t s_m,size_t s_p, matrix_block encrypted_A,matrix_block encrypted_B,double scale, seal::CKKSEncoder& encoder,
                                             seal::Evaluator &evaluator,seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys)
{
    size_t block_rows_A=encrypted_A.get_block_rows();
    size_t block_cols_A=encrypted_A.get_block_cols();
    size_t block_rows_B=encrypted_B.get_block_rows();
    size_t block_cols_B=encrypted_B.get_block_cols();


    if(block_cols_A!=block_rows_B){
        cerr<<"Matrix blocks not equal"<<endl;
    }

    vector<vector<seal::Ciphertext>> cipher_matrix(block_rows_A,vector<seal::Ciphertext>(block_cols_B));
    seal::Ciphertext cipher_tmp;
    for(size_t i=0;i<block_rows_A;i++){
        for(size_t j=0;j<block_cols_B;j++){
            vector<seal::Ciphertext> cipher_vector;
            for(size_t k=0;k<block_cols_A;k++){
                // cout<<"i:"<<i<<"  j:"<<j<<"  k:"<<k<<endl;
                cipher_matrix_coprime matrix_A(encrypted_A.get_cipher_matrix()[i][k],s_n,s_m);
                cipher_matrix_coprime matrix_B(encrypted_B.get_cipher_matrix()[k][j],s_m,s_p);
                cipher_matrix_coprime matrix_C;

                encrypted_coprime_matrix_multiplication_imporve_log(encoder,scale,evaluator,relin_keys,galois_keys,matrix_A,matrix_B,matrix_C);
                cipher_vector.push_back(matrix_C.get_cipher_matrix());
            }
            evaluator.add_many(cipher_vector,cipher_tmp);
            cipher_matrix[i][j]=cipher_tmp;
        }
    }
    return matrix_block(cipher_matrix,encrypted_A.get_rows(),encrypted_B.get_cols(),block_rows_A,block_cols_B);
}