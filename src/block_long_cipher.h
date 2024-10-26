#include <iostream>
#include <seal/seal.h>
#include "matrix.h"
#include "utils.h"
#include "long_enc_cipher.h"

using namespace std;
using namespace seal;


class long_matrix_block{
private:
    size_t n;//rows
    size_t m;//cols
    size_t b_rows;//block rows
    size_t b_cols;//block cols
    vector<vector<long_cipher>> block_matrix_cipher;//block matrix
public:
    long_matrix_block();
    long_matrix_block(vector<vector<long_cipher>> data,size_t n,size_t m,size_t b_rows,size_t b_cols);
    long_matrix_block(Matrix<double> data,size_t s_n,size_t s_m,double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    long_matrix_block extract_matrix(size_t start_row,size_t end_row,size_t start_col,size_t end_col);
    void dec_block_long_cipher(Matrix<double> &destination,size_t s_n,size_t s_m,seal::CKKSEncoder &encoder,seal::Decryptor &decryptor);
    ~long_matrix_block();
    size_t get_rows();
    size_t get_cols();
    size_t get_block_rows();
    size_t get_block_cols();
    vector<vector<long_cipher>> get_cipher_matrix();
    vector<vector<seal::Ciphertext>> get_sub_cipher_matrix(seal::Evaluator &evaluator,seal::CKKSEncoder &encoder);
};

long_matrix_block::long_matrix_block()
{
    ;
}

long_matrix_block::long_matrix_block(vector<vector<long_cipher>> data, size_t n, size_t m, size_t b_rows, size_t b_cols)
{
    block_matrix_cipher=data;
    this->n=n;
    this->m=m;
    this->b_rows=b_rows;
    this->b_cols=b_cols;
}

inline long_matrix_block::long_matrix_block(Matrix<double> data, size_t s_n, size_t s_m, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
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
            
            long_cipher tmp(s_n,s_m,da,scale,encoder,encryptor);
            block_matrix_cipher[i].push_back(tmp);
        }
    }
}

inline void long_matrix_block::dec_block_long_cipher(Matrix<double> &destination, size_t s_n, size_t s_m, seal::CKKSEncoder &encoder, seal::Decryptor &decryptor)
{
    destination.resize(b_rows*s_n,b_cols*s_m);
    Matrix<double> tmp_matrix;
    for(size_t i=0;i<b_rows;i++){
        for(size_t j=0;j<b_cols;j++){
            vector<double> dec_vector;
            block_matrix_cipher[i][j].dec_diag_vector(dec_vector,encoder,decryptor);
            
            tmp_matrix.clear();
            tmp_matrix.resize(s_n,s_m);
            
            tmp_matrix=diagonal_decoding(dec_vector,s_n,s_m);
            // tmp_matrix.print(6,6);
            destination.Merage_matrix(i*s_n,(i+1)*s_n,j*s_m,(j+1)*s_m,tmp_matrix);
        }
    }
}

long_matrix_block::~long_matrix_block()
{
    ;
}

inline size_t long_matrix_block::get_rows()
{
    return n;
}

inline size_t long_matrix_block::get_cols()
{
    return m;
}

inline size_t long_matrix_block::get_block_rows()
{
    return b_rows;
}

inline size_t long_matrix_block::get_block_cols()
{
    return b_cols;
}

inline vector<vector<long_cipher>> long_matrix_block::get_cipher_matrix()
{
    return block_matrix_cipher;
}

// Naive block
long_matrix_block Naive_long_block_coprime(size_t s_n,size_t s_m,size_t s_p, long_matrix_block encrypted_A,long_matrix_block encrypted_B,double scale, seal::CKKSEncoder& encoder,
                                           seal::Evaluator &evaluator,seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys,seal::Encryptor &encryptor)
{
    size_t block_rows_A=encrypted_A.get_block_rows();
    size_t block_cols_A=encrypted_A.get_block_cols();
    size_t block_rows_B=encrypted_B.get_block_rows();
    size_t block_cols_B=encrypted_B.get_block_cols();


    if(block_cols_A!=block_rows_B){
        cerr<<"Matrix blocks not equal"<<endl;
    }

    vector<vector<long_cipher>> cipher_matrix(block_rows_A,vector<long_cipher>(block_cols_B));
    for(size_t i=0;i<block_rows_A;i++){
        for(size_t j=0;j<block_cols_B;j++){
            long_cipher long_tmp(s_n,s_p,scale,encoder,encryptor);
            for(size_t k=0;k<block_cols_A;k++){
                // cout<<"i:"<<i<<"  j:"<<j<<"  k:"<<k<<endl;
                long_cipher matrix_A=encrypted_A.get_cipher_matrix()[i][k];
                long_cipher matrix_B=encrypted_B.get_cipher_matrix()[k][j];
                long_cipher matrix_C(s_n,s_p,scale,encoder,encryptor);

                encrypted_long_matrix_multiplication(encoder, scale, evaluator, relin_keys, galois_keys, matrix_A, matrix_B, matrix_C);
                vector<seal::Ciphertext> cipher_tmp=matrix_C.get_cipher();
                long_tmp.LongAdd(cipher_tmp,evaluator);
            }
            cipher_matrix[i][j]=long_tmp;
        }
    }
    return long_matrix_block(cipher_matrix,encrypted_A.get_rows(),encrypted_B.get_cols(),block_rows_A,block_cols_B);
}