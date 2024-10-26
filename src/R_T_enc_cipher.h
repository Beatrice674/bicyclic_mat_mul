#include <iostream>
#include <seal/seal.h>
#include "matrix.h"
#include "utils.h"

using namespace std;
using namespace seal;

class cipher_matrix_r_t{
private:
    size_t n;//rows
    size_t m;//cols
    size_t d;
    seal::Ciphertext cipher_matrix;
public:
    cipher_matrix_r_t();
    cipher_matrix_r_t(seal::Ciphertext &cipher,size_t n,size_t m,size_t d);
    cipher_matrix_r_t(Matrix<double> data,double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    void enc_matrix_cipher(vector<double> &data,double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    void dec_matrix_cipher(Matrix<double> &destination,seal::CKKSEncoder &encoder,seal::Decryptor &decryptor);
    virtual ~cipher_matrix_r_t();

    // change to matrix A
    void change_cipher_to_matrix_A(double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator,
                                   seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys);
    // change to matrix B
    void change_cipher_to_matrix_B(double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator,
                                   seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys);
    // get i-th rows by ciphertext 
    void get_cipher_row(size_t index,seal::Ciphertext& destination, double scale, seal::CKKSEncoder& encoder,
                        seal::Evaluator &evaluator,seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys);
    // get i-th cols by ciphertext 
    void get_cipher_col(size_t index,seal::Ciphertext& destination, double scale, seal::CKKSEncoder& encoder,
                        seal::Evaluator &evaluator,seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys);
    size_t get_rows();
    size_t get_cols();
    size_t get_matrix_mul_size();
    Ciphertext get_matrix_cipher();                        
};

cipher_matrix_r_t::cipher_matrix_r_t()
{
    ;
}

inline cipher_matrix_r_t::cipher_matrix_r_t(seal::Ciphertext &cipher, size_t n, size_t m, size_t d)
{
    cipher_matrix=cipher;
    this->n=n;
    this->m=m;
    this->d=d;
}

inline cipher_matrix_r_t::cipher_matrix_r_t(Matrix<double> data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    this->n=data.get_rows();
    this->m=data.get_cols();
    d=round(pow(encoder.slot_count(),1.0/3.0));
    data.resize(d,d);
    vector<double> flatten_data=data.flatten_matrix_to_rows_vector();
    enc_matrix_cipher(flatten_data,scale,encoder,encryptor);
}

inline void cipher_matrix_r_t::enc_matrix_cipher(vector<double> &data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    if(data.size()>d*d){
        cerr << "  !!! the number of slot is not enough for the Matrix" << endl;
    }

    seal::Plaintext plain_tmp;
    encoder.encode(data,scale,plain_tmp);
    encryptor.encrypt(plain_tmp,cipher_matrix);
}

inline void cipher_matrix_r_t::dec_matrix_cipher(Matrix<double> &destination, seal::CKKSEncoder &encoder, seal::Decryptor &decryptor)
{
    //decrypte ciphertext
    Plaintext plain_tmp;
    vector<double> vec_tmp;
    decryptor.decrypt(cipher_matrix,plain_tmp);
    encoder.decode(plain_tmp,vec_tmp);

    destination.resize(d,d);
    for(size_t i=0;i<d;i++){
        for(size_t j=0;j<d;j++){
            destination.set(i,j,vec_tmp[i*d+j]);
        }
    }

    destination.resize(n,m);
}

inline cipher_matrix_r_t::~cipher_matrix_r_t()
{
    ;
}

inline void cipher_matrix_r_t::change_cipher_to_matrix_A(double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys)
{
    seal::Ciphertext cipher_tmp,cipher_a;
    vector<seal::Ciphertext> cipher_vector;
    // cout<<"strat change matrix a"<<endl;
    // generate matrix a 
    size_t rotate_number=0;
    for(size_t i=0;i<d;i++){
        get_cipher_row(i,cipher_tmp,scale,encoder,evaluator,galois_keys,relin_keys);
        cipher_vector.push_back(cipher_tmp);
        rotate_number++;
    }
    evaluator.add_many(cipher_vector,cipher_a);
    evaluator.relinearize_inplace(cipher_a,relin_keys);
    evaluator.rescale_to_next_inplace(cipher_a);

    // cout<<"rotate cipher a"<<endl;
    // rotate and add
    for(size_t i=1;i<d;i=i*2){
        evaluator.rotate_vector(cipher_a,-i,galois_keys,cipher_tmp);
        evaluator.add_inplace(cipher_a,cipher_tmp);
        rotate_number++;
    }

    cipher_matrix=cipher_a;
    // cout<<rotate_number<<endl;
}

inline void cipher_matrix_r_t::change_cipher_to_matrix_B(double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys)
{
    seal::Ciphertext cipher_tmp,cipher_b;
    vector<seal::Ciphertext> cipher_vector;
    // generate matrix a 
    for(size_t i=0;i<d;i++){
        get_cipher_col(i,cipher_tmp,scale,encoder,evaluator,galois_keys,relin_keys);
        cipher_vector.push_back(cipher_tmp);
    }
    evaluator.add_many(cipher_vector,cipher_b);
    evaluator.relinearize_inplace(cipher_b,relin_keys);
    evaluator.rescale_to_next_inplace(cipher_b);

    // rotate and add
    for(size_t i=d;i<d*d;i=i*2){
        evaluator.rotate_vector(cipher_b,-i,galois_keys,cipher_tmp);
        evaluator.add_inplace(cipher_b,cipher_tmp);
    }

    cipher_matrix=cipher_b;
}

inline void cipher_matrix_r_t::get_cipher_row(size_t index, seal::Ciphertext &destination, double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys)
{
    vector<double> vec_tmp(d*d,0);
    for(size_t i=0;i<d;i++){
        vec_tmp[i*d+index]=1;
    }

    seal::Plaintext plain_tmp;
    encoder.encode(vec_tmp,scale,plain_tmp);
    evaluator.mod_switch_to_inplace(plain_tmp,cipher_matrix.parms_id());

    evaluator.multiply_plain(cipher_matrix,plain_tmp,destination);
    evaluator.rotate_vector_inplace(destination,-(d*d*index-index),galois_keys);
}

inline void cipher_matrix_r_t::get_cipher_col(size_t index, seal::Ciphertext &destination, double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys)
{
    vector<double> vec_tmp(d*d,0);
    for(size_t i=0;i<d;i++){
        vec_tmp[i+index*d]=1;
    }

    seal::Plaintext plain_tmp;
    encoder.encode(vec_tmp,scale,plain_tmp);
    evaluator.mod_switch_to_inplace(plain_tmp,cipher_matrix.parms_id());

    evaluator.multiply_plain(cipher_matrix,plain_tmp,destination);
    evaluator.rotate_vector_inplace(destination,-(d*d*index-index*d),galois_keys);
}

inline size_t cipher_matrix_r_t::get_rows()
{
    return n;
}

inline size_t cipher_matrix_r_t::get_cols()
{
    return m;
}

inline size_t cipher_matrix_r_t::get_matrix_mul_size()
{
    return d;
}

inline Ciphertext cipher_matrix_r_t::get_matrix_cipher()
{
    return cipher_matrix;
}


void encrypted_r_t_matrix_multiplication(seal::CKKSEncoder &encoder, // seal::Decryptor &decryptor,
                                          double scale, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys,
                                          seal::GaloisKeys &galois_keys, cipher_matrix_r_t &encrypted_A,
                                          cipher_matrix_r_t &encrypted_B, cipher_matrix_r_t &destination)
{
    if(encrypted_A.get_cols()!=encrypted_B.get_rows()){
        cerr << "Matrix dimensions are not equal to multiply.";
    }
    size_t n=encrypted_A.get_rows();
    size_t p=encrypted_B.get_cols();
    size_t d=encrypted_A.get_matrix_mul_size();

    seal::Ciphertext cipher_a,cipher_b,cipher_c,cipher_tmp;
    cipher_a=encrypted_A.get_matrix_cipher();
    cipher_b=encrypted_B.get_matrix_cipher();

    evaluator.multiply(cipher_a,cipher_b,cipher_c);
    evaluator.relinearize_inplace(cipher_c,relin_keys);
    evaluator.rescale_to_next_inplace(cipher_c);

    for(size_t i=d*d;i<d*d*d;i=i*2){
        evaluator.rotate_vector(cipher_c,i,galois_keys,cipher_tmp);
        evaluator.add_inplace(cipher_c,cipher_tmp);
    }

    cipher_matrix_r_t ct_C(cipher_c,n,p,d);
    destination=ct_C;
}