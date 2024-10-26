#include <iostream>
#include <seal/seal.h>
#include "matrix.h"
#include "utils.h"

using namespace std;
using namespace seal;

class cipher_matrix_jiang{
private:
    size_t n;//rows
    size_t m;//cols
    size_t d;
    seal::Ciphertext cipher_matrix;
public:
    cipher_matrix_jiang();
    cipher_matrix_jiang(seal::Ciphertext &cipher,size_t n,size_t m,size_t d);
    cipher_matrix_jiang(Matrix<double> data,double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    void enc_matrix_cipher(vector<double> &data,double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    void dec_matrix_cipher(Matrix<double> &destination,seal::CKKSEncoder &encoder,seal::Decryptor &decryptor);
    virtual ~cipher_matrix_jiang();

    // change to matrix A
    void change_cipher_to_matrix_A(double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator,
                                   seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys);
    // change to matrix B
    void change_cipher_to_matrix_B(double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator,
                                   seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys);
    //rotate matrix A
    void rotate_ctA(int i ,seal::Ciphertext& destination,seal::Evaluator &evaluator,
                    seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys, seal::CKKSEncoder& encoder);
    //rotate matrix B
    void rotate_ctB(int i ,seal::Ciphertext& destination,seal::Evaluator &evaluator,
                    seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys, seal::CKKSEncoder& encoder);
    size_t get_rows();
    size_t get_cols();
    size_t get_matrix_mul_size();   
    seal::Ciphertext get_cipher_matrix();                      
};

inline cipher_matrix_jiang::cipher_matrix_jiang()
{
    n=0;
    m=0;
}

inline cipher_matrix_jiang::cipher_matrix_jiang(seal::Ciphertext &cipher,size_t n,size_t m,size_t d)
{
    cipher_matrix=cipher;
    this->n=n;
    this->m=m;
    this->d=d;
}

inline cipher_matrix_jiang::cipher_matrix_jiang(Matrix<double> data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    this->n=data.get_rows();
    this->m=data.get_cols();
    d=sqrt(encoder.slot_count());
    data.resize(d,d);
    vector<double> flatten_data=data.flatten_matrix_to_rows_vector();
    enc_matrix_cipher(flatten_data,scale,encoder,encryptor);
}

inline void cipher_matrix_jiang::enc_matrix_cipher(vector<double> &data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    if(data.size()>encoder.slot_count()){
        cerr << "  !!! the number of slot is not enough for the Matrix" << endl;
    }

    seal::Plaintext plain_tmp;
    encoder.encode(data,scale,plain_tmp);
    encryptor.encrypt(plain_tmp,cipher_matrix);
}

inline void cipher_matrix_jiang::dec_matrix_cipher(Matrix<double> &destination, seal::CKKSEncoder &encoder, seal::Decryptor &decryptor)
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

inline cipher_matrix_jiang::~cipher_matrix_jiang()
{
    ;
}

inline void cipher_matrix_jiang::change_cipher_to_matrix_A(double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys)
{
    // cout<<"start change cipher to matrix A"<<endl;
    Matrix<double> u_sigma;
    u_sigma.generate_u_sigma(d, d);
    Ciphertext cipher_tmp;//save intermediate variables.
    vector<Ciphertext> cipher_rotate,cipher_result;//save rotate result vector
    Plaintext plain_tmp;
    vector<double> vec_tmp;
    Ciphertext destination;


    //start change matrixA
    // auto start = std::chrono::high_resolution_clock::now();
    //step1. determining the dimensions
    int rotate_size = 2 * d - 1;
    int rotate_outside_size = ceil(sqrt(rotate_size));
    int rotate_inside_size = ceil(double(rotate_size) / double(rotate_outside_size));
    // step2. pre-rotate ciphertext
    for (int i = 0; i < rotate_inside_size; i++) {
        evaluator.rotate_vector(cipher_matrix, -d+d*d+i+1, galois_keys, cipher_tmp);
        cipher_rotate.push_back(cipher_tmp);
    }

    //step3. starting step function 
    for (int i = 0; i < rotate_outside_size; i++) {
        vector<Ciphertext> cipher_vector;
        for (int j = 0; j < rotate_inside_size; j++) {
            vec_tmp = u_sigma.diag_vector(i * rotate_inside_size + j-d+1,encoder.slot_count());
            if (std::all_of(vec_tmp.begin(), vec_tmp.end(), [](double num) { return num == 0; }) == 1) {
                continue;
            }
            else {
                std::rotate(vec_tmp.rbegin(), vec_tmp.rbegin() + i * rotate_inside_size, vec_tmp.rend());
                encoder.encode(vec_tmp, cipher_matrix.scale(), plain_tmp);
                evaluator.mod_switch_to_inplace(plain_tmp, cipher_matrix.parms_id());
                evaluator.multiply_plain(cipher_rotate[j], plain_tmp, cipher_tmp);
                cipher_vector.push_back(cipher_tmp);
            }
        }
        // if(cipher_vector.size()==0){
        //     continue;
        // }
        evaluator.add_many(cipher_vector, cipher_tmp);
        evaluator.relinearize_inplace(cipher_tmp, relin_keys);
        evaluator.rescale_to_next_inplace(cipher_tmp);
        evaluator.rotate_vector_inplace(cipher_tmp, i * rotate_inside_size, galois_keys);
        cipher_result.push_back(cipher_tmp);
    }
    evaluator.add_many(cipher_result, destination);
    cipher_matrix=destination;
}

inline void cipher_matrix_jiang::change_cipher_to_matrix_B(double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys)
{
    // cout<<"start change cipher to matrix B"<<endl;
    // generate u_tau matrix
    Matrix<double> u_tau;
    u_tau.generate_u_tau(d, d);
    Ciphertext cipher_tmp;//save intermediate variables
    vector<Ciphertext> cipher_rotate, cipher_result;//save rotate result vector
    Plaintext plain_tmp;
    vector<double> vec_tmp;
    Ciphertext destination;

    //start change matrixB
    //step1. determining the dimensions
    int rotate_size = d;
    int rotate_outside_size = ceil(sqrt(rotate_size));
    int rotate_inside_size = ceil(double(rotate_size) / double(rotate_outside_size));
    // step2. pre-rotate ciphertext 
    for (int i = 0; i < rotate_inside_size; i++) {
        evaluator.rotate_vector(cipher_matrix, i*d, galois_keys, cipher_tmp);
        cipher_rotate.push_back(cipher_tmp);
    }
    //step3. starting step function 
    for (int i = 0; i < rotate_outside_size; i++) {
        // cout<<i<<endl;
        vector<Ciphertext> cipher_vector;
        for (int j = 0; j < rotate_inside_size; j++) {
            if(size_t(i*rotate_inside_size+j)>=d){
                continue;
            }
            vec_tmp = u_tau.diag_vector((i*rotate_inside_size+j)*d,encoder.slot_count());
            // for(size_t k=0;k<512;k++){
            //     if(vec_tmp[k]==1){
            //         cout<<k<<"  ";
            //     }
            // }
            // cout<<endl;
            if (std::all_of(vec_tmp.begin(), vec_tmp.end(), [](double num) { return num == 0; }) == 1) {
                continue;
            }
            else {
                std::rotate(vec_tmp.rbegin(), vec_tmp.rbegin() + i * rotate_inside_size*d, vec_tmp.rend());
                encoder.encode(vec_tmp, cipher_matrix.scale(), plain_tmp);
                evaluator.mod_switch_to_inplace(plain_tmp, cipher_matrix.parms_id());
                evaluator.multiply_plain(cipher_rotate[j], plain_tmp, cipher_tmp);
                cipher_vector.push_back(cipher_tmp);
            }
        }
        // if(cipher_vector.size()==0){
        //     continue;
        // }
        evaluator.add_many(cipher_vector, cipher_tmp);
        evaluator.relinearize_inplace(cipher_tmp, relin_keys);
        evaluator.rescale_to_next_inplace(cipher_tmp);
        evaluator.rotate_vector_inplace(cipher_tmp, i * rotate_inside_size*d, galois_keys);
        cipher_result.push_back(cipher_tmp);
    }
    evaluator.add_many(cipher_result, destination);
    cipher_matrix=destination;
}

inline void cipher_matrix_jiang::rotate_ctA(int i, seal::Ciphertext &destination, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys, seal::CKKSEncoder &encoder)
{
    // step1.1.1 generate vector
    if (i != 0) {
        vector<double> rotate_right(d * d), rotate_left(d * d);
        seal::Plaintext right_plain, left_plain;
        seal::Ciphertext right_cipher, left_cipher, result_cipher;
        for (size_t j = 0; j < d; j++) {
            size_t k = 0;
            while (k < d) {
                if (k >= d-i) {
                    rotate_right[j * d + k] = 1;
                }
                else {
                    rotate_left[j * d + k] = 1;
                }
                k++;
            }
        }

        //step1.1.2 encode vector
        encoder.encode(rotate_right, cipher_matrix.scale(), right_plain);
        evaluator.mod_switch_to_inplace(right_plain, cipher_matrix.parms_id());
        encoder.encode(rotate_left, cipher_matrix.scale(), left_plain);
        evaluator.mod_switch_to_inplace(left_plain, cipher_matrix.parms_id());


        //step1.1.3 rotate and multiply
        evaluator.rotate_vector(cipher_matrix, d * d - d + i, galois_keys, right_cipher);
        evaluator.rotate_vector(cipher_matrix, i, galois_keys, left_cipher);
        evaluator.multiply_plain_inplace(right_cipher, right_plain);
        evaluator.multiply_plain_inplace(left_cipher, left_plain);
        evaluator.add(right_cipher, left_cipher, result_cipher);
        evaluator.relinearize_inplace(result_cipher, relin_keys);
        evaluator.rescale_to_next_inplace(result_cipher);
        destination = result_cipher;
    }
    else {
        vector<double> rotate_all(d * d, 1);
        Plaintext plain_tmp;
        Ciphertext result_cipher;
        encoder.encode(rotate_all, cipher_matrix.scale(), plain_tmp);
        evaluator.mod_switch_to_inplace(plain_tmp, cipher_matrix.parms_id());
        evaluator.multiply_plain(cipher_matrix, plain_tmp, result_cipher);
        evaluator.relinearize_inplace(result_cipher, relin_keys);
        evaluator.rescale_to_next_inplace(result_cipher);
        destination = result_cipher;
    }
}

inline void cipher_matrix_jiang::rotate_ctB(int i, seal::Ciphertext &destination, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys, seal::CKKSEncoder &encoder)
{
    evaluator.rotate_vector(cipher_matrix,i,galois_keys,destination);
}

inline size_t cipher_matrix_jiang::get_rows()
{
    return n;
}

inline size_t cipher_matrix_jiang::get_cols()
{
    return m;
}

inline size_t cipher_matrix_jiang::get_matrix_mul_size()
{
    return d;
}

inline seal::Ciphertext cipher_matrix_jiang::get_cipher_matrix()
{
    return cipher_matrix;
}

void encrypted_jiang_matrix_multiplication(seal::CKKSEncoder &encoder, // seal::Decryptor &decryptor,
                                          double scale, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys,
                                          seal::GaloisKeys &galois_keys, cipher_matrix_jiang &encrypted_A,
                                          cipher_matrix_jiang &encrypted_B, cipher_matrix_jiang &destination)
{
    if(encrypted_A.get_cols()!=encrypted_B.get_rows()){
        cerr << "Matrix dimensions are not equal to multiply.";
    }
    size_t n=encrypted_A.get_rows();
    size_t p=encrypted_B.get_cols();
    size_t d=encrypted_A.get_matrix_mul_size();

    vector<Ciphertext>  cipher_vector;
    Ciphertext cipher_tmpA,cipher_tmpB,cipher_tmp,encrypted_C;//save intermediate variables
    Plaintext plain_tmp;
    vector<double> vec_tmp;

    for (size_t i = 0; i < d; i++) {
        //step1.1 move ct.A(0)-->ct.A(i)
        encrypted_A.rotate_ctA(i, cipher_tmpA, evaluator, galois_keys, relin_keys, encoder);
        encrypted_B.rotate_ctB(d*i, cipher_tmpB, evaluator, galois_keys, relin_keys, encoder);
        evaluator.mod_switch_to_inplace(cipher_tmpB, cipher_tmpA.parms_id());
        evaluator.multiply(cipher_tmpA, cipher_tmpB, cipher_tmp);
        cipher_vector.push_back(cipher_tmp);
    }
    evaluator.add_many(cipher_vector, encrypted_C);
    evaluator.relinearize_inplace(encrypted_C, relin_keys);
    evaluator.rescale_to_next_inplace(encrypted_C);

    cipher_matrix_jiang ct_C(encrypted_C,n,p,d);
    destination=ct_C;
}