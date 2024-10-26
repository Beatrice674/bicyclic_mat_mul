#include <iostream>
#include <seal/seal.h>
#include "matrix.h"
#include "utils.h"

using namespace std;
using namespace seal;

class lu_cipher_B
{
private:
    size_t n;                                       // n: rows
    size_t m;                                       // m: cols
    vector<vector<seal::Ciphertext>> cipher_matrix; // Ciphertext
public:
    lu_cipher_B();
    lu_cipher_B(vector<vector<double>> &data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    size_t get_rows();
    size_t get_cols();
    lu_cipher_B &operator=(const lu_cipher_B &other);
    vector<seal::Ciphertext> get_rows_element_cipher(size_t i);
    vector<seal::Ciphertext> mul_cipher(size_t i,seal::Ciphertext &cipher_tmp,seal::Evaluator &evaluator);
    void dec_lu_matrix(vector<vector<double>> &B, seal::CKKSEncoder &encoder, seal::Decryptor &decryptor);
    void resize(vector<vector<Ciphertext>> &data,size_t n, size_t m);
};

lu_cipher_B::lu_cipher_B()
{
    n = 0;
    m = 0;
}

inline lu_cipher_B::lu_cipher_B(vector<vector<double>> &data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    size_t ciphertext_count = data.size();
    size_t slots_count = encoder.slot_count();

    seal::Ciphertext cipher_tmp;
    seal::Plaintext plain_tmp;

    for (size_t i = 0; i < ciphertext_count; i++)
    {
        vector<seal::Ciphertext> split_vector;
        for (size_t j = 0; j < data[i].size(); j += slots_count)
        {
            vector<double> sub_vector(data[i].begin() + j, data[i].begin() + min(j + slots_count, data[i].size()));
            encoder.encode(sub_vector, scale, plain_tmp);
            encryptor.encrypt(plain_tmp, cipher_tmp);
            split_vector.push_back(cipher_tmp);
        }
        cipher_matrix.push_back(split_vector);
    }
    this->n = ciphertext_count;
    this->m = data[0].size(); 
}

inline size_t lu_cipher_B::get_rows()
{
    return n;
}

inline size_t lu_cipher_B::get_cols()
{
    return m;
}

inline lu_cipher_B &lu_cipher_B::operator=(const lu_cipher_B &other)
{
    if (this == &other)
    {
        return *this; // Avoid self-assignment
    }

    // Copy data members from the other object
    n = other.n;
    m = other.m;
    cipher_matrix = other.cipher_matrix;

    // Return the reference to the updated object
    return *this;
}

inline vector<seal::Ciphertext> lu_cipher_B::get_rows_element_cipher(size_t i)
{
    return cipher_matrix[i];
}

inline vector<seal::Ciphertext> lu_cipher_B::mul_cipher(size_t index, seal::Ciphertext &cipher, seal::Evaluator &evaluator)
{
    vector<seal::Ciphertext> destination=cipher_matrix[index];
    for(size_t i=0;i<destination.size();i++){
        evaluator.mod_switch_to_inplace(destination[i],cipher.parms_id());
    }
    // cout<<index<<endl;
    for(size_t i=0;i<destination.size();i++){
        evaluator.multiply_inplace(destination[i],cipher);
    }
    return destination;
}

inline void lu_cipher_B::dec_lu_matrix(vector<vector<double>> &B, seal::CKKSEncoder &encoder, seal::Decryptor &decryptor)
{
    
    Plaintext plain_tmp;
    size_t length=encoder.slot_count();
    size_t cols=(m+length-1)/length;
    for(size_t i=0;i<n;i++){
        vector<double> rows_element,rows_tmp;
        for(size_t j=0;j<cols;j++){
            // cout<<i<<"  "<<j<<m-j<<endl;
            decryptor.decrypt(cipher_matrix[i][j],plain_tmp);
            encoder.decode(plain_tmp,rows_tmp);
            rows_element.insert(rows_element.end(),rows_tmp.begin(),rows_tmp.begin()+min(length,m-j*length));
        }
        B.push_back(rows_element);
    }
}

inline void lu_cipher_B::resize(vector<vector<Ciphertext>> &data,size_t n,size_t m)
{
    this->n=n;
    this->m=m;
    cipher_matrix=data;
}

class lu_cipher_A
{
private:
    size_t n;
    size_t m;
    size_t batch_size;
    vector<seal::Ciphertext> cipher_matrix;

public:
    lu_cipher_A();
    lu_cipher_A(vector<vector<double>> &data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    lu_cipher_A(vector<vector<double>> &data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor, size_t batch_size);
    vector<seal::Ciphertext> replicate(size_t i, double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator,
                                       seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys);
    seal::Plaintext gen_plain_vec(size_t start_index, double scale, seal::CKKSEncoder &encoder);
    seal::Plaintext gen_replicate_plain(size_t replicate_index, double scale, seal::CKKSEncoder &encoder);
    size_t get_rows();
    size_t get_cols();
};

lu_cipher_A::lu_cipher_A()
{
    n = 0;
    m = 0;
    batch_size=0;
}

inline lu_cipher_A::lu_cipher_A(vector<vector<double>> &data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    n = data.size();
    m = data[0].size();
    this->batch_size=8;

    seal::Ciphertext cipher_tmp;
    seal::Plaintext plain_tmp;

    for (size_t i = 0; i < n; i++)
    {
        encoder.encode(data[i], scale, plain_tmp);
        encryptor.encrypt(plain_tmp, cipher_tmp);
        cipher_matrix.push_back(cipher_tmp);
    }
}

inline lu_cipher_A::lu_cipher_A(vector<vector<double>> &data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor, size_t batch_size)
{
    if(batch_size>0&&((batch_size&(batch_size-1))!=0)){
        cerr<<" Batch size not power of 2."<<endl;
    }

    n = data.size();
    m = data[0].size();
    this->batch_size=batch_size;

    seal::Ciphertext cipher_tmp;
    seal::Plaintext plain_tmp;

    for (size_t i = 0; i < n; i++)
    {
        encoder.encode(data[i], scale, plain_tmp);
        encryptor.encrypt(plain_tmp, cipher_tmp);
        cipher_matrix.push_back(cipher_tmp);
    }
}

inline vector<seal::Ciphertext> lu_cipher_A::replicate(size_t i, double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys)
{
    Ciphertext rotate_cipher=cipher_matrix[i];
    vector<seal::Ciphertext> split_rotate_cipher;
    size_t l=encoder.slot_count();
    seal::Plaintext plain_tmp;
    seal::Ciphertext cipher_tmp;

    vector<seal::Ciphertext> destination;
    
    // Calculate the number of batches
    size_t num_batches = (m + batch_size - 1) / batch_size;
    if(m>batch_size)
    {
        // Split the rotate_cipher into batches
        for (size_t j = 0; j < num_batches; j++)
        {
            plain_tmp=gen_plain_vec(j*batch_size,scale,encoder);
            evaluator.mod_switch_to_inplace(plain_tmp,rotate_cipher.parms_id());
            evaluator.multiply_plain(rotate_cipher,plain_tmp,cipher_tmp);
            evaluator.relinearize_inplace(cipher_tmp,relin_keys);
            evaluator.rescale_to_next_inplace(cipher_tmp);
            split_rotate_cipher.push_back(cipher_tmp);
        }
    }
    else
    {
        split_rotate_cipher.push_back(rotate_cipher);
    }

    // batch rotate
    for(size_t i=0;i<num_batches;i++){
        seal::Ciphertext out_rotate_cipher=split_rotate_cipher[i];
        for(size_t j=batch_size;j<l;j=j*2){
            evaluator.rotate_vector(out_rotate_cipher,j,galois_keys,cipher_tmp);
            evaluator.add_inplace(out_rotate_cipher,cipher_tmp);
        }

        for(size_t j=0;j<min(batch_size,m-i*batch_size);j++){
            seal::Ciphertext inline_rotate_cipher;
            plain_tmp=gen_replicate_plain(j,scale,encoder);
            evaluator.mod_switch_to_inplace(plain_tmp,out_rotate_cipher.parms_id());

            evaluator.multiply_plain(out_rotate_cipher,plain_tmp,inline_rotate_cipher);
            evaluator.relinearize_inplace(inline_rotate_cipher,relin_keys);
            evaluator.rescale_to_next_inplace(inline_rotate_cipher);
            
            for(size_t k=1;k<batch_size;k=k*2){
                evaluator.rotate_vector(inline_rotate_cipher,k,galois_keys,cipher_tmp);
                evaluator.add_inplace(inline_rotate_cipher,cipher_tmp);
            }
            destination.push_back(inline_rotate_cipher);
        }
    }

    return destination;
}

inline seal::Plaintext lu_cipher_A::gen_plain_vec(size_t start_index, double scale, seal::CKKSEncoder &encoder)
{
    size_t length=encoder.slot_count();
    vector<double> vec_tmp(length,0);
    for(size_t i=0;i<batch_size;i++){
        vec_tmp[(i+start_index)%length]=1;
    }
    // for(size_t i=0;i<length;i++){
    //     cout<<vec_tmp[i]<<",";
    // }
    // cout<<endl;
    seal::Plaintext destination;
    encoder.encode(vec_tmp,scale,destination);
    return destination;
}

inline seal::Plaintext lu_cipher_A::gen_replicate_plain(size_t replicate_index, double scale, seal::CKKSEncoder &encoder)
{
    size_t length=encoder.slot_count();
    vector<double> vec_tmp(length,0);
    for(size_t i=replicate_index;i<length;i+=batch_size){
        vec_tmp[i]=1;
    }
    seal::Plaintext destination;
    encoder.encode(vec_tmp,scale,destination);
    return destination;
}

inline size_t lu_cipher_A::get_rows()
{
    return n;
}

inline size_t lu_cipher_A::get_cols()
{
    return m;
}

void encrypted_lu_matrix_multiplication(seal::CKKSEncoder &encoder, // seal::Decryptor &decryptor,
                                          double scale, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys,
                                          seal::GaloisKeys &galois_keys, lu_cipher_A &encrypted_A,
                                          lu_cipher_B &encrypted_B, lu_cipher_B &destination)
{
    if (encrypted_A.get_cols() != encrypted_B.get_rows())
    {
        cerr << "Matrix dimensions are not equal to multiply.";
    }
    size_t n, p;
    n = encrypted_A.get_rows();
    p = encrypted_B.get_cols();
    vector<vector<seal::Ciphertext>> encrypted_C;

    for (std::size_t i = 0; i < n; i++)
    {
        vector<seal::Ciphertext> repliacate_cipher=encrypted_A.replicate(i,scale,encoder,evaluator,galois_keys,relin_keys);
        vector<seal::Ciphertext> cipher_vector,rows_element_tmp;
        // cout<<i<<endl;
        for(size_t j=0;j<repliacate_cipher.size();j++){
            rows_element_tmp = encrypted_B.mul_cipher(j,repliacate_cipher[j],evaluator);
            if(j==0){
                cipher_vector=rows_element_tmp;
            }
            else{
                for(size_t k=0;k<cipher_vector.size();k++){
                    evaluator.add_inplace(cipher_vector[k],rows_element_tmp[k]);
                }
            }
        }

        //relinearize and rescale
        for(size_t k=0;k<cipher_vector.size();k++){
            evaluator.relinearize_inplace(cipher_vector[k],relin_keys);
            evaluator.rescale_to_next_inplace(cipher_vector[k]);
        }
        encrypted_C.push_back(cipher_vector);
    }
    destination.resize(encrypted_C,n,p);
    // destination=encrypted_C;
}
