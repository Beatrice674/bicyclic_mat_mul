#include <iostream>
#include <seal/seal.h>
#include "matrix.h"
#include "utils.h"

using namespace std;
using namespace seal;

class cipher_matrix_coprime{
private:
    size_t n;//rows
    size_t m;//cols
    seal::Ciphertext cipher_matrix;
public:
    cipher_matrix_coprime();
    cipher_matrix_coprime(seal::Ciphertext &cipher,size_t n,size_t m);
    cipher_matrix_coprime(Matrix<double> data,size_t repeate,double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    void enc_matrix_cipher(vector<double> data,size_t repeate,double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    void dec_matrix_cipher(Matrix<double> &destination,seal::CKKSEncoder &encoder,seal::Decryptor &decryptor);
    virtual ~cipher_matrix_coprime();

    //rotate matrix
    void rotate_ct(int i ,seal::Ciphertext& destination,seal::Evaluator &evaluator,
                   seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys, seal::CKKSEncoder& encoder);
    size_t get_rows();
    size_t get_cols(); 
    seal::Ciphertext get_cipher_matrix();                      
};

inline cipher_matrix_coprime::cipher_matrix_coprime()
{
    n=0;
    m=0;
}

inline cipher_matrix_coprime::cipher_matrix_coprime(seal::Ciphertext &cipher,size_t n,size_t m)
{
    cipher_matrix=cipher;
    this->n=n;
    this->m=m;
}

inline cipher_matrix_coprime::cipher_matrix_coprime(Matrix<double> data, size_t repeate, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    this->n=data.get_rows();
    this->m=data.get_cols();

    vector<double> vector_data;
    vector_data = data.diagonal_encoding();
    enc_matrix_cipher(vector_data,repeate,scale,encoder,encryptor);
}

inline void cipher_matrix_coprime::enc_matrix_cipher(vector<double> data,size_t repeate, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    if(data.size()>encoder.slot_count()/2){
        cerr << "  !!! the number of slot is not enough for the Matrix" << endl;
    }

    vector<double> dA;
    dA.resize(encoder.slot_count());

    for (size_t i = 0; i < n * m; i++)
    {
        for (size_t j = 0; j < repeate; j++)
        {
            dA[i + n * m * j] = data[i];
        }
    }

    seal::Plaintext plain_tmp;
    // cout << "the size of dA is: " << dA.size() << endl;
    encoder.encode(dA,scale,plain_tmp);
    // cout<<dA.size()<<endl;
    encryptor.encrypt(plain_tmp,cipher_matrix);
    
}

inline void cipher_matrix_coprime::dec_matrix_cipher(Matrix<double> &destination, seal::CKKSEncoder &encoder, seal::Decryptor &decryptor)
{
    //decrypte ciphertext
    Plaintext plain_tmp;
    vector<double> vec_tmp;
    decryptor.decrypt(cipher_matrix,plain_tmp);
    encoder.decode(plain_tmp,vec_tmp);

    destination=diagonal_decoding(vec_tmp,n,m);
}

inline cipher_matrix_coprime::~cipher_matrix_coprime()
{
    ;
}

inline void cipher_matrix_coprime::rotate_ct(int i, seal::Ciphertext &destination, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys, seal::CKKSEncoder &encoder)
{
    evaluator.rotate_vector(cipher_matrix,i,galois_keys,destination);
}

inline size_t cipher_matrix_coprime::get_rows()
{
    return n;
}

inline size_t cipher_matrix_coprime::get_cols()
{
    return m;
}

inline seal::Ciphertext cipher_matrix_coprime::get_cipher_matrix()
{
    return cipher_matrix;
}

inline void reorder_hoisted_B(vector<Ciphertext> &hoisted_B, size_t n, size_t m, size_t p, size_t r) 
{
    vector<Ciphertext> T;
    vector<int> order_list;
    for(size_t i = 0; i < m; i++) {
        order_list.push_back(int (mod(i * (r*m - n),  (m*p)) / p));
    }
    for (size_t i = 0; i < m; i++) {
        T.push_back(hoisted_B[order_list[i]]);
    }
    hoisted_B = T; // Update hoisted_B with the reordered elements     
}

inline void seg_hoisting(std::vector<Ciphertext> &destination, const Ciphertext &encrypt, std::vector<int> steps, std::size_t m,std::size_t t, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys) 
{
    vector<Ciphertext> hoisted_tmp;
    Ciphertext encrypt_cpoy = encrypt; 
    // size_t t = sqrt(m - 1); 
    size_t t_for = (m - 1) / t;
    // cout << "t_for: " << t_for << endl;
    int t_0 = m - 1 - t * t_for;
    // cout << "t_0: " << t_0 << endl;

    // cout << " hoisting steps: " << t << endl;
    // cout << " for times " << t_for << endl;
    // cout << " t_0:  " << t_0 << endl;

    destination.push_back(encrypt_cpoy);

    for (size_t i = 0; i < t_for; i++) {
        // cout << "i = " << i << endl;
        evaluator.rotate_vector_hoist(encrypt_cpoy, steps, galois_keys, hoisted_tmp);
        destination.insert(destination.end(), hoisted_tmp.begin(), hoisted_tmp.end());
        // cout << destination.size() << endl; 
        encrypt_cpoy = hoisted_tmp.back();; 
        hoisted_tmp.clear();
    }
    if (t_0 > 0) 
    {
        vector<int> step_s (steps.begin(), steps.begin() + t_0);
        evaluator.rotate_vector_hoist(encrypt_cpoy, step_s, galois_keys, hoisted_tmp);
        destination.insert(destination.end(), hoisted_tmp.begin(), hoisted_tmp.end());
    }
    // cout << destination.size() << endl;
}

void encrypted_coprime_matrix_multiplication(seal::CKKSEncoder &encoder, // seal::Decryptor &decryptor,
                                          double scale, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys,
                                          seal::GaloisKeys &galois_keys, cipher_matrix_coprime &encrypted_A,
                                          cipher_matrix_coprime &encrypted_B, cipher_matrix_coprime &destination, int method_choice)
{
    uint64_t r;
    size_t n=encrypted_A.get_rows();
    size_t m=encrypted_A.get_cols();
    size_t p=encrypted_B.get_cols();
    size_t slot_count=encoder.slot_count();
    // cout << "slot_count: " << slot_count << endl;
    // size_t t = sqrt(m-1);
    size_t t = (log2(slot_count))/2 - 1;
    // cout << "t: " << t << endl;

    if ((slot_count < n * m * (ceil(double(p) / m) + 1)) || (slot_count < m * p * (ceil(double(n) / m) + 1)))
    {
        cout << "slot_count:" << slot_count << endl;
        cout << n << m << ceil(double(p) / m) + 1 << endl;
        cout << "n * m * (ceil(double(p) / m) + 1):" << n * m * (ceil(double(p) / m) + 1) << endl;

        cerr << "  !!! the number of slot is not enough for the  (" << n << ", "
             << m << ", " << p << ") matrix multiplication" << endl;
        return;
    }
    r = smallest_r(n, m, p);
    // cout << "the smallest r is: " << r << endl;

    vector<Ciphertext>  cipher_vector;
    Ciphertext cipher_tmpA,cipher_tmpB,cipher_tmp,encrypted_C;//save intermediate variables
    Plaintext plain_tmp;
    vector<double> vec_tmp;
    // cout<<n<<"  "<<m<<"  "<<p<<endl;

    if (method_choice == 1)
    {
        // xyy add : seg-hoisting
        vector<int> A_steps, B_steps;
        vector<Ciphertext> hoisted_A, hoisted_B, hoisted_C;

        // auto start_time = chrono::high_resolution_clock::now();
        for (size_t i = 1; i < t + 1; i++) 
        {
            A_steps.push_back(mod(i * n, n*m));
            B_steps.push_back(mod(i * p, m*p));
        }
        // cout << "test: 1" << endl;
        cipher_tmpA = encrypted_A.get_cipher_matrix();
        cipher_tmpB = encrypted_B.get_cipher_matrix();
        // cout << "test: 2" << endl;
        // seg_hoisting(hoisted_A, cipher_tmpA, A_steps, m, evaluator, galois_keys);
        // seg_hoisting(hoisted_B, cipher_tmpB, B_steps, m, evaluator, galois_keys);
        seg_hoisting(hoisted_A, cipher_tmpA, A_steps, m, t, evaluator, galois_keys);
        seg_hoisting(hoisted_B, cipher_tmpB, B_steps, m, t, evaluator, galois_keys);
        // reorder 
        reverse(hoisted_A.begin() + 1, hoisted_A.end());
        reorder_hoisted_B(hoisted_B, n, m, p, r);

        // auto end_time = chrono::high_resolution_clock::now();
        // auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
        // cout << "The Rot timing cost is: " << time_diff.count() / 1e6 << " s" << endl;
        // cout << "hoisted_B size: " << hoisted_B.size() << endl;
        // cout << "test: 3" << endl;
        // for (size_t i = 0; i < hoisted_A.size(); i++) {
        //     evaluator.multiply(hoisted_A[i], hoisted_B[i], cipher_tmp);
        //     hoisted_C.push_back(cipher_tmp);
        // }
        // // cout << "test: 4" << endl;
        // evaluator.add_many(hoisted_C, encrypted_C);
        // xyy add : lazy reduction
        // auto start_time1 = chrono::high_resolution_clock::now();
        evaluator.multiply_lazy(hoisted_A, hoisted_B, encrypted_C);
        
        evaluator.relinearize_inplace(encrypted_C, relin_keys);
        evaluator.rescale_to_next_inplace(encrypted_C);
        // auto end_time1 = chrono::high_resolution_clock::now();
        // auto time_diff1 = chrono::duration_cast<chrono::microseconds>(end_time1 - start_time1);
        // cout << "The Mult+Add timing cost is: " << time_diff1.count() / 1e6 << " s" << endl;
    }

    // xyy add: for hoisting
    // vector<int> A_steps, B_steps;
    // vector<Ciphertext> hoisted_A, hoisted_B, hoisted_C;

    // for (size_t i = 0; i < m; i++) {
    //     A_steps.push_back(mod(-i * n, n*m));
    //     B_steps.push_back(mod(i * (r*m - n),  (m*p)));
    // }
    
    // cipher_tmpA = encrypted_A.get_cipher_matrix();
    // cipher_tmpB = encrypted_B.get_cipher_matrix();

    // evaluator.rotate_vector_hoist(cipher_tmpA, A_steps, galois_keys, hoisted_A);
    // evaluator.rotate_vector_hoist(cipher_tmpB, B_steps, galois_keys, hoisted_B);

    // for (size_t i = 0; i < A_steps.size(); i++) {
    //     evaluator.multiply(hoisted_A[i], hoisted_B[i], cipher_tmp);
    //     hoisted_C.push_back(cipher_tmp);
    // }
    // evaluator.add_many(hoisted_C, encrypted_C);
    // evaluator.relinearize_inplace(encrypted_C, relin_keys);
    // evaluator.rescale_to_next_inplace(encrypted_C);

    if (method_choice == 2)
    {
        // chen : native method
        // auto start_time = chrono::high_resolution_clock::now();

        vector<Ciphertext> test_A, test_B;
        for(size_t i = 0; i < m; i++) {
            seal::Ciphertext ct_a, ct_b, ct_tmp;
            encrypted_A.rotate_ct(mod(-i * n, n*m),ct_a,evaluator,galois_keys,relin_keys,encoder);
            // cout << mod(-i * n, n*m) << ", ";
            encrypted_B.rotate_ct(mod(i * (r*m - n),  (m*p)),ct_b,evaluator,galois_keys,relin_keys,encoder);
            // cout << mod(i * (r*m - n),  (m*p)) << ", ";
            test_A.push_back(ct_a);
            test_B.push_back(ct_b);

            evaluator.multiply(ct_a, ct_b, ct_tmp);      
            cipher_vector.push_back(ct_tmp);
        }
        // auto end_time = chrono::high_resolution_clock::now();
        // auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
        // cout << "The Rot timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

        // auto start_time1 = chrono::high_resolution_clock::now();
        // for(size_t i = 0; i < m; i++) {
        //     evaluator.multiply(test_A[i], test_B[i], cipher_tmp);
        //     cipher_vector.push_back(cipher_tmp);
        // }
        evaluator.add_many(cipher_vector, encrypted_C);
        evaluator.relinearize_inplace(encrypted_C, relin_keys);
        evaluator.rescale_to_next_inplace(encrypted_C);

        // auto end_time1 = chrono::high_resolution_clock::now();
        // auto time_diff1 = chrono::duration_cast<chrono::microseconds>(end_time1 - start_time1);
        // cout << "The Mult+Add timing cost is: " << time_diff1.count() / 1e6 << " s" << endl;
    }
    cipher_matrix_coprime ct_C(encrypted_C,n,p);
    destination=ct_C;
}


// Function: Returns the largest power of two less than m
size_t largestPowerOfTwoLessThan(size_t m) {
    int result = m - 1; // Start with m-1 to ensure the result is always less than m

    // Gradually extend the 1s in the binary representation
    result |= result >> 1;
    result |= result >> 2;
    result |= result >> 4;
    result |= result >> 8;
    result |= result >> 16;

    // For systems where int is larger than 32 bits, include this line to cover all bits
    // result |= result >> 32;

    return result - (result >> 1); // Subtract itself after a shift operation to get the result
}


void encrypted_coprime_matrix_multiplication_imporve_log(seal::CKKSEncoder &encoder, // seal::Decryptor &decryptor,
                                                         double scale, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys,
                                                         seal::GaloisKeys &galois_keys, cipher_matrix_coprime &encrypted_A,
                                                         cipher_matrix_coprime &encrypted_B, cipher_matrix_coprime &destination)
{
    size_t n=encrypted_A.get_rows();
    size_t m=encrypted_A.get_cols();
    size_t p=encrypted_B.get_cols();
    size_t slot_count=encoder.slot_count();

    if (slot_count < n * m * p )
    {

        cerr << "  !!! the number of slot is not enough for the  (" << n << ", "
             << m << ", " << p << ") matrix multiplication" << endl;
        return;
    }


    vector<Ciphertext>  cipher_vector;
    Ciphertext ct_a,ct_b,ct_tmp,encrypted_C;//save intermediate variables
    Plaintext plain_tmp;
    vector<double> vec_tmp;
    ct_a=encrypted_A.get_cipher_matrix();
    ct_b=encrypted_B.get_cipher_matrix();
    // xyy add:
    vector<int> steps;
    for (size_t i = 1; i < m; i = i * 2)
    {
        steps.push_back(i * n * p);
    }

    // If m is a power of two, directly add using the algorithm
    if(m > 0 && (m & (m - 1)) == 0){
        evaluator.multiply(ct_a,ct_b,encrypted_C);
        evaluator.relinearize_inplace(encrypted_C,relin_keys);
        evaluator.rescale_to_next_inplace(encrypted_C);
        for(size_t i=1;i<m;i=i*2){
            evaluator.rotate_vector(encrypted_C,i*n*p,galois_keys,ct_tmp);
            evaluator.add_inplace(encrypted_C,ct_tmp);
        }

    }
    // SegSum
    else{
        size_t tmp_m=largestPowerOfTwoLessThan(m);
        // cout<<tmp_m<<endl;
        evaluator.multiply(ct_a,ct_b,encrypted_C);
        evaluator.relinearize_inplace(encrypted_C,relin_keys);
        evaluator.rescale_to_next_inplace(encrypted_C);

        vector<double> vec_1(tmp_m*n*p,1),vec_2((m-tmp_m)*n*p,1);
        Plaintext plain_1,plain_2;
        encoder.encode(vec_1,scale,plain_1);
        encoder.encode(vec_2,scale,plain_2);
        evaluator.mod_switch_to_inplace(plain_1,encrypted_C.parms_id());
        evaluator.mod_switch_to_inplace(plain_2,encrypted_C.parms_id());
        Ciphertext cipher_1,cipher_2;
        evaluator.multiply_plain(encrypted_C,plain_1,cipher_1);
        evaluator.rotate_vector(encrypted_C,tmp_m*n*p,galois_keys,cipher_2);
        evaluator.multiply_plain_inplace(cipher_2,plain_2);
        evaluator.add(cipher_1,cipher_2,encrypted_C);
        evaluator.rescale_to_next_inplace(encrypted_C);

        for(size_t i=1;i<tmp_m;i=i*2){
            evaluator.rotate_vector(encrypted_C,i*n*p,galois_keys,ct_tmp);
            evaluator.add_inplace(encrypted_C,ct_tmp);
        }

    }

    
    cipher_matrix_coprime ct_C(encrypted_C,n,p);
    destination=ct_C;
}


// void encrypted_coprime_matrix_multiplication_improve_case1(seal::CKKSEncoder &encoder, // seal::Decryptor &decryptor,
//                                           double scale, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys,
//                                           seal::GaloisKeys &galois_keys, cipher_matrix_coprime &encrypted_A,
//                                           cipher_matrix_coprime &encrypted_B, cipher_matrix_coprime &destination)
// {
//     size_t n=encrypted_A.get_rows();
//     size_t m=encrypted_A.get_cols();
//     size_t p=encrypted_B.get_cols();
//     size_t slot_count=encoder.slot_count();

//     if ((slot_count < n * m * (ceil(double(p) / m) + 1)) || (slot_count < m * p * (ceil(double(n) / m) + 1)))
//     {

//         cerr << "  !!! the number of slot is not enough for the  (" << n << ", "
//              << m << ", " << p << ") matrix multiplication" << endl;
//         return;
//     }

//     vector<Ciphertext>  cipher_vector;
//     Ciphertext cipher_tmpA,cipher_tmpB,cipher_tmp,encrypted_C,ct_last;//save intermediate variables
//     Plaintext plain_tmp;
//     vector<double> vec_tmp;
//     // cout<<n<<"  "<<m<<"  "<<p<<endl;

//     // for (size_t i = 0; i < m; i+=2) {
//     //     // cout<<"i:"<<i<<"  ( "<<mod(i * n * p, n*m)<<" , "<<mod(i * n * p, m*p)<<" )"<<endl;
//     //     seal::Ciphertext ct_a, ct_b, ct_tmp;
//     //     encrypted_A.rotate_ct(mod(i * n * p, n*m),ct_a,evaluator,galois_keys,relin_keys,encoder);
//     //     encrypted_B.rotate_ct(mod(i * n * p, m*p),ct_b,evaluator,galois_keys,relin_keys,encoder);
//     //     evaluator.multiply(ct_a, ct_b, ct_tmp);
//     //     if(i==m-1 && m%2==1){
//     //         evaluator.relinearize_inplace(ct_tmp,relin_keys);
//     //         evaluator.rescale_to_next_inplace(ct_tmp);
//     //         ct_last=ct_tmp;
//     //         break;
//     //     }      
//     //     cipher_vector.push_back(ct_tmp);
//     // }
    
//     // evaluator.add_many(cipher_vector, encrypted_C);
//     // evaluator.relinearize_inplace(encrypted_C, relin_keys);
//     // evaluator.rescale_to_next_inplace(encrypted_C);

//     // evaluator.rotate_vector(encrypted_C,n*p,galois_keys,cipher_tmp);
//     // evaluator.add_inplace(encrypted_C,cipher_tmp);

//     // if(m%2==1){
//     //     ct_last.scale()=encrypted_C.scale();
//     //     evaluator.add_inplace(encrypted_C,ct_last);
//     // }

//     for (size_t i = 0; i < m; i++) {
//         if(2*n*p+mod(i * n * p, n*m)>2*n*m||2*n*p+mod(i * n * p, m*p)>2*m*p){
//             cout<<"i:"<<i<<"  ( "<<mod(i * n * p, n*m)<<" , "<<mod(i * n * p, m*p)<<" )"<<endl;
//         }
        
//             // cout<<"i:"<<i<<"  ( "<<mod(i * n * p, n*m)<<" , "<<mod(i * n * p, m*p)<<" )"<<endl;
        
//         seal::Ciphertext ct_a, ct_b, ct_tmp;
//         encrypted_A.rotate_ct(mod(i * n * p, n*m),ct_a,evaluator,galois_keys,relin_keys,encoder);
//         encrypted_B.rotate_ct(mod(i * n * p, m*p),ct_b,evaluator,galois_keys,relin_keys,encoder);
//         evaluator.multiply(ct_a, ct_b, ct_tmp);     
//         cipher_vector.push_back(ct_tmp);
//     }
    
//     evaluator.add_many(cipher_vector, encrypted_C);
//     evaluator.relinearize_inplace(encrypted_C, relin_keys);
//     evaluator.rescale_to_next_inplace(encrypted_C);

//     cipher_matrix_coprime ct_C(encrypted_C,n,p);
//     destination=ct_C;
// }