#include <iostream>
#include <seal/seal.h>
#include "matrix.h"
#include "utils.h"

using namespace std;
using namespace seal;

class long_cipher
{
private:
    /* data */
    size_t n; // matrix rows
    size_t m; // matrix cols
    vector<seal::Ciphertext> LongCipher;
    size_t l; // Cipher parameter slot length

protected:
    // split plain vector
    std::vector<std::vector<double>> splitVector(const std::vector<double> &input, size_t length);
    // generate select vector
    std::vector<double> generateVector(int n, int m);
    // genetate plain vector to select number
    void gen_plain_vec(int start_index, int end_index, double scale,
                       seal::CKKSEncoder &encoder, seal::Plaintext &destination);

public:
    long_cipher(/* args */);
    long_cipher(int n, int m,double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    long_cipher(int n, int m, vector<double> &plain_value, double scale,
                seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    long_cipher(int n, int m, vector<seal::Ciphertext> &cipher_value, int l);
    ~long_cipher();

    // display data
    void print_data();

    // 1. encrypte plain data
    void enc_diag_vector(int n, int m, vector<double> &plain_value, double scale,
                         seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    // 2. decrypt cipher data to plain data
    void dec_diag_vector(vector<double> &plain_value, seal::CKKSEncoder &encoder, seal::Decryptor &decryptor);

    /*********************************************************************
     * 3. LongRot
     * In fact, the LongRot operation needs to be divided into two situations,
     * namely the size relationship between r and v.(Divided into three types
     * in the experiment, greater than, less than, and equal to).
     * Each situation needs to be divided into 0 to w-u-1 and w-u-1 to w-1 and
     * considered separately. Special attention needs to be paid to the opera-
     * -tion of item w-1.
     **********************************************************************/
    void LongRot(int k, int tau, vector<seal::Ciphertext> &destination,
                 double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator,
                 seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys);
    vector<seal::Ciphertext> LongRot(int rotSize, int length, double scale,
                                     seal::CKKSEncoder &encoder, seal::Evaluator &evaluator,
                                     seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys);

     /*********************************************************************
     * 3. AutoLongRot
     * In practice, various factors need to be considered when performing rotations.
     * Depending on the length of 'length', adjustments may be made to the rotation result.
     * Therefore, the algorithm cannot provide a specific number of rotations.
     **********************************************************************/
    void AutoLongRot(int rotSize, int length, vector<seal::Ciphertext> &destination,
                     double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator,
                     seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys);
    

    // 4.LongAdd: add all LongRot result
    void LongAdd(vector<seal::Ciphertext> &cipher_vec, seal::Evaluator &evaluator);

    // plain select 
    // In fact, when v_tmp==0 or v_tmp==l, not need multiply and rotate, just mod_switch_to_next. 
    void gen_select_plain(size_t v_tmp,double scale,seal::CKKSEncoder &encoder,Plaintext &destination1,Plaintext &destination2, int &flag_mul);
    void gen_select_plain(size_t start_index, size_t v_tmp,size_t end_index,double scale,seal::CKKSEncoder &encoder,Plaintext &destination1,Plaintext &destination2, int &flag_mul); 
    void select_mul(int &flag_mul,seal::Plaintext &plain1,seal::Ciphertext &cipher1,
                    seal::Plaintext &plain2,seal::Ciphertext cipher2,seal::Ciphertext &destination,
                    seal::Evaluator &evaluator,seal::RelinKeys &relin_keys);

    long_cipher &operator=(const long_cipher &other);

    size_t get_rows() const;
    size_t get_cols() const;
    size_t get_cipher_length() const;
    vector<seal::Ciphertext> get_cipher() const;
    vector<seal::Ciphertext> get_cipher();
    size_t get_slot() const;
    size_t g_func(int x, int y) const;
};

inline std::vector<std::vector<double>> long_cipher::splitVector(const std::vector<double> &input, size_t length)
{
    std::vector<std::vector<double>> result;
    size_t size = input.size();

    for (size_t i = 0; i < size; i += length)
    {
        std::vector<double> subVector;
        subVector.insert(subVector.end(), input.begin() + i, input.begin() + std::min(i + length, size));
        result.push_back(subVector);
    }
    return result;
}

inline std::vector<double> long_cipher::generateVector(int n, int m)
{
    std::vector<double> result(l, 0);
    // generate n-m number of 1
    for (int i = n; i < m; i++)
    {
        result[i] = 1;
    }
    return result;
}

inline void long_cipher::gen_plain_vec(int start_index, int end_index, double scale, seal::CKKSEncoder &encoder, seal::Plaintext &destination)
{
    vector<double> select = generateVector(start_index, end_index);
    encoder.encode(select, scale, destination);
}

long_cipher::long_cipher(/* args */)
{
    n = 0;
    m = 0;
    l = 0;
    LongCipher.resize(0);
}

inline long_cipher::long_cipher(int n, int m,double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    this->n = n;
    this->m = m;
    this->l = encoder.slot_count();
    vector<double> vec(l,0);
    // print_vector(vec);
    Plaintext plain_tmp;
    Ciphertext cipher_tmp;
    encoder.encode(vec,scale,plain_tmp);
    encryptor.encrypt(plain_tmp,cipher_tmp);
    for(int i=0;i<ceil(double(n*m)/double(l));i++){
        this->LongCipher.push_back(cipher_tmp);
    }
    // this->LongCipher.resize(ceil(double(n*m)/double(l)));
}

inline long_cipher::long_cipher(int n, int m, vector<double> &plain_value, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    enc_diag_vector(n, m, plain_value, scale, encoder, encryptor);
}

inline long_cipher::long_cipher(int n, int m, vector<seal::Ciphertext> &cipher_value, int l)
{
    this->n = n;
    this->m = m;
    this->LongCipher = cipher_value;
    this->l = l;
}

long_cipher::~long_cipher()
{
}

inline void long_cipher::print_data()
{
    std::cout << "Matrix dimension:(" << n << "," << m << ")" << std::endl;
    std::cout << "Number of slots: " << l << std::endl;
    std::cout << "Number of cipher: " << LongCipher.size() << std::endl;
}

inline void long_cipher::enc_diag_vector(int n, int m, vector<double> &plain_value, double scale,
                                         seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    this->n = n;
    this->m = m;
    l = encoder.slot_count();
    vector<vector<double>> splitvec = splitVector(plain_value, l);
    seal::Plaintext plain_tmp;
    seal::Ciphertext cipher_tmp;
    for (size_t i = 0; i < splitvec.size(); i++)
    {
        encoder.encode(splitvec[i], scale, plain_tmp);
        encryptor.encrypt(plain_tmp, cipher_tmp);
        LongCipher.push_back(cipher_tmp);
    }
}

inline void long_cipher::dec_diag_vector(vector<double> &plain_value, seal::CKKSEncoder &encoder, seal::Decryptor &decryptor)
{
    int plain_size = n * m;
    // cout<<plain_size<<endl;
    vector<double> vec_tmp;
    seal::Plaintext plain_tmp;
    for (size_t i = 0; i < LongCipher.size(); i++)
    {
        decryptor.decrypt(LongCipher[i], plain_tmp);
        encoder.decode(plain_tmp, vec_tmp);
        if (i == LongCipher.size() - 1)
        {
            plain_value.insert(plain_value.end(), vec_tmp.begin(), vec_tmp.begin() + plain_size - i * l);
        }
        else
        {
            plain_value.insert(plain_value.end(), vec_tmp.begin(), vec_tmp.end());
        }
    }
}

inline void long_cipher::LongRot(int k, int tau, vector<seal::Ciphertext> &destination,
                                 double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator,
                                 seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys)
{
    // In fact, Each time only the elements of length need to be considered for rotation
    // int rot_num=ceil(double(length)/double(length));
    int d = n * m;
    int w = LongCipher.size();
    int z = d % l;

    int u = k / l;
    int v = k % l;

    int r = (tau - (d - k)) / d;
    int delta = (tau - (d - k)) % d;

    int alpha = delta / l;
    int beta = delta % l;

    int iota = tau / l;
    int kappa = tau % l;

    vector<int> v_vector(1,v);
    int h=0;

    // cout<<"  d:"<<d<<"  u:"<<u<<"  v:"<<v<<"  r:"<<r<<"  tau:"<<tau<<"  z:"<<z
    //     <<"  delta:"<<delta<<"  alpha:"<<alpha<<"  beta:"<<beta<<"  iota:"<<iota<<"  kappa:"<<kappa<<endl;
    // cout<<"  k:"<<k<<"  w:"<<w<<"  w - u + g(z, v0) - 2:"<<w - u + g_func(z, v_vector[0]) - 2<<endl;


    // cout<<"r:"<<r<<"  v:"<<v<<"  u:"<<u<<"  w:"<<w<<endl;

    //+ : left | {1,2,3,4,5,6,0...0}  (2)---> {3,4,5,6,0...0,1,2}
    //- : right| {1,2,3,4,5,6,0...0}  (-2)---> {0,0,1,2,3,4,5,6,0...0}

    seal::Ciphertext cipher_tmp,cipher0,cipher1,cipher2,cipher3,cipher_a_0;
    seal::Plaintext plain_tmp,plain0,plain1,plain2,plain3;
    vector<seal::Ciphertext> cipher_vector,cipher_result;
    int flag_mul=0;

    // case tau<= d-k
    if ((tau - (d - k)) <= 0)
    {
        // a)
        if (iota == 0 && tau <= int(l - v))
        {
            evaluator.rotate_vector(LongCipher[u], v, galois_keys, cipher0);
            gen_plain_vec(0, tau, scale, encoder, plain_tmp);
            evaluator.multiply_plain(cipher0, plain_tmp, cipher_tmp);
            evaluator.relinearize_inplace(cipher_tmp, relin_keys);
            evaluator.rescale_to_next_inplace(cipher_tmp);
            cipher_result.push_back(cipher_tmp);
        }
        // b)
        else if (iota == 0 && tau > int(l - v))
        {
            evaluator.rotate_vector(LongCipher[u], v, galois_keys, cipher0);
            evaluator.rotate_vector(LongCipher[u + 1], -(l - v), galois_keys, cipher1);
            gen_plain_vec(0, l - v, scale, encoder, plain0);
            gen_plain_vec(l - v, tau, scale, encoder, plain1);
            evaluator.multiply_plain_inplace(cipher0, plain0);
            evaluator.multiply_plain_inplace(cipher1, plain1);
            evaluator.add(cipher0, cipher1, cipher_tmp);
            evaluator.relinearize_inplace(cipher_tmp, relin_keys);
            evaluator.rescale_to_next_inplace(cipher_tmp);
            cipher_result.push_back(cipher_tmp);
        }
        // c)
        else if(iota>0){
            // I.
            for(int i=0;i<iota+1;i++){
                evaluator.rotate_vector(LongCipher[u+i],v,galois_keys,cipher_tmp);
                cipher_vector.push_back(cipher_tmp);
            }
            // destination=cipher_vector;
            // II.
            // gen_plain_vec(0,l-v,scale,encoder,plain0);
            // gen_plain_vec(l-v,l,scale,encoder,plain1);
            gen_select_plain(v,scale,encoder,plain0,plain1,flag_mul);
            for(int i=0;i<iota;i++){
                // evaluator.multiply_plain(cipher_vector[i],plain0,cipher0);
                // evaluator.multiply_plain(cipher_vector[i+1],plain1,cipher1);
                // evaluator.add(cipher0,cipher1,cipher_tmp);
                // evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                // evaluator.rescale_to_next_inplace(cipher_tmp);
                select_mul(flag_mul,plain0,cipher_vector[i],plain1,cipher_vector[i+1],cipher_tmp,evaluator,relin_keys);
                cipher_result.push_back(cipher_tmp);
            }
            // III.
            h=h+iota;
            if(kappa<=int(l-v)){
                gen_plain_vec(0,kappa,scale,encoder,plain0);
                evaluator.multiply_plain(cipher_vector[iota],plain0,cipher_tmp);
                evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                evaluator.rescale_to_next_inplace(cipher_tmp);
                cipher_result.push_back(cipher_tmp);
            }
            else{
                evaluator.rotate_vector(LongCipher[u+tau+1],-(l-v),galois_keys,cipher_tmp);
                cipher_vector.push_back(cipher_tmp);
                gen_plain_vec(0,l-v,scale,encoder,plain0);
                gen_plain_vec(l-v,kappa,scale,encoder,plain1);
                evaluator.multiply_plain(cipher_vector[tau],plain0,cipher0);
                evaluator.multiply_plain(cipher_vector[tau+1],plain1,cipher1);
                evaluator.add(cipher0,cipher1,cipher_tmp);
                evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                evaluator.rescale_to_next_inplace(cipher_tmp);
                cipher_result.push_back(cipher_tmp);
            }
        }
    }
    else if((tau - (d - k)) > 0){
        // a)
        vector<int> i_vector(1, w - u + g_func(z, v_vector[0]) - 2);
        // b)
        for (int j = 1; j <= r + 1; j++)
        {
            if (z > v_vector[j - 1])
            {
                v_vector.push_back(l - (z - v_vector[j - 1]));
            }
            else
            {
                v_vector.push_back(v_vector[j - 1] - z);
            }
            i_vector.push_back(w + g_func(z, v_vector[j]) - 2);
        }
        for(size_t i=0;i<v_vector.size();i++){
            cout<<v_vector[i]<<"  ";
        }
        cout<<endl;
        for(size_t i=0;i<i_vector.size();i++){
            cout<<i_vector[i]<<"  ";
        }
        cout<<endl;
        // c)
        // cout<<"case c"<<endl;
        for(int i=0;i<=i_vector[0];i++){
            evaluator.rotate_vector(LongCipher[u+i],v_vector[0],galois_keys,cipher_tmp);
            cipher_vector.push_back(cipher_tmp);
        }
        // gen_plain_vec(0,l-v_vector[0],scale,encoder,plain0);
        // gen_plain_vec(l-v_vector[0],l,scale,encoder,plain1);
        gen_select_plain(v_vector[0],scale,encoder,plain0,plain1,flag_mul);
        for(int i=0;i<i_vector[0];i++){
            // evaluator.multiply_plain(cipher_vector[i],plain0,cipher0);
            // evaluator.multiply_plain(cipher_vector[i+1],plain1,cipher1);
            // evaluator.add(cipher0,cipher1,cipher_tmp);
            // evaluator.relinearize_inplace(cipher_tmp,relin_keys);
            // evaluator.rescale_to_next_inplace(cipher_tmp);
            select_mul(flag_mul,plain0,cipher_vector[i],plain1,cipher_vector[i+1],cipher_tmp,evaluator,relin_keys);
            cipher_result.push_back(cipher_tmp);
        }
        // evaluator.rotate_vector(LongCipher[0],v_vector[0]-z,galois_keys,cipher_tmp);
        // cipher_vector[0]=cipher_tmp;
        evaluator.rotate_vector(LongCipher[0],v_vector[0]-z,galois_keys,cipher_a_0);
        h=h+i_vector[0];
        // d)
        if(z>v_vector[0]){
            if(int((i_vector[0]+1)*l)>tau && beta<v_vector[1] && r==0){
                // cout<<" case d-1"<<endl;
                gen_plain_vec(0,z-v_vector[0],scale,encoder,plain0);
                gen_plain_vec(z-v_vector[0],z-v_vector[0]+beta,scale,encoder,plain1);
                evaluator.multiply_plain(cipher_vector[i_vector[0]],plain0,cipher0);
                evaluator.multiply_plain(cipher_a_0,plain1,cipher1);
                evaluator.add(cipher0,cipher1,cipher_tmp);
                evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                evaluator.rescale_to_next_inplace(cipher_tmp);
                cipher_result.push_back(cipher_tmp);
            }
            else{
                // cout<<" case d-2"<<endl;
                gen_plain_vec(0,z-v_vector[0],scale,encoder,plain0);
                gen_plain_vec(z-v_vector[0],z-v_vector[0]+v_vector[1],scale,encoder,plain1);
                evaluator.multiply_plain(cipher_vector[i_vector[0]],plain0,cipher0);
                evaluator.multiply_plain(cipher_a_0,plain1,cipher1);
                evaluator.add(cipher0,cipher1,cipher_tmp);
                evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                evaluator.rescale_to_next_inplace(cipher_tmp);
                cipher_result.push_back(cipher_tmp);
            }
        }
        // e)
        else if(z<=v_vector[0]){
            // cout<<" case e"<<endl;
            evaluator.rotate_vector(LongCipher[w-1],-(l-v_vector[0]),galois_keys,cipher_tmp);
            cipher_vector.push_back(cipher_tmp);

            gen_plain_vec(0,l-v_vector[0],scale,encoder,plain0);
            gen_plain_vec(l-v_vector[0],l-v_vector[0]+z,scale,encoder,plain1);
            evaluator.multiply_plain(cipher_vector[i_vector[0]],plain0,cipher0);
            evaluator.multiply_plain(cipher_vector[i_vector[0]+1],plain1,cipher1);
            evaluator.add_inplace(cipher0,cipher1);
            if(int((h+1)*l)>tau && beta<v_vector[1] && r==0){
                // cout<<" case e-1"<<endl;
                if(l-v_vector[0]!=l-v_vector[0]+z+beta){
                    gen_plain_vec(l-v_vector[0]+z,l-v_vector[0]+z+beta,scale,encoder,plain2);
                    evaluator.multiply_plain(cipher_a_0,plain2,cipher2);
                    evaluator.add(cipher0,cipher2,cipher_tmp);
                    evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                    evaluator.rescale_to_next_inplace(cipher_tmp);
                    cipher_result.push_back(cipher_tmp);    
                }
                else{
                    evaluator.relinearize_inplace(cipher0,relin_keys);
                    evaluator.rescale_to_next_inplace(cipher0);
                    cipher_result.push_back(cipher0);
                }
            }
            else{
                // cout<<" case e-2"<<endl;
                if(l-v_vector[0]+z!=l-v_vector[0]+z+v_vector[1]){
                    gen_plain_vec(l-v_vector[0]+z,l-v_vector[0]+z+v_vector[1],scale,encoder,plain2);
                    evaluator.multiply_plain(cipher_a_0,plain2,cipher2);
                    evaluator.add(cipher0,cipher2,cipher_tmp);
                    evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                    evaluator.rescale_to_next_inplace(cipher_tmp);
                    cipher_result.push_back(cipher_tmp);
                }
                else{
                    evaluator.relinearize_inplace(cipher0,relin_keys);
                    evaluator.rescale_to_next_inplace(cipher0);
                    cipher_result.push_back(cipher0);
                }
            }
        }
        cipher_vector[0]=cipher_a_0;
        // f)
        for(int j=1;j<=r;j++){
            // cout<<"case f"<<endl;
            // I.
            cipher_vector.resize(i_vector[j]);
            for(int i=1;i<=i_vector[j];i++){
                evaluator.rotate_vector(LongCipher[i],v_vector[j],galois_keys,cipher_tmp);
                cipher_vector[i]=cipher_tmp;
            }
            for(int i=0;i<i_vector[j];i++){
                // gen_plain_vec(0,l-v_vector[j],scale,encoder,plain0);
                // gen_plain_vec(l-v_vector[j],l,scale,encoder,plain1);
                gen_select_plain(v_vector[j],scale,encoder,plain0,plain1,flag_mul);
                // evaluator.multiply_plain(cipher_vector[i],plain0,cipher0);
                // evaluator.multiply_plain(cipher_vector[i+1],plain1,cipher1);
                // evaluator.add(cipher0,cipher1,cipher_tmp);
                // evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                // evaluator.rescale_to_next_inplace(cipher_tmp);
                select_mul(flag_mul,plain0,cipher_vector[i],plain1,cipher_vector[i+1],cipher_tmp,evaluator,relin_keys);
                cipher_result.push_back(cipher_tmp);
            }
            // evaluator.rotate_vector(LongCipher[0],v_vector[j]-z,galois_keys,cipher_tmp);
            // cipher_vector[0]=cipher_tmp;
            evaluator.rotate_vector(LongCipher[0],v_vector[j]-z,galois_keys,cipher_a_0);
            h=h+i_vector[j]+1;
            // II.
            if(z>v_vector[j]){
                if(int((h+1)*l)>tau &&beta<=v_vector[1] && j==r){
                    gen_plain_vec(0,z-v_vector[0],scale,encoder,plain0);
                    gen_plain_vec(z-v_vector[0],z-v_vector[0]+beta,scale,encoder,plain1);
                    evaluator.multiply_plain(cipher_vector[i_vector[j]],plain0,cipher0);
                    evaluator.multiply_plain(cipher_a_0,plain1,cipher1);
                    evaluator.add(cipher0,cipher1,cipher_tmp);
                    evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                    evaluator.rescale_to_next_inplace(cipher_tmp);
                    cipher_result.push_back(cipher_tmp);
                }
                else{
                    gen_plain_vec(0,z-v_vector[j],scale,encoder,plain0);
                    gen_plain_vec(z-v_vector[j],z-v_vector[j]+v_vector[j+1],scale,encoder,plain1);
                    evaluator.multiply_plain(cipher_vector[i_vector[j]],plain0,cipher0);
                    evaluator.multiply_plain(cipher_a_0,plain1,cipher1);
                    evaluator.add(cipher0,cipher1,cipher_tmp);
                    evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                    evaluator.rescale_to_next_inplace(cipher_tmp);
                    cipher_result.push_back(cipher_tmp);
                }
            }
            // III.
            else if(z<=v_vector[j]){
                evaluator.rotate_vector(LongCipher[w-1],-(l-v_vector[j]),galois_keys,cipher_tmp);
                cipher_vector.push_back(cipher_tmp);

                gen_plain_vec(0,l-v_vector[j],scale,encoder,plain0);
                gen_plain_vec(l-v_vector[j],l-v_vector[j]+z,scale,encoder,plain1);
                evaluator.multiply_plain(cipher_vector[i_vector[j]],plain0,cipher0);
                evaluator.multiply_plain(cipher_vector[i_vector[j]+1],plain1,cipher1);
                evaluator.add_inplace(cipher0,cipher1);
                if(int((h+1)*l)>tau && beta<v_vector[1] && j==r){
                    gen_plain_vec(l-v_vector[j]+z,l-v_vector[j]+z+beta,scale,encoder,plain2);
                    evaluator.multiply_plain(cipher_a_0,plain2,cipher2);
                    evaluator.add(cipher0,cipher2,cipher_tmp);
                    evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                    evaluator.rescale_to_next_inplace(cipher_tmp);
                    cipher_result.push_back(cipher_tmp);
                }
                else{
                    // cout<<" case e-2"<<endl;
                    gen_plain_vec(l-v_vector[j]+z,l-v_vector[j]+z+v_vector[j+1],scale,encoder,plain2);
                    evaluator.multiply_plain(cipher_a_0,plain2,cipher2);
                    evaluator.add(cipher0,cipher2,cipher_tmp);
                    evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                    evaluator.rescale_to_next_inplace(cipher_tmp);
                    cipher_result.push_back(cipher_tmp);
                }
            }
            cipher_vector[0]=cipher_a_0;
        }
        // g)
        // cout<<"case g"<<endl;
        cipher_vector.resize(alpha);
        for(int i=1;i<alpha;i++){
            evaluator.rotate_vector(LongCipher[i],v_vector[r+1],galois_keys,cipher_tmp);
            cipher_vector[i]=cipher_tmp;
        }
        // cout<<"case g-1"<<endl;
        // gen_plain_vec(0,l-v_vector[r+1],scale,encoder,plain0);
        // gen_plain_vec(l-v_vector[r+1],l,scale,encoder,plain1);
        gen_select_plain(v_vector[r+1],scale,encoder,plain0,plain1,flag_mul);
        // cout<<"alpha-1:"<<alpha-1<<"  size:"<<cipher_vector.size()<<endl;
        for(int i=0;i<alpha-1;i++){
            // evaluator.multiply_plain(cipher_vector[i],plain0,cipher0);
            // evaluator.multiply_plain(cipher_vector[i+1],plain1,cipher1);
            // evaluator.add(cipher0,cipher1,cipher_tmp);
            // evaluator.relinearize_inplace(cipher_tmp,relin_keys);
            // evaluator.rescale_to_next_inplace(cipher_tmp);
            select_mul(flag_mul,plain0,cipher_vector[i],plain1,cipher_vector[i+1],cipher_tmp,evaluator,relin_keys);
            cipher_result.push_back(cipher_tmp);
        }
        // cout<<"case g-2"<<endl;
        h=h+std::max(0,alpha-1);

        // h)
        if(alpha==0 && beta>v_vector[r+1]){
            // cout<<"case h"<<endl;
            gen_plain_vec(0,beta-v_vector[r+1],scale,encoder,plain0);
            evaluator.multiply_plain(cipher_vector[0],plain0,cipher_tmp);
            evaluator.relinearize_inplace(cipher_tmp,relin_keys);
            evaluator.rescale_to_next_inplace(cipher_tmp);
            cipher_result.push_back(cipher_tmp);
            h=h+1;
        }
        // i)
        else if(alpha>0){
            // cout<<"case i"<<endl;
            // cout<<"step:"<<-(int(l)-v_vector[r+1])%l<<"  v_vector[r+1]:"<<v_vector[r+1]<<endl;
            evaluator.rotate_vector(LongCipher[alpha],-(int(l)-v_vector[r+1])%l,galois_keys,cipher_tmp);
            cipher_vector.push_back(cipher_tmp);
            if(beta>v_vector[r+1]){
                // gen_plain_vec(0,l-v_vector[r+1],scale,encoder,plain0);
                // gen_plain_vec(l-v_vector[r+1],l,scale,encoder,plain1);
                gen_select_plain(v_vector[r+1],scale,encoder,plain0,plain1,flag_mul);
                // evaluator.multiply_plain(cipher_vector[alpha-1],plain0,cipher0);
                // evaluator.multiply_plain(cipher_vector[alpha],plain1,cipher1);
                // evaluator.add(cipher0,cipher1,cipher_tmp);
                // evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                // evaluator.rescale_to_next_inplace(cipher_tmp);
                select_mul(flag_mul,plain0,cipher_vector[alpha-1],plain1,cipher_vector[alpha],cipher_tmp,evaluator,relin_keys);
                cipher_result.push_back(cipher_tmp);

                gen_plain_vec(0,beta-v_vector[r+1],scale,encoder,plain2);
                evaluator.multiply_plain(cipher_vector[alpha],plain2,cipher_tmp);
                evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                evaluator.rescale_to_next_inplace(cipher_tmp);
                cipher_result.push_back(cipher_tmp);
                h=h+2;
            }
            else if(beta<=v_vector[r+1]){
                // cout<<"case i-2"<<endl;
                // gen_plain_vec(0,l-v_vector[r+1],scale,encoder,plain0);
                // gen_plain_vec(l-v_vector[r+1],l-v_vector[r+1]+beta,scale,encoder,plain1);
                // cout<<"l-v_vector[r+1]:"<<l-v_vector[r+1]<<"  l-v_vector[r+1]+beta:"<<l-v_vector[r+1]+beta<<endl;
                gen_select_plain(0,l-v_vector[r+1],l-v_vector[r+1]+beta,scale,encoder,plain0,plain1,flag_mul);
                // evaluator.multiply_plain(cipher_vector[alpha-1],plain0,cipher0);
                // evaluator.multiply_plain(cipher_vector[alpha],plain1,cipher1);
                // evaluator.add(cipher0,cipher1,cipher_tmp);
                // evaluator.relinearize_inplace(cipher_tmp,relin_keys);
                // evaluator.rescale_to_next_inplace(cipher_tmp);
                select_mul(flag_mul,plain0,cipher_vector[alpha-1],plain1,cipher_vector[alpha],cipher_tmp,evaluator,relin_keys);
                cipher_result.push_back(cipher_tmp);
                h=h+1;
            }
        }
    }
    destination=cipher_result;
}

inline vector<seal::Ciphertext> long_cipher::LongRot(int rotSize, int length, double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys)
{
    vector<seal::Ciphertext> destination;
    LongRot(rotSize, length, destination, scale, encoder, evaluator, galois_keys, relin_keys);
    return destination;
}

inline void long_cipher::AutoLongRot(int rotSize, int length, vector<seal::Ciphertext> &destination, double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys)
{
    int r = (n * m) % l;    // Size of the ciphertext residue 
    int w = LongCipher.size();  // Size of the ciphertext 

    rotSize=rotSize%(n*m);
    int v = rotSize % l;    // Ciphertext residue 
    int u = rotSize / l;    // Ciphertext start index 

    int stop_length = ceil(double(length)/double(l));   //Ciphertext number of destination
    size_t last_r=length%l;

    // cout<<"r:"<<r<<"  v:"<<v<<"  u:"<<u<<"  w:"<<w<<"  stop_length:"<<stop_length<<endl;

    seal::Ciphertext cipher1,cipher2,cipher3;
    vector<seal::Ciphertext> rotate_cipher;

    int v_tmp=v;
    int k=1;
    for(int i=0;i<stop_length+k;i++){
        // cout<<v_tmp<<endl;
        if((u+i)%w!=w-1){
            // cout<<"case1"<<endl;
            if(v_tmp!=0){
                evaluator.rotate_vector(LongCipher[(u+i)%w],v_tmp,galois_keys,cipher1);
                rotate_cipher.push_back(cipher1);
            }
            else{
                rotate_cipher.push_back(LongCipher[(u+i)%w]);
            }
            
        }
        else{
            if(v_tmp!=0){
                // cout<<"case2_1"<<endl;
                evaluator.rotate_vector(LongCipher[(u+i)%w],v_tmp,galois_keys,cipher1);
                rotate_cipher.push_back(cipher1);
            }
            else{
                // cout<<"case2_2"<<endl;
                rotate_cipher.push_back(LongCipher[(u+i)%w]);
            }
            if(v_tmp<=r){
                v_tmp=(l-r+v_tmp)%l;
                if(v_tmp==0){
                    k++;
                }
            }
            else{
                v_tmp=(v_tmp-r)%l;
                k++;
            }
        }
    }
    // cout<<"rotate_length:"<<rotate_cipher.size()<<endl;
    // destination=rotate_cipher;

    Plaintext plain_tmp1,plain_tmp2,plain_tmp3;
    v_tmp=v;
    int flag_mul=0;
    // gen_plain_vec(0,l-v_tmp,scale,encoder,plain_tmp1);
    // gen_plain_vec(l-v_tmp,l,scale,encoder,plain_tmp2);
    gen_select_plain(v_tmp,scale,encoder,plain_tmp1,plain_tmp2,flag_mul);
    int i=0;
    // For the case where the destination is not the last one,
    // each ciphertext in l slots contains all the values.
    while(int(destination.size())<int(stop_length)-1){
        // When (u + i) % w == w - 1 and v_tmp<=r, special handling is required
        // In this case, we need to update v_tmp and reconstruct the plaintext vector
        if ((u + i) % w == w - 1 && v_tmp <= r)
        {
            // cout << "Update v_tmp" << endl;
            v_tmp = (l - r + v_tmp) % l;
            gen_select_plain(v_tmp,scale,encoder,plain_tmp1,plain_tmp2,flag_mul);
            // gen_plain_vec(0, l - v_tmp, scale, encoder, plain_tmp1);
            // gen_plain_vec(l - v_tmp, l, scale, encoder, plain_tmp2);
            if(v_tmp==0){
                i++;
                continue;
            }
        }

        // When (u + i) % w == w - 2 and v_tmp > r, we need to merge three ciphertexts into one.
        // Therefore, if w - 1, which is the content of the last ciphertext, is not enough,
        // we need to add the content from the first ciphertext.
        // After that, we update v_tmp and reconstruct the plaintext vector.
        if((u+i)%w==w-2 && v_tmp>r){
            // std::cout<<"Merage three ciphertext"<<endl;
            gen_plain_vec(l - v_tmp + r, l, scale, encoder, plain_tmp3);
            evaluator.multiply_plain(rotate_cipher[i% rotate_cipher.size()], plain_tmp1, cipher1);
            evaluator.multiply_plain(rotate_cipher[(i + 1) % rotate_cipher.size()], plain_tmp2, cipher2);
            evaluator.multiply_plain(rotate_cipher[(i + 2) % rotate_cipher.size()], plain_tmp3, cipher3);
            evaluator.add_inplace(cipher1, cipher2);
            evaluator.add_inplace(cipher1, cipher3);
            evaluator.relinearize_inplace(cipher1, relin_keys);
            evaluator.rescale_to_next_inplace(cipher1);
            destination.push_back(cipher1);
            v_tmp = (v_tmp - r) % l;
            gen_select_plain(v_tmp,scale,encoder,plain_tmp1,plain_tmp2,flag_mul);
            // gen_plain_vec(0, l - v_tmp, scale, encoder, plain_tmp1);
            // gen_plain_vec(l - v_tmp, l, scale, encoder, plain_tmp2);
            i=i+2; 
            continue;
        }

        // Regular operation: Simply concatenate the first l - v_tmp elements of the current ciphertext
        // with the last v_tmp elements of the next ciphertext to form a new ciphertext.
        // std::cout<<"  v_tmp:"<<v_tmp<<"  i:"<<i% rotate_cipher.size()<<endl;
        // evaluator.multiply_plain(rotate_cipher[i% rotate_cipher.size()],plain_tmp1,cipher1);
        // evaluator.multiply_plain(rotate_cipher[(i+1)% rotate_cipher.size()],plain_tmp2,cipher2);
        // evaluator.add_inplace(cipher1,cipher2);
        // evaluator.relinearize_inplace(cipher1, relin_keys);
        // evaluator.rescale_to_next_inplace(cipher1);
        select_mul(flag_mul,plain_tmp1,rotate_cipher[i% rotate_cipher.size()],plain_tmp2,rotate_cipher[(i+1)% rotate_cipher.size()],cipher1,evaluator,relin_keys);
        destination.push_back(cipher1);
        i++;
    }

    // For the case where the destination is the last one,
    // ciphertext in l slots contains length % l values.
    if(int(destination.size())==stop_length-1){
        // cout<<"v_tmp:"<<v_tmp<<endl;
        // cout<<"last_r:"<<last_r<<"  v_tmp:"<<v_tmp<<"  i:"<<i<<endl;
        if ((u+i)%w==w-2 && v_tmp>r)
        {
            if(last_r<=l-v_tmp)
            {
                // cout<<"last_r case1_1"<<endl;
                gen_plain_vec(0, last_r, scale, encoder, plain_tmp1);
                evaluator.multiply_plain(rotate_cipher[i% rotate_cipher.size()],plain_tmp1,cipher1);
                evaluator.relinearize_inplace(cipher1, relin_keys);
                evaluator.rescale_to_next_inplace(cipher1);
                destination.push_back(cipher1);
            }
            else if (last_r<=l-v_tmp+r)
            {
                // cout<<"last_r case1_2"<<endl;
                gen_plain_vec(0, l - v_tmp, scale, encoder, plain_tmp1);
                gen_plain_vec(l - v_tmp, last_r, scale, encoder, plain_tmp2);
                evaluator.multiply_plain(rotate_cipher[i% rotate_cipher.size()],plain_tmp1,cipher1);
                evaluator.multiply_plain(rotate_cipher[(i+1)% rotate_cipher.size()],plain_tmp2,cipher2);
                evaluator.add_inplace(cipher1,cipher2);
                evaluator.relinearize_inplace(cipher1, relin_keys);
                evaluator.rescale_to_next_inplace(cipher1);
                destination.push_back(cipher1);
            }
            else
            {
                // cout<<"last_r case1_3"<<endl;
                gen_plain_vec(0, l - v_tmp, scale, encoder, plain_tmp1);
                gen_plain_vec(l - v_tmp, l - v_tmp + r, scale, encoder, plain_tmp2);
                gen_plain_vec(l - v_tmp + r, last_r, scale, encoder, plain_tmp3);
                evaluator.multiply_plain(rotate_cipher[i% rotate_cipher.size()], plain_tmp1, cipher1);
                evaluator.multiply_plain(rotate_cipher[(i + 1) % rotate_cipher.size()], plain_tmp2, cipher2);
                evaluator.multiply_plain(rotate_cipher[(i + 2) % rotate_cipher.size()], plain_tmp3, cipher3);
                evaluator.add_inplace(cipher1, cipher2);
                evaluator.add_inplace(cipher1, cipher3);
                evaluator.relinearize_inplace(cipher1, relin_keys);
                evaluator.rescale_to_next_inplace(cipher1);
                destination.push_back(cipher1);
            }
            
        }
        else
        {
            if ((u + i) % w == w - 1 && v_tmp <= r){
                v_tmp = (l - r + v_tmp) % l;
                // cout<<"last_r change v_tmp"<<endl;
                if(v_tmp==0){
                i++;
            }
            }
            if(last_r<=l - v_tmp){
                // cout<<"last_r case2_1"<<endl;
                gen_plain_vec(0, last_r, scale, encoder, plain_tmp1);
                evaluator.multiply_plain(rotate_cipher[i% rotate_cipher.size()],plain_tmp1,cipher1);
                evaluator.relinearize_inplace(cipher1, relin_keys);
                evaluator.rescale_to_next_inplace(cipher1);
                destination.push_back(cipher1);
            }
            else{
                // cout<<"last_r case2_2"<<endl;
                gen_plain_vec(0, l - v_tmp, scale, encoder, plain_tmp1);
                gen_plain_vec(l - v_tmp, last_r, scale, encoder, plain_tmp2);
                evaluator.multiply_plain(rotate_cipher[i% rotate_cipher.size()],plain_tmp1,cipher1);
                evaluator.multiply_plain(rotate_cipher[(i+1)% rotate_cipher.size()],plain_tmp2,cipher2);
                evaluator.add_inplace(cipher1,cipher2);
                evaluator.relinearize_inplace(cipher1, relin_keys);
                evaluator.rescale_to_next_inplace(cipher1);
                destination.push_back(cipher1);
            }
        }
        
    }

}

inline void long_cipher::LongAdd(vector<seal::Ciphertext> &cipher_vec, seal::Evaluator &evaluator)
{
    if (cipher_vec.size() != LongCipher.size())
    {
        cout << "Two ciphertext length are not equal." << endl;
    }

    // Check if modulus switching is needed
    bool parms_id_switch = false;
    if (LongCipher[0].parms_id() != cipher_vec[0].parms_id())
    {
        parms_id_switch = true;
    }

    // Check if scale is equal
    bool scale_change = false;
    if (LongCipher[0].scale() != cipher_vec[0].scale())
    {
        // cout<<"LongCipher[0].scale():"<<LongCipher[0].scale()<<endl;
        // cout<<"cipher_vec[0].scale():"<<cipher_vec[0].scale()<<endl;
        // cerr<<"scale is not equal."<<endl;
        scale_change = true;
    }

    for (size_t i = 0; i < LongCipher.size(); i++)
    {
        if (parms_id_switch)
        {
            evaluator.mod_switch_to_inplace(LongCipher[i], cipher_vec[i].parms_id());
        }
        if (scale_change)
        {
            cipher_vec[i].scale() = LongCipher[i].scale();
        }
        evaluator.add_inplace(LongCipher[i], cipher_vec[i]);
    }
}


inline void long_cipher::gen_select_plain(size_t v_tmp, double scale, seal::CKKSEncoder &encoder, Plaintext &destination1, Plaintext &destination2, int &flag_mul)
{
    if(v_tmp==0){
        flag_mul=1;
    }
    else if (v_tmp==l)
    {
        flag_mul=2;
    }
    else{
        flag_mul=0;
        gen_plain_vec(0,l-v_tmp,scale,encoder,destination1);
        gen_plain_vec(l-v_tmp,l,scale,encoder,destination2);
    }
    
}

inline void long_cipher::gen_select_plain(size_t start_index, size_t v_tmp, size_t end_index, double scale, seal::CKKSEncoder &encoder, Plaintext &destination1, Plaintext &destination2, int &flag_mul)
{
    if(v_tmp==start_index){
        flag_mul=1;
    }
    else if (v_tmp==end_index)
    {
        flag_mul=2;
    }
    else{
        flag_mul=0;
        gen_plain_vec(0,l-v_tmp,scale,encoder,destination1);
        gen_plain_vec(l-v_tmp,l,scale,encoder,destination2);
    }
}

inline void long_cipher::select_mul(int &flag_mul, seal::Plaintext &plain1, seal::Ciphertext &cipher1, seal::Plaintext &plain2, seal::Ciphertext cipher2, seal::Ciphertext &destination, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys)
{
    seal::Ciphertext cipher_tmp1,cipher_tmp2;
    if(flag_mul==0){
        evaluator.multiply_plain(cipher1,plain1,cipher_tmp1);
        evaluator.multiply_plain(cipher2,plain2,cipher_tmp2);
        evaluator.add(cipher_tmp1,cipher_tmp2,destination);
        evaluator.relinearize_inplace(destination, relin_keys);
        evaluator.rescale_to_next_inplace(destination);
    }else if (flag_mul==1)
    {
        evaluator.mod_switch_to_next(cipher1,destination);
    }else if (flag_mul==2)
    {
        evaluator.mod_switch_to_next(cipher2,destination);
    }
    
    
}

inline long_cipher &long_cipher::operator=(const long_cipher &other)
{
    if (this == &other)
    {
        return *this; // Avoid self-assignment
    }

    // Copy data members from the other object
    n = other.n;
    m = other.m;
    LongCipher = other.LongCipher;
    l = other.l;

    // Return the reference to the updated object
    return *this;
}

inline size_t long_cipher::get_cols() const
{
    return m;
}

inline size_t long_cipher::get_cipher_length() const
{
    return LongCipher.size();
}

inline vector<seal::Ciphertext> long_cipher::get_cipher() const
{
    return LongCipher;
}

inline vector<seal::Ciphertext> long_cipher::get_cipher()
{
    return LongCipher;
}

inline size_t long_cipher::get_slot() const
{
    return l;
}

inline size_t long_cipher::get_rows() const
{
    return n;
}

long_cipher LongMul(long_cipher &A, long_cipher &B,
                    seal::Evaluator &evaluator, seal::RelinKeys &relin_keys)
{
    vector<seal::Ciphertext> destination,cipher_A,cipher_B;
    seal::Ciphertext cipher_tmp;    
    for (size_t i = 0; i < std::min(A.get_cipher_length(), B.get_cipher_length()); i++)
    {
        evaluator.multiply(A.get_cipher()[i], B.get_cipher()[i], cipher_tmp);
        evaluator.relinearize_inplace(cipher_tmp, relin_keys);
        evaluator.rescale_to_next_inplace(cipher_tmp);
        destination.push_back(cipher_tmp);
    }
    long_cipher C(A.get_rows(), B.get_cols(), destination, A.get_slot());
    return C;
}

vector<seal::Ciphertext> LongMul(vector<seal::Ciphertext> &A, vector<seal::Ciphertext> &B,
                                 seal::Evaluator &evaluator, seal::RelinKeys &relin_keys)
{
    vector<seal::Ciphertext> destination;
    seal::Ciphertext cipher_tmp;
    for (size_t i = 0; i < std::min(A.size(), B.size()); i++)
    {
        evaluator.multiply(A[i], B[i], cipher_tmp);
        evaluator.relinearize_inplace(cipher_tmp, relin_keys);
        evaluator.rescale_to_next_inplace(cipher_tmp);
        destination.push_back(cipher_tmp);
    }
    return destination;
}

// ciphertext multiply
void encrypted_long_matrix_multiplication(seal::CKKSEncoder &encoder, // seal::Decryptor &decryptor,
                                          double scale, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys,
                                          seal::GaloisKeys &galois_keys, long_cipher &encrypted_A,
                                          long_cipher &encrypted_B, long_cipher &destination)
{
    size_t r;
    size_t n, m, p;
    n = encrypted_A.get_rows();
    m = encrypted_A.get_cols();
    p = encrypted_B.get_cols();
    if (m != encrypted_B.get_rows())
    {
        cerr << "Matrix dimensions are not equal to multiply.";
    }
    r = smallest_r(n, m, p);

    // multiply
    // long_cipher encrypted_C =LongMul(encrypted_A,encrypted_B,evaluator,relin_keys);

    for (std::size_t i = 0; i < m; i++)
    {
        vector<seal::Ciphertext> ct_a, ct_b, ct_c;
        // cout<<"A("<<mod(-i * n, n * m)<<"): ";
        encrypted_A.AutoLongRot(mod(-i * n, n * m), n * p, ct_a, scale, encoder, evaluator, galois_keys, relin_keys);
        // cout<<"B("<<mod(i * (r * m - n), (m * p))<<"):";
        encrypted_B.AutoLongRot(mod(i * (r * m - n), (m * p)), n * p, ct_b, scale, encoder, evaluator, galois_keys, relin_keys);
        ct_c = LongMul(ct_a, ct_b, evaluator, relin_keys);
        destination.LongAdd(ct_c, evaluator);
    }
    // destination=encrypted_C;
}

// inline void long_cipher::LongRot(int rotSize, int length, vector<seal::Ciphertext> &destination,
//                                  double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator,
//                                  seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys)
// {
//     // In fact, Each time only the elements of length need to be considered for rotation
//     // int rot_num=ceil(double(length)/double(length));

//     int r = (n * m) % l;
//     int v = rotSize % l;
//     int u = rotSize / l;
//     int w = LongCipher.size();

//     // cout<<"r:"<<r<<"  v:"<<v<<"  u:"<<u<<"  w:"<<w<<endl;

//     //+ : left | {1,2,3,4,5,6,0...0}  (2)---> {3,4,5,6,0...0,1,2}
//     //- : right| {1,2,3,4,5,6,0...0}  (-2)---> {0,0,1,2,3,4,5,6,0...0}

//     // generate rotate vector
//     vector<Ciphertext> rotate_part;
//     seal::Ciphertext cipher_tmp;

//     // a)
//     for (int i = 0; i < w - u; i++)
//     {
//         evaluator.rotate_vector(LongCipher[(u + i) % w], v, galois_keys, cipher_tmp);
//         rotate_part.push_back(cipher_tmp);
//     }

//     // b)
//     for (int i = w - u; i < w + 1; i++)
//     {
//         evaluator.rotate_vector(LongCipher[(i - w + u) % w], v - r, galois_keys, cipher_tmp);
//         rotate_part.push_back(cipher_tmp);
//     }

//     if (r < v)
//     {
//         Plaintext select_plain_1, select_plain_2, select_plain_3;
//         Ciphertext cipher1, cipher2, cipher3;
//         // c)
//         gen_plain_vec(0, l - v, scale, encoder, select_plain_1);
//         gen_plain_vec(l - v, l, scale, encoder, select_plain_2);
//         for (int i = 0; i < w - u - 2; i++)
//         {
//             // cout<<"d:"<<i<<endl;
//             evaluator.multiply_plain(rotate_part[i], select_plain_1, cipher1);
//             evaluator.multiply_plain(rotate_part[i + 1], select_plain_2, cipher2);
//             evaluator.add_inplace(cipher1, cipher2);

//             evaluator.relinearize_inplace(cipher1, relin_keys);
//             evaluator.rescale_to_next_inplace(cipher1);
//             destination.push_back(cipher1);
//         }

//         // e
//         gen_plain_vec(0, l - v, scale, encoder, select_plain_1);
//         gen_plain_vec(l - v, l - v + r, scale, encoder, select_plain_2);
//         gen_plain_vec(l - v + r, l, scale, encoder, select_plain_3);
//         evaluator.multiply_plain(rotate_part[w - u - 2], select_plain_1, cipher1);
//         evaluator.multiply_plain(rotate_part[w - u - 1], select_plain_2, cipher2);
//         evaluator.multiply_plain(rotate_part[w - u], select_plain_3, cipher3);
//         evaluator.add_inplace(cipher1, cipher2);
//         evaluator.add_inplace(cipher1, cipher3);
//         evaluator.relinearize_inplace(cipher1, relin_keys);
//         evaluator.rescale_to_next_inplace(cipher1);
//         destination.push_back(cipher1);

//         // f
//         gen_plain_vec(0, l - v + r, scale, encoder, select_plain_1);
//         gen_plain_vec(l - v + r, l, scale, encoder, select_plain_2);
//         for (int i = w - u; i < w; i++)
//         {
//             // cout<<"f:"<<i<<endl;
//             evaluator.multiply_plain(rotate_part[i], select_plain_1, cipher1);
//             evaluator.multiply_plain(rotate_part[i + 1], select_plain_2, cipher2);
//             evaluator.add_inplace(cipher1, cipher2);

//             evaluator.relinearize_inplace(cipher1, relin_keys);
//             evaluator.rescale_to_next_inplace(cipher1);
//             destination.push_back(cipher1);
//         }

//         // g
//         gen_plain_vec(0, r, scale, encoder, select_plain_1);
//         evaluator.multiply_plain(rotate_part[w], select_plain_1, cipher_tmp);
//         evaluator.relinearize_inplace(cipher_tmp, relin_keys);
//         evaluator.rescale_to_next_inplace(cipher_tmp);
//         destination.push_back(cipher_tmp);
//     }
//     else if (r >= v)
//     {
//         Plaintext select_plain_1, select_plain_2;
//         Ciphertext cipher1, cipher2;
//         // c)
//         gen_plain_vec(0, l - v, scale, encoder, select_plain_1);
//         gen_plain_vec(l - v, l, scale, encoder, select_plain_2);
//         for (int i = 0; i < w - u - 1; i++)
//         {
//             evaluator.multiply_plain(rotate_part[i], select_plain_1, cipher1);
//             evaluator.multiply_plain(rotate_part[i + 1], select_plain_2, cipher2);
//             evaluator.add_inplace(cipher1, cipher2);

//             evaluator.relinearize_inplace(cipher1, relin_keys);
//             evaluator.rescale_to_next_inplace(cipher1);
//             destination.push_back(cipher1);
//         }

//         // d)
//         gen_plain_vec(0, r - v, scale, encoder, select_plain_1);
//         gen_plain_vec(r - v, l, scale, encoder, select_plain_2);
//         for (int i = w - u - 1; i < w - 1; i++)
//         {
//             evaluator.multiply_plain(rotate_part[i], select_plain_1, cipher1);
//             evaluator.multiply_plain(rotate_part[i + 1], select_plain_2, cipher2);
//             evaluator.add_inplace(cipher1, cipher2);

//             evaluator.relinearize_inplace(cipher1, relin_keys);
//             evaluator.rescale_to_next_inplace(cipher1);
//             destination.push_back(cipher1);
//         }

//         // e)
//         gen_plain_vec(0, r - v, scale, encoder, select_plain_1);
//         gen_plain_vec(r - v, r, scale, encoder, select_plain_2);
//         evaluator.multiply_plain(rotate_part[w - 1], select_plain_1, cipher1);
//         evaluator.multiply_plain(rotate_part[w], select_plain_2, cipher2);
//         evaluator.add_inplace(cipher1, cipher2);
//         evaluator.relinearize_inplace(cipher1, relin_keys);
//         evaluator.rescale_to_next_inplace(cipher1);
//         destination.push_back(cipher1);
//     }
// }


