#include "utils.h"

using namespace seal::util;

int64_t mod(int64_t x, int64_t q){
    if (q < 0 || q == 0){
        std::cerr << "mod: modulus must be an positive integer " << std::endl;
        return 0;
    }
    else{
        if (x > 0 || x == 0){
            return  x % q;
        }
        else{
            if (abs(x) % q == 0){
                return 0;
            }
            else{
                return x + ((-x)/q+1)*q;
            }
        }
    }
}


int64_t mods(int64_t x, int64_t q){
    int r = mod(x, q);
    if (r > q/2)
        r -= q;
    return r;
}






uint64_t smallest_r(uint64_t n, uint64_t m, uint64_t p){
    uint64_t r = 1;
    int64_t s = r * m - n;
    while (( s % p != 0) || (s < 0)){
        r += 1;
        s += m;
    }
    return r;
}




// decrypt and then decode the plaintext into a vector
void decrypt_decode(seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, 
                    seal::Ciphertext &x_encrypted, std::vector<double> &x)
{
    seal::Plaintext x_decrypted;
    decryptor.decrypt(x_encrypted, x_decrypted);
    encoder.decode(x_decrypted, x);
}





// void decrypt_decode(seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, 
//                     std::vector<seal::Ciphertext> &x_encrypted, std::vector<std::vector<double>> &x)
// {
//     std::size_t n = x_encrypted.size();
//     x.resize(n);
//     for (std::size_t i = 0; i < n; i++){
//         seal::Plaintext x_decrypted;
//         decryptor.decrypt(x_encrypted[i], x_decrypted);
//         encoder.decode(x_decrypted, x[i]);
//     }
    
// }

void decrypt_decode(seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, std::vector<std::vector<seal::Ciphertext>> &x_encrypted, std::vector<std::vector<std::vector<double>>> &x)
{
    std::size_t n = x_encrypted.size();
    std::size_t m = x_encrypted[0].size();
    // std::cout << "decrypt_decode: (n, p): (" << n << ", " << m << ")" << std::endl;
    x.resize(n);
    for (std::size_t i = 0; i < n; i++){
        x[i].resize(m);
        for (std::size_t j = 0; j < m; j++){
            seal::Plaintext x_decrypted;
            decryptor.decrypt(x_encrypted[i][j], x_decrypted);
            encoder.decode(x_decrypted, x[i][j]);
        }
    }   
}






// return 1 if n is a power of 2, 0 otherwise
bool is_power_of_two(std::size_t n)
{
    if (n - pow(2, round(log2(double(n)))) == 0)
        return true;
    else
        return false;
}


// return 1 if n is a square number, 0 otherwise
bool is_square(std::size_t n)
{
    if (n - std::pow(std::round(std::sqrt(double(n))), 2) == 0)
        return true;
    else
        return false;
}





// return a factor of n closest to sqrt(n)
uint64_t factor_around_sqrt(uint64_t n)
{
    uint64_t k =  round(sqrt(double(n)));
    uint64_t r = 0;
    if (is_square(n)){
        r = k;
    }
    else{
        // for (uint64_t i = k; i > 0; i--){
        for (uint64_t i = k; ; i--){
            if (mod(n, i) == 0) {
                r = i;
                break;
            }
        }
    }
    return r;
}












// return the set of indices of the nonzero diagonal vectors 
// of the transformation matrix for submatrix extracting
void nonzero_index_set(std::size_t n, std::size_t m, std::size_t k, std::size_t l, std::vector<std::size_t> &ind)
{
    auto gcd = seal::util::xgcd(n, m);
    uint64_t g = get<0>(gcd);
    if ((g != 1) || (std::gcd(k, l) != 1)){
        std::cerr << "submatrix_transformation: the dimension (" << n << ", " << m 
             << ") or (" << k << ", " << l << ") are not coprime " << std::endl;
        return;
    }
    int64_t s = get<1>(gcd);
    int64_t t = get<2>(gcd);
    ind.clear();
    ind.push_back(0);
    for(std::size_t i = 1; i < k*l; i++){
        std::size_t r = mod( (mod(i, k) - i)*t*m + (mod(i, l) - i)*s*n, n*m);
        bool repeated = false;
        for (std::size_t j = 0; j < ind.size(); j++){
            if ( r == ind[j]){
                repeated = true;
                break;
            }
        }
        if (!repeated){
            ind.push_back(r);
        }
    }
}


// return the nearest square number at least n
uint64_t nearest_square_number(uint64_t n)
{
    uint64_t r = std::pow(std::ceil(std::sqrt(double(n))), 2);
    return r;
}




// split a vector into 'num' sections.
void split(std::vector<double> &x, std::size_t num, std::vector<std::vector<double>> &x_data)
{
    std::size_t d = std::ceil((double)x.size() / (double)num);
    std::size_t n = x.size();
    if (n / d + 1 < num)
    {
        x_data.resize(n / d + 1);
    }
    else
    {
        x_data.resize(num);
    }
    if (d == 1)
    {
        x_data.resize(n);
    }
    std::size_t dim = x_data.size();
    // cout << "   d:" << d << endl;
    // cout << "   n/d+1:" << n/d + 1 << endl;
    // cout << "   size of data is:" << dim << endl;
    for (std::size_t i = 0; i < dim - 1; i++)
    {
        x_data[i].clear();
        for (std::size_t j = 0; j < d; j++)
        {
            x_data[i].push_back(x[i * d + j]);
        }
    }
    x_data[dim - 1].clear();
    for (std::size_t i = d * (dim - 1); i < n; i++)
    {
        x_data[dim - 1].push_back(x[i]);
    }
}//end split


// encode and encrypt a vector<double>, return a vector of Ciphertext
void encode_encrypt(seal::CKKSEncoder &encoder, seal::Encryptor &encryptor, std::vector<double> &x,
                    double scale, std::vector<seal::Ciphertext> &x_encrypted)
{
    std::size_t slot_count = encoder.slot_count();
    std::size_t dim = x.size();
    if (dim < slot_count)
    {
        seal::Plaintext x_plain;
        encoder.encode(x, scale, x_plain);
        x_encrypted.resize(1);
        encryptor.encrypt(x_plain, x_encrypted[0]);
    }
    else
    { // split x into num_sect sections, each section forms a new vector
        std::size_t num_sect = std::ceil((double)dim / (double)slot_count);
        std::vector<std::vector<double>> x_data;
        split(x, num_sect, x_data);
        x_encrypted.resize(num_sect);
        for (std::size_t i = 0; i < num_sect; i++)
        {
            seal::Plaintext x_plain;
            encoder.encode(x_data[i], scale, x_plain);
            // cout << "encode_encrypt: x_data[" << i << "]:" << endl;
            // print_vector(x_data[i]);
            encryptor.encrypt(x_plain, x_encrypted[i]);
        }
    }
}


void encode_encrypt(seal::CKKSEncoder &encoder, seal::Encryptor &encryptor, std::vector<double> &x,
                    double scale, seal::Ciphertext &x_encrypted)
{
    std::size_t slot_count = encoder.slot_count();


    std::size_t dim = x.size();
    
    std::vector<double> y;

    if (dim > slot_count){    
        for (std::size_t i = 0; i < slot_count; i++){
            y.push_back(x[i]);
        }    
    }
    else{
        for (std::size_t i = 0; i < dim; i++){
            y.push_back(x[i]);
        }
    }
    seal::Plaintext x_plain;
    encoder.encode(y, scale, x_plain);
    encryptor.encrypt(x_plain, x_encrypted);
}


// encode and encrypt a double number
void encode_encrypt(seal::CKKSEncoder &encoder, seal::Encryptor &encryptor, double x, double scale, seal::Ciphertext &x_encrypted)
{
    std::vector<double> v;
    for (std::size_t i = 0; i < encoder.slot_count(); i++)
    {
        v.push_back(x);
    }
    seal::Plaintext x_plain;
    encoder.encode(v, scale, x_plain);
    encryptor.encrypt(x_plain, x_encrypted);
}








