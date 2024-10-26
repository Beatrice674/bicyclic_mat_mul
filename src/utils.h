#ifndef UTILS_H
#define UTILS_H

#include "seal/seal.h"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// using namespace seal;
// using namespace std;

// return the nonzero mod
int64_t mod(int64_t x, int64_t q);


// return the symmetric mod
int64_t mods(int64_t x, int64_t q);




/*
Helper function: Prints the parameters in a SEALContext.
*/
inline void print_parameters(const seal::SEALContext &context)
{
    auto &context_data = *context.key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::bfv:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::ckks:
        scheme_name = "CKKS";
        break;
    case seal::scheme_type::bgv:
        scheme_name = "BGV";
        break;
    default:
        throw std::invalid_argument("unsupported scheme");
    }
    std::cout << "/" << std::endl;
    std::cout << "| Encryption parameters :" << std::endl;
    std::cout << "|   scheme: " << scheme_name << std::endl;
    std::cout << "|   poly_modulus_degree: " << context_data.parms().poly_modulus_degree() << std::endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_modulus_size = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_modulus_size - 1; i++)
    {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }
    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;

         /*
    Print the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus: ";
    // std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    // auto coeff_modulus = context_data.parms().coeff_modulus();



    for (std::size_t i = 0; i < coeff_modulus_size-1 ; i++){
        if (context_data.total_coeff_modulus()[i] != 0){
            std::cout << context_data.total_coeff_modulus()[i] << " * ";
        }
        else{
            std::cout <<  "1 * ";
        }
    }
    if (context_data.total_coeff_modulus()[coeff_modulus_size] != 0){
        std::cout << context_data.total_coeff_modulus()[coeff_modulus_size - 1];
    }
    else{
        std::cout <<  "1 ";
    }
    std::cout << std::endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::bfv)
    {
        std::cout << "|   plain_modulus: " << context_data.parms().plain_modulus().value() << std::endl;
    }

    std::cout << "\\" << std::endl;
}

/*
Helper function: Prints the `parms_id' to std::ostream.
*/
inline std::ostream &operator<<(std::ostream &stream, seal::parms_id_type parms_id)
{
    /*
    Save the formatting information for std::cout.
    */
    std::ios old_fmt(nullptr);
    old_fmt.copyfmt(std::cout);

    stream << std::hex << std::setfill('0') << std::setw(16) << parms_id[0] << " " << std::setw(16) << parms_id[1]
           << " " << std::setw(16) << parms_id[2] << " " << std::setw(16) << parms_id[3] << " ";

    /*
    Restore the old std::cout formatting.
    */
    std::cout.copyfmt(old_fmt);

    return stream;
}

/*
Helper function: Prints a vector of floating-point values.
*/
template <typename T>
inline void print_vector(std::vector<T> vec, std::size_t print_size = 4, int prec = 3)
{
    /*
    Save the formatting information for std::cout.
    */
    std::ios old_fmt(nullptr);
    old_fmt.copyfmt(std::cout);

    std::size_t slot_count = vec.size();

    std::cout << std::fixed << std::setprecision(prec);
    std::cout << std::endl;
    if (slot_count <= 2 * print_size)
    {
        std::cout << "    [";
        for (std::size_t i = 0; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    else
    {
        vec.resize(std::max(vec.size(), 2 * print_size));
        std::cout << "    [";
        for (std::size_t i = 0; i < print_size; i++)
        {
            std::cout << " " << vec[i] << ",";
        }
        if (vec.size() > 2 * print_size)
        {
            std::cout << " ...,";
        }
        for (std::size_t i = slot_count - print_size; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    std::cout << std::endl;

    /*
    Restore the old std::cout formatting.
    */
    std::cout.copyfmt(old_fmt);
}

/*
Helper function: Print line number.
*/
inline void print_line(int line_number)
{
    std::cout << "Line " << std::setw(3) << line_number << " --> ";
}

/*
Helper function: Convert a value into a hexadecimal string, e.g., uint64_t(17) --> "11".
*/
inline std::string uint64_to_hex_string(std::uint64_t value)
{
    return seal::util::uint_to_hex_string(&value, std::size_t(1));
}


uint64_t smallest_r(uint64_t n, uint64_t m, uint64_t p);



// decrypt and then decode the plaintext into a vector
void decrypt_decode(seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, seal::Ciphertext &x_encrypted, std::vector<double> &x);

// void decrypt_decode(seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, std::vector<seal::Ciphertext> &x_encrypted, std::vector<std::vector<double>> &x);

void decrypt_decode(seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, std::vector<std::vector<seal::Ciphertext>> &x_encrypted, std::vector<std::vector<std::vector<double>>> &x);

// rotate the vector v by k step and extend dim of v from k to d
template <typename T>
void rotate_vector(std::vector<T> & v, int64_t k, size_t d, std::vector<T> &t)
{
    
    std::vector<T> x;
    x.resize(d, 0);
    t.resize(d);
    std::size_t n = v.size();
    
    // int64_t r = mod(k, n);
    for (std::size_t i = 0; i < n; i++){
        x[i] = v[i];
    }
    
    for (std::size_t i = 0; i < d; i++){
        t[i] = x[mod(i + k, n)];
    }
}




// rotate the vector by k step
template <typename T>
void rotate_vector(std::vector<T> & v, int64_t k, std::vector<T> &t)
{
    std::size_t n = v.size();
    t.resize(n);
    // int64_t r = mod(k, n);
    for (std::size_t i = 0; i < n; i++){
        t[i] = v[mod(i + k, n)];
    }
}

// return 1 if n is a power of 2, 0 otherwise
bool is_power_of_two(std::size_t n);

// return 1 if n is a square number, 0 otherwise
bool is_square(std::size_t n);


// compute the Hadamard product of two vectors, 
// the result is the last parameter
template <typename T>
void hadamard_product(std::vector<T> &a, std::vector<T> &b, std::vector<T> &c)
{
    std::size_t n = a.size();
    if (n != b.size()){
        std::cerr << "hadamard_product: the size of the two input vectors ("<< n << ", " << b.size() 
                  <<") does not match " << std::endl;
        return;
    }
    c.clear();
    for(std::size_t i = 0; i < n; i++){
        c.push_back(a[i] * b[i]);
    }
}

// addtion of two vectors, 
// the result is the last parameter
template <typename T>
void add_vector(std::vector<T> &a, std::vector<T> &b, std::vector<T> &c)
{
    std::size_t n = a.size();
    if (n != b.size()){
        std::cerr << "add_vector: the size of the two input vectors ("<< n << ", " << b.size() 
                <<") does not match " << std::endl;
        return;
    }
    c.clear();
    for(std::size_t i = 0; i < n; i++){
        c.push_back(a[i] + b[i]);
    }
}

// return a factor of n closest to sqrt(n)
uint64_t factor_around_sqrt(uint64_t n);













// return the set of indices of the nonzero diagonal vectors 
// of the transformation matrix for submatrix extracting
void nonzero_index_set(std::size_t n, std::size_t m, std::size_t k, std::size_t l, std::vector<size_t> &ind);


// return the nearest square number at least n
uint64_t nearest_square_number(uint64_t n);



                                            

// split a vector into 'num' sections.
void split(std::vector<double> &x, std::size_t num, std::vector<std::vector<double>> &x_data);



// encode and encrypt a vector<double>, return a vector of Ciphertext
void encode_encrypt(seal::CKKSEncoder &encoder, seal::Encryptor &encryptor, std::vector<double> &x,
                    double scale, std::vector<seal::Ciphertext> &x_encrypted);
                    

// encode and encrypt a vector<double>, return a Ciphertext
void encode_encrypt(seal::CKKSEncoder &encoder, seal::Encryptor &encryptor, std::vector<double> &x,
                    double scale, seal::Ciphertext &x_encrypted);

// encode and encrypt a double number
void encode_encrypt(seal::CKKSEncoder &encoder, seal::Encryptor &encryptor, double x, double scale, 
                    seal::Ciphertext &x_encrypted);



// encode and encrypt a Matrix<double> using diagonal encoding
#endif // utils.h
                                 






