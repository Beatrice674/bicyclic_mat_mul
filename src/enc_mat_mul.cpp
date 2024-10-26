#include "enc_mat_mul.h"






//! note that this is not the 2*sqrt() algorithm, 
//! since it  may not hold that 
//! ind[k*i + j] = ind[k*i] + ind[j]
void submatrix_extracting(seal::Decryptor &decryptor, 
                          seal::CKKSEncoder &encoder, seal::Encryptor &encryptor, seal::Evaluator &evaluator, 
                          seal::RelinKeys &relin_keys, seal::GaloisKeys &galois_keys, 
                          seal::Ciphertext &ct_mat, std::size_t n, std::size_t m, std::size_t p, std::size_t q, 
                          seal::Ciphertext &ct_r)
{
    std::size_t slot_count = encoder.slot_count();
    double scale = ct_mat.scale();
    
    Matrix<double> T;
    submatrix_transformation(n, m, p, q, T);
    
    std::vector<std::size_t> index;
    nonzero_index_set(n, m, p, q, index);
    

    
    std::vector<std::vector<double>> rotated_diagonal_vectors;
    rotated_diagonal_vectors_for_submatrix(n, m, p, q, slot_count, T, index, rotated_diagonal_vectors);

    std::size_t n_r = index.size();



    

    std::cout << "the number of nonzero diagonal vectors: " << n_r << std::endl;
    std::cout << "the indices of nonzero diagonal vectors: " << std::endl;
    print_vector(index, index.size());


    std::size_t r = ceil(sqrt(double(n_r)));
    std::size_t l, k;
    k = r;
    if (is_square(n_r)){
        l = r;
    }
    else{
        l = n_r/r + 1;
    }

    std::vector<int64_t> ind(l*k);
    for(size_t i = 0; i < l * k; i++){
        if (i < n_r)
            ind[i] = index[i];
        else
            ind[i] = 0;
    }
    // ind[2] = 13;

    std::cout << "l = " << l << ", k = " << k << std::endl;
    


    std::vector<double> giant_sum(slot_count, 0);
    std::vector<seal::Ciphertext> giant_sum_encrypted;
    encode_encrypt(encoder, encryptor, giant_sum, scale, giant_sum_encrypted);

    ct_r = giant_sum_encrypted[0];

    for(size_t i = 0; i < l; i++){
        std::cout << "i: " << i << std::endl;
        std::vector<double> baby_sum(slot_count, 0);
        std::vector<seal::Ciphertext> baby_sum_encrypted;
        encode_encrypt(encoder, encryptor, baby_sum, scale, baby_sum_encrypted);

        for(std::size_t j = 0; j < k; j++){
            std::cout << "j: " << j << std::endl;
            std::vector<seal::Ciphertext> diagonal_vector_encrypted;
            std::cout << "rotated_diagonal_vectors: " << std::endl;
            print_vector(rotated_diagonal_vectors[k*i + j], 15);
            seal::Plaintext rot_diag_vec;
            encoder.encode(rotated_diagonal_vectors[k*i + j], scale,  rot_diag_vec);

            // encode_encrypt(encoder, encryptor, rotated_diagonal_vectors[k*i + j], scale, diagonal_vector_encrypted);

            seal::Ciphertext rotated_ct_mat;
            // cout << "ind[k*i+j] - ind[k*i] = " << ind[k*i+j] << " - " <<  ind[k*i] <<" = " << ind[k*i+j] - ind[k*i] << endl;
            // evaluator.rotate_vector(ct_mat, mod(ind[k*i+j] - ind[k*i], n*m), galois_keys, rotated_ct_mat);
            evaluator.rotate_vector(ct_mat, ind[k*i+j] - ind[k*i], galois_keys, rotated_ct_mat);
            
            std::cout << "rotated_diagonal_encoding of the matrix: " << std::endl;
            std::vector<double> mat;
            decrypt_decode(encoder, decryptor, rotated_ct_mat, mat);
            print_vector(mat, 15);

            std::cout << "    + Scale of rot_diag_vec before multiply: " << rot_diag_vec.scale()  << std::endl;
            std::cout << "    + Scale of rotated_ct_mat before multiply: " << rotated_ct_mat.scale() << std::endl;

            seal::parms_id_type last_parms_id_b = rotated_ct_mat.parms_id();
            evaluator.mod_switch_to_inplace(rot_diag_vec, last_parms_id_b);

            // evaluator.mod_switch_to_inplace(diagonal_vector_encrypted[0], last_parms_id_b);
            // rotated_ct_mat.scale() = scale;
            // diagonal_vector_encrypted[0].scale() = scale;


            evaluator.multiply_plain_inplace(rotated_ct_mat, rot_diag_vec);

            std::cout << "    + Scale of rotated_ct_mat after multiply / before add: " 
                      << log2(rotated_ct_mat.scale()) << " bits" << std::endl;

            evaluator.rescale_to_next_inplace(rotated_ct_mat);

            
            last_parms_id_b = rotated_ct_mat.parms_id();
            evaluator.mod_switch_to_inplace(baby_sum_encrypted[0], last_parms_id_b);
            baby_sum_encrypted[0].scale() = scale;
            rotated_ct_mat.scale() = scale;
            evaluator.add_inplace(baby_sum_encrypted[0], rotated_ct_mat);

            // cout << "    + Scale of baby_sum_encrypted after add baby steps: " << log2(baby_sum_encrypted[0].scale()) << " bits" << endl;
        }
        
        
        seal::Ciphertext rotated_baby_sum_encypted;
        std::cout << "before rotate ind[k*i] " << ind[k*i] <<": ";
        std::vector<double> y;
        decrypt_decode(encoder, decryptor, baby_sum_encrypted[0], y);
        print_vector(y, 15);
        evaluator.rotate_vector(baby_sum_encrypted[0], ind[k*i], galois_keys, rotated_baby_sum_encypted);
        evaluator.relinearize_inplace(rotated_baby_sum_encypted, relin_keys);
        std::cout << "after rotate ind[k*i] " << ind[k*i] <<": ";
        decrypt_decode(encoder, decryptor, rotated_baby_sum_encypted, y);
        print_vector(y, 15);

        seal::parms_id_type last_parms_id = rotated_baby_sum_encypted.parms_id();
        evaluator.mod_switch_to_inplace(ct_r, last_parms_id);
        ct_r.scale() = scale;
        rotated_baby_sum_encypted.scale() = scale;

        // cout << "    + Scale of baby_sum_encrypted before add giant steps: " << log2(baby_sum_encrypted[0].scale()) << " bits" << endl;
        evaluator.add_inplace(ct_r, rotated_baby_sum_encypted);
        // cout << "    + Scale of baby_sum_encrypted before add baby steps: " << log2(baby_sum_encrypted[0].scale()) << " bits" << endl;
    }
}




// The following only works for m = max(n, m, p) for the moment
void encrypted_matrix_multiplication(//seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, 
                                     seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, 
                                     seal::GaloisKeys &galois_keys, uint64_t n, uint64_t m, uint64_t p,
                                     seal::Ciphertext &encrypted_A, seal::Ciphertext &encrypted_B, 
                                     seal::Ciphertext &encrypted_C)
{
    uint64_t r;
    r = smallest_r(n, m, p);                                    
    evaluator.multiply(encrypted_A, encrypted_B, encrypted_C);
    evaluator.relinearize_inplace(encrypted_C, relin_keys);
    evaluator.rescale_to_next_inplace(encrypted_C);

    for(std::size_t i = 1; i < m; i++){
        seal::Ciphertext ct_a, ct_b, ct_tmp;
        std::vector<double> a, b;

        evaluator.rotate_vector(encrypted_A, mod(-i * n, n*m), galois_keys, ct_a);
        // decrypt_decode(encoder, decryptor, ct_a, a);
        // cout << "rot(a, -" << i*n << "): " << endl;
        // print_vector(a);
        evaluator.rotate_vector(encrypted_B, mod(i * (r*m - n),  (m*p)), galois_keys, ct_b);
        // decrypt_decode(encoder, decryptor, ct_b, b);
        // cout << "rot(b, " << i*(r*m - n) << "): " << endl;
        // print_vector(b);
        evaluator.multiply(ct_a, ct_b, ct_tmp);
        evaluator.relinearize_inplace(ct_tmp, relin_keys);
        evaluator.rescale_to_next_inplace(ct_tmp);       
        evaluator.add_inplace(encrypted_C, ct_tmp);
    }
    // encrypted_C.scale() = pow(2.0, 40);
}

void rotate_cipher_vector(seal::CKKSEncoder &encoder, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys,seal::RelinKeys &relin_keys,
                          int64_t move_step, int64_t slot_all,
                          seal::Ciphertext &cipher, seal::Ciphertext &destination)
{
    seal::Ciphertext left_cipher,right_cipher;
    seal::Plaintext plain_left,plain_right;
    std::vector<double> left(slot_all-move_step,1),right(move_step,1);
    double scale=cipher.scale();
    // double scale=pow(2,20);
    // std::cout<<scale<<"   "<<cipher.scale()<<std::endl;

    encoder.encode(left,scale,plain_left);
    encoder.encode(right,scale,plain_right);
    evaluator.rotate_vector(cipher, move_step, galois_keys, left_cipher);
    evaluator.multiply_plain_inplace(left_cipher,plain_left);
    evaluator.multiply_plain(cipher,plain_right,right_cipher);
    evaluator.rotate_vector_inplace(right_cipher ,move_step-slot_all, galois_keys);
    evaluator.add(left_cipher,right_cipher,destination);
    evaluator.relinearize_inplace(destination,relin_keys);
    evaluator.rescale_to_next_inplace(destination);
    // evaluator.mod_switch_to_next_inplace(destination);
}

void encrypted_matrix_multiplication_v2(seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, 
                                     seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, 
                                     seal::GaloisKeys &galois_keys, uint64_t n, uint64_t m, uint64_t p,
                                     seal::Ciphertext &encrypted_A, seal::Ciphertext &encrypted_B, 
                                     seal::Ciphertext &encrypted_C)
{
    uint64_t r;
    r = smallest_r(n, m, p);                                    
    evaluator.multiply(encrypted_A, encrypted_B, encrypted_C);
    evaluator.relinearize_inplace(encrypted_C, relin_keys);
    evaluator.rescale_to_next_inplace(encrypted_C);
    evaluator.mod_switch_to_next_inplace(encrypted_C);
    // double scale=encrypted_C.scale();
    std::vector<double> c;
    // std::cout<<"encrypted_C scale:"<<log2(encrypted_C.scale())<<std::endl;

    for(std::size_t i = 1; i < m; i++){
        seal::Ciphertext ct_a, ct_b, ct_tmp;
        std::vector<double> a, b;
        
        
        rotate_cipher_vector(encoder,evaluator,galois_keys,relin_keys,mod(-i * n, n*m),n*m,encrypted_A,ct_a);
        // decrypt_decode(encoder, decryptor, ct_a, a);
        // std::cout << "rot(a, -" << i<< "): " << std::endl;
        // print_vector(a);
        // std::cout<<"ct_a scale:"<<log2(ct_a.scale())<<std::endl;
        
        rotate_cipher_vector(encoder,evaluator,galois_keys,relin_keys,mod(i * (r*m - n),  (m*p)),m*p,encrypted_B,ct_b);
        // decrypt_decode(encoder, decryptor, ct_b, b);
        // std::cout << "rot(b, " << i*(r*m - n) << "): " << std::endl;
        // print_vector(b);
        // std::cout<<"ct_b scale:"<<log2(ct_b.scale())<<std::endl;
        
        // std::cout<<"---------1----------"<<std::endl;
        evaluator.multiply(ct_a, ct_b, ct_tmp);
        // std::cout<<"---------2----------"<<std::endl;
        // std::cout<<"tmp scale:"<<log2(ct_tmp.scale())<<std::endl;
        evaluator.relinearize_inplace(ct_tmp, relin_keys);
        evaluator.rescale_to_next_inplace(ct_tmp);
        // std::cout<<"tmp scale:"<<log2(ct_tmp.scale())<<std::endl;
        if(encrypted_C.scale()!=ct_tmp.scale()){
            encrypted_C.scale()=ct_tmp.scale();
        }
        evaluator.add_inplace(encrypted_C, ct_tmp);
    }
    // decrypt_decode(encoder, decryptor, encrypted_C, c);
    // std::cout << "ct_c: " << std::endl;
    // print_vector(c);
    // encrypted_C.scale() = pow(2.0, 50);
}



void encode_encrypt_block_matrix(seal::CKKSEncoder &encoder, seal::Encryptor &encryptor,  double scale, 
                                Matrix<double> &A, Matrix<double> &B, std::size_t b_n, std::size_t b_m, std::size_t b_p, 
                                 std::vector<std::vector<seal::Ciphertext>> &A_encrypted, 
                                 std::vector<std::vector<seal::Ciphertext>> &B_encrypted)
{
    std::size_t m = A.get_cols();
    std::size_t p = B.get_rows();

    if (p != m){
        std::cerr << "the number of cols of A is not equal to the number of rows of B" << std::endl;
        return;
    }

    // std::size_t n = A.get_rows();
    p = B.get_cols();
    

    if ((std::gcd(b_n, b_m) != 1) ||(std::gcd(b_p, b_m) != 1) || (std::gcd(b_n, b_p) != 1)) {
        std::cerr << "encode_encrypt: the dimensions (n, m, p) of the blocked matrices must be coprime" << std::endl;
        return;
    }

    // if (b_p != b_m){
    //     std::cerr << "the number of cols of the blocks of A is not equal to the number of rows of the blocks of B" << std::endl;
    //     return;
    // }

    std::size_t slot_count = encoder.slot_count();
    
    std::vector<std::vector<Matrix<double>>> SA;
    std::vector<std::vector<Matrix<double>>> SB;
    split_matrix(A, b_n, b_m, SA);
    split_matrix(B, b_m, b_p, SB);
    
    std::size_t n_rows_a = SA.size();
    std::size_t n_rows_b = SB.size();
    std::size_t n_cols_a = SA[0].size();
    std::size_t n_cols_b = SB[0].size();

    A_encrypted.resize(n_rows_a);
    B_encrypted.resize(n_rows_b);

    // std::cout << "number of blocked rows: " << n_rows_a << std::endl;
    // std::cout << "number of blocked cols: " << n_cols_a << std::endl;
    // std::cout << "number of resulting blocked cols: " << n_cols_b << std::endl;

    // double scale = encoder.scale();

    for(std::size_t i = 0; i < n_rows_a; i++){
        // A_encrypted[i].resize(n_cols_a);
        for(std::size_t j = 0; j < n_cols_a; j++){
            //! to be written!
            std::vector<double> dA;
            packing_left_matrix(SA[i][j], dA, b_p, slot_count);
            seal::Ciphertext tmp;
            encode_encrypt(encoder, encryptor, dA, scale, tmp);
            A_encrypted[i].push_back(tmp);
        }
    }

    for(std::size_t i = 0; i < n_rows_b; i++){
        // B_encrypted[i].resize(n_cols_b);
        for(std::size_t j = 0; j < n_cols_b; j++){
            //! to be written!
            std::vector<double> dB;
            packing_right_matrix(SB[i][j], dB, b_n, slot_count);
            seal::Ciphertext tmp;
            encode_encrypt(encoder, encryptor, dB, scale, tmp);
            B_encrypted[i].push_back(tmp);
        }
    }
    // std::cout << "number of blocked cols: " << n_rows_b << std::endl;
}

void encode_encrypt_block_matrix_v2(seal::CKKSEncoder &encoder, seal::Encryptor &encryptor,  double scale, 
                                Matrix<double> &A, Matrix<double> &B, std::size_t b_n, std::size_t b_m, std::size_t b_p, 
                                 std::vector<std::vector<seal::Ciphertext>> &A_encrypted, 
                                 std::vector<std::vector<seal::Ciphertext>> &B_encrypted)
{
    std::size_t m = A.get_cols();
    std::size_t p = B.get_rows();

    if (p != m){
        std::cerr << "the number of cols of A is not equal to the number of rows of B" << std::endl;
        return;
    }

    // std::size_t n = A.get_rows();
    p = B.get_cols();
    

    if ((std::gcd(b_n, b_m) != 1) ||(std::gcd(b_p, b_m) != 1) || (std::gcd(b_n, b_p) != 1)) {
        std::cerr << "encode_encrypt: the dimensions (n, m, p) of the blocked matrices must be coprime" << std::endl;
        return;
    }

    // if (b_p != b_m){
    //     std::cerr << "the number of cols of the blocks of A is not equal to the number of rows of the blocks of B" << std::endl;
    //     return;
    // }

    std::size_t slot_count = encoder.slot_count();
    
    std::vector<std::vector<Matrix<double>>> SA;
    std::vector<std::vector<Matrix<double>>> SB;
    split_matrix(A, b_n, b_m, SA);
    split_matrix(B, b_m, b_p, SB);
    
    std::size_t n_rows_a = SA.size();
    std::size_t n_rows_b = SB.size();
    std::size_t n_cols_a = SA[0].size();
    std::size_t n_cols_b = SB[0].size();

    A_encrypted.resize(n_rows_a);
    B_encrypted.resize(n_rows_b);

    // std::cout << "number of blocked rows: " << n_rows_a << std::endl;
    // std::cout << "number of blocked cols: " << n_cols_a << std::endl;
    // std::cout << "number of resulting blocked cols: " << n_cols_b << std::endl;

    // double scale = encoder.scale();

    for(std::size_t i = 0; i < n_rows_a; i++){
        // A_encrypted[i].resize(n_cols_a);
        for(std::size_t j = 0; j < n_cols_a; j++){
            //! to be written!
            std::vector<double> dA;
            packing_left_matrix_v2(SA[i][j], dA, b_p, slot_count);
            seal::Ciphertext tmp;
            encode_encrypt(encoder, encryptor, dA, scale, tmp);
            A_encrypted[i].push_back(tmp);
        }
    }

    for(std::size_t i = 0; i < n_rows_b; i++){
        // B_encrypted[i].resize(n_cols_b);
        for(std::size_t j = 0; j < n_cols_b; j++){
            //! to be written!
            std::vector<double> dB;
            packing_right_matrix_v2(SB[i][j], dB, b_n, slot_count);
            seal::Ciphertext tmp;
            encode_encrypt(encoder, encryptor, dB, scale, tmp);
            B_encrypted[i].push_back(tmp);
        }
    }
    // std::cout << "number of blocked cols: " << n_rows_b << std::endl;
}


void encrypted_block_matrix_multiplication(//seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, 
                                     seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, 
                                     seal::GaloisKeys &galois_keys, uint64_t b_n, uint64_t b_m, uint64_t b_p,
                                     std::vector<std::vector<seal::Ciphertext>> &A_encrypted, 
                                     std::vector<std::vector<seal::Ciphertext>> &B_encrypted,
                                     std::vector<std::vector<seal::Ciphertext>> &C_encrypted)
{
    std::size_t n_rows_a = A_encrypted.size();
    std::size_t n_rows_b = B_encrypted.size();
    std::size_t n_cols_a = A_encrypted[0].size();
    std::size_t n_cols_b = B_encrypted[0].size();

    if (n_cols_a != n_rows_b){
        std::cerr << "the number of blocked cols A is not equal to the number of blocked rows of B" << std::endl;
        return;
    }

    

    C_encrypted.resize(n_rows_a);

    for(std::size_t i = 0; i < n_rows_a; i++){
        seal::Ciphertext tmp1;
        for(std::size_t j = 0; j < n_cols_b; j++){
            std::vector<seal::Ciphertext> t;
            seal::Ciphertext tmp;
            for(std::size_t k = 0; k < n_cols_a; k++){
                encrypted_matrix_multiplication(evaluator, relin_keys, galois_keys, 
                                                b_n, b_m, b_p, A_encrypted[i][k], 
                                                B_encrypted[k][j], tmp);
                t.push_back(tmp);
            }
            evaluator.add_many(t, tmp1);
            C_encrypted[i].push_back(tmp1);
        }
    }

    // std::cout << "encrypted_block_matrix_multiplication: (n, p): (" <<  C_encrypted.size() << ", " << C_encrypted[0].size() << ")" << std::endl;
}
void encrypted_block_matrix_multiplication_v2(seal::CKKSEncoder &encoder, seal::Decryptor &decryptor, 
                                     seal::Evaluator &evaluator, seal::RelinKeys &relin_keys, 
                                     seal::GaloisKeys &galois_keys, uint64_t b_n, uint64_t b_m, uint64_t b_p,
                                     std::vector<std::vector<seal::Ciphertext>> &A_encrypted, 
                                     std::vector<std::vector<seal::Ciphertext>> &B_encrypted,
                                     std::vector<std::vector<seal::Ciphertext>> &C_encrypted)
{
    std::size_t n_rows_a = A_encrypted.size();
    std::size_t n_rows_b = B_encrypted.size();
    std::size_t n_cols_a = A_encrypted[0].size();
    std::size_t n_cols_b = B_encrypted[0].size();

    if (n_cols_a != n_rows_b){
        std::cerr << "the number of blocked cols A is not equal to the number of blocked rows of B" << std::endl;
        return;
    }

    

    C_encrypted.resize(n_rows_a);

    for(std::size_t i = 0; i < n_rows_a; i++){
        seal::Ciphertext tmp1;
        for(std::size_t j = 0; j < n_cols_b; j++){
            std::vector<seal::Ciphertext> t;
            seal::Ciphertext tmp;
            for(std::size_t k = 0; k < n_cols_a; k++){
                encrypted_matrix_multiplication_v2(encoder,decryptor,evaluator, relin_keys, galois_keys, 
                                                b_n, b_m, b_p, A_encrypted[i][k], 
                                                B_encrypted[k][j], tmp);
                t.push_back(tmp);
            }
            evaluator.add_many(t, tmp1);
            C_encrypted[i].push_back(tmp1);
        }
    }

    // std::cout << "encrypted_block_matrix_multiplication: (n, p): (" <<  C_encrypted.size() << ", " << C_encrypted[0].size() << ")" << std::endl;
}


