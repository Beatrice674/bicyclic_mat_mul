#include "../src/matrix.h"
// #include "../src/long_enc_cipher.h"
#include "../src/lu_enc_cipher.h"
#include "../src/jiang_enc_cipher.h"
#include "../src/R_T_enc_cipher.h"
#include "../src/strassen_block.h"
// #include "../src/coprime_enc_cipher.h"
#include "../src/block_long_cipher.h"
#include <cstdlib>
#include <chrono>

using namespace std;
using namespace seal;

void run_coprime_mat_mul_I(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void run_coprime_mat_mul_II(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_coprime_mat_mul(size_t n, size_t m, size_t p,
                          CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator,
                          RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor, int method_choice);

void run_long_mat_mul_I(size_t n, size_t m, size_t p,size_t poly_modulus_degree);
void run_long_mat_mul_II(size_t n, size_t m, size_t p,size_t poly_modulus_degree);
void test_long_mat_mul(size_t n, size_t m, size_t p,
                       CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator,
                       RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);

void run_naive_block_long_mat_mul(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree);
void test_naive_block_long_mat_mul(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, 
                                      CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                                      RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);

void run_lu_mat_mul(size_t n, size_t m, size_t p, size_t batch_size);
void test_lu_mat_mul(size_t n, size_t m, size_t p, size_t batch_size,
                     CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator,
                     RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);

void run_jiang_mat_mul(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_jiang_mat_mul(size_t n, size_t m, size_t p, 
                        CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                        RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);

void run_r_t_mat_mul(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_r_t_mat_mul(size_t n, size_t m, size_t p, 
                      CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                      RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);

void run_strassen_block_jiang_mat_mul(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void run_strassen_block_jiang_mat_mul_open(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_strassen_block_jiang_mat_mul(size_t n, size_t m, size_t p, 
                                       CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                                       RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);
void test_strassen_block_jiang_mat_mul_open(size_t n, size_t m, size_t p, 
                                       CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                                       RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);

void run_strassen_block_coprime_mat_mul(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree);
void run_strassen_block_coprime_mat_mul_open(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree);
void run_strassen_block_coprime_mat_mul_I(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree);
void test_strassen_block_coprime_mat_mul(size_t n, size_t m, size_t p, size_t s_n,size_t s_m,size_t s_p,
                                         CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                                         RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor, int method_choice);

void test_strassen_block_coprime_mat_mul_open(size_t n, size_t m, size_t p, size_t s_n,size_t s_m,size_t s_p,
                                         CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                                         RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor, int method_choice);


void run_naive_block_jiang_mat_mul(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_naive_block_jiang_mat_mul(size_t n, size_t m, size_t p, 
                                    CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                                    RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);

void run_naive_block_coprime_mat_mul(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree);
void test_naive_block_coprime_mat_mul(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, 
                                      CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                                      RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor, int method_choice);


void run_coprime_mat_mul_imporve_log(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void run_coprime_mat_mul_imporve_log_I(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_coprime_mat_mul_imporve_log(size_t n, size_t m, size_t p, CKKSEncoder &encoder, 
                                        Encryptor &encryptor, double scale, Evaluator &evaluator, 
                                        RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);


void run_naive_block_coprime_mat_mul_imporve_log(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree);
void run_naive_block_coprime_mat_mul_imporve_log_I(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree);
void test_naive_block_coprime_mat_mul_imporve_log(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);


void run_strassen_block_coprime_mat_mul_imporve_log(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree);
void run_strassen_block_coprime_mat_mul_imporve_log_I(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree);
void test_strassen_block_coprime_mat_mul_imporve_log(size_t n, size_t m, size_t p, size_t s_n, size_t s_m, size_t s_p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);

int main()
{
    while (true)
    {
        cout << "\nExamples:" << endl
             << endl;
        cout << "  1.  Bicyclic matrix multiplication for small dimensions" << endl;
        cout << "  2.  Bicyclic matrix multiplication for large dimensions" << endl;
        cout << "  3.  OpenMP optimization" << endl;
        // cout << "  4.  Lu matrix multiplication base" << endl;
        // cout << "  5.  R-T matrix multiplication base" << endl;
        // cout << "  6.  Block matrix multiplication by LongRot action" << endl;
        // cout << "  7.  Block matrix multiplication strassen and naive(jiang)" << endl;
        // cout << "  8.  Block matrix multiplication strassen and naive(Bicyclic)" << endl;
        // cout << "  9.  Base matrix multiplication in Table 9" << endl;
        // cout << " 10.  Base matrix multiplication in Table 10" << endl;
        // cout << " 11.  Block matrix multiplication in Table 11(128,128,128)" << endl;
        // cout << " 12.  Block matrix multiplication in Table 12(256,256,256)" << endl;
        // cout << " 13.  Block matrix multiplication in Table 13(512,512,512)" << endl;
        // cout << " 14.  Block matrix multiplication in Table 14" << endl;
        // cout << " 15.  Block matrix multiplication in Table 15" << endl;
        // cout << " 16.  Block matrix multiplication in Table 16" << endl;
        cout << " 0. Exit" << endl;

        /*
        Print how much memory we have allocated from the current memory pool.
        By default the memory pool will be a static global pool and the
        MemoryManager class can be used to change it. Most users should have
        little or no reason to touch the memory allocation system.
        */
        cout << "\nTotal memory allocated from the current memory pool: "
             << (MemoryManager::GetPool().alloc_byte_count() >> 20) << " MB" << endl;

        int selection = 0;
        cout << endl
             << "Run example: ";
        if (!(cin >> selection))
        {
            cout << "Invalid option." << endl;
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            continue;
        }

        switch (selection)
        {
        case 1:
            // small dimensions
            cout << "Running coprime matrix multiplication with small dimensions..." << endl;
            run_coprime_mat_mul_I(64, 31, 29, 8192);
            run_coprime_mat_mul_I(44, 45, 43, 8192);
            run_coprime_mat_mul_I(31, 64, 29, 8192);

            run_coprime_mat_mul_I(128, 31, 29, 16384);
            run_coprime_mat_mul_I(61, 64, 63, 16384);
            run_coprime_mat_mul_I(29, 128, 31, 16384);

            run_coprime_mat_mul_I(128, 63, 61, 32768);
            run_coprime_mat_mul_I(89, 91, 90, 32768);
            run_coprime_mat_mul_I(63, 128, 61, 32768); 

            // run_coprime_mat_mul_I(8, 8, 8, 8192);
            // run_coprime_mat_mul_I(44, 45, 43, 8192); // 140: 0.338726 s
            // run_coprime_mat_mul_I(31, 64, 29, 8192); 
            // run_coprime_mat_mul_I(61, 64, 63, 16384);// 140: 1.16957 s
            // run_coprime_mat_mul_I(29, 128, 31, 16384);// 140: 1.16957 s
            // run_coprime_mat_mul_I(89, 91, 90, 32768);// 140: 3.70669 s
            // run_coprime_mat_mul_I(63, 128, 61, 32768); 
            // test for different dims
            // run_coprime_mat_mul_I(64, 31, 29, 8192); 
            // run_coprime_mat_mul_I(128, 31, 29, 16384);
            // run_coprime_mat_mul_I(128, 63, 61, 32768);
            // with zheng
            // run_coprime_mat_mul_I(16, 19, 17, 8192);
            // run_coprime_mat_mul_I(32, 35, 33, 8192);
            // run_strassen_block_coprime_mat_mul(2*32,2*35,2*33,32,35,33,8192);

            // run_coprime_mat_mul_I(16, 19, 17, 16384);
            // run_coprime_mat_mul_I(32, 35, 33, 16384);
            // run_strassen_block_coprime_mat_mul(2*32,2*35,2*33,32,35,33,16384);
            break;
        // case 2:
        //     run_coprime_mat_mul_imporve_log_I(15, 16, 17, 8192);// 140: 0.012744 s
        //     run_coprime_mat_mul_imporve_log_I(21, 16, 23, 16384);// 140: 0.031736 s
        //     run_coprime_mat_mul_imporve_log_I(31, 16, 33, 32768);// 140: 0.051467 s
        //     break;
        // case 3:
        //     run_jiang_mat_mul(64,64,64,8192);// 1.64031 s
        //     run_jiang_mat_mul(128,128,128,32768);//17.4258 s
        //     break;
        // case 4:
        //     run_lu_mat_mul(16,16,5*4096,8);// 2.66471 s
        //     run_lu_mat_mul(16,16,5*4096,1);// 6.77389 s
        //     run_lu_mat_mul(64,64,5*4096,8);// 36.1684 s
        //     run_lu_mat_mul(64,64,5*4096,1);// 102.865 s
        //     break;
        // case 5:
        //     run_r_t_mat_mul(16,16,16,8192); // 0.229158 s
        //     break;
        // case 6:
        //     run_long_mat_mul_I(256, 257, 17,8192); // keygen: 19.6271 s , Matmul: 11.5917 s
        //     run_long_mat_mul_I(256, 17, 257,8192); // keygen: 19.3821 s , Matmul: 5.90978 s

        //     run_long_mat_mul_II(256, 257, 17,8192); // 16.7642 s
        //     run_long_mat_mul_II(256, 17, 257,8192); // 9.10757 s
        //     break;
        // case 7:
        //     run_naive_block_jiang_mat_mul(4*64,4*64,4*64,8192);// 105.528 s
        //     run_strassen_block_jiang_mat_mul(4*64,4*64,4*64,8192);// 81.4588s

        //     run_naive_block_jiang_mat_mul(8*64,8*64,8*64,8192);// 868.31 s
        //     run_strassen_block_jiang_mat_mul(8*64,8*64,8*64,8192);// 572.986 s
        //     break;
        // case 8:
        //     run_naive_block_coprime_mat_mul(4*44,4*45,4*43,44,45,43,8192);// 21.9163 s
        //     run_strassen_block_coprime_mat_mul(4*44,4*45,4*43,44,45,43,8192);// 17.0068 s

        //     run_naive_block_coprime_mat_mul(8*44,8*45,8*43,44,45,43,8192);// 178.584 s
        //     run_strassen_block_coprime_mat_mul(8*44,8*45,8*43,44,45,43,8192);// 121.481 s
        //     break;
        // case 9:
        //     run_r_t_mat_mul(16,16,16,8192);// 0.259183 s
        //     run_coprime_mat_mul_II(16,19,17,8192);// 0.169384 s
        //     run_coprime_mat_mul_I(16,19,17,8192);// 0.105432 s
        //     // run_coprime_mat_mul_imporve_log(15, 17, 16, 8192);// 170: 0.018019 s
        //     run_coprime_mat_mul_imporve_log_I(15, 16, 17, 8192);// 140: 0.012744 s
        //     break;
        // case 10:
        //     run_jiang_mat_mul(64,64,64,8192);// 1.83336 s
        //     run_jiang_mat_mul(64,64,64,32768);// 18.8702 s

        //     run_coprime_mat_mul_I(43,45,44,8192);// 0.349032 s
        //     run_coprime_mat_mul_I(61,64,63,16384);// 1.22242 s
        //     run_coprime_mat_mul_I(89,91,90,32768);// 3.88849 s

        //     // run_coprime_mat_mul_imporve_log(15, 17, 16, 8192);// 170: 0.018019 s
        //     // run_coprime_mat_mul_imporve_log(19, 21, 20, 16384);// 170: 0.047813 s
        //     // run_coprime_mat_mul_imporve_log(23, 27, 25, 32768);// 170: 0.099075 s

        //     run_coprime_mat_mul_imporve_log_I(15, 16, 17, 8192);// 140: 0.012744 s
        //     run_coprime_mat_mul_imporve_log_I(21, 16, 23, 16384);// 140: 0.031736 s
        //     run_coprime_mat_mul_imporve_log_I(31, 16, 33, 32768);// 140: 0.051467 s
        //     break;
        case 2:
            // (64,64,64)
            cout << "-----------------------(64,64,64)-----------------------" << endl;
            run_jiang_mat_mul(64,64,64,8192);// 1.83336 s
    
            run_naive_block_coprime_mat_mul(2*32,2*35,2*33,32,35,33,8192);
            run_strassen_block_coprime_mat_mul(2*32,2*35,2*33,32,35,33,8192);

            run_coprime_mat_mul_I(64, 65, 67, 32768); 

            //(128,128,128)
            cout << "-------------------------(128,128,128)---------------------" << endl;
            run_naive_block_jiang_mat_mul(2*64,2*64,2*64,8192);// 12.347194 s
            run_naive_block_jiang_mat_mul(128,128,128,32768);// 16.251344 s

            run_strassen_block_jiang_mat_mul(2*64,2*64,2*64,8192);// 10.848159 s

            run_naive_block_coprime_mat_mul(3*43,3*45,3*44,43,45,44,8192);// 8.131559 s
            run_naive_block_coprime_mat_mul(2*64,2*67,2*65,64,67,65,32768);// 14.546507 s

            run_strassen_block_coprime_mat_mul(4*32,4*35,4*33,32,35,33,8192);// 8.95519 s
            run_strassen_block_coprime_mat_mul(2*64,2*67,2*65,64,67,65,32768);// 13.010182 s

            //(256,256,256)
            cout << "-------------------------(256,256,256)---------------------" << endl;
            run_naive_block_jiang_mat_mul(4*64,4*64,4*64,8192);// 89.33 s
            run_naive_block_jiang_mat_mul(2*128,2*128,2*128,32768);// 112.01 s

            run_strassen_block_jiang_mat_mul(4*64,4*64,4*64,8192);// 71.54 s
            run_strassen_block_jiang_mat_mul(2*128,2*128,2*128,32768);// 97.70 s

            run_naive_block_coprime_mat_mul(6*44,6*45,6*43,44,45,43,8192);// 59.57 s
            run_naive_block_coprime_mat_mul(3*88,3*89,3*87,88,89,87,32768);// 73.35 s

            run_strassen_block_coprime_mat_mul(8*32,8*35,8*33,32,35,33,8192);// 56.98 s
            run_strassen_block_coprime_mat_mul(4*64,4*67,4*65,64,67,65,32768);// 80.99 s

            //(512,512,512)
            cout << "-------------------------(512,512,512)---------------------" << endl;
            // 8192
            run_naive_block_jiang_mat_mul(8*64,8*64,8*64,8192);// 727.73 s
            run_strassen_block_jiang_mat_mul(8*64,8*64,8*64,8192);// 479.19 s

            run_naive_block_coprime_mat_mul(12*43,12*45,12*44,43,45,44,8192);// 489.83 s
            run_strassen_block_coprime_mat_mul(16*32,16*35,16*33,32,35,33,8192);// 390.26 s

            // 32768
            run_naive_block_jiang_mat_mul(4*128,4*128,4*128,32768);
            run_strassen_block_jiang_mat_mul(4*128,4*128,4*128,32768);
            
            run_naive_block_coprime_mat_mul(6*88,6*89,6*87,88,89,87,32768);
            run_strassen_block_coprime_mat_mul(8*64,8*67,8*65,64,67,65,32768);

            break;
        // case 12:
        //     //(256,256,256)
        //     run_naive_block_jiang_mat_mul(4*64,4*64,4*64,8192);// 89.33 s
        //     run_naive_block_jiang_mat_mul(2*128,2*128,2*128,32768);// 112.01 s

        //     run_strassen_block_jiang_mat_mul(4*64,4*64,4*64,8192);// 71.54 s
        //     run_strassen_block_jiang_mat_mul(2*128,2*128,2*128,32768);// 97.70 s

        //     run_naive_block_coprime_mat_mul(6*44,6*45,6*43,44,45,43,8192);// 59.57 s
        //     run_naive_block_coprime_mat_mul(3*88,3*89,3*87,88,89,87,32768);// 73.35 s

        //     run_strassen_block_coprime_mat_mul(8*32,8*35,8*33,32,35,33,8192);// 56.98 s
        //     run_strassen_block_coprime_mat_mul(4*64,4*67,4*65,64,67,65,32768);// 80.99 s

        //     // run_naive_block_coprime_mat_mul_imporve_log_I(18*15,16*16,16*17,15,16,17,8192);// 81.6179 s
        //     // run_naive_block_coprime_mat_mul_imporve_log_I(13*21,8*32,12*23,21,32,23,32768);//90.0473 s

        //     // run_long_mat_mul_I(256, 259, 257,8192); // 42.98 s
        //     // // run_long_mat_mul_I(256, 259, 257,32768); // 160 G memory, 281 for keygen + 86s 
            
        //     // run_long_mat_mul_II(256, 259, 257,32768); // 110.99 s
        //     // // run_long_mat_mul_II(256, 259, 257,8192); // 144.362 s
        //     break;
        // case 13:
            // (512,512,512)
            // 8192
            // run_naive_block_jiang_mat_mul(8*64,8*64,8*64,8192);// 727.73 s
            // run_strassen_block_jiang_mat_mul(8*64,8*64,8*64,8192);// 479.19 s

            // run_naive_block_coprime_mat_mul(12*43,12*45,12*44,43,45,44,8192);// 489.83 s
            // run_strassen_block_coprime_mat_mul(16*32,16*35,16*33,32,35,33,8192);// 390.26 s

            // // 32768
            // run_naive_block_jiang_mat_mul(4*128,4*128,4*128,32768);
            // run_strassen_block_jiang_mat_mul(4*128,4*128,4*128,32768);
            
            // run_naive_block_coprime_mat_mul(6*88,6*89,6*87,88,89,87,32768);
            // run_strassen_block_coprime_mat_mul(8*64,8*67,8*65,64,67,65,32768);

            // run_long_mat_mul_I(512,515,513,8192);// 181.03 s
            
            // (1024,1024,1024)
            // run_naive_block_jiang_mat_mul(16*64,16*64,16*64,8192);// 6028.15 s
            // run_naive_block_coprime_mat_mul(24*43,24*45,24*44,43,45,44,8192);// 4239.80 s

            // run_strassen_block_jiang_mat_mul(16*64,16*64,16*64,8192);// 3514.14 s
            // run_strassen_block_coprime_mat_mul_I(32*32,32*35,32*33,32,35,33,8192);// log\Delta = 40. 2757.16 s
            // run_long_mat_mul_I(1024,1027,1025,8192);// 1200.76 s
        //     break;
        // case 14:
        //     // Algo.5
        //     run_long_mat_mul_I(256, 257, 17,8192); // keygen: 16.15 s , total: 23.60 s
        //     run_long_mat_mul_I(256, 17, 257,8192); // keygen: 14.69 s , total: 19.37 s

        //     // run_long_mat_mul_II(256, 257, 17,8192); // 16.7642 s
        //     // run_long_mat_mul_II(256, 17, 257,8192); // 9.10757 s

        //     run_long_mat_mul_I(1024, 1025, 17,8192); // keygen: 14.68 s , total: 55.79 s
        //     run_long_mat_mul_I(1024, 17, 1025,8192); // keygen: 14.47 s , total: 43.64 s

        //     run_long_mat_mul_I(2048, 2049, 11,8192); // keygen: 14.49 s , total: 98.01 s
        //     run_long_mat_mul_I(2049, 8, 2051,8192); // keygen: 14.63 s , total: 81.09 s

        //     // run_long_mat_mul_I(8191, 8192, 15, 8192); // keygen: 27.6571 s , total: 139.218 s
        //     // run_long_mat_mul_I(8191, 15, 8192, 8192); // keygen: 20.5851 s , total: 112.152 s
        //     break;
        // case 15:
        //     // Naive block(Jiang)
        //     run_naive_block_jiang_mat_mul(4,1636,5,8192);//36.92 s
        //     run_naive_block_jiang_mat_mul(8,3405,9,8192);//78.58 s
        //     run_naive_block_jiang_mat_mul(16,6903,17,8192);//157.00 s
        //     run_naive_block_jiang_mat_mul(32,13847,33,8192);//317.31 s
        //     // run_naive_block_jiang_mat_mul(32*2,13847*2,33*2,8192);

        //     // Naive block(Algo.5)
        //     run_naive_block_coprime_mat_mul(4, 4*409, 5, 4,409,5,8192);//9.76 s
        //     run_naive_block_coprime_mat_mul(8, 15*227, 9, 8,227,9,8192);//19.92 s
        //     run_naive_block_coprime_mat_mul(16,59*117,17,16,117,17,8192);//38.57 s
        //     run_naive_block_coprime_mat_mul(32,227*61,33,32,61,33,8192);//76.73 s
            

        //     break;
        // case 16:
        //     // Jiang
        //     run_naive_block_jiang_mat_mul(4,5,1636,8192);//36.31 s
        //     run_naive_block_jiang_mat_mul(8,9,3405,8192);//77.56 s
        //     run_naive_block_jiang_mat_mul(16,17,6903,8192);//157.32 s
        //     run_naive_block_jiang_mat_mul(32,33,13847,8192);//316.76 s
        //     // run_naive_block_jiang_mat_mul(63,64,16383,8192);//490.316 s

        //     // algo.3
        //     run_naive_block_coprime_mat_mul(4,5,4*409,4,5,409,8192);//0.10 s
        //     run_naive_block_coprime_mat_mul(8,9,15*227,8,9,227,8192);//0.65 s
        //     run_naive_block_coprime_mat_mul(16,17,59*117,16,17,117,8192);//4.96 s
        //     run_naive_block_coprime_mat_mul(32,33,227*61,32,33,61,8192);//38.82 s
        //     // run_naive_block_coprime_mat_mul(63,64,43*381,63,64,43,32768);//929.822 s

        //     // algo.7
        //     run_lu_mat_mul(4,5,1636,8);// 0.302308 s
        //     run_lu_mat_mul(8,9,3405,8);// 0.774307 s
        //     run_lu_mat_mul(16,17,6903,8);// 2.64838 s
        //     run_lu_mat_mul(32,33,8225,8);// 9.78487 s
        //     // run_lu_mat_mul(63,64,13847,8);// 36.022 s
        //     break;

        case 3:
            // openmp
            // (64,64,64)
            cout << "-----------------------(64,64,64)-----------------------" << endl;
            cout << "-----------------------8192 part-----------------------" << endl;
            run_jiang_mat_mul(64,64,64,8192);// 1.83336 s
            run_strassen_block_coprime_mat_mul_open(2*32,2*35,2*33,32,35,33,8192);

            cout << "-----------------------32768 part-----------------------" << endl;            
            run_coprime_mat_mul_I(64, 65, 67, 32768); 

            //(128,128,128)
            cout << "-------------------------(128,128,128)---------------------" << endl;
            cout << "-----------------------8192 part-----------------------" << endl;
            run_strassen_block_jiang_mat_mul(2*64,2*64,2*64,8192);// 10.848159 s
            run_strassen_block_coprime_mat_mul_open(4*32,4*35,4*33,32,35,33,8192);// 8.95519 s

            cout << "-----------------------32768 part-----------------------" << endl;  
            run_naive_block_jiang_mat_mul(128,128,128,32768);// 16.251344 s
            run_strassen_block_coprime_mat_mul_open(2*64,2*67,2*65,64,67,65,32768);// 13.010182 s

            //(256,256,256)
            cout << "-------------------------(256,256,256)---------------------" << endl;
            cout << "-----------------------8192 part-----------------------" << endl;
            run_strassen_block_jiang_mat_mul(4*64,4*64,4*64,8192);// 71.54 s
            run_strassen_block_coprime_mat_mul_open(8*32,8*35,8*33,32,35,33,8192);// 56.98 s
            
            cout << "-----------------------32768 part-----------------------" << endl;
            run_strassen_block_jiang_mat_mul(2*128,2*128,2*128,32768);// 97.70 s
            run_strassen_block_coprime_mat_mul_open(4*64,4*67,4*65,64,67,65,32768);// 80.99 s

            // (512,512,512)
            cout << "-------------------------(512,512,512)---------------------" << endl;
            cout << "-----------------------8192 part-----------------------" << endl;
            // 8192
            run_strassen_block_jiang_mat_mul(8*64,8*64,8*64,8192);// 479.19 s
            run_strassen_block_coprime_mat_mul_open(16*32,16*35,16*33,32,35,33,8192);// 390.26 s

            // 32768
            cout << "-----------------------32768 part-----------------------" << endl;
            run_strassen_block_jiang_mat_mul(4*128,4*128,4*128,32768);
            run_strassen_block_coprime_mat_mul_open(8*64,8*67,8*65,64,67,65,32768);

            break;
        //  case 17:
        //     // run_naive_block_coprime_mat_mul_imporve_log_I(98*21,8,90*23,21,8,23,8192);//95.4541 s

        //     //Table VIII and IX
        //     run_coprime_mat_mul_imporve_log_I(15, 16, 17, 8192);// 140: 0.012744 s
        //     run_coprime_mat_mul_imporve_log_I(21, 16, 23, 16384);// 140: 0.031736 s
        //     run_coprime_mat_mul_imporve_log_I(31, 16, 33, 32768);// 140: 0.051467 s

        //     //Table XI(256,256,256)
        //     run_naive_block_coprime_mat_mul_imporve_log_I(18*15,16*16,16*17,15,16,17,8192);// 81.6179 s
        //     run_naive_block_coprime_mat_mul_imporve_log_I(13*21,8*32,12*23,21,32,23,32768);//90.0473 s

        //     //Table XII(512,512,512) and (1024,1024,1024)
        //     run_naive_block_coprime_mat_mul_imporve_log_I(35*15,32*16,31*17,15,16,17,8192);
        //     run_naive_block_coprime_mat_mul_imporve_log_I(69*15,64*16,61*17,15,16,17,8192);

        //     // Table  Naive block(Algo.6)
        //     run_naive_block_coprime_mat_mul_imporve_log(4,9*203,5,4,203,5,8192);//0.200701 s
        //     run_naive_block_coprime_mat_mul_imporve_log(8,62*55,9,8,55,9,8192);//1.516552 s
        //     run_naive_block_coprime_mat_mul_imporve_log(16,461*15,17,16,15,17,8192);// 23.020219 s
        //     run_naive_block_coprime_mat_mul_imporve_log(2*16,924*15,2*17,16,15,17,8192);//513.067469 s
        //     break;

        // case 18:
        //     run_naive_block_long_mat_mul(2*128,2*131,2*129,128,131,129,8192);// 30.8473 s + 61.9948 s
        //     run_naive_block_long_mat_mul(2*256,2*259,2*257,256,259,257,8192);//33.0755 s + 256.217 s
        //     run_naive_block_long_mat_mul(2*512,2*515,2*513,512,515,513,8192);// 25.6464 s +1547.26 s
        //     break;
        
        case 0:
            return 0;
        default:
            cout << "Invalid option." << endl;
        }
    }
    return 0;
}


void run_coprime_mat_mul_I(size_t n, size_t m, size_t p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    GaloisKeys gal_keys;

    cout << "Please choose the method: 1 for our method, 2 for chen et al. method" << endl;
    int method_choice;
    cin >> method_choice;
    if (method_choice == 1)
    {
        // xyy add: seg-hoisting
        vector<int> step_s;
        size_t size_m = (log2(poly_modulus_degree/2))/2 - 1;
        // cout << "size_m = " << size_m << endl;
        // for (size_t i = 1; i < size_t(sqrt(m)) + 1; i++) 
        for (size_t i = 1; i < size_m + 1; i++)
        {
            step_s.push_back(mod(i * n, n*m));
            step_s.push_back(mod(i * p, m*p));
            // cout << mod(i * n, n*m) << ", ";
            // cout << mod(i * p, m*p) << ", ";
        }
        cout << endl;
        keygen.create_galois_keys(step_s, gal_keys);
    }
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;
    if (method_choice == 2)
    {
        // chen: native method
        keygen.create_galois_keys(gal_keys);
    }

    test_coprime_mat_mul(n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor, method_choice);
}

void run_coprime_mat_mul_II(size_t n, size_t m, size_t p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30,30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);
    cout << endl;

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_coprime_mat_mul(n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor, 2);
}

void test_coprime_mat_mul(size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor, int method_choice)
{
    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    // double total_time_us = 0.0;
    // int repeat = 1000;
    // for (int i = 0; i < repeat; i++) {
    auto start_time = chrono::high_resolution_clock::now();

    cipher_matrix_coprime cipher_A(A,ceil(double(p)/m)+1,scale,encoder,encryptor);
    cipher_matrix_coprime cipher_B(B,ceil(double(n)/m)+1,scale,encoder,encryptor);
    cipher_matrix_coprime cipher_C;

    // auto end_time = chrono::high_resolution_clock::now();
    // auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    // cout << "The encode timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    encrypted_coprime_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A,cipher_B,cipher_C,method_choice);

    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
        << "----------------- Diag encode matrix multiplication --------------------" << endl;
    cout << endl
        << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // total_time_us += time_diff.count();

    Matrix<double> computed_C;
    cipher_C.dec_matrix_cipher(computed_C,encoder,decryptor);

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl <<endl;
    // }
    // double avg_time_us = total_time_us / repeat;         // 平均微秒
    // double avg_time_ms = avg_time_us / 1000.0;           // 平均毫秒
    // cout << endl << "----------------- Diag encode matrix multiplication --------------------" << endl;
    // cout << "Average time over " << repeat << " runs: " << avg_time_ms << " ms" << endl;
}

void run_long_mat_mul_I(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 40, 40, 60}));
    double scale = pow(2.0, 40);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    auto start_time = chrono::high_resolution_clock::now();
    vector<int> steps;
    for (size_t i = 0; i < slot_count; i++)
    {
        steps.push_back(i);
    }
    GaloisKeys gal_keys;
    keygen.create_galois_keys(steps, gal_keys);

    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- keygen galois_keys --------------------" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    test_long_mat_mul(n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void run_long_mat_mul_II(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 40, 40, 60}));
    double scale = pow(2.0, 40);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_long_mat_mul(n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void test_long_mat_mul(size_t n, size_t m, size_t p,
                       CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator,
                       RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{

    if ((gcd(n, m) != 1) || (gcd(m, p) != 1) || (gcd(n, p) != 1))
    {
        cerr << "the dimensions (" << n << ", " << m << ", " << p << ") are not coprime" << endl;
        return;
    }

    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print();
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();

    vector<double> da, dA, db, dc;

    da = A.diagonal_encoding();
    // cout << "diagonal encoding of A: " << endl;
    // print_vector(da);

    db = B.diagonal_encoding();
    // cout << "diagonal encoding of B: " << endl;
    // print_vector(db);

    long_cipher encrypted_A(n, m, da, scale, encoder, encryptor);
    long_cipher encrypted_B(m, p, db, scale, encoder, encryptor);
    long_cipher encrypted_C(n, p, scale, encoder, encryptor);

    encrypted_long_matrix_multiplication(encoder, scale, evaluator, relin_keys, gal_keys, encrypted_A, encrypted_B, encrypted_C);

    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- Diag LongRot matrix multiplication --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result
    vector<double> decrypted_C;
    encrypted_C.dec_diag_vector(decrypted_C, encoder, decryptor);

    Matrix<double> computed_C;
    computed_C = diagonal_decoding(decrypted_C, n, p);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "-------------------------------------------------------------------------" << endl;
    cout << endl <<endl;
}

void run_naive_block_long_mat_mul(size_t n, size_t m, size_t p, size_t s_n, size_t s_m, size_t s_p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 40, 40, 60}));
    double scale = pow(2.0, 40);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);

    size_t slot_count = encoder.slot_count();
    auto start_time = chrono::high_resolution_clock::now();
    vector<int> steps;
    for (size_t i = 0; i < slot_count; i++)
    {
        steps.push_back(i);
    }
    GaloisKeys gal_keys;
    keygen.create_galois_keys(steps, gal_keys);

    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- keygen galois_keys --------------------" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;


    // GaloisKeys gal_keys;
    // keygen.create_galois_keys(gal_keys);
    test_naive_block_long_mat_mul(n, m, p, s_n, s_m, s_p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void test_naive_block_long_mat_mul(size_t n, size_t m, size_t p, size_t s_n, size_t s_m, size_t s_p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();

    long_matrix_block encrypted_A(A,s_n,s_m,scale,encoder,encryptor);
    long_matrix_block encrypted_B(B,s_m,s_p,scale,encoder,encryptor);

    long_matrix_block encrypted_C;
    encrypted_C=Naive_long_block_coprime(s_n,s_m,s_p ,encrypted_A,encrypted_B,scale,encoder,evaluator,gal_keys,relin_keys,encryptor);
    
   
    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- Long cipher naive block matrix multiplication --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result

    Matrix<double> computed_C;
    encrypted_C.dec_block_long_cipher(computed_C,s_n,s_p,encoder,decryptor);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "-----------------------------------------------------------------------------" << endl;
    cout << endl <<endl;
}

void run_lu_mat_mul(size_t n, size_t m, size_t p, size_t batch_size)
{
    EncryptionParameters parms(scheme_type::ckks);
    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30, 30, 30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_lu_mat_mul(n, m, p, batch_size, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void test_lu_mat_mul(size_t n, size_t m, size_t p, size_t batch_size,
                     CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator,
                     RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{

    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();

    vector<vector<double>> da, db, result_c;
    da = A.get_matrix_by_rows();
    db = B.get_matrix_by_rows();

    lu_cipher_A cipher_mat_a(da, scale, encoder, encryptor, batch_size);
    lu_cipher_B cipher_mat_b(db, scale, encoder, encryptor);
    lu_cipher_B cipher_mat_c;

    encrypted_lu_matrix_multiplication(encoder, scale, evaluator, relin_keys, gal_keys, cipher_mat_a, cipher_mat_b, cipher_mat_c);

    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- Lu matrix multiplication --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "Batch size: " << batch_size << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result
    Matrix<double> computed_C;
    cipher_mat_c.dec_lu_matrix(result_c, encoder, decryptor);
    computed_C.gen_matrix_by_rows(result_c);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "---------------------------------------------------------------" << endl;
    cout << endl <<endl;
}

void run_jiang_mat_mul(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30, 30, 30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_jiang_mat_mul(n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void test_jiang_mat_mul(size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    Matrix<double> A(n, m);
    Matrix<double> B(m, p);
    size_t d=sqrt(encoder.slot_count());

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();

    cipher_matrix_jiang cipher_A(A,scale,encoder,encryptor);
    cipher_matrix_jiang cipher_B(B,scale,encoder,encryptor);
    cipher_matrix_jiang cipher_C;

    cipher_A.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A,cipher_B,cipher_C);

    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- Jiang matrix multiplication --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "Actual Multiply:"
         << " (" << d << "x" << d << ")X(" << d << "x" << d << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result

    Matrix<double> computed_C;
    cipher_C.dec_matrix_cipher(computed_C,encoder,decryptor);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "------------------------------------------------------------------" << endl;
    cout << endl <<endl;
}

void run_r_t_mat_mul(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30, 30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_r_t_mat_mul(n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void test_r_t_mat_mul(size_t n, size_t m, size_t p, 
                      CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                      RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();

    cipher_matrix_r_t cipher_A(A,scale,encoder,encryptor);
    cipher_matrix_r_t cipher_B(B,scale,encoder,encryptor);
    cipher_matrix_r_t cipher_C;

    cipher_A.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);

    encrypted_r_t_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A,cipher_B,cipher_C);

    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- R-T matrix multiplication --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result
    Matrix<double> computed_C;
    cipher_C.dec_matrix_cipher(computed_C,encoder,decryptor);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "----------------------------------------------------------------" << endl;
    cout << endl <<endl;
}


void run_strassen_block_jiang_mat_mul(size_t n, size_t m, size_t p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30,30,30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;
    
    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_strassen_block_jiang_mat_mul(n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void run_strassen_block_jiang_mat_mul_open(size_t n, size_t m, size_t p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30,30,30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;
    
    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_strassen_block_jiang_mat_mul_open(n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void test_strassen_block_jiang_mat_mul(size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();

    matrix_block encrypted_A(A,scale,encoder,encryptor);
    matrix_block encrypted_B(B,scale,encoder,encryptor);
    matrix_block encrypted_C;
    encrypted_C = Strassen_block_jiang(encrypted_A,encrypted_B,scale,encoder,evaluator,gal_keys,relin_keys);
    
   
    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- Jiang strassen block matrix multiplication --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result

    Matrix<double> computed_C;
    encrypted_C.dec_block_matrix_jiang(computed_C,encoder,decryptor);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "---------------------------------------------------------------------------------" << endl;
    cout << endl <<endl;
}

void test_strassen_block_jiang_mat_mul_open(size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();

    matrix_block encrypted_A(A,scale,encoder,encryptor);
    matrix_block encrypted_B(B,scale,encoder,encryptor);
    matrix_block encrypted_C;
    encrypted_C = Strassen_block_jiang_open(encrypted_A,encrypted_B,scale,encoder,evaluator,gal_keys,relin_keys);
    
   
    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- Jiang strassen block matrix multiplication --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result

    Matrix<double> computed_C;
    encrypted_C.dec_block_matrix_jiang(computed_C,encoder,decryptor);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "---------------------------------------------------------------------------------" << endl;
    cout << endl <<endl;
}


void run_strassen_block_coprime_mat_mul(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    GaloisKeys gal_keys;
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    cout << "Please choose the method: 1 for our method, 2 for chen et al. method" << endl;
    int method_choice;
    cin >> method_choice;
    if (method_choice == 1)
    {
        // xyy add: seg-hoisting
        vector<int> step_s;
        size_t size_m = (log2(poly_modulus_degree/2))/2 - 1;
        // cout << "size_m = " << size_m << endl;
        // for (size_t i = 1; i < size_t(sqrt(m)) + 1; i++) 
        for (size_t i = 1; i < size_m + 1; i++)
        {
            step_s.push_back(mod(i * s_n, s_n*s_m));
            step_s.push_back(mod(i * s_p, s_m*s_p));
            // cout << mod(i * n, n*m) << ", ";
            // cout << mod(i * p, m*p) << ", ";
        }
        cout << endl;
        keygen.create_galois_keys(step_s, gal_keys);
    }
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;
    if (method_choice == 2)
    {
        // chen: native method
        keygen.create_galois_keys(gal_keys);
    }

    test_strassen_block_coprime_mat_mul(n,m,p,s_n,s_m,s_p,encoder,encryptor,scale,evaluator,relin_keys,gal_keys,decryptor, method_choice);
}

void run_strassen_block_coprime_mat_mul_open(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    GaloisKeys gal_keys;
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    cout << "Please choose the method: 1 for our method, 2 for chen et al. method" << endl;
    int method_choice;
    cin >> method_choice;
    if (method_choice == 1)
    {
        // xyy add: seg-hoisting
        vector<int> step_s;
        size_t size_m = (log2(poly_modulus_degree/2))/2 - 1;
        // cout << "size_m = " << size_m << endl;
        // for (size_t i = 1; i < size_t(sqrt(m)) + 1; i++) 
        for (size_t i = 1; i < size_m + 1; i++)
        {
            step_s.push_back(mod(i * s_n, s_n*s_m));
            step_s.push_back(mod(i * s_p, s_m*s_p));
            // cout << mod(i * n, n*m) << ", ";
            // cout << mod(i * p, m*p) << ", ";
        }
        cout << endl;
        keygen.create_galois_keys(step_s, gal_keys);
    }
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;
    if (method_choice == 2)
    {
        // chen: native method
        keygen.create_galois_keys(gal_keys);
    }

    test_strassen_block_coprime_mat_mul(n,m,p,s_n,s_m,s_p,encoder,encryptor,scale,evaluator,relin_keys,gal_keys,decryptor, method_choice);
}


void run_strassen_block_coprime_mat_mul_I(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 40, 60}));
    double scale = pow(2.0, 40);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;
    
    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_strassen_block_coprime_mat_mul(n,m,p,s_n,s_m,s_p,encoder,encryptor,scale,evaluator,relin_keys,gal_keys,decryptor,2);
}

void test_strassen_block_coprime_mat_mul(size_t n, size_t m, size_t p, size_t s_n, size_t s_m, size_t s_p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor, int method_choice)

{
    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();

    matrix_block encrypted_A(A,s_n,s_m,ceil(double(s_p)/s_m)+1,scale,encoder,encryptor);
    matrix_block encrypted_B(B,s_m,s_p,ceil(double(s_n)/s_m)+1,scale,encoder,encryptor);
    matrix_block encrypted_C;
    encrypted_C = Strassen_block_coprime(s_n,s_m,s_p,encrypted_A,encrypted_B,scale,encoder,evaluator,gal_keys,relin_keys, method_choice);
    
   
    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "-----------------Diag strassen block matrix multiplication --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result

    Matrix<double> computed_C;
    encrypted_C.dec_block_matrix_coprime(computed_C,s_n,s_p,encoder,decryptor);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "-------------------------------------------------------------------------------" << endl;
    cout << endl <<endl;
}

void test_strassen_block_coprime_mat_mul_open(size_t n, size_t m, size_t p, size_t s_n, size_t s_m, size_t s_p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor, int method_choice)

{
    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();

    matrix_block encrypted_A(A,s_n,s_m,ceil(double(s_p)/s_m)+1,scale,encoder,encryptor);
    matrix_block encrypted_B(B,s_m,s_p,ceil(double(s_n)/s_m)+1,scale,encoder,encryptor);
    matrix_block encrypted_C;
    encrypted_C = Strassen_block_coprime_open(s_n,s_m,s_p,encrypted_A,encrypted_B,scale,encoder,evaluator,gal_keys,relin_keys, method_choice);
    
   
    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "-----------------Diag strassen block matrix multiplication --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result

    Matrix<double> computed_C;
    encrypted_C.dec_block_matrix_coprime(computed_C,s_n,s_p,encoder,decryptor);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "-------------------------------------------------------------------------------" << endl;
    cout << endl <<endl;
}

void run_naive_block_jiang_mat_mul(size_t n, size_t m, size_t p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30,30,30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;
    
    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_naive_block_jiang_mat_mul(n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void test_naive_block_jiang_mat_mul(size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();

    matrix_block encrypted_A(A,scale,encoder,encryptor);
    matrix_block encrypted_B(B,scale,encoder,encryptor);
    matrix_block encrypted_C;
    encrypted_C = Naive_block_jiang(encrypted_A,encrypted_B,scale,encoder,evaluator,gal_keys,relin_keys);
    
   
    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- Jiang naive block matrix multiplication --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result

    Matrix<double> computed_C;
    encrypted_C.dec_block_matrix_jiang(computed_C,encoder,decryptor);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "------------------------------------------------------------------------------" << endl;
    cout << endl <<endl;
}

void run_naive_block_coprime_mat_mul(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    GaloisKeys gal_keys;
    cout << "Please choose the method: 1 for our method, 2 for chen et al. method" << endl;
    int method_choice;
    cin >> method_choice;
    if (method_choice == 1)
    {
        // xyy add: seg-hoisting
        vector<int> step_s;
        size_t size_m = (log2(poly_modulus_degree/2))/2 - 1;
        // cout << "size_m = " << size_m << endl;
        // for (size_t i = 1; i < size_t(sqrt(m)) + 1; i++) 
        for (size_t i = 1; i < size_m + 1; i++)
        {
            step_s.push_back(mod(i * s_n, s_n*s_m));
            step_s.push_back(mod(i * s_p, s_m*s_p));
        }
        cout << endl;
        keygen.create_galois_keys(step_s, gal_keys);
    }
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;
    if (method_choice == 2)
    {
        // chen: native method
        keygen.create_galois_keys(gal_keys);
    }

    test_naive_block_coprime_mat_mul(n, m, p,s_n,s_m,s_p,encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor, method_choice);
}

void test_naive_block_coprime_mat_mul(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor, int method_choice)
{
    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();

    matrix_block encrypted_A(A,s_n,s_m,ceil(double(s_p)/s_m)+1,scale,encoder,encryptor);
    matrix_block encrypted_B(B,s_m,s_p,ceil(double(s_n)/s_m)+1,scale,encoder,encryptor);
    matrix_block encrypted_C;
    encrypted_C = Naive_block_coprime(s_n,s_m,s_p ,encrypted_A,encrypted_B,scale,encoder,evaluator,gal_keys,relin_keys,method_choice);
    
   
    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- Diag naive block matrix multiplication --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result

    Matrix<double> computed_C;
    encrypted_C.dec_block_matrix_coprime(computed_C,s_n,s_p,encoder,decryptor);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "-----------------------------------------------------------------------------" << endl;
    cout << endl <<endl;
}

void run_coprime_mat_mul_imporve_log(size_t n, size_t m, size_t p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30,30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_coprime_mat_mul_imporve_log(n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void run_coprime_mat_mul_imporve_log_I(size_t n, size_t m, size_t p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_coprime_mat_mul_imporve_log(n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void test_coprime_mat_mul_imporve_log(size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();
    if(n*m*p>encoder.slot_count()){
        cerr<<"Not support O(logm) algorithm."<<endl;
        exit(1);
    }

    cipher_matrix_coprime cipher_A(A,p,scale,encoder,encryptor);
    cipher_matrix_coprime cipher_B(B,n,scale,encoder,encryptor);
    cipher_matrix_coprime cipher_C;

    encrypted_coprime_matrix_multiplication_imporve_log(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A,cipher_B,cipher_C);

    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- Diag encode matrix multiplication with O(logm) --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result

    Matrix<double> computed_C;
    cipher_C.dec_matrix_cipher(computed_C,encoder,decryptor);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "------------------------------------------------------------------------" << endl;
    cout << endl <<endl;
}


void run_naive_block_coprime_mat_mul_imporve_log(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30, 30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;
    
    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_naive_block_coprime_mat_mul_imporve_log(n, m, p,s_n,s_m,s_p,encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void run_naive_block_coprime_mat_mul_imporve_log_I(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;
    
    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_naive_block_coprime_mat_mul_imporve_log(n, m, p,s_n,s_m,s_p,encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void test_naive_block_coprime_mat_mul_imporve_log(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();

    matrix_block encrypted_A(A,s_n,s_m,s_p,scale,encoder,encryptor);
    matrix_block encrypted_B(B,s_m,s_p,s_n,scale,encoder,encryptor);
    matrix_block encrypted_C;
    encrypted_C = Naive_block_coprime_imporve_log(s_n,s_m,s_p ,encrypted_A,encrypted_B,scale,encoder,evaluator,gal_keys,relin_keys);
    
   
    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- Diag naive block faster matrix multiplication --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result

    Matrix<double> computed_C;
    encrypted_C.dec_block_matrix_coprime(computed_C,s_n,s_p,encoder,decryptor);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "-----------------------------------------------------------------------------" << endl;
    cout << endl <<endl;
}


void run_strassen_block_coprime_mat_mul_imporve_log(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30, 30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;
    
    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_strassen_block_coprime_mat_mul_imporve_log(n,m,p,s_n,s_m,s_p,encoder,encryptor,scale,evaluator,relin_keys,gal_keys,decryptor);
}

void run_strassen_block_coprime_mat_mul_imporve_log_I(size_t n, size_t m, size_t p,size_t s_n,size_t s_m,size_t s_p, size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30, 60}));
    double scale = pow(2.0, 30);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;
    
    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);

    test_strassen_block_coprime_mat_mul_imporve_log(n,m,p,s_n,s_m,s_p,encoder,encryptor,scale,evaluator,relin_keys,gal_keys,decryptor);
}

void test_strassen_block_coprime_mat_mul_imporve_log(size_t n, size_t m, size_t p, size_t s_n, size_t s_m, size_t s_p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    Matrix<double> A(n, m);
    Matrix<double> B(m, p);

    random_matrix_generator(n, m, A);
    random_matrix_generator(m, p, B);

    // cout << "the " << n << " * " << m << " matrix A: " << endl;
    // A.print();
    // cout << "the " << m << " * " << p << " matrix B: " << endl;
    // B.print(6,6);
    Matrix<double> C;
    C = A.multiply(B);
    // cout << "the resulting matrix C = A * B: " << endl;
    // C.print();

    auto start_time = chrono::high_resolution_clock::now();

    matrix_block encrypted_A(A,s_n,s_m,s_p,scale,encoder,encryptor);
    matrix_block encrypted_B(B,s_m,s_p,s_n,scale,encoder,encryptor);
    matrix_block encrypted_C;
    encrypted_C = Strassen_block_coprime_imporve_log(s_n,s_m,s_p,encrypted_A,encrypted_B,scale,encoder,evaluator,gal_keys,relin_keys);
    
   
    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "-----------------Diag strassen block faster matrix multiplication --------------------" << endl;
    cout << endl
         << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    // decrypte result

    Matrix<double> computed_C;
    encrypted_C.dec_block_matrix_coprime(computed_C,s_n,s_p,encoder,decryptor);

    // cout << "the computed C = A * B: " << endl;
    // computed_C.print();

    // cout << "the true C = A * B: " << endl;
    // C.print();

    computed_C.scaling_inplace(-1.);
    computed_C = computed_C.add(C);
    computed_C.is_zero(-2);
    cout << "-------------------------------------------------------------------------------" << endl;
    cout << endl <<endl;
}

