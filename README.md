# Homomorphic Matrix Operations under Bicyclic Encoding

## Compile and Run (Ubuntu 22.04.2 LTS)

### Dependencies

- cmake 3.13 or higher
- clang
- [Microsoft SEAL](https://github.com/microsoft/seal)

### Compile

Download the source code from GitHub, unzip and enter the repository root directory, and then run the following:

    cd bicyclic_mat_mul
    mkdir build
    cd build
    cmake ../test/
    make

### How to use 

- Open a terminal, enter the `build` directory and run 

      ./test_enc_mat_mul


### Test Instructions

Once the file is executed, it will output the following information:

```bash
Examples:

  1.  Bicyclic matrix multiplication base
  2.  Faster Bicyclic matrix multiplication base
  3.  Jiang matrix multiplication base
  4.  Lu matrix multiplication base
  5.  R-T matrix multiplication base
  6.  Block matrix multiplication by LongRot action
  7.  Block matrix multiplication strassen and naive(jiang)
  8.  Block matrix multiplication strassen and naive(Bicyclic)
  9.  Base matrix multiplication in Table 9
 10.  Base matrix multiplication in Table 10
 11.  Block matrix multiplication in Table 11(128,128,128)
 12.  Block matrix multiplication in Table 12(256,256,256)
 13.  Block matrix multiplication in Table 13(512,512,512)
 14.  Block matrix multiplication in Table 14
 15.  Block matrix multiplication in Table 15
 16.  Block matrix multiplication in Table 16
 0. Exit
```

In this case, `case 1-7` refer to the basic test cases, while `case 8-14` include all the scenarios from Table 8 to Table 14. Here, $(n,m,p)$ represents matrices $A\in \mathcal{R}^{n \times m}$ and $B\in \mathcal{R}^{m \times p}$. The specific test scenarios are as follows:

1. `case 1`:

   - `Bicyclic encode matrix multiplication`: $(44, 45, 43)$
   - `Bicyclic encode matrix multiplication`: $(61, 64, 63)$
   - `Bicyclic encode matrix multiplication`: $(29, 128, 31)$
   - `Bicyclic encode matrix multiplication`: $(89, 91, 90)$

2. `case 2`:
   * `Faster Bicyclic encode matrix multiplication`: $(15, 16, 17)$
   * `Faster Bicyclic encode matrix multiplication`: $(21, 16, 23)$
   * `Faster Bicyclic encode matrix multiplication`: $(31, 16, 33)$
   
3. `case 3`:

   - `Jiang matrix multiplication`: $(64,64,64)$
   - `Jiang matrix multiplication`: $(128,128,128)$ | N=32768

4. `case 4`:

   - `Lu matrix multiplication segment 8`: $(16,16,5 \times 4096)$
   - `Lu matrix multiplication not segment`: $(16,16,5 \times 4096)$
   - `Lu matrix multiplication segment 8`: $(64,64,5 \times 4096)$
   - `Lu matrix multiplication not segment`: $(64,64,5 \times 4096)$

5. `case 5`:

   - `R-T matrix multiplication`: $(16,16,16)$

6. `case 6`:

   - `Bicyclic LongRot matrix multiplication with pre-generate`: $(256, 257, 17)$
   - `Bicyclic LongRot matrix multiplication with pre-generate`: $(256, 17, 257)$
   - `Bicyclic LongRot matrix multiplication`: $(256, 257, 17)$
   - `Bicyclic LongRot matrix multiplication`: $(256, 17, 257)$

7. `case 7`:

   - `Jiang naive block matrix multiplication`: $(64,64,64) \times 4$
   - `Jiang strassen block matrix multiplication`: $(64,64,64) \times 4$
   - `Jiang naive block matrix multiplication`: $(64,64,64) \times 8$
   - `Jiang strassen block matrix multiplication`: $(64,64,64) \times 8$

8. `case 7`:

   - `Bicyclic encode naive block matrix multiplication`: $(44, 45, 43) \times 4$
   - `Bicyclic encode strassen block matrix multiplication`: $(44, 45, 43) \times 4$
   - `Bicyclic encode naive block matrix multiplication`: $(44, 45, 43) \times 8$
   - `Bicyclic encode strassen block matrix multiplication`: $(44, 45, 43) \times 8$

9. `case 9`:
   - `R-T matrix multiplication`: $(16,16,16)$
   - `Bicyclic encode matrix multiplication`: $(16, 19, 17)  | \log q=170$
   - `Bicyclic encode matrix multiplication`: $(16, 19, 17)  | \log q=140$
   - `Faster Bicyclic encode matrix multiplication`: $(15, 16, 17)$
   
10. `case 10`:
   - `Jiang matrix multiplication`: $(64,64,64)$
   - `Jiang matrix multiplication`: $(128,128,128) | N=32768$
   - `Bicyclic encode matrix multiplication`: $(44, 45, 43)$
   - `Bicyclic encode matrix multiplication`: $(61, 64, 63) | N=16384$
   - `Bicyclic encode matrix multiplication`: $(89, 91, 90) | N=32768$
   - `Faster Bicyclic encode matrix multiplication`: $(15, 16, 17)$
   - `Faster Bicyclic encode matrix multiplication`: $(21, 16, 23)$
   - `Faster Bicyclic encode matrix multiplication`: $(31, 16, 33)$

11. `case 11`:
    * `Jiang naive block matrix multiplication`: $(64,64,64) \times 2$

    * `Jiang naive block matrix multiplication`: $(128,128,128) \times 1 |N=32768$

    * `Jiang strassen block matrix multiplication`: $(64,64,64) \times 2$

      

    * `Bicyclic encode naive block matrix multiplication`: $(43, 45, 44) \times 3$

    * `Bicyclic encode naive block matrix multiplication`: $(64, 67, 65) \times 2 | N=32768$

    * `Bicyclic encode strassen block matrix multiplication`: $(32, 35, 33) \times 4$

    * `Bicyclic encode strassen block matrix multiplication`: $(64, 67, 65) \times 2 | N=32768$

      

    * `Faster Bicyclic encode naive block matrix multiplication`: $(9\times 15, 8\times 16, 8\times 17)$ 

    * `Faster Bicyclic encode naive block matrix multiplication`: $(7\times 21, 4\times 32, 6\times 23 ) | N=32768$

    * `Faster Bicyclic encode strassen block matrix multiplication`: $(11, 8, 9) \times 16$

    * `Faster Bicyclic encode strassen block matrix multiplication`: $(19, 16, 17) \times 8 | N=32768$

      

    * `Bicyclic LongRot matrix multiplication `: $(128, 131, 129)$

    * `Bicyclic LongRot matrix multiplication`: $(128, 131, 129) | N=32768$

12. `case 12`:

    - `Jiang naive block matrix multiplication`: $(64,64,64) \times 4$

    - `Jiang naive block matrix multiplication`: $(128,128,128) \times 2 | N=32768$

    - `Jiang strassen block matrix multiplication`: $(64,64,64) \times 4$

    - `Jiang strassen block matrix multiplication`:  $(128,128,128) \times 2 | N=32768$

      

    - `Bicyclic encode naive block matrix multiplication`: $(44, 45, 43)\times 6$

    - `Bicyclic encode naive block matrix multiplication`: $(88, 89, 87) \times 3 | N=32768$

    - `Bicyclic encode strassen block matrix multiplication`: $(44, 45, 43) \times 8$

    - `Bicyclic encode strassen block matrix multiplication`: $(64, 65, 67) \times 4 | N=32768$

      

    - `Faster Bicyclic encode naive block matrix multiplication`: $(18\times 15, 16\times 16, 16\times 17)$ 

    - `Faster Bicyclic encode naive block matrix multiplication`: $(13\times 21, 8\times 32, 12\times 23 ) | N=32768$

      

    - `Bicyclic LongRot matrix multiplication with pre-generate`: $(256, 259, 257)$

    - `Bicyclic LongRot matrix multiplication`: $(256, 259, 257)| N=32768$

13. `case 13`:

    - `Jiang naive block matrix multiplication`: $(64,64,64) \times 8$

    - `Bicyclic encode naive block matrix multiplication`: $(44, 45, 43) \times 12$

    - `Jiang strassen block matrix multiplication`: $(64,64,64) \times 8$

    - `Bicyclic encode strassen block matrix multiplication`: $(44, 45, 43) \times 16$

    - `Bicyclic LongRot matrix multiplication with pre-generate`: $(512, 515, 513)$

      

    - `Jiang naive block matrix multiplication`: $(64,64,64) \times 16$

    - `Bicyclic encode naive block matrix multiplication`: $(44, 45, 43) \times 24$

    - `Jiang strassen block matrix multiplication`: $(64,64,64) \times 16$

    - `Bicyclic encode strassen block matrix multiplication`: $(32, 35, 33) \times 32 | \log q=150$

    - `Bicyclic LongRot matrix multiplication with pre-generate`: $(1024, 1027, 1025)$

14. `case 14`:

    - `Bicyclic LongRot matrix multiplication with pre-generate`: $(256, 257, 17)$
    - `Bicyclic LongRot matrix multiplication with pre-generate`: $(256, 17, 257)$
    - `Bicyclic LongRot matrix multiplication with pre-generate`: $(1024, 1025, 17)$
    - `Bicyclic LongRot matrix multiplication with pre-generate`: $(1024, 17, 1025)$
    - `Bicyclic LongRot matrix multiplication with pre-generate`: $(2048, 2049, 11)$
    - `Bicyclic LongRot matrix multiplication with pre-generate`: $(2049, 8, 2051)$

15. `case 15`:

    - `Jiang naive block matrix multiplication`: $(4,1636,5)$

    - `Jiang naive block matrix multiplication`: $(8,3045,9)$

    - `Jiang naive block matrix multiplication`: $(16,6093,17)$

    - `Jiang naive block matrix multiplication`: $(32,13847,33)$

      

    - `Bicyclic encode naive block matrix multiplication`: $(4,1636,5)$

    - `Bicyclic encode naive block matrix multiplication`: $(8,3045,9)$

    - `Bicyclic encode naivee block matrix multiplication`: $(16,6093,17)$

    - `Bicyclic encode naive block matrix multiplication`: $(32,13847,33)$

16. `case 16`:

    - `Jiang naive block matrix multiplication`: $(4,5,1636)$

    - `Jiang naive block matrix multiplication`: $(8,9,3405)$

    - `Jiang naive block matrix multiplication`: $(16,17,6093)$

    - `Jiang naive block matrix multiplication`: $(32,33,13847)$

      

    - `Bicyclic encode naive block matrix multiplication`: $(4,5,1636)$

    - `Bicyclic encode naive block matrix multiplication`: $(8,9,3405)$

    - `Bicyclic encode naivee block matrix multiplication`: $(16,17,6093)$

    - `Bicyclic encode naive block matrix multiplication`: $(32,33,13847)$

      

    - `Lu matrix multiplication segment 8`:  $(4,5,1636)$

    - `Lu matrix multiplication segment 8`:  $(8,9,3405)$

    - `Lu matrix multiplication segment 8`:  $(16,17,6093)$

    - `Lu matrix multiplication segment 8`:  $(32,33,13847)$
