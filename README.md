# Improved Homomorphic Matrix Operations under Bicyclic Encoding

### Dependencies

- cmake 3.13 or higher
- clang
- [Microsoft SEAL](https://github.com/microsoft/seal)

### Compile the code

Download the source code from GitHub, unzip and enter the repository root directory, and then run the following:

```
cd bicyclic_mat_mul
```

- Install the SEAL library

```
cd seal_install
cmake -S . -B build 
cmake --build build
sudo cmake --install build
cd native/examples
cmake -S . -B build
cmake --build build
```

- Execute the experimental code

```
mkdir build
cd build
cmake ../test/
make
./test_enc_mat_mul
```


### Test Instructions

Once the file is executed, it will output the following information:

```bash
Examples:

  1.  Bicyclic matrix multiplication for small dimensions
  2.  Bicyclic matrix multiplication for large dimensions
  3.  OpenMP optimization
 	0. Exit
```

In this case, `case 1`  refer to matrix multiplication for small dimensions, `case 2`  refer to matrix multiplication for large dimensions, `case 3`  refer to OpenMP optimization for Strassen algorithm, respectively.Here $N$ is the polynomial degree and $(n,m,p)$ represents matrices $A\in \mathbb{R}^{n \times m}$ and $B\in \mathbb{R}^{m \times p}$. The specific test scenarios are as follows:

1. `case 1`:

   💡 $N$ = 8192

   - `Improved bicyclic encode matrix multiplication`: $(64, 31, 29)$
   - `Improved bicyclic encode matrix multiplication`: $(44, 45, 43)$
   - `Improved bicyclic encode matrix multiplication`: $(31, 64, 29)$

   💡 $N$ = 16384

   - `Improved bicyclic encode matrix multiplication`: $(128, 31, 29)$
   - `Improved bicyclic encode matrix multiplication`: $(61, 64, 63)$
   - `Improved bicyclic encode matrix multiplication`: $(29, 128, 31)$

   💡 $N$ = 32768

   - `Improved bicyclic encode matrix multiplication`: $(128, 63, 61)$

   - `Improved bicyclic encode matrix multiplication`: $(89, 91, 90)$

   - `Improved bicyclic encode matrix multiplication`: $(63, 128, 61)$

     

2. `case 2`:

   ▪️Matrix dimensions $(64, 64, 64)$

   💡 $N$ = 8192

   * `Jiang matrix multiplication`: $(64,64,64)$
   * `Bicyclic encode naive block matrix multiplication`: $(32, 35, 33) \times 2$ 
   * `Bicyclic encode strassen block matrix multiplication`: $(32, 35, 33) \times 2$ 

   💡 $N$ = 32768

   * `Improved bicyclic encode matrix multiplication`: $(64,67,65)$

   ▪️Matrix dimensions $(128, 128, 128)$

   💡 $N$ = 8192

   * `Jiang naive block matrix multiplication`: $(64,64,64) \times 2$
   * `Jiang strassen block matrix multiplication`: $(64,64,64) \times 2$
   * `Bicyclic encode naive block matrix multiplication`: $(43, 45, 44) \times 3$ 
   * `Bicyclic encode strassen block matrix multiplication`: $(32, 35, 33) \times 4$ 

   💡 $N$ = 32768

   * `Jiang matrix multiplication`: $(128,128,128)$ 
   * `Bicyclic encode naive block matrix multiplication`: $(64, 67, 65) \times 2 $
   * `Bicyclic encode strassen block matrix multiplication`: $(64, 67, 65) \times 2 $

   ▪️Matrix dimensions $(256, 256, 256)$

   💡 $N$ = 8192

   * `Jiang naive block matrix multiplication`: $(64,64,64) \times 4$
   * `Jiang strassen block matrix multiplication`: $(64,64,64) \times 4$
   * `Bicyclic encode naive block matrix multiplication`: $(43, 45, 44) \times 6$ 
   * `Bicyclic encode strassen block matrix multiplication`: $(32, 35, 33) \times 8$ 

   💡 $N$ = 32768

   * `Jiang naive block matrix multiplication`: $(128,128,128) \times 2$
   * `Jiang strassen block matrix multiplication`:  $(128,128,128) \times 2$
   * `Bicyclic encode naive block matrix multiplication`: $(88, 89, 87) \times 3 $
   * `Bicyclic encode strassen block matrix multiplication`: $(64, 67, 65) \times 4 $

   ▪️Matrix dimensions $(512, 512, 512)$

   💡 $N$ = 8192

   * `Jiang naive block matrix multiplication`: $(64,64,64) \times 8$
   * `Jiang strassen block matrix multiplication`: $(64,64,64) \times 8$
   * `Bicyclic encode naive block matrix multiplication`: $(43, 45, 44) \times 12$ 
   * `Bicyclic encode strassen block matrix multiplication`: $(32, 35, 33) \times 16$ 

   💡 $N$ = 32768

   * `Jiang naive block matrix multiplication`: $(128,128,128) \times 4$
   * `Jiang strassen block matrix multiplication`:  $(128,128,128) \times 4$
   * `Bicyclic encode naive block matrix multiplication`: $(88, 89, 87) \times 6 $
   * `Bicyclic encode strassen block matrix multiplication`: $(64, 67, 65) \times 8 $

   

3. `case 3`:

   ▪️Matrix dimensions $(64, 64, 64)$

   💡 $N$ = 8192

   * `Bicyclic encode strassen block matrix multiplication with OpenMP`: $(32, 35, 33) \times 2$ 

   ▪️Matrix dimensions $(128, 128, 128)$

   💡 $N$ = 8192

   * `Jiang strassen block matrix multiplication with OpenMP`: $(64,64,64) \times 2$
   * `Bicyclic encode strassen block matrix multiplication with OpenMP`: $(32, 35, 33) \times 4$ 

   💡 $N$ = 32768

   * `Bicyclic encode strassen block matrix multiplication with OpenMP`: $(64, 67, 65) \times 2 $

   ▪️Matrix dimensions $(256, 256, 256)$

   💡 $N$ = 8192

   * `Jiang strassen block matrix multiplication with OpenMP`: $(64,64,64) \times 4$
   * `Bicyclic encode strassen block matrix multiplication with OpenMP`: $(32, 35, 33) \times 8$ 

   💡 $N$ = 32768

   * `Jiang strassen block matrix multiplication with OpenMP`:  $(128,128,128) \times 2$
   * `Bicyclic encode strassen block matrix multiplication with OpenMP`: $(64, 67, 65) \times 4 $

   ▪️Matrix dimensions $(512, 512, 512)$

   💡 $N$ = 8192

   * `Jiang strassen block matrix multiplication with OpenMP`: $(64,64,64) \times 8$
   * `Bicyclic encode strassen block matrix multiplication with OpenMP`: $(32, 35, 33) \times 16$ 

   💡 $N$ = 32768

   * `Jiang strassen block matrix multiplication with OpenMP`:  $(128,128,128) \times 4$
   * `Bicyclic encode strassen block matrix multiplication with OpenMP`: $(64, 67, 65) \times 8 $

   
