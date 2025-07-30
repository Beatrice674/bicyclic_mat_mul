cmake_minimum_required(VERSION 3.13)

## We must specify the C and C++ compilers before the project command
# set(CMAKE_C_COMPILER "/usr/bin/clang")
# set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_STANDARD 17)

project(test_bicyclic_mat_mul)

# Import Microsoft SEAL
find_package(SEAL 4.1.2 EXACT REQUIRED
              PATHS /usr/local) 

# Import OpenMP
FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


set(REAL_SOURCE_DIR ${PROJECT_SOURCE_DIR}/../src)
set(TEST_SOURCE_DIR ${PROJECT_SOURCE_DIR})

add_executable(test_enc_mat_mul
               ${TEST_SOURCE_DIR}/test_enc_mat_mul.cpp
               ${TEST_SOURCE_DIR}/../src/utils.cpp
               ${TEST_SOURCE_DIR}/../src/matrix.cpp
            #    ${TEST_SOURCE_DIR}/../src/enc_mat_mul.cpp
               ${TEST_SOURCE_DIR}/../src/long_enc_cipher.h
               ${TEST_SOURCE_DIR}/../src/jiang_enc_cipher.h
               ${TEST_SOURCE_DIR}/../src/strassen_block.h
               ${TEST_SOURCE_DIR}/../src/R_T_enc_cipher.h
               ${TEST_SOURCE_DIR}/../src/coprime_enc_cipher.h
               ${TEST_SOURCE_DIR}/../src/block_long_cipher.h
               )
                
set_target_properties(test_enc_mat_mul PROPERTIES LINK_FLAGS "-static")
                    
if(TARGET SEAL::seal)
    target_link_libraries(test_enc_mat_mul SEAL::seal)
elseif(TARGET SEAL::seal_shared)
    target_link_libraries(test_enc_mat_mul SEAL::seal_shared)
else()
    message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
endif()
