# Install script for directory: /home/delld/xyy/seall/native/src/seal

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SEAL-4.1/seal" TYPE FILE FILES
    "/home/delld/xyy/seall/native/src/seal/batchencoder.h"
    "/home/delld/xyy/seall/native/src/seal/ciphertext.h"
    "/home/delld/xyy/seall/native/src/seal/ckks.h"
    "/home/delld/xyy/seall/native/src/seal/modulus.h"
    "/home/delld/xyy/seall/native/src/seal/context.h"
    "/home/delld/xyy/seall/native/src/seal/decryptor.h"
    "/home/delld/xyy/seall/native/src/seal/dynarray.h"
    "/home/delld/xyy/seall/native/src/seal/encryptionparams.h"
    "/home/delld/xyy/seall/native/src/seal/encryptor.h"
    "/home/delld/xyy/seall/native/src/seal/evaluator.h"
    "/home/delld/xyy/seall/native/src/seal/galoiskeys.h"
    "/home/delld/xyy/seall/native/src/seal/keygenerator.h"
    "/home/delld/xyy/seall/native/src/seal/kswitchkeys.h"
    "/home/delld/xyy/seall/native/src/seal/memorymanager.h"
    "/home/delld/xyy/seall/native/src/seal/plaintext.h"
    "/home/delld/xyy/seall/native/src/seal/publickey.h"
    "/home/delld/xyy/seall/native/src/seal/randomgen.h"
    "/home/delld/xyy/seall/native/src/seal/randomtostd.h"
    "/home/delld/xyy/seall/native/src/seal/relinkeys.h"
    "/home/delld/xyy/seall/native/src/seal/seal.h"
    "/home/delld/xyy/seall/native/src/seal/secretkey.h"
    "/home/delld/xyy/seall/native/src/seal/serializable.h"
    "/home/delld/xyy/seall/native/src/seal/serialization.h"
    "/home/delld/xyy/seall/native/src/seal/valcheck.h"
    "/home/delld/xyy/seall/native/src/seal/version.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/delld/xyy/seall/build/native/src/seal/util/cmake_install.cmake")

endif()

