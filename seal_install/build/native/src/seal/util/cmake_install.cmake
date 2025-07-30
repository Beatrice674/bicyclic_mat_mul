# Install script for directory: /home/delld/xyy/seall/native/src/seal/util

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SEAL-4.1/seal/util" TYPE FILE FILES
    "/home/delld/xyy/seall/native/src/seal/util/blake2.h"
    "/home/delld/xyy/seall/native/src/seal/util/blake2-impl.h"
    "/home/delld/xyy/seall/native/src/seal/util/clang.h"
    "/home/delld/xyy/seall/native/src/seal/util/clipnormal.h"
    "/home/delld/xyy/seall/native/src/seal/util/common.h"
    "/home/delld/xyy/seall/native/src/seal/util/croots.h"
    "/home/delld/xyy/seall/native/src/seal/util/defines.h"
    "/home/delld/xyy/seall/native/src/seal/util/dwthandler.h"
    "/home/delld/xyy/seall/native/src/seal/util/fips202.h"
    "/home/delld/xyy/seall/native/src/seal/util/galois.h"
    "/home/delld/xyy/seall/native/src/seal/util/gcc.h"
    "/home/delld/xyy/seall/native/src/seal/util/globals.h"
    "/home/delld/xyy/seall/native/src/seal/util/hash.h"
    "/home/delld/xyy/seall/native/src/seal/util/hestdparms.h"
    "/home/delld/xyy/seall/native/src/seal/util/iterator.h"
    "/home/delld/xyy/seall/native/src/seal/util/locks.h"
    "/home/delld/xyy/seall/native/src/seal/util/mempool.h"
    "/home/delld/xyy/seall/native/src/seal/util/msvc.h"
    "/home/delld/xyy/seall/native/src/seal/util/numth.h"
    "/home/delld/xyy/seall/native/src/seal/util/pointer.h"
    "/home/delld/xyy/seall/native/src/seal/util/polyarithsmallmod.h"
    "/home/delld/xyy/seall/native/src/seal/util/polycore.h"
    "/home/delld/xyy/seall/native/src/seal/util/rlwe.h"
    "/home/delld/xyy/seall/native/src/seal/util/rns.h"
    "/home/delld/xyy/seall/native/src/seal/util/scalingvariant.h"
    "/home/delld/xyy/seall/native/src/seal/util/ntt.h"
    "/home/delld/xyy/seall/native/src/seal/util/streambuf.h"
    "/home/delld/xyy/seall/native/src/seal/util/uintarith.h"
    "/home/delld/xyy/seall/native/src/seal/util/uintarithmod.h"
    "/home/delld/xyy/seall/native/src/seal/util/uintarithsmallmod.h"
    "/home/delld/xyy/seall/native/src/seal/util/uintcore.h"
    "/home/delld/xyy/seall/native/src/seal/util/ztools.h"
    )
endif()

