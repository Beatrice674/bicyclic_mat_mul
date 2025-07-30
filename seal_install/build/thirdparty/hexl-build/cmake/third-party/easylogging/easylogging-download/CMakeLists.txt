# Copyright (C) 2020 Intel Corporation
# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.13)

project(easylogging-download NONE)

include(ExternalProject)
ExternalProject_Add(easylogging
  GIT_REPOSITORY    https://github.com/amrayn/easyloggingpp.git
  GIT_TAG           8489989 # master as of Oct 26, 2021
  SOURCE_DIR        "/home/delld/xyy/seall/build/thirdparty/hexl-build/cmake/third-party/easylogging/easylogging-src"
  BINARY_DIR        "/home/delld/xyy/seall/build/thirdparty/hexl-build/cmake/third-party/easylogging/easylogging-build"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
  UPDATE_COMMAND    ""
)
