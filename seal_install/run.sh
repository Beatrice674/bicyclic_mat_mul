#!/bin/bash

# 顶层构建
# 删除 build 文件夹
if [ -d "build" ]; then
    echo "检测到 build 文件夹，正在删除..."
    rm -rf build
else
    echo "未检测到 build 文件夹，跳过删除步骤。"
fi  # ✅ 正确结束 if
# cmake -S . -B build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
sudo cmake --install build

# 进入 native/examples 目录
cd native/examples

# 删除 build 文件夹
if [ -d "build" ]; then
    echo "检测到 build 文件夹，正在删除..."
    rm -rf build
else
    echo "未检测到 build 文件夹，跳过删除步骤。"
fi  # ✅ 正确结束 if

# 重新生成构建系统
echo "正在运行 cmake -S . -B build ..."
# cmake -S . -B build
# 开启debug模式
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug

# 关闭debug模式
# cmake -S . -B build -DCMAKE_BUILD_TYPE=Release


# 编译项目
echo "正在编译项目..."
cmake --build build

# 进入执行目录并运行示例
cd build/bin
./sealexamples

# 为脚本添加执行权限
# chmod +x run.sh
