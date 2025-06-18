#!/bin/bash

set -e # Exit on error

# Folder names
DEPS_DIR="build-deps"
BUILD_DIR="build"

# Conan settings
BUILD_TYPE=${1:-Release}

echo ">>> Installing Conan dependencies in $BUILD_TYPE mode..."
conan install . \
  --output-folder=$DEPS_DIR \
  --build=missing \
  -s build_type=$BUILD_TYPE

echo ">>> Configuring CMake with Conan toolchain..."
cmake -S . -B $BUILD_DIR \
  -DCMAKE_TOOLCHAIN_FILE=$DEPS_DIR/conan_toolchain.cmake \
  -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

echo ">>> Building the project..."
cmake --build $BUILD_DIR

echo ">>> Linking compile_commands.json for clangd..."
ln -sf $BUILD_DIR/compile_commands.json .

echo "✅ Done!"
