#!/bin/bash

# Main build script

set -e

# Configuration
BUILD_TYPE=${1:-Release}
BUILD_DIR="build"
CONAN_PROFILE=${2:-default}

echo "Building AeroacousticBEM in $BUILD_TYPE mode..."

# Create build directory
mkdir -p $BUILD_DIR
cd $BUILD_DIR

# Install dependencies with Conan
echo "Installing dependencies..."
conan install .. \
  --build=missing \
  --profile=$CONAN_PROFILE \
  -s build_type=$BUILD_TYPE \
  -o enable_testing=True \
  -o enable_examples=True

# Configure with CMake
echo "Configuring CMake..."
cmake .. \
  -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake \
  -DCMAKE_BUILD_TYPE=$BUILD_TYPE

# Build
echo "Building..."
cmake --build . --config $BUILD_TYPE -j$(nproc)

echo "Build completed successfully!"
