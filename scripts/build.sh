#!/bin/bash

set -e

# Default to Release if no argument is provided
BUILD_TYPE=${1:-Release}
CONAN_PROFILE=${2:-default}
BUILD_DIR="build/$BUILD_TYPE"

echo "🔧 Building AeroacousticBEM in $BUILD_TYPE mode..."

# Clean up previous build if it exists
if [[ -d "build" ]]; then
  echo "🧹 Cleaning previous build directory..."
  rm -rf "build"
fi

# Install dependencies with Conan from the project root
echo "📦 Installing dependencies with Conan..."
conan install . \
  --build=missing \
  --profile="$CONAN_PROFILE" \
  -s build_type="$BUILD_TYPE" \
  --output-folder=.

# The cmake_layout() in conanfile.py will create the appropriate structure
# cd "$BUILD_DIR"

# Find toolchain
TOOLCHAIN_FILE=""
if [[ -f "conan_toolchain.cmake" ]]; then
  TOOLCHAIN_FILE="conan_toolchain.cmake"
elif [[ -f "build/$BUILD_TYPE/generators/conan_toolchain.cmake" ]]; then
  TOOLCHAIN_FILE="build/$BUILD_TYPE/generators/conan_toolchain.cmake"
elif [[ -f "generators/conan_toolchain.cmake" ]]; then
  TOOLCHAIN_FILE="generators/conan_toolchain.cmake"
fi

echo "✅ Found toolchain file: $TOOLCHAIN_FILE"
echo "📁 Working directory: $(pwd)"

# Configure with CMake - pointing back to the source directory
echo "🛠️ Configuring with CMake..."
cmake -S . -B "$BUILD_DIR" \
  -DCMAKE_TOOLCHAIN_FILE="$TOOLCHAIN_FILE" \
  -DCMAKE_BUILD_TYPE="$BUILD_TYPE"

# Build
echo "🚀 Building..."
cmake --build "$BUILD_DIR" --config "$BUILD_TYPE" -j$(nproc)

echo "✅ Build completed successfully!"

# Link compile_commands.json to project root for clangd
echo "🔗 Linking compile_commands.json to project root..."
ln -sf "$BUILD_DIR/compile_commands.json" ./compile_commands.json

# Optional: Run tests if they exist and were built
if [[ -d "tests" ]] || [[ -f "CTestTestfile.cmake" ]]; then
  echo "🧪 Running tests..."
  pushd "$BUILD_DIR" >/dev/null
  ctest --build-config "$BUILD_TYPE" --output-on-failure
  popd >/dev/null
fi
