#!/bin/bash

set -e

# Help support
if [[ "$1" == "--help" ]]; then
  echo "Usage: ./build.sh [BuildType] [ConanProfile] [clean]"
  echo "Defaults: BuildType=Release, Profile=default"
  exit 0
fi

# Default to Release if no argument is provided
BUILD_TYPE=${1:-Release}
CONAN_PROFILE=${2:-default}
BUILD_DIR="build/$BUILD_TYPE"

# Clean up previous build if it exists
if [[ "$1" == "clean" ]]; then
  echo "🧹 Cleaning build directory..."
  rm -rf "$BUILD_DIR"
  echo "✅ Build directory was removed successfully!"
  exit 0
fi

echo "🔧 Building AeroacousticBEM in $BUILD_TYPE mode..."

# Only configure once (unless CMakeLists.txt or dependencies change)
if [[ ! -f "$BUILD_DIR/CMakeCache.txt" ||
  CMakeLists.txt -nt "$BUILD_DIR/CMakeCache.txt" ||
  conanfile.* -nt "$BUILD_DIR/CMakeCache.txt" ]]; then

  # Install dependencies with Conan from the project root
  echo "📦 Installing dependencies with Conan..."
  conan install . \
    --build=missing \
    --profile="$CONAN_PROFILE" \
    -s build_type="$BUILD_TYPE" \
    --output-folder=.

  # Configure with CMake - pointing back to the source directory
  echo "🛠️ Configuring with CMake..."
  cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_TOOLCHAIN_FILE="$BUILD_DIR/generators/conan_toolchain.cmake" \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
fi

# Build (incremental!)
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
  # ctest --build-config "$BUILD_TYPE" --output-on-failure
  ctest --build-config "$BUILD_TYPE" --output-on-failure -V --timeout 10
  popd >/dev/null
fi
