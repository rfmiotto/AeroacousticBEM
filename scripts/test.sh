#!/bin/bash
# Run tests

set -e

BUILD_DIR="build"
TEST_RESULTS_DIR="test_results"

if [ ! -d "$BUILD_DIR" ]; then
  echo "Build directory not found. Please run build.sh first."
  exit 1
fi

cd $BUILD_DIR

echo "Running unit tests..."
mkdir -p $TEST_RESULTS_DIR

# Run tests with verbose output
ctest --verbose --output-on-failure

# Generate coverage report if available
if command -v gcov &>/dev/null && [ "$1" = "--coverage" ]; then
  echo "Generating coverage report..."
  gcov -r ../src/**/*.cpp

  if command -v lcov &>/dev/null; then
    lcov --capture --directory . --output-file coverage.info
    lcov --remove coverage.info '/usr/*' --output-file coverage.info
    lcov --remove coverage.info '*/tests/*' --output-file coverage.info

    if command -v genhtml &>/dev/null; then
      genhtml coverage.info --output-directory coverage_html
      echo "Coverage report generated in coverage_html/"
    fi
  fi
fi

echo "Tests completed!"
