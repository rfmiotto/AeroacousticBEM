#!/bin/bash
# Format code using clang-format

set -e

# Find all C++ source files
find src tests examples -name "*.hpp" -o -name "*.cpp" | while read -r file; do
  echo "Formatting $file"
  clang-format -i "$file"
done

echo "Code formatting completed!"
