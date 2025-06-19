cmake_minimum_required(VERSION 3.20)

project(AeroacousticBEM
    VERSION 1.0.0
    DESCRIPTION "Boundary Element Method solver for aeroacoustic problems"
    LANGUAGES CXX
)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler-specific options
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-Wall -Wextra -Wpedantic -Werror)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-g -O0 -fsanitize=address,undefined)
    add_link_options(-fsanitize=address,undefined)
  else()
    add_compile_options(-O3 -DNDEBUG -march=native)
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  add_compile_options(/W4 /WX)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(/Od /Zi)
  else()
    add_compile_options(/O2 /DNDEBUG)
  endif()
endif()

# Find Conan dependencies
find_package(Eigen3 CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(GTest CONFIG REQUIRED)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/foundation
    ${CMAKE_CURRENT_SOURCE_DIR}/src/infrastructure
    ${CMAKE_CURRENT_SOURCE_DIR}/src/domain
    ${CMAKE_CURRENT_SOURCE_DIR}/src/application
)

# Foundation library
file(GLOB_RECURSE FOUNDATION_SOURCES "src/foundation/*.cpp")
add_library(bem_foundation STATIC ${FOUNDATION_SOURCES})
target_link_libraries(bem_foundation PUBLIC
    fmt::fmt
    spdlog::spdlog
    Eigen3::Eigen
)

# Infrastructure library
# file(GLOB_RECURSE INFRASTRUCTURE_SOURCES "src/infrastructure/*.cpp")
# add_library(bem_infrastructure STATIC ${INFRASTRUCTURE_SOURCES})
# target_link_libraries(bem_infrastructure PUBLIC
#     bem_foundation
#     Eigen3::Eigen
#     spdlog::spdlog
# )

# Domain library
# file(GLOB_RECURSE DOMAIN_SOURCES "src/domain/*.cpp")
# add_library(bem_domain STATIC ${DOMAIN_SOURCES})
# target_link_libraries(bem_domain PUBLIC
#     bem_foundation
#     bem_infrastructure
#     Eigen3::Eigen
# )

# Application library
# file(GLOB_RECURSE APPLICATION_SOURCES "src/application/*.cpp")
# add_library(bem_application STATIC ${APPLICATION_SOURCES})
# target_link_libraries(bem_application PUBLIC
#     bem_domain
#     spdlog::spdlog
# )

# Main executable
add_executable(aeroacoustic_bem src/main.cpp)
# target_link_libraries(aeroacoustic_bem PRIVATE
#     bem_application
# )

# Enable testing
enable_testing()
add_subdirectory(tests)

# Install targets
# install(TARGETS aeroacoustic_bem
#     RUNTIME DESTINATION bin
# )

# Copy configuration files
# install(FILES
#     ${CMAKE_CURRENT_SOURCE_DIR}/examples/
#     DESTINATION share/aeroacoustic_bem/examples
# )

# Package configuration
include(CPack)
set(CPACK_PACKAGE_NAME "AeroacousticBEM")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
