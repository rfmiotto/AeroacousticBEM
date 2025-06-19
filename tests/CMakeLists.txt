# Find required packages
find_package(GTest REQUIRED)

# Test executable
file(GLOB_RECURSE TEST_SOURCES
    "unit/*.cpp"
    "integration/*.cpp"
)

add_executable(run_tests
    main_test.cpp
    ${TEST_SOURCES}
)

target_link_libraries(run_tests PRIVATE
    # bem_application
    # bem_domain
    # bem_infrastructure
    bem_foundation
    GTest::gtest
    GTest::gtest_main
    GTest::gmock
    Eigen3::Eigen
)

# Add tests to CTest
include(GoogleTest)
gtest_discover_tests(run_tests
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Coverage target (if available)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  option(ENABLE_COVERAGE "Enable coverage reporting" OFF)
  if(ENABLE_COVERAGE)
    target_compile_options(run_tests PRIVATE --coverage)
    target_link_options(run_tests PRIVATE --coverage)
  endif()
endif()
