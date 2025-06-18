from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps
from conan.tools.files import save, load
import os


class AeroacousticBEMConan(ConanFile):
    name = "aeroacoustic-bem"
    version = "1.0.0"
    description = "Boundary Element Method solver for aeroacoustic problems"
    author = "Your Name <your.email@example.com>"
    url = "https://github.com/yourname/aeroacoustic-bem"
    license = "MIT"

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "with_openmp": [True, False],
        "with_mkl": [True, False],
        "enable_testing": [True, False],
        "enable_examples": [True, False],
    }
    default_options = {
        "shared": False,
        "fPIC": True,
        "with_openmp": True,
        "with_mkl": False,
        "enable_testing": True,
        "enable_examples": True,
    }

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = "CMakeLists.txt", "src/*", "tests/*", "examples/*"

    def config_options(self):
        if self.settings.os == "Windows":
            self.options.rm_safe("fPIC")

    def configure(self):
        if self.options.shared:
            self.options.rm_safe("fPIC")

    def layout(self):
        cmake_layout(self)

    def requirements(self):
        # Linear algebra
        self.requires("armadillo/12.6.4")

        # BLAS backend - choose one
        if self.options.with_mkl:
            self.requires("intel-oneapi-mkl/2024.0")
        else:
            self.requires("openblas/0.3.25")

        # Logging and formatting
        self.requires("spdlog/1.12.0")
        self.requires("fmt/10.1.1")

        # JSON for configuration
        self.requires("nlohmann_json/3.11.3")

        # File I/O and utilities
        self.requires("tinyxml2/10.0.0")  # For XML mesh files

        # Optional: Parallel computing
        if self.options.with_openmp:
            if self.settings.compiler == "gcc":
                self.requires("openmp/system")

    def build_requirements(self):
        if self.options.enable_testing:
            self.requires("gtest/1.14.0")

        # Build tools
        self.tool_requires("cmake/[>=3.20]")

    def generate(self):
        deps = CMakeDeps(self)
        deps.generate()

        tc = CMakeToolchain(self)

        # Configure Armadillo
        tc.variables["ARMA_DONT_USE_WRAPPER"] = True
        tc.variables["ARMA_USE_LAPACK"] = True
        tc.variables["ARMA_USE_BLAS"] = True

        if self.options.with_mkl:
            tc.variables["BLA_VENDOR"] = "Intel10_64lp"
            tc.variables["ARMA_USE_MKL_ALLOC"] = True
        else:
            tc.variables["BLA_VENDOR"] = "OpenBLAS"

        # OpenMP support
        if self.options.with_openmp:
            tc.variables["CMAKE_CXX_FLAGS"] = "-fopenmp"
            tc.variables["CMAKE_EXE_LINKER_FLAGS"] = "-fopenmp"

        # Testing
        tc.variables["ENABLE_TESTING"] = self.options.enable_testing
        tc.variables["ENABLE_EXAMPLES"] = self.options.enable_examples

        # Build type specific settings
        if self.settings.build_type == "Debug":
            tc.variables["CMAKE_CXX_FLAGS_DEBUG"] = "-g -O0 -DDEBUG"
        else:
            tc.variables["CMAKE_CXX_FLAGS_RELEASE"] = "-O3 -DNDEBUG -march=native"

        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

        if self.options.enable_testing:
            cmake.test()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        # Main library
        self.cpp_info.libs = [
            "bem_application",
            "bem_domain",
            "bem_infrastructure",
            "bem_foundation",
        ]

        # Include directories
        self.cpp_info.includedirs = ["include"]

        # System libraries
        if self.settings.os in ["Linux", "FreeBSD"]:
            self.cpp_info.system_libs = ["m", "pthread"]
        elif self.settings.os == "Windows":
            self.cpp_info.system_libs = ["ws2_32"]

        # Preprocessor definitions
        if self.options.with_openmp:
            self.cpp_info.defines = ["BEM_USE_OPENMP"]

        if self.options.with_mkl:
            self.cpp_info.defines.append("BEM_USE_MKL")

        # Compiler flags
        if self.settings.compiler == "gcc":
            if self.options.with_openmp:
                self.cpp_info.cxxflags = ["-fopenmp"]
                self.cpp_info.sharedlinkflags = ["-fopenmp"]
                self.cpp_info.exelinkflags = ["-fopenmp"]

    def deploy(self):
        # Copy binaries to a bin folder
        self.copy("*.exe", dst="bin", keep_path=False)
        self.copy("aeroacoustic_bem", dst="bin", keep_path=False)

        # Copy examples and documentation
        if self.options.enable_examples:
            self.copy("examples/*", dst="examples", keep_path=True)


# Additional configuration files
def create_profile():
    """Create a default Conan profile for the project"""
    profile_content = """
[settings]
os=Linux
arch=x86_64
compiler=gcc
compiler.version=11
compiler.libcxx=libstdc++11
build_type=Release

[options]
*:shared=False
armadillo*:shared=False
openblas*:shared=False

[conf]
tools.system.package_manager:mode=install
tools.system.package_manager:sudo=True
"""
    return profile_content


def create_conandata():
    """Create conandata.yml for version management"""
    conandata_content = """
sources:
  "1.0.0":
    url: "https://github.com/yourname/aeroacoustic-bem/archive/v1.0.0.tar.gz"
    sha256: "your_sha256_here"
"""
    return conandata_content
