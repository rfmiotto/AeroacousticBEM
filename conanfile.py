# pyright: reportOptionalCall=false

from conan import ConanFile
from conan.tools.cmake import cmake_layout


class BEMConan(ConanFile):
    name = "aeroacoustic-bem"
    version = "1.0"
    settings = "os", "arch", "compiler", "build_type"
    generators = "CMakeDeps", "CMakeToolchain"
    exports_sources = "CMakeLists.txt", "src/*", "include/*", "tests/*"

    def requirements(self):
        self.requires("eigen/3.4.0")
        self.requires("openblas/0.3.25")
        self.requires("fmt/10.2.1")
        self.requires("spdlog/1.12.0")
        self.requires("nlohmann_json/3.11.3")
        self.requires("tinyxml2/10.0.0")
        self.requires("gtest/1.14.0")  # Added GTest dependency

    def build_requirements(self):
        self.tool_requires("cmake/[>=3.20]")

    def layout(self):
        cmake_layout(self)
