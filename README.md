# Aeroacoustic BEM Solver

A high-performance C++ implementation of the Boundary Element Method (BEM) for solving 2D Helmholtz equations in aeroacoustic problems, inspired by the methodology in "Numerical computation of aeroacoustic transfer functions for realistic airfoils" by Miotto et al (2017).

## Features

- **Modern C++20** implementation with SOLID principles
- **High-performance linear algebra** using Eigen3
- **Flexible boundary conditions** (Dirichlet, Neumann, Robin)
- **Multiple integration strategies** (Gauss quadrature, adaptive integration)
- **Comprehensive test suite** with GoogleTest
- **Professional development tools** (clang-format, clang-tidy, coverage)

## Quick Start

### Prerequisites

- C++20 compatible compiler (GCC 11+, Clang 12+, MSVC 2019+)  
- CMake 3.20+
- Conan 2.0+
- Python 3.7+

### Installation

1. **Clone the repository:**
```bash
git clone https://github.com/yourname/aeroacoustic-bem.git
cd aeroacoustic-bem
```

2. **Install dependencies and build:**
```bash
# Install Conan if not already installed
pip install conan

# Create default profile
conan profile detect --force

# Build the project
./scripts/build.sh Release
```

3. **Run tests:**
```bash
./scripts/test.sh
```

4. **Run example:**
```bash
./build/aeroacoustic_bem examples/simple_scattering/config.json
```

### Build Options

The build system supports several configuration options:

```bash
# Debug build with sanitizers
./scripts/build.sh Debug

# Release build with Intel MKL
conan install . --build=missing -o with_mkl=True
./scripts/build.sh Release

# Build with coverage reporting
./scripts/test.sh --coverage
```

## Usage Example

```cpp
#include "application/App.hpp"
#include "domain/physics/HelmholtzProblem.hpp"
#include "domain/bem/HelmholtzBEMSolver.hpp"

int main() {
    // Create problem configuration
    auto config = bem::Configuration::loadFromFile("config.json");
    
    // Set up the Helmholtz problem
    bem::HelmholtzProblem problem(config.frequency, config.wave_speed);
    
    // Load mesh and boundary conditions
    auto mesh = bem::MeshFactory::createFromFile(config.mesh_file);
    auto boundary_conditions = bem::BoundaryConditionFactory::createFromConfig(config);
    
    // Create and configure solver
    bem::HelmholtzBEMSolver solver(std::make_unique<bem::DirectSolver>());
    solver.setIntegrationStrategy(std::make_unique<bem::GaussQuadrature>(8));
    
    // Solve the problem
    auto solution = solver.solve(problem, *mesh, boundary_conditions);
    
    // Export results
    bem::VTKWriter writer("results.vtk");
    writer.write(*mesh, solution);
    
    return 0;
}
```

## Configuration File Format

The solver uses JSON configuration files:

```json
{
    "problem": {
        "frequency": 1000.0,
        "wave_speed": 343.0,
        "domain_type": "exterior"
    },
    "mesh": {
        "file": "airfoil.msh",
        "format": "gmsh"
    },
    "boundary_conditions": [
        {
            "element_ids": [1, 2, 3],
            "type": "dirichlet",
            "value": {
                "real": 1.0,
                "imag": 0.0
            }
        }
    ],
    "solver": {
        "type": "direct",
        "integration_order": 8,
        "singularity_treatment": "analytical"
    },
    "output": {
        "format": "vtk",
        "file": "solution.vtk"
    }
}
```

## Architecture Overview

The software follows a layered architecture with clear separation of concerns:

```
Application Layer    → CLI, Configuration, Main Program
     ↓
Domain Layer         → BEM Logic, Physics, Boundary Conditions  
     ↓
Infrastructure Layer → Linear Algebra, I/O, Logging
     ↓
Foundation Layer     → Utilities, Types, Error Handling
```

Key design patterns used:
- **Strategy Pattern** for numerical methods and solvers
- **Factory Pattern** for object creation
- **Observer Pattern** for progress monitoring
- **Template Method Pattern** for BEM algorithm structure

## Development

### Code Style

The project uses automated formatting and static analysis:

```bash
# Format code
./scripts/format.sh

# Run static analysis
clang-tidy src/**/*.cpp
```

### Testing

Tests are organized into unit and integration tests:

```bash
# Run all tests
./scripts/test.sh

# Run specific test suite
cd build && ctest -R "test_BEMSolver"

# Generate coverage report
./scripts/test.sh --coverage
```

### Adding New Features

1. **Write tests first** (TDD approach)
2. **Follow SOLID principles** 
3. **Update documentation**
4. **Run full test suite**

Example of adding a new boundary condition:

```cpp
// 1. Define interface (already exists)
class BoundaryCondition { /* ... */ };

// 2. Implement new condition
class RobinBC : public BoundaryCondition {
public:
    Complex evaluate(const Point2D& point) const override;
    // ...
};

// 3. Add to factory
class BoundaryConditionFactory {
    static std::unique_ptr<BoundaryCondition> create(const std::string& type) {
        if (type == "robin") return std::make_unique<RobinBC>();
        // ...
    }
};

// 4. Write tests
TEST(RobinBCTest, EvaluatesCorrectly) {
    // Test implementation
}
```

## Performance Notes

- **BLAS Backend**: OpenBLAS provides good performance for most cases. Use Intel MKL for maximum performance on Intel processors.
- **OpenMP**: Enable with `-o with_openmp=True` for parallel matrix assembly.
- **Memory**: Large problems may require careful memory management. The solver supports iterative methods for memory-constrained systems.

## Examples

The `examples/` directory contains several demonstration cases:

- **simple_scattering/**: Basic scattering from a circular cylinder
- **airfoil_analysis/**: NACA airfoil with realistic boundary conditions  
- **validation_cases/**: Comparison with analytical solutions

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## References

- Miotto, R., Wolf, W. R., & Azevedo, J. L. F. (2017). "Numerical computation of aeroacoustic transfer functions for realistic airfoils." *Journal of Sound and Vibration*, 407, 253-277.
- Kirkup, S. (2007). *The Boundary Element Method in Acoustics*. Integrated Sound Software.
- Wu, T. W. (Ed.). (2000). *Boundary element acoustics: fundamentals and computer codes*. WIT press.

## Support

For questions and support:
- Create an issue on GitHub
- Check the [documentation](docs/)
- Contact: renato.fmiotto@gmail.com
