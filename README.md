# AKATOSH - Neutron Transport Code

AKATOSH is a high-performance neutron transport and diffusion solver developed in Fortran, designed for nuclear reactor physics calculations. The code implements both finite element and spectral element methods for solving the neutron transport equation in 1D, 2D, and 3D geometries.

## Features

- **Multiple Solution Methods**: Supports both diffusion and transport equation solvers
- **Advanced Numerical Methods**: Implements finite element and spectral element discretizations
- **NURBS Support**: Non-Uniform Rational B-Splines for geometry representation
- **High Performance**: Parallel execution with OpenMP and PETSc integration
- **Multiple Geometries**: 1D, 2D, and 3D problem support
- **Material Libraries**: Extensive nuclear data libraries for various reactor configurations

## Project Structure

```
├── src/               # Fortran source code modules
├── input_files/       # Example input files for various problems
├── materials/         # Nuclear material data files
├── py_scripts/        # Python analysis and plotting utilities
├── auto_scripts/      # Automated problem setup scripts
├── main.f90          # Main program entry point
├── CMakeLists.txt    # CMake build configuration
└── build.sh          # Build script
```

## Dependencies

### Required
- **Fortran Compiler**: GNU Fortran (gfortran) or Intel Fortran (ifx)
- **CMake**: Version 3.6 or higher
- **PETSc**: Portable, Extensible Toolkit for Scientific Computation
- **OpenBLAS**: Optimized BLAS library
- **OpenMP**: For parallel execution

### Optional
- **Python 3**: For analysis scripts and visualization
- **Matplotlib**: For plotting results
- **NumPy/SciPy**: For numerical analysis

## Installation

### 1. Install PETSc

Choose one of the following configurations:

#### For Intel Fortran Compiler:
```bash
./configure --prefix=/path/to/petsc-ifx \
            --with-cc=gcc \
            --with-cxx=g++ \
            --with-fc=ifx \
            --download-mpich \
            --download-openblas
```

#### For GNU Fortran Compiler:
```bash
./configure --prefix=/path/to/petsc-openmp \
            --with-debugging=0 \
            --download-openblas \
            --with-threadsafety=1 \
            --with-openmp=1 \
            --download-mpich=1 \
            COPTFLAGS='-O2 -march=native -mtune=native' \
            CXXOPTFLAGS='-O2 -march=native -mtune=native' \
            FOPTFLAGS='-O2 -march=native -mtune=native'
```

### 2. Build AKATOSH

```bash
# Set environment variables
export PETSC_DIR=/path/to/your/petsc/installation
export OPENBLAS_NUM_THREADS=1

# Build the project
./build.sh
# or manually:
mkdir build && cd build
cmake ..
make -j4
```

## Usage

### Basic Execution

```bash
# Run with default input file
./akatosh

# Run with specific input file
./akatosh input_files/pincell.in

# Run with OpenMP parallelization
export OMP_NUM_THREADS=4
./akatosh input_files/c5g7_pin.in
```

### Example Problems

The `input_files/` directory contains several example problems:

- `pincell.in` - Simple pin cell calculation
- `c5g7_pin.in` - C5G7 benchmark pin cell
- `c5g7_assy.in` - C5G7 assembly calculation
- `htr2D.in` / `htr3D.in` - High-temperature reactor problems
- `supercell.in` - Supercell calculations

### Python Analysis Scripts

```bash
# Analyze results
python py_scripts/read_outdata.py

# Plot flux distributions
python py_scripts/analytic_analysis.py

# Generate animations
python py_scripts/animate.py
```

## Performance Tuning

For optimal performance:

```bash
# Set number of OpenBLAS threads
export OPENBLAS_NUM_THREADS=1

# Set number of OpenMP threads
export OMP_NUM_THREADS=8

# Profile execution
gprof ./akatosh gmon.out > prof.txt
```

## Contributing

This is a research code developed as part of PhD studies. For questions or collaboration opportunities, please contact the author.

## License

This project is part of academic research. Please contact the author for usage permissions.

## Author

Ciaran Jones  
PhD Candidate  
Imperial College London
ciaran.jones19@imperial.ac.uk

 --------------------------------------------------
 | Reached maximum iterations:    1               |
 | keff  =  0.8412483958E+00                      |
 | Error =  1.8870954767E-01                      |
 --------------------------------------------------
