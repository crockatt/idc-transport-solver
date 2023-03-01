#-------------------------------------------------------------------------------------------------------------
# Makefile variable configurations for using:
#   - MPI compiler wrappers (compiler agnostic: tested with gcc, intel, clang)
#   - MPI enabled
#   - petsc enabled
#   - hwloc enabled
#   - openblas
#-------------------------------------------------------------------------------------------------------------

SHELL := /bin/bash
MAKE  := make

ENVIRONMENT_MODULES := gcc/9.4.0 hwloc openblas openmpi petsc

# Compilers
CC  = mpicc
CXX = mpic++

# Flags
CFLAGS   := -std=gnu11
CXXFLAGS := -std=gnu++14

CONFIGFLAGS_1D := -DSPACE_DIMS=1
CONFIGFLAGS_2D := -DSPACE_DIMS=2

ALLFLAGS := -Wall -Wextra -Wshadow -Wpedantic \
            -DCOLOR_TERM -DLOGLEVEL=1 \
            -DENABLE_HWLOC -DUSE_HWLOC_ALLOC -DUSE_ALIGNED_ALLOC \
            -DENABLE_PETSC \
            -DENABLE_SIMD_BLOCKING -DSIMD_LEN=16
OPTFLAGS := -O3 -march=native -fopenmp -fopenmp-simd -fdiagnostics-color

# Library paths.
LIB = -lpetsc -lmpi -lhwloc \
      -lopenblas -lgomp -lpthread \
      -lm -ldl

# Include paths.
INCLUDE = -I ./src/ \
          -I ../ \
          -I ../../ \
          -I ../../../
