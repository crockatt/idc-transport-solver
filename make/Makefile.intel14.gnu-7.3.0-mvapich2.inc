#-------------------------------------------------------------------------------------------------------------
# Make options for using GNU-4.8.3 compiler with BLCR restart-capable executable on MSU HPCC.
#
# Author:    Michael M. Crockatt
# Modified:  February 2018
#-------------------------------------------------------------------------------------------------------------

SHELL := /bin/bash
MAKE  := make

# Compilers
CC  = mpicc
CXX = mpic++

# Flags
CFLAGS   := -std=gnu11
CXXFLAGS := -std=gnu++14
ALLFLAGS := -Wall -Wextra -Wpedantic \
            -DLOGLEVEL=1 -DSPACE_DIMS=1 \
            -DENABLE_INIT_PROJ_LIMITER \
            -DENABLE_PETSC -DENABLE_AZTEC_EPETRA -DENABLE_BELOS_TPETRA \
            -DENABLE_HWLOC -DUSE_HWLOC_ALLOC \
            -DENABLE_SIMD_BLOCKING -DSIMD_LEN=4
OPTFLAGS := -O3 -march=native -fopenmp -fopenmp-simd -fdiagnostics-color

# Library paths.
LIB	= -lhwloc \
      -lpetsc \
      -lepetra -laztecoo \
      -lbelos -lbelostpetra -ltpetra \
      -lteuchoskokkoscompat -lteuchosnumerics -lkokkoscore \
      -lteuchoscore -lteuchoscomm -lteuchosparameterlist \
      -lmpi \
      -Wl,--start-group \
      /opt/software/ClusterStudio/2017.4/mkl/lib/intel64/libmkl_intel_lp64.a \
      /opt/software/ClusterStudio/2017.4/mkl/lib/intel64/libmkl_sequential.a \
      /opt/software/ClusterStudio/2017.4/mkl/lib/intel64/libmkl_core.a \
      -Wl,--end-group \
      -ldl -lm -lrt


# Include paths.
INCLUDE = -I ./src/ \
          -I ../ \
          -I ../../ \
          -I ../../../ \
