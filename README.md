
This software is licensed under the BSD 3-clause license: see [LICENSE](LICENSE).

# Build Instructions

1. Copy or link to desired build configuration to `make.inc` in the top level directory, e.g.:
```
ln -s make/Makefile.workstation.inc make.inc
```
2. Edit module, compiler, and build options as appropriate.
1. Execute make command: `make -j <N>`.
1. Executables '*.x' will be created in the top level directory for 1D and 2D versions of the code.
1. A bash environment can be setup to run the executables through the generated `setenv.sh` file. E.g., run
```
source setenv.sh
```
