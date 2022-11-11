#-------------------------------------------------------------------------------------------------------------
# Makefile
#       Master Makefile for transport code.
#-------------------------------------------------------------------------------------------------------------

export

BUILDDIRS1 := $(wildcard src/*/.)
BUILDDIRS2 := $(wildcard src/*/*/.)
BUILDDIRS3 := $(wildcard src/*/*/*/.)
CLEANDIRS := $(addsuffix .clean,$(BUILDDIRS1) $(BUILDDIRS2) $(BUILDDIRS3))

INC := make.inc
include $(INC)

# --- Build commands. --------------------------------------------------------------------------------- #

# Compile solver, postproc, and deckmaker.
.PHONY: all
all: solver postproc deckmaker
	@echo "compilation finished."

.PHONY: $(BUILDDIRS1) $(BUILDDIRS2) $(BUILDDIRS3)

# Compile solver executables.
.PHONY: solver solver_1d solver_2d
solver: solver_1d solver_2d
	@echo "."

solver_1d: $(BUILDDIRS1) $(BUILDDIRS2) $(BUILDDIRS3) sha setenv
	$(CXX) -o solver_1d.x ./src/main_solver.cpp src/*/*.o1 src/*/*/*.o1 src/*/*/*/*.o1 \
	$(ALLFLAGS) $(CONFIGFLAGS_1D) $(OPTFLAGS) $(CXXFLAGS) \
	$(LIB) $(INCLUDE)

solver_2d: $(BUILDDIRS1) $(BUILDDIRS2) $(BUILDDIRS3) sha setenv
	$(CXX) -o solver_2d.x ./src/main_solver.cpp src/*/*.o2 src/*/*/*.o2 src/*/*/*/*.o2 \
	$(ALLFLAGS) $(CONFIGFLAGS_2D) $(OPTFLAGS) $(CXXFLAGS) \
	$(LIB) $(INCLUDE)

# Compile postproc executables.
.PHONY: postproc postproc_1d postproc_2d
postproc: postproc_1d postproc_2d
	@echo "."

postproc_1d: $(BUILDDIRS1) $(BUILDDIRS2) $(BUILDDIRS3) sha setenv
	$(CXX) -o postproc_1d.x ./src/main_postproc.cpp src/*/*.o1 src/*/*/*.o1 src/*/*/*/*.o1 \
	$(ALLFLAGS) $(CONFIGFLAGS_1D) $(OPTFLAGS) $(CXXFLAGS) \
	$(LIB) $(INCLUDE)

postproc_2d: $(BUILDDIRS1) $(BUILDDIRS2) $(BUILDDIRS3) sha setenv
	$(CXX) -o postproc_2d.x ./src/main_postproc.cpp src/*/*.o2 src/*/*/*.o2 src/*/*/*/*.o2 \
	$(ALLFLAGS) $(CONFIGFLAGS_2D) $(OPTFLAGS) $(CXXFLAGS) \
	$(LIB) $(INCLUDE)

# Compile deckmaker executables.
.PHONY: deckmaker deckmaker_1d deckmaker_2d
deckmaker: deckmaker_1d deckmaker_2d
	@echo "."

deckmaker_1d: $(BUILDDIRS1) $(BUILDDIRS2) $(BUILDDIRS3) sha setenv
	$(CXX) -o deckmaker_1d.x ./src/main_deckmaker.cpp src/*/*.o1 src/*/*/*.o1 src/*/*/*/*.o1 \
	$(ALLFLAGS) $(CONFIGFLAGS_1D) $(OPTFLAGS) $(CXXFLAGS) \
	$(LIB) $(INCLUDE)

deckmaker_2d: $(BUILDDIRS1) $(BUILDDIRS2) $(BUILDDIRS3) sha setenv
	$(CXX) -o deckmaker_2d.x ./src/main_deckmaker.cpp src/*/*.o2 src/*/*/*.o2 src/*/*/*/*.o2 \
	$(ALLFLAGS) $(CONFIGFLAGS_2D) $(OPTFLAGS) $(CXXFLAGS) \
	$(LIB) $(INCLUDE)

.PHONY: sha
sha:
	git log --pretty=format:"%H" -n 1 > src/sha
	xxd -i src/sha > src/sha.h

.PHONY: setenv
setenv:
	echo "#!/bin/bash" > setenv.sh
	echo "module load $(ENVIRONMENT_MODULES)" >> setenv.sh
	echo -n "PATH=$$" >> setenv.sh
	echo "{PATH}:$$(pwd)/" >> setenv.sh


# Compile object files.
$(BUILDDIRS1):
	cp make/Makefile.sub1 $@/Makefile
	$(MAKE) -C $@

$(BUILDDIRS2):
	cp make/Makefile.sub2 $@/Makefile
	$(MAKE) -C $@

$(BUILDDIRS3):
	cp make/Makefile.sub3 $@/Makefile
	$(MAKE) -C $@


# Compile documentation.
.PHONY: doxygen
doxygen:
	doxygen Doxyfile


# Plotting commands.
.PHONY: plot1
plot1:
	gnuplot idc.gnuplot
	xdg-open idc.png 2> /dev/null

.PHONY: plot2
plot2:
	gnuplot heatmap.gnuplot
	xdg-open heatmap.png 2> /dev/null

.PHONY: plotdiff2
plotdiff2:
	gnuplot heatmap_diff.gnuplot
	xdg-open heatmap_diff.png 2> /dev/null


# Package code for distribution.
.PHONY: archive
archive: $(eval SHA=$(shell git log --pretty=format:"%H" -n 1))
archive:
	git archive --format=tar --prefix=solver/ HEAD | bzip2 > solver_$(SHA).tar.bz2


# Clean commands.
.PHONY: clean $(CLEANDIRS)
clean: $(CLEANDIRS)
	-rm *.x
	-rm setenv.sh

$(CLEANDIRS): %.clean:
	-$(MAKE) -C $* clean

.PHONY: cleanall clean
cleanall: clean
	-rm -r *.out *.png *.af? *.sd? *.cs? *.log *.chk html latex C_* H[0-9]* U_* src/sha*

