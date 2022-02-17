# LAMMPS multiple-machine -*- Makefile -*-

SHELL = /bin/bash
PYTHON = python

#.IGNORE:

# Definitions

ROOT =	 lmp
EXE =	 lmp_$@
ARLIB =  liblammps_$@.a
SHLIB =	 liblammps_$@.so
ARLINK = liblammps.a
SHLINK = liblammps.so

OBJDIR =   Obj_$@
OBJSHDIR = Obj_shared_$@

SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)

SRCLIB = $(filter-out main.cpp,$(SRC))
OBJLIB = $(filter-out main.o,$(OBJ))

# Command-line options for mode: exe (default), shexe, lib, shlib

mode = exe
objdir = $(OBJDIR)

ifeq ($(mode),shexe)
objdir = $(OBJSHDIR)
endif

ifeq ($(mode),lib)
objdir = $(OBJDIR)
endif

ifeq ($(mode),shlib)
objdir = $(OBJSHDIR)
endif

# Package variables

PACKAGE = asphere body class2 colloid compress coreshell dipole gpu \
	  granular kim \
	  kokkos kspace manybody mc meam misc molecule mpiio opt peri poems \
	  python qeq reax replica rigid shock snap srd voronoi

PACKUSER = user-atc user-awpmd user-cg-cmm user-colvars \
	   user-diffraction user-dpd user-drude user-eff user-fep user-h5md \
	   user-intel user-lb user-manifold user-mgpt user-misc user-molfile \
	   user-nc-dump user-omp user-phonon user-qmmm user-qtb \
	   user-quip user-reaxc user-smd user-smtbq user-sph user-tally \
	   user-vtk user-nufeb user-cfddem

PACKLIB = compress gpu kim kokkos meam mpiio poems python reax voronoi \
	  user-atc user-awpmd user-colvars user-h5md user-lb user-molfile \
	  user-nc-dump user-qmmm user-quip user-smd user-vtk

PACKALL = $(PACKAGE) $(PACKUSER)

PACKAGEUC = $(shell echo $(PACKAGE) | tr a-z A-Z)
PACKUSERUC = $(shell echo $(PACKUSER) | tr a-z A-Z)

YESDIR = $(shell echo $(@:yes-%=%) | tr a-z A-Z)
NODIR  = $(shell echo $(@:no-%=%) | tr a-z A-Z)

# List of all targets

help:
	@echo ''
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'
	@echo 'make mpi-stubs           build dummy MPI library in STUBS'
	@echo 'make install-python      install LAMMPS wrapper in Python'
	@echo 'make tar                 create lmp_src.tar.gz of src dir and packages'
	@echo ''
	@echo 'make package                 list available packages'
	@echo 'make package-status (ps)     status of all packages'
	@echo 'make yes-package             install a single pgk in src dir'
	@echo 'make no-package              remove a single pkg from src dir'
	@echo 'make yes-all                 install all pgks in src dir'
	@echo 'make no-all                  remove all pkgs from src dir'
	@echo 'make yes-standard (yes-std)  install all standard pkgs'
	@echo 'make no-standard (no-std)    remove all standard pkgs'
	@echo 'make yes-user                install all user pkgs'
	@echo 'make no-user                 remove all user pkgs'
	@echo 'make no-lib                  remove all pkgs with external libs'
	@echo ''
	@echo 'make package-update (pu) replace src files with updated package files'
	@echo 'make package-overwrite   replace package files with src files'
	@echo 'make package-diff (pd)   diff src files against package files'
	@echo 'make purge               purge obsolete copies of source files'
	@echo ''
	@echo 'make machine             build LAMMPS for machine'
	@echo 'make mode=lib machine    build LAMMPS as static lib for machine'
	@echo 'make mode=shlib machine  build LAMMPS as shared lib for machine'
	@echo 'make mode=shexe machine  build LAMMPS as shared exe for machine'
	@echo 'make makelist            create Makefile.list used by old makes'
	@echo 'make -f Makefile.list machine     build LAMMPS for machine (old)'
	@echo ''
	@echo 'machine is one of these from src/MAKE:'
	@echo ''
	@files="`ls MAKE/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''
	@echo '... or one of these from src/MAKE/OPTIONS:'
	@echo ''
	@files="`ls MAKE/OPTIONS/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''
	@echo '... or one of these from src/MAKE/MACHINES:'
	@echo ''
	@files="`ls MAKE/MACHINES/Makefile.*`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''
	@echo '... or one of these from src/MAKE/MINE:'
	@echo ''
	@files="`ls MAKE/MINE/Makefile.* 2>/dev/null`"; \
	  for file in $$files; do head -1 $$file; done
	@echo ''

# Build LAMMPS in one of 4 modes
# exe =   exe with static compile in Obj_machine (default)
# shexe = exe with shared compile in Obj_shared_machine
# lib =   static lib in Obj_machine
# shlib = shared lib in Obj_shared_machine

.DEFAULT:
	@if [ $@ = "serial" -a ! -f STUBS/libmpi_stubs.a ]; \
	  then $(MAKE) mpi-stubs; fi
	@test -f MAKE/Makefile.$@ -o -f MAKE/OPTIONS/Makefile.$@ -o \
	  -f MAKE/MACHINES/Makefile.$@ -o -f MAKE/MINE/Makefile.$@
	@if [ ! -d $(objdir) ]; then mkdir $(objdir); fi
	@$(SHELL) Make.sh style
	@if [ -f MAKE/MACHINES/Makefile.$@ ]; \
	  then cp MAKE/MACHINES/Makefile.$@ $(objdir)/Makefile; fi
	@if [ -f MAKE/OPTIONS/Makefile.$@ ]; \
	  then cp MAKE/OPTIONS/Makefile.$@ $(objdir)/Makefile; fi
	@if [ -f MAKE/Makefile.$@ ]; \
	  then cp MAKE/Makefile.$@ $(objdir)/Makefile; fi
	@if [ -f MAKE/MINE/Makefile.$@ ]; \
	  then cp MAKE/MINE/Makefile.$@ $(objdir)/Makefile; fi
	@if [ ! -e Makefile.package ]; \
	  then cp Makefile.package.empty Makefile.package; fi
	@if [ ! -e Makefile.package.settings ]; \
	  then cp Makefile.package.settings.empty Makefile.package.settings; fi
	@cp Makefile.package Makefile.package.settings $(objdir)
	@cd $(objdir); rm -f .depend; \
	$(MAKE) $(MFLAGS) "SRC = $(SRC)" "INC = $(INC)" depend || :
ifeq ($(mode),exe)
	@cd $(objdir); \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJ)" "INC = $(INC)" "SHFLAGS =" \
	  "EXE = ../$(EXE)" ../$(EXE)
endif
ifeq ($(mode),shexe)
	@cd $(objdir); \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJ)" "INC = $(INC)" \
	  "EXE = ../$(EXE)" ../$(EXE)
endif
ifeq ($(mode),lib)
	@cd $(objdir); \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJLIB)" "INC = $(INC)" "SHFLAGS =" \
	  "EXE = ../$(ARLIB)" lib
	@rm -f $(ARLINK)
	@ln -s $(ARLIB) $(ARLINK)
endif
ifeq ($(mode),shlib)
	@cd $(objdir); \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJLIB)" "INC = $(INC)" \
	  "EXE = ../$(SHLIB)" shlib
	@rm -f $(SHLINK)
	@ln -s $(SHLIB) $(SHLINK)
endif

# Remove machine-specific object files

clean:
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'

clean-all:
	rm -rf Obj_*
clean-%:
	rm -rf Obj_$(@:clean-%=%) Obj_shared_$(@:clean-%=%)

# Create Makefile.list

makelist:
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh Makefile.list

# Make MPI STUBS library

mpi-stubs:
	@cd STUBS; $(MAKE) clean; $(MAKE)

# install LAMMPS shared lib and Python wrapper for Python usage

install-python:
	@$(PYTHON) ../python/install.py

# Create a tarball of src dir and packages

tar:
	@cd STUBS; $(MAKE) clean
	@cd ..; tar cvzf src/$(ROOT)_src.tar.gz \
	  src/Make* src/Package.sh src/Depend.sh src/Install.sh \
	  src/MAKE src/DEPEND src/*.cpp src/*.h src/STUBS \
	  $(patsubst %,src/%,$(PACKAGEUC)) $(patsubst %,src/%,$(PACKUSERUC)) \
          --exclude=*/.svn
	@cd STUBS; $(MAKE)
	@echo "Created $(ROOT)_src.tar.gz"

# Package management

package:
	@echo 'Standard packages:' $(PACKAGE)
	@echo ''
	@echo 'User-contributed packages:' $(PACKUSER)
	@echo ''
	@echo 'make package                 list available packages'
	@echo 'make package-status (ps)     status of all packages'
	@echo 'make yes-package             install a single pgk in src dir'
	@echo 'make no-package              remove a single pkg from src dir'
	@echo 'make yes-all                 install all pgks in src dir'
	@echo 'make no-all                  remove all pkgs from src dir'
	@echo 'make yes-standard (yes-std)  install all standard pkgs'
	@echo 'make no-standard (no-srd)    remove all standard pkgs'
	@echo ''
	@echo 'make yes-user                install all user pkgs'
	@echo 'make no-user                 remove all user pkgs'
	@echo 'make no-lib                  remove all pkgs with external libs'
	@echo 'make package-update (pu)  replace src files with package files'
	@echo 'make package-overwrite    replace package files with src files'
	@echo 'make package-diff (pd)    diff src files against package file'

yes-all:
	@for p in $(PACKALL); do $(MAKE) yes-$$p; done

no-all:
	@for p in $(PACKALL); do $(MAKE) no-$$p; done

yes-standard yes-std:
	@for p in $(PACKAGE); do $(MAKE) yes-$$p; done

no-standard no-std:
	@for p in $(PACKAGE); do $(MAKE) no-$$p; done

yes-user:
	@for p in $(PACKUSER); do $(MAKE) yes-$$p; done

no-user:
	@for p in $(PACKUSER); do $(MAKE) no-$$p; done

no-lib:
	@for p in $(PACKLIB); do $(MAKE) no-$$p; done

yes-%:
	@if [ ! -e Makefile.package ]; \
	  then cp Makefile.package.empty Makefile.package; fi
	@if [ ! -e Makefile.package.settings ]; \
	  then cp Makefile.package.settings.empty Makefile.package.settings; fi
	@if [ ! -e $(YESDIR) ]; then \
	  echo "Package $(@:yes-%=%) does not exist"; \
	elif [ -e $(YESDIR)/Install.sh ]; then \
	  echo "Installing package $(@:yes-%=%)"; \
	  cd $(YESDIR); $(SHELL) Install.sh 1; cd ..; \
		$(SHELL) Depend.sh $(YESDIR) 1; \
	else \
	  echo "Installing package $(@:yes-%=%)"; \
	  cd $(YESDIR); $(SHELL) ../Install.sh 1; cd ..; \
		$(SHELL) Depend.sh $(YESDIR) 1; \
	fi;

no-%:
	@if [ ! -e $(NODIR) ]; then \
	  echo "Package $(@:no-%=%) does not exist"; \
	elif [ -e $(NODIR)/Install.sh ]; then \
	  echo "Uninstalling package $(@:no-%=%)"; \
	  cd $(NODIR); $(SHELL) Install.sh 0; cd ..; \
		$(SHELL) Depend.sh $(NODIR) 0; \
	else \
	  echo "Uninstalling package $(@:no-%=%)"; \
	  cd $(NODIR); $(SHELL) ../Install.sh 0; cd ..; \
		$(SHELL) Depend.sh $(NODIR) 0; \
        fi;

# status = list src files that differ from package files
# update = replace src files with newer package files
# overwrite = overwrite package files with newer src files
# diff = show differences between src and package files
# purge = delete obsolete and auto-generated package files

package-status ps:
	@for p in $(PACKAGEUC); do $(SHELL) Package.sh $$p status; done
	@echo ''
	@for p in $(PACKUSERUC); do $(SHELL) Package.sh $$p status; done

package-update pu:
	@for p in $(PACKAGEUC); do $(SHELL) Package.sh $$p update; done
	@echo ''
	@for p in $(PACKUSERUC); do $(SHELL) Package.sh $$p update; done

package-overwrite:
	@for p in $(PACKAGEUC); do $(SHELL) Package.sh $$p overwrite; done
	@echo ''
	@for p in $(PACKUSERUC); do $(SHELL) Package.sh $$p overwrite; done

package-diff pd:
	@for p in $(PACKAGEUC); do $(SHELL) Package.sh $$p diff; done
	@echo ''
	@for p in $(PACKUSERUC); do $(SHELL) Package.sh $$p diff; done

purge: Purge.list
	@echo 'Purging obsolete and auto-generated source files'
	@for f in `grep -v '#' Purge.list` ;		\
	    do test -f $$f && rm $$f && echo $$f || : ;		\
	done
