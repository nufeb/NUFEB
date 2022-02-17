Build LAMMPS with make
======================

Building LAMMPS with traditional makefiles requires that you have a
Makefile."machine" file appropriate for your system in the src/MAKE,
src/MAKE/MACHINES, src/MAKE/OPTIONS, or src/MAKE/MINE directory (see
below).  It can include various options for customizing your LAMMPS
build with a number of global compilation options and features.

To include LAMMPS packages (i.e. optional commands and styles) you
must install them first, as discussed on the :doc:`Build package <Build_package>` doc page.  If the packages require
provided or external libraries, you must build those libraries before
building LAMMPS.  Building :doc:`LAMMPS with CMake <Build_cmake>` can
automate all of this for many types of machines, especially
workstations, desktops and laptops, so we suggest you try it first.

These commands perform a default LAMMPS build, producing the LAMMPS
executable lmp\_serial or lmp\_mpi in lammps/src:


.. parsed-literal::

   cd lammps/src
   make serial     # build a serial LAMMPS executable
   make mpi        # build a parallel LAMMPS executable with MPI
   make            # see a variety of make options

This initial compilation can take a long time, since LAMMPS is a large
project with many features. If your machine has multiple CPU cores
(most do these days), using a command like "make -jN mpi" (with N =
the number of available CPU cores) can be much faster.  If you plan to
do development on LAMMPS or need to re-compile LAMMPS repeatedly, the
installation of the ccache (= Compiler Cache) software may speed up
compilation even more.

After the initial build, whenever you edit LAMMPS source files, or add
or remove new files to the source directory (e.g. by installing or
uninstalling packages), you must re-compile and relink the LAMMPS
executable with the same "make" command.  This makefiles dependencies
should insure that only the subset of files that need to be are
re-compiled.

.. note::

   When you build LAMMPS for the first time, a long list of \*.d
   files will be printed out rapidly.  This is not an error; it is the
   Makefile doing its normal creation of dependencies.


----------


The lammps/src/MAKE tree contains all the Makefile.machine files
included in the LAMMPS distribution.  Typing "make machine" uses
Makefile.machine.  Thus the "make serial" or "make mpi" lines above
use Makefile.serial and Makefile.mpi.  Others are in these dirs:


.. parsed-literal::

   OPTIONS      # Makefiles which enable specific options
   MACHINES     # Makefiles for specific machines
   MINE         # customized Makefiles you create (you may need to create this folder)

Typing "make" lists all the available Makefile.machine files.  A file
with the same name can appear in multiple folders (not a good idea).
The order the dirs are searched is as follows: src/MAKE/MINE,
src/MAKE, src/MAKE/OPTIONS, src/MAKE/MACHINES.  This gives preference
to a customized file you put in src/MAKE/MINE.

Makefiles you may wish to try include these (some require a package
first be installed).  Many of these include specific compiler flags
for optimized performance.  Please note, however, that some of these
customized machine Makefile are contributed by users.  Since both
compilers, OS configurations, and LAMMPS itself keep changing, their
settings may become outdated:


.. parsed-literal::

   make mac             # build serial LAMMPS on a Mac
   make mac_mpi         # build parallel LAMMPS on a Mac
   make intel_cpu       # build with the USER-INTEL package optimized for CPUs
   make knl             # build with the USER-INTEL package optimized for KNLs
   make opt             # build with the OPT package optimized for CPUs
   make omp             # build with the USER-OMP package optimized for OpenMP
   make kokkos_omp      # build with the KOKKOS package for OpenMP
   make kokkos_cuda_mpi # build with the KOKKOS package for GPUs
   make kokkos_phi      # build with the KOKKOS package for KNLs


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
