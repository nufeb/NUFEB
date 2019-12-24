Source and script files for building NUFEB third-party packages.

---------------------------------------------------------------------------
### HDF5 (v1.10.5)
Enable NUFEB to output simulation results in HDF5 format. NUFEB command that requires this library:

- *dump bio/hdf5*

### VTK (v8.0.0)
Enable NUFEB to output simulation results in VTK format. NUFEB commands that require this library:

- *dump custom/vtk*
- *dump grid*

### OpenFoam (v2.4.0)
A computational fluid dynamics (CFD) software package for the simulation of 3D fluid flow. The package 
is required by SediFOAM. We recommand you refer the website 
https://openfoamwiki.net/index.php/Installation/Linux/OpenFOAM-2.4.0 for the installation.

### SediFoam
A hybrid CFD-DEM solver for simulating microbial community in fluid dynamics. 

> **NOTE**: A 'mpi.h no such file' error may occur in some Linux systems due to the incorrected openmpi link provided in 
*/thirdparty/sediFoam/lammpsFoam/Make/options-* * files. To resolve the problem, you can either correct the openmpi path in the files, or
use the following way to build the solver:

> <pre>
> ./install-openmpi.sh
> ./install-sedifoam.sh --local-ompi
> </pre>
