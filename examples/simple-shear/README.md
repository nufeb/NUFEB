Biofilm growth model and removal. 

The biofilm is grown from 40 initial HETs inoculated on the substratum for ~4.5 days (4x10^5s).
Then shear force is applied with the direction along x axis to remove the particle from the domain.

<pre>       
Inputscript-vtk.lammps           Simulation inputscript (with VTK outputs) 
Inputscript-jpeg.lammps          Simulation inputscript (with jpeg outputs) 
atom.in                          Data file defining initial microbes, species, nutrients, their kinetic parameters, etc
Allclean.sh                      Script for cleanup
</pre>

To run the simulation with VTK outputs:
<pre>
mpirun -np 4 lmp_mpi -in Inputscript-vtk.lammps
</pre>
or 
<pre>
mpirun -np 4 ../../lammps/src/./lmp_mpi -in Inputscript-vtk.lammps
</pre>

To cleanup output files
<pre>
./Allclean.sh
</pre>


Simulation result:

(↓ Youtube video)

[![youtube video](http://i3.ytimg.com/vi/DLCRtAS9xvA/maxresdefault.jpg)](https://www.youtube.com/watch?v=DLCRtAS9xvA)
