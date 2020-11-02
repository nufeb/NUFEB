A simple biofilm growth model that consists of heterotrophs (HETs) and their EPS production. 
The biofilm is grown from 40 initial HETs inoculated on the substratum for ~9.25 days (8x10^5s).

<pre>       
Inputscript-vtk.lammps           Simulation inputscript (with VTK outputs) 
Inputscript.lammps               Simulation inputscript (without VTK outputs) 
Inputscript-comments.lammps      Simulation inputscript (with command explanation) 
atom.in                          Data file defining initial microbes, species, nutrients, etc 
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

[![youtube video](https://img.youtube.com/vi/DrDD7_OZNQ4/0.jpg)](https://www.youtube.com/watch?v=DrDD7_OZNQ4)
