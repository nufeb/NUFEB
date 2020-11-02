A multi-species biofilm model that consists of 
ammonia oxidizing bacteria (AOB), nitrite oxidizing bacteria (NOB),
heterotrophs (HET) and their EPS production. 
The biofilm is grown from 100 initial microbes of each species that
are randomly placed on the substratum.

<pre>       
Inputscript.lammps               Simulation inputscript
atom.in                          Data file defining species, nutrients, their kinetic parameters, etc
Allclean.sh                      Script for cleanup
</pre>

To run the simulation:
<pre>
mpirun -np 4 lmp_mpi -in Inputscript.lammps
</pre>
or 
<pre>
mpirun -np 4 ../../lammps/src/./lmp_mpi -in Inputscript.lammps
</pre>

To cleanup output files
<pre>
./Allclean.sh
</pre>
