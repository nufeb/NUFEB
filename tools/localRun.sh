cd runs
for f in *.lammps
do
mpirun -np 8 lmp_png -in "$f"
done
cd ..
