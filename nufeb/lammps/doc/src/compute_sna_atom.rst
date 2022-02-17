.. index:: compute sna/atom

compute sna/atom command
========================

compute snad/atom command
=========================

compute snav/atom command
=========================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID sna/atom rcutfac rfac0 twojmax R_1 R_2 ... w_1 w_2 ... keyword values ...
   compute ID group-ID snad/atom rcutfac rfac0 twojmax R_1 R_2 ... w_1 w_2 ... keyword values ...
   compute ID group-ID snav/atom rcutfac rfac0 twojmax R_1 R_2 ... w_1 w_2 ... keyword values ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* sna/atom = style name of this compute command
* rcutfac = scale factor applied to all cutoff radii (positive real)
* rfac0 = parameter in distance to angle conversion (0 < rcutfac < 1)
* twojmax = band limit for bispectrum components (non-negative integer)
* R\_1, R\_2,... = list of cutoff radii, one for each type (distance units)
* w\_1, w\_2,... = list of neighbor weights, one for each type
* zero or more keyword/value pairs may be appended
* keyword = *rmin0* or *switchflag* or *bzeroflag* or *quadraticflag*
  
  .. parsed-literal::
  
       *rmin0* value = parameter in distance to angle conversion (distance units)
       *switchflag* value = *0* or *1*
          *0* = do not use switching function
          *1* = use switching function
       *bzeroflag* value = *0* or *1*
          *0* = do not subtract B0
          *1* = subtract B0
       *quadraticflag* value = *0* or *1*
          *0* = do not generate quadratic terms
          *1* = generate quadratic terms



Examples
""""""""


.. parsed-literal::

   compute b all sna/atom 1.4 0.99363 6 2.0 2.4 0.75 1.0 rmin0 0.0
   compute db all sna/atom 1.4 0.95 6 2.0 1.0
   compute vb all sna/atom 1.4 0.95 6 2.0 1.0

Description
"""""""""""

Define a computation that calculates a set of bispectrum components
for each atom in a group.

Bispectrum components of an atom are order parameters characterizing
the radial and angular distribution of neighbor atoms. The detailed
mathematical definition is given in the paper by Thompson et
al. :ref:`(Thompson) <Thompson20141>`

The position of a neighbor atom *i'* relative to a central atom *i* is
a point within the 3D ball of radius *R\_ii' = rcutfac\*(R\_i + R\_i')*

Bartok et al. :ref:`(Bartok) <Bartok20101>`, proposed mapping this 3D ball
onto the 3-sphere, the surface of the unit ball in a four-dimensional
space.  The radial distance *r* within *R\_ii'* is mapped on to a third
polar angle *theta0* defined by,

.. image:: Eqs/compute_sna_atom1.jpg
   :align: center

In this way, all possible neighbor positions are mapped on to a subset
of the 3-sphere.  Points south of the latitude *theta0max=rfac0\*Pi*
are excluded.

The natural basis for functions on the 3-sphere is formed by the 4D
hyperspherical harmonics *U\^j\_m,m'(theta, phi, theta0).*  These
functions are better known as *D\^j\_m,m',* the elements of the Wigner
*D*\ -matrices :ref:`(Meremianin <Meremianin2006>`,
:ref:`Varshalovich) <Varshalovich1987>`.

The density of neighbors on the 3-sphere can be written as a sum of
Dirac-delta functions, one for each neighbor, weighted by species and
radial distance. Expanding this density function as a generalized
Fourier series in the basis functions, we can write each Fourier
coefficient as

.. image:: Eqs/compute_sna_atom2.jpg
   :align: center

The *w\_i'* neighbor weights are dimensionless numbers that are chosen
to distinguish atoms of different types, while the central atom is
arbitrarily assigned a unit weight.  The function *fc(r)* ensures that
the contribution of each neighbor atom goes smoothly to zero at
*R\_ii'*:

.. image:: Eqs/compute_sna_atom4.jpg
   :align: center

The expansion coefficients *u\^j\_m,m'* are complex-valued and they are
not directly useful as descriptors, because they are not invariant
under rotation of the polar coordinate frame. However, the following
scalar triple products of expansion coefficients can be shown to be
real-valued and invariant under rotation :ref:`(Bartok) <Bartok20101>`.

.. image:: Eqs/compute_sna_atom3.jpg
   :align: center

The constants *H\^jmm'\_j1m1m1'\_j2m2m2'* are coupling coefficients,
analogous to Clebsch-Gordan coefficients for rotations on the
2-sphere. These invariants are the components of the bispectrum and
these are the quantities calculated by the compute *sna/atom*\ . They
characterize the strength of density correlations at three points on
the 3-sphere. The j2=0 subset form the power spectrum, which
characterizes the correlations of two points. The lowest-order
components describe the coarsest features of the density function,
while higher-order components reflect finer detail.  Note that the
central atom is included in the expansion, so three point-correlations
can be either due to three neighbors, or two neighbors and the central
atom.

Compute *snad/atom* calculates the derivative of the bispectrum components
summed separately for each atom type:

.. image:: Eqs/compute_sna_atom5.jpg
   :align: center

The sum is over all atoms *i'* of atom type *I*\ .  For each atom *i*\ ,
this compute evaluates the above expression for each direction, each
atom type, and each bispectrum component.  See section below on output
for a detailed explanation.

Compute *snav/atom* calculates the virial contribution due to the
derivatives:

.. image:: Eqs/compute_sna_atom6.jpg
   :align: center

Again, the sum is over all atoms *i'* of atom type *I*\ .  For each atom
*i*\ , this compute evaluates the above expression for each of the six
virial components, each atom type, and each bispectrum component.  See
section below on output for a detailed explanation.

The value of all bispectrum components will be zero for atoms not in
the group. Neighbor atoms not in the group do not contribute to the
bispectrum of atoms in the group.

The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (i.e. each time a snapshot of atoms
is dumped).  Thus it can be inefficient to compute/dump this quantity
too frequently.

The argument *rcutfac* is a scale factor that controls the ratio of
atomic radius to radial cutoff distance.

The argument *rfac0* and the optional keyword *rmin0* define the
linear mapping from radial distance to polar angle *theta0* on the
3-sphere.

The argument *twojmax* defines which
bispectrum components are generated. See section below on output for a
detailed explanation of the number of bispectrum components and the
ordered in which they are listed.

The keyword *switchflag* can be used to turn off the switching
function.

The keyword *bzeroflag* determines whether or not *B0*\ , the bispectrum
components of an atom with no neighbors, are subtracted from
the calculated bispectrum components. This optional keyword
normally only affects compute *sna/atom*\ . However, when
*quadraticflag* is on, it also affects *snad/atom* and *snav/atom*\ .

The keyword *quadraticflag* determines whether or not the
quadratic analogs to the bispectrum quantities are generated.
These are formed by taking the outer product of the vector
of bispectrum components with itself.
See section below on output for a
detailed explanation of the number of quadratic terms and the
ordered in which they are listed.

.. note::

   If you have a bonded system, then the settings of
   :doc:`special\_bonds <special_bonds>` command can remove pairwise
   interactions between atoms in the same bond, angle, or dihedral.  This
   is the default setting for the :doc:`special\_bonds <special_bonds>`
   command, and means those pairwise interactions do not appear in the
   neighbor list.  Because this fix uses the neighbor list, it also means
   those pairs will not be included in the calculation.  One way to get
   around this, is to write a dump file, and use the :doc:`rerun <rerun>`
   command to compute the bispectrum components for snapshots in the dump
   file.  The rerun script can use a :doc:`special\_bonds <special_bonds>`
   command that includes all pairs in the neighbor list.

;line

**Output info:**

Compute *sna/atom* calculates a per-atom array, each column
corresponding to a particular bispectrum component.  The total number
of columns and the identity of the bispectrum component contained in
each column depend of the value of *twojmax*\ , as
described by the following piece of python code:


.. parsed-literal::

   for j1 in range(0,twojmax+1):
       for j2 in range(0,j1+1):
           for j in range(j1-j2,min(twojmax,j1+j2)+1,2):
               if (j>=j1): print j1/2.,j2/2.,j/2.

.. note::

   the *diagonal* keyword allowing other possible choices
   for the number of bispectrum components was removed in 2019,
   since all potentials use the value of 3, corresponding to the
   above set of bispectrum components.

Compute *snad/atom* evaluates a per-atom array. The columns are
arranged into *ntypes* blocks, listed in order of atom type *I*\ .  Each
block contains three sub-blocks corresponding to the *x*\ , *y*\ , and *z*
components of the atom position.  Each of these sub-blocks contains
one column for each bispectrum component, the same as for compute
*sna/atom*

Compute *snav/atom* evaluates a per-atom array. The columns are
arranged into *ntypes* blocks, listed in order of atom type *I*\ .  Each
block contains six sub-blocks corresponding to the *xx*\ , *yy*\ , *zz*\ ,
*yz*\ , *xz*\ , and *xy* components of the virial tensor in Voigt
notation.  Each of these sub-blocks contains one column for each
bispectrum component, the same as for compute *sna/atom*

For example, if *K* =30 and ntypes=1, the number of columns in the per-atom
arrays generated by *sna/atom*\ , *snad/atom*\ , and *snav/atom*
are 30, 90, and 180, respectively. With *quadratic* value=1,
the numbers of columns are 930, 2790, and 5580, respectively.

If the *quadratic* keyword value is set to 1, then additional
columns are generated, corresponding to
the products of all distinct pairs of  bispectrum components. If the
number of bispectrum components is *K*\ , then the number of distinct pairs
is  *K*\ (\ *K*\ +1)/2.
For compute *sna/atom* these columns are appended to existing *K* columns.
The ordering of quadratic terms is upper-triangular,
(1,1),(1,2)...(1,\ *K*\ ),(2,1)...(\ *K*\ -1,\ *K*\ -1),(\ *K*\ -1,\ *K*\ ),(\ *K*\ ,\ *K*\ ).
For computes *snad/atom* and *snav/atom* each set of *K*\ (\ *K*\ +1)/2
additional columns is inserted directly after each of sub-block
of linear terms i.e. linear and quadratic terms are contiguous.
So the nesting order from inside to outside is bispectrum component,
linear then quadratic, vector/tensor component, type.

These values can be accessed by any command that uses per-atom values
from a compute as input.  See the :doc:`Howto output <Howto_output>` doc
page for an overview of LAMMPS output options.

Restrictions
""""""""""""


These computes are part of the SNAP package.  They are only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair\_style snap <pair_snap>`

Default
"""""""

The optional keyword defaults are *rmin0* = 0,
*switchflag* = 1, *bzeroflag* = 1, *quadraticflag* = 0,


----------


.. _Thompson20141:



**(Thompson)** Thompson, Swiler, Trott, Foiles, Tucker, under review, preprint
available at `arXiv:1409.3880 <http://arxiv.org/abs/1409.3880>`_

.. _Bartok20101:



**(Bartok)** Bartok, Payne, Risi, Csanyi, Phys Rev Lett, 104, 136403 (2010).

.. _Meremianin2006:



**(Meremianin)** Meremianin, J. Phys. A,  39, 3099 (2006).

.. _Varshalovich1987:



**(Varshalovich)** Varshalovich, Moskalev, Khersonskii, Quantum Theory
of Angular Momentum, World Scientific, Singapore (1987).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
