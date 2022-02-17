.. index:: improper\_style harmonic

improper\_style harmonic command
================================

improper\_style harmonic/intel command
======================================

improper\_style harmonic/kk command
===================================

improper\_style harmonic/omp command
====================================

Syntax
""""""


.. parsed-literal::

   improper_style harmonic

Examples
""""""""


.. parsed-literal::

   improper_style harmonic
   improper_coeff 1 100.0 0

Description
"""""""""""

The *harmonic* improper style uses the potential

.. image:: Eqs/improper_harmonic.jpg
   :align: center

where X is the improper angle, X0 is its equilibrium value, and K is a
prefactor.  Note that the usual 1/2 factor is included in K.

If the 4 atoms in an improper quadruplet (listed in the data file read
by the :doc:`read\_data <read_data>` command) are ordered I,J,K,L then X
is the angle between the plane of I,J,K and the plane of J,K,L.
Alternatively, you can think of atoms J,K,L as being in a plane, and
atom I above the plane, and X as a measure of how far out-of-plane I
is with respect to the other 3 atoms.

Note that defining 4 atoms to interact in this way, does not mean that
bonds necessarily exist between I-J, J-K, or K-L, as they would in a
linear dihedral.  Normally, the bonds I-J, I-K, I-L would exist for an
improper to be defined between the 4 atoms.

The following coefficients must be defined for each improper type via
the :doc:`improper\_coeff <improper_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands:

* K (energy/radian\^2)
* X0 (degrees)

X0 is specified in degrees, but LAMMPS converts it to radians
internally; hence the units of K are in energy/radian\^2.


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


Restrictions
""""""""""""


This improper style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`improper\_coeff <improper_coeff>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
