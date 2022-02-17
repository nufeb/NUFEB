.. index:: pair\_style kolmogorov/crespi/z

pair\_style kolmogorov/crespi/z command
=======================================

Syntax
""""""


.. parsed-literal::

   pair_style [hybrid/overlay ...] kolmogorov/crespi/z cutoff

Examples
""""""""


.. parsed-literal::

   pair_style hybrid/overlay kolmogorov/crespi/z 20.0
   pair_coeff \* \* none
   pair_coeff 1 2 kolmogorov/crespi/z  CC.KC   C C

   pair_style hybrid/overlay rebo kolmogorov/crespi/z 14.0
   pair_coeff \* \* rebo                 CH.rebo    C C
   pair_coeff 1 2 kolmogorov/crespi/z  CC.KC      C C

Description
"""""""""""

The *kolmogorov/crespi/z* style computes the Kolmogorov-Crespi interaction
potential as described in :ref:`(Kolmogorov) <KC05>`. An important simplification is made,
which is to take all normals along the z-axis.

.. image:: Eqs/pair_kolmogorov_crespi_z.jpg
   :align: center

It is important to have a sufficiently large cutoff to ensure smooth forces.
Energies are shifted so that they go continuously to zero at the cutoff assuming
that the exponential part of *Vij* (first term) decays sufficiently fast.
This shift is achieved by the last term in the equation for *Vij* above.

This potential is intended for interactions between two layers of graphene.
Therefore, to avoid interaction between layers in multi-layered materials,
each layer should have a separate atom type and interactions should only
be computed between atom types of neighboring layers.

The parameter file (e.g. CC.KC), is intended for use with metal
:doc:`units <units>`, with energies in meV. An additional parameter, *S*\ ,
is available to facilitate scaling of energies in accordance with
:ref:`(vanWijk) <vanWijk>`.

This potential must be used in combination with hybrid/overlay.
Other interactions can be set to zero using pair\_style *none*\ .

Restrictions
""""""""""""


This fix is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`,
:doc:`pair\_none <pair_none>`,
:doc:`pair\_style hybrid/overlay <pair_hybrid>`,
:doc:`pair\_style drip <pair_drip>`,
:doc:`pair\_style ilp/graphene/hbn <pair_ilp_graphene_hbn>`.
:doc:`pair\_style kolmogorov/crespi/full <pair_kolmogorov_crespi_full>`,
:doc:`pair\_style lebedeva/z <pair_lebedeva_z>`

**Default:** none


----------


.. _KC05:



**(Kolmogorov)** A. N. Kolmogorov, V. H. Crespi, Phys. Rev. B 71, 235415 (2005)

.. _vanWijk:



**(vanWijk)** M. M. van Wijk, A. Schuring, M. I. Katsnelson, and A. Fasolino,
Physical Review Letters, 113, 135504 (2014)


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
