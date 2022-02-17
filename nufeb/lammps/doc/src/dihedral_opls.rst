.. index:: dihedral\_style opls

dihedral\_style opls command
============================

dihedral\_style opls/intel command
==================================

dihedral\_style opls/kk command
===============================

dihedral\_style opls/omp command
================================

Syntax
""""""


.. parsed-literal::

   dihedral_style opls

Examples
""""""""


.. parsed-literal::

   dihedral_style opls
   dihedral_coeff 1 1.740 -0.157 0.279 0.00   # CT-CT-CT-CT
   dihedral_coeff 2 0.000 0.000 0.366 0.000   # CT-CT-CT-HC
   dihedral_coeff 3 0.000 0.000 0.318 0.000   # HC-CT-CT-HC

Description
"""""""""""

The *opls* dihedral style uses the potential

.. image:: Eqs/dihedral_opls.jpg
   :align: center

Note that the usual 1/2 factor is not included in the K values.

This dihedral potential is used in the OPLS force field and is
described in :ref:`(Watkins) <Watkins>`.

The following coefficients must be defined for each dihedral type via the
:doc:`dihedral\_coeff <dihedral_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* K1 (energy)
* K2 (energy)
* K3 (energy)
* K4 (energy)


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


This dihedral style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`dihedral\_coeff <dihedral_coeff>`

**Default:** none


----------


.. _Watkins:



**(Watkins)** Watkins and Jorgensen, J Phys Chem A, 105, 4118-4125 (2001).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
