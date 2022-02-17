.. index:: angle_style cosine/periodic

angle_style cosine/periodic command
===================================

angle_style cosine/periodic/omp command
=======================================

Syntax
""""""


.. code-block:: LAMMPS

   angle_style cosine/periodic

Examples
""""""""


.. code-block:: LAMMPS

   angle_style cosine/periodic
   angle_coeff * 75.0 1 6

Description
"""""""""""

The *cosine/periodic* angle style uses the following potential, which
is commonly used in the :doc:`DREIDING <Howto_bioFF>` force field,
particularly for organometallic systems where :math:`n` = 4 might be used
for an octahedral complex and :math:`n` = 3 might be used for a trigonal
center:

.. math::

   E = C \left[ 1 - B(-1)^n\cos\left( n\theta\right) \right]


where :math:`C`, :math:`B` and :math:`n` are coefficients defined for each angle type.

See :ref:`(Mayo) <cosine-Mayo>` for a description of the DREIDING force field

The following coefficients must be defined for each angle type via the
:doc:`angle\_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* :math:`C` (energy)
* :math:`B` = 1 or -1
* :math:`n` = 1, 2, 3, 4, 5 or 6 for periodicity

Note that the prefactor :math:`C` is specified and not the overall force
constant :math:`K = \frac{C}{n^2}`.  When :math:`B = 1`, it leads to a minimum for the
linear geometry.  When :math:`B = -1`, it leads to a maximum for the linear
geometry.


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


This angle style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`angle\_coeff <angle_coeff>`

**Default:** none


----------


.. _cosine-Mayo:



**(Mayo)** Mayo, Olfason, Goddard III, J Phys Chem, 94, 8897-8909
(1990).
