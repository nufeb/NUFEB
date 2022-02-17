.. index:: pair\_style lj/smooth/linear

pair\_style lj/smooth/linear command
====================================

pair\_style lj/smooth/linear/omp command
========================================

Syntax
""""""


.. parsed-literal::

   pair_style lj/smooth/linear cutoff

* cutoff = global cutoff for Lennard-Jones interactions (distance units)

Examples
""""""""


.. parsed-literal::

   pair_style lj/smooth/linear 2.5
   pair_coeff \* \* 1.0 1.0
   pair_coeff 1 1 0.3 3.0 9.0

Description
"""""""""""

Style *lj/smooth/linear* computes a truncated and force-shifted LJ
interaction (aka Shifted Force Lennard-Jones) that combines the
standard 12/6 Lennard-Jones function and subtracts a linear term based
on the cutoff distance, so that both, the potential and the force, go
continuously to zero at the cutoff Rc :ref:`(Toxvaerd) <Toxvaerd>`:

.. image:: Eqs/pair_lj_smooth_linear.jpg
   :align: center

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair\_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands, or by mixing as described below:

* epsilon (energy units)
* sigma (distance units)
* cutoff (distance units)

The last coefficient is optional. If not specified, the global
LJ cutoff specified in the pair\_style command is used.


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


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance can be mixed. The default mix value is geometric.
See the "pair\_modify" command for details.

This pair style does not support the :doc:`pair\_modify <pair_modify>`
shift option for the energy of the pair interaction, since it goes
to 0.0 at the cutoff by construction.

The :doc:`pair\_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style does not support the :doc:`pair\_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure, since the energy of the pair interaction is smoothed to 0.0
at the cutoff.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run\_style respa <run_style>` command.  It does not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`, :doc:`pair lj/smooth <pair_lj_smooth>`

**Default:** none


----------


.. _Toxvaerd:



**(Toxvaerd)** Toxvaerd, Dyre, J Chem Phys, 134, 081102 (2011).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
