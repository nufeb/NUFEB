.. index:: pair\_style lj/expand

pair\_style lj/expand command
=============================

pair\_style lj/expand/gpu command
=================================

pair\_style lj/expand/kk command
================================

pair\_style lj/expand/omp command
=================================

pair\_style lj/expand/coul/long command
=======================================

pair\_style lj/expand/coul/long/gpu command
===========================================

Syntax
""""""


.. parsed-literal::

   pair_style lj/expand cutoff

* cutoff = global cutoff for lj/expand interactions (distance units)

Examples
""""""""


.. parsed-literal::

   pair_style lj/expand 2.5
   pair_coeff \* \* 1.0 1.0 0.5
   pair_coeff 1 1 1.0 1.0 -0.2 2.0

   pair_style lj/expand/coul/long 2.5
   pair_style lj/expand/coul/long 2.5 4.0
   pair_coeff \* \* 1.0 1.0 0.5
   pair_coeff 1 1 1.0 1.0 -0.2 3.0

Description
"""""""""""

Style *lj/expand* computes a LJ interaction with a distance shifted by
delta which can be useful when particles are of different sizes, since
it is different that using different sigma values in a standard LJ
formula:

.. image:: Eqs/pair_lj_expand.jpg
   :align: center

Rc is the cutoff which does not include the delta distance.  I.e. the
actual force cutoff is the sum of cutoff + delta.

For all of the *lj/expand* pair styles, the following coefficients must
be defined for each pair of atoms types via the
:doc:`pair\_coeff <pair_coeff>` command as in the examples above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands, or by mixing as
described below:

* epsilon (energy units)
* sigma (distance units)
* delta (distance units)
* cutoff (distance units)

The delta values can be positive or negative.  The last coefficient is
optional.  If not specified, the global LJ cutoff is used.

For *lj/expand/coul/long* only the LJ cutoff can be specified since a
Coulombic cutoff cannot be specified for an individual I,J type pair.
All type pairs use the same global Coulombic cutoff specified in the
pair\_style command.


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

For atom type pairs I,J and I != J, the epsilon, sigma, and shift
coefficients and cutoff distance for this pair style can be mixed.
Shift is always mixed via an *arithmetic* rule.  The other
coefficients are mixed according to the pair\_modify mix value.  The
default mix value is *geometric*\ .  See the "pair\_modify" command for
details.

This pair style supports the :doc:`pair\_modify <pair_modify>` shift
option for the energy of the pair interaction.

The :doc:`pair\_modify <pair_modify>` table option is not relevant
for this pair style.

This pair style supports the :doc:`pair\_modify <pair_modify>` tail
option for adding a long-range tail correction to the energy and
pressure of the pair interaction.

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

:doc:`pair\_coeff <pair_coeff>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
