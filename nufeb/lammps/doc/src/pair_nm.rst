.. index:: pair\_style nm/cut

pair\_style nm/cut command
==========================

pair\_style nm/cut/coul/cut command
===================================

pair\_style nm/cut/coul/long command
====================================

pair\_style nm/cut/omp command
==============================

pair\_style nm/cut/coul/cut/omp command
=======================================

pair\_style nm/cut/coul/long/omp command
========================================

Syntax
""""""


.. parsed-literal::

   pair_style style args

* style = *nm/cut* or *nm/cut/coul/cut* or *nm/cut/coul/long*
* args = list of arguments for a particular style
  
  .. parsed-literal::
  
       *nm/cut* args = cutoff
         cutoff = global cutoff for Pair interactions (distance units)
       *nm/cut/coul/cut* args = cutoff (cutoff2)
         cutoff = global cutoff for Pair (and Coulombic if only 1 arg) (distance units)
         cutoff2 = global cutoff for Coulombic (optional) (distance units)
       *nm/cut/coul/long* args = cutoff (cutoff2)
         cutoff = global cutoff for Pair (and Coulombic if only 1 arg) (distance units)
         cutoff2 = global cutoff for Coulombic (optional) (distance units)



Examples
""""""""


.. parsed-literal::

   pair_style nm/cut 12.0
   pair_coeff \* \* 0.01 5.4 8.0 7.0
   pair_coeff 1 1 0.01 4.4 7.0 6.0

   pair_style nm/cut/coul/cut 12.0 15.0
   pair_coeff \* \* 0.01 5.4 8.0 7.0
   pair_coeff 1 1 0.01 4.4 7.0 6.0

   pair_style nm/cut/coul/long 12.0 15.0
   pair_coeff \* \* 0.01 5.4 8.0 7.0
   pair_coeff 1 1 0.01 4.4 7.0 6.0

Description
"""""""""""

Style *nm* computes site-site interactions based on the N-M potential
by :ref:`Clarke <Clarke>`, mainly used for ionic liquids.  A site can
represent a single atom or a united-atom site.  The energy of an
interaction has the following form:

.. image:: Eqs/pair_nm.jpg
   :align: center

Rc is the cutoff.

Style *nm/cut/coul/cut* adds a Coulombic pairwise interaction given by

.. image:: Eqs/pair_coulomb.jpg
   :align: center

where C is an energy-conversion constant, Qi and Qj are the charges on
the 2 atoms, and epsilon is the dielectric constant which can be set
by the :doc:`dielectric <dielectric>` command.  If one cutoff is
specified in the pair\_style command, it is used for both the NM and
Coulombic terms.  If two cutoffs are specified, they are used as
cutoffs for the NM and Coulombic terms respectively.

Styles *nm/cut/coul/long* compute the same
Coulombic interactions as style *nm/cut/coul/cut* except that an
additional damping factor is applied to the Coulombic term so it can
be used in conjunction with the :doc:`kspace\_style <kspace_style>`
command and its *ewald* or *pppm* option.  The Coulombic cutoff
specified for this style means that pairwise interactions within this
distance are computed directly; interactions outside that distance are
computed in reciprocal space.

For all of the *nm* pair styles, the following coefficients must
be defined for each pair of atoms types
via the :doc:`pair\_coeff <pair_coeff>` command as in the
examples above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands.

* E0 (energy units)
* r0 (distance units)
* n (unitless)
* m (unitless)
* cutoff1 (distance units)
* cutoff2 (distance units)

The latter 2 coefficients are optional.  If not specified, the global
NM and Coulombic cutoffs specified in the pair\_style command are used.
If only one cutoff is specified, it is used as the cutoff for both NM
and Coulombic interactions for this type pair.  If both coefficients
are specified, they are used as the NM and Coulombic cutoffs for this
type pair.  You cannot specify 2 cutoffs for style *nm*\ , since it
has no Coulombic terms.

For *nm/cut/coul/long* only the NM cutoff can be specified since a
Coulombic cutoff cannot be specified for an individual I,J type pair.
All type pairs use the same global Coulombic cutoff specified in the
pair\_style command.


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

These pair styles do not support mixing. Thus, coefficients for all
I,J pairs must be specified explicitly.

All of the *nm* pair styles supports the
:doc:`pair\_modify <pair_modify>` shift option for the energy of the pair
interaction.

The *nm/cut/coul/long* pair styles support the
:doc:`pair\_modify <pair_modify>` table option since they can tabulate
the short-range portion of the long-range Coulombic interaction.

All of the *nm* pair styles support the :doc:`pair\_modify <pair_modify>`
tail option for adding a long-range tail correction to the energy and
pressure for the NM portion of the pair interaction.

All of the *nm* pair styles write their information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

All of the *nm* pair styles can only be used via the *pair* keyword of
the :doc:`run\_style respa <run_style>` command.  They do not support the
*inner*\ , *middle*\ , *outer* keywords.


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

Restrictions
""""""""""""


These pair styles are part of the MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`

**Default:** none


----------


.. _Clarke:



**(Clarke)** Clarke and Smith, J Chem Phys, 84, 2290 (1986).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
