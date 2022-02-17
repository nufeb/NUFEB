.. index:: pair\_style morse

pair\_style morse command
=========================

pair\_style morse/gpu command
=============================

pair\_style morse/omp command
=============================

pair\_style morse/opt command
=============================

pair\_style morse/smooth/linear command
=======================================

pair\_style morse/smooth/linear/omp command
===========================================

pair\_style morse/kk command
============================

Syntax
""""""


.. parsed-literal::

   pair_style style args

* style = *morse* or *morse/smooth/linear* or *morse/soft*
* args = list of arguments for a particular style


.. parsed-literal::

    *morse* args = cutoff
      cutoff = global cutoff for Morse interactions (distance units)
    *morse/smooth/linear* args = cutoff
      cutoff = global cutoff for Morse interactions (distance units)

Examples
""""""""


.. parsed-literal::

   pair_style morse 2.5
   pair_style morse/smooth/linear 2.5
   pair_coeff \* \* 100.0 2.0 1.5
   pair_coeff 1 1 100.0 2.0 1.5 3.0

Description
"""""""""""

Style *morse* computes pairwise interactions with the formula

.. image:: Eqs/pair_morse.jpg
   :align: center

Rc is the cutoff.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair\_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands:

* D0 (energy units)
* alpha (1/distance units)
* r0 (distance units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global morse
cutoff is used.


----------


The *morse/smooth/linear* variant is similar to the lj/smooth/linear
variant in that it adds to the potential a shift and a linear term
so that both, potential energy and force, go to zero at the cut-off:

.. image:: Eqs/pair_morse_smooth_linear.jpg
   :align: center

The syntax of the pair\_style and pair\_coeff commands are the same for
the *morse* and *morse/smooth/linear* styles.


----------


A version of the *morse* style with a soft core, *morse/soft*\ , suitable for use in
free energy calculations, is part of the USER-FEP package and is documented with
the :doc:`pair\_fep\_soft <pair_fep_soft>` styles. The version with soft core is only
available if LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.


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

None of these pair styles support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

All of these pair styles support the :doc:`pair\_modify <pair_modify>`
shift option for the energy of the pair interaction.

The :doc:`pair\_modify <pair_modify>` table options is not relevant for
the Morse pair styles.

None of these pair styles support the :doc:`pair\_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

All of these pair styles write their information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
:doc:`run\_style respa <run_style>` command.  They do not support the
*inner*\ , *middle*\ , *outer* keywords.


----------


Restrictions
""""""""""""


The *morse/smooth/linear* pair style is only enabled if LAMMPS was
built with the USER-MISC package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`, :doc:`pair\_fep\_soft <pair_fep_soft>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
