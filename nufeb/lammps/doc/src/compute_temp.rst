.. index:: compute temp

compute temp command
====================

compute temp/kk command
=======================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID temp

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp = style name of this compute command

Examples
""""""""


.. parsed-literal::

   compute 1 all temp
   compute myTemp mobile temp

Description
"""""""""""

Define a computation that calculates the temperature of a group of
atoms.  A compute of this style can be used by any command that
computes a temperature, e.g. :doc:`thermo\_modify <thermo_modify>`, :doc:`fix temp/rescale <fix_temp_rescale>`, :doc:`fix npt <fix_nh>`, etc.

The temperature is calculated by the formula KE = dim/2 N k T, where
KE = total kinetic energy of the group of atoms (sum of 1/2 m v\^2),
dim = 2 or 3 = dimensionality of the simulation, N = number of atoms
in the group, k = Boltzmann constant, and T = temperature.

A kinetic energy tensor, stored as a 6-element vector, is also
calculated by this compute for use in the computation of a pressure
tensor.  The formula for the components of the tensor is the same as
the above formula, except that v\^2 is replaced by vx\*vy for the xy
component, etc.  The 6 components of the vector are ordered xx, yy,
zz, xy, xz, yz.

The number of atoms contributing to the temperature is assumed to be
constant for the duration of the run; use the *dynamic* option of the
:doc:`compute\_modify <compute_modify>` command if this is not the case.

This compute subtracts out degrees-of-freedom due to fixes that
constrain molecular motion, such as :doc:`fix shake <fix_shake>` and
:doc:`fix rigid <fix_rigid>`.  This means the temperature of groups of
atoms that include these constraints will be computed correctly.  If
needed, the subtracted degrees-of-freedom can be altered using the
*extra* option of the :doc:`compute\_modify <compute_modify>` command.

A compute of this style with the ID of "thermo\_temp" is created when
LAMMPS starts up, as if this command were in the input script:


.. parsed-literal::

   compute thermo_temp all temp

See the "thermo\_style" command for more details.

See the :doc:`Howto thermostat <Howto_thermostat>` doc page for a
discussion of different ways to compute temperature and perform
thermostatting.


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


**Output info:**

This compute calculates a global scalar (the temperature) and a global
vector of length 6 (KE tensor), which can be accessed by indices 1-6.
These values can be used by any command that uses global scalar or
vector values from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options.

The scalar value calculated by this compute is "intensive".  The
vector values are "extensive".

The scalar value will be in temperature :doc:`units <units>`.  The
vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute temp/partial <compute_temp_partial>`, :doc:`compute temp/region <compute_temp_region>`, :doc:`compute pressure <compute_pressure>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
