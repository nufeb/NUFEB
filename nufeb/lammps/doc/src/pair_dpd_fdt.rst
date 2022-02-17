.. index:: pair\_style dpd/fdt

pair\_style dpd/fdt command
===========================

pair\_style dpd/fdt/energy command
==================================

pair\_style dpd/fdt/energy/kk command
=====================================

Syntax
""""""


.. parsed-literal::

   pair_style style args

* style = *dpd/fdt* or *dpd/fdt/energy*
* args = list of arguments for a particular style


.. parsed-literal::

     *dpd/fdt* args = T cutoff seed
       T = temperature (temperature units)
       cutoff = global cutoff for DPD interactions (distance units)
       seed = random # seed (positive integer)
     *dpd/fdt/energy* args = cutoff seed
       cutoff = global cutoff for DPD interactions (distance units)
       seed = random # seed (positive integer)

Examples
""""""""


.. parsed-literal::

   pair_style dpd/fdt 300.0 2.5 34387
   pair_coeff \* \* 3.0 1.0 2.5

   pair_style dpd/fdt/energy 2.5 34387
   pair_coeff \* \* 3.0 1.0 0.1 2.5

Description
"""""""""""

Styles *dpd/fdt* and *dpd/fdt/energy* compute the force for dissipative
particle dynamics (DPD) simulations.  The *dpd/fdt* style is used to
perform DPD simulations under isothermal and isobaric conditions,
while the *dpd/fdt/energy* style is used to perform DPD simulations
under isoenergetic and isoenthalpic conditions (see :ref:`(Lisal) <Lisal3>`).
For DPD simulations in general, the force on atom I due to atom J is
given as a sum of 3 terms

.. image:: Eqs/pair_dpd.jpg
   :align: center

where Fc is a conservative force, Fd is a dissipative force, and Fr is
a random force.  Rij is a unit vector in the direction Ri - Rj, Vij is
the vector difference in velocities of the two atoms = Vi - Vj, alpha
is a Gaussian random number with zero mean and unit variance, dt is
the timestep size, and w(r) is a weighting factor that varies between
0 and 1.  Rc is the cutoff.  The weighting factor, omega\_ij, varies
between 0 and 1, and is chosen to have the following functional form:

.. image:: Eqs/pair_dpd_omega.jpg
   :align: center

Note that alternative definitions of the weighting function exist, but
would have to be implemented as a separate pair style command.

For style *dpd/fdt*\ , the fluctuation-dissipation theorem defines gamma
to be set equal to sigma\*sigma/(2 T), where T is the set point
temperature specified as a pair style parameter in the above examples.
The following coefficients must be defined for each pair of atoms types
via the :doc:`pair\_coeff <pair_coeff>` command as in the examples above,
or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>` commands:

* A (force units)
* sigma (force\*time\^(1/2) units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global DPD
cutoff is used.

Style *dpd/fdt/energy* is used to perform DPD simulations
under isoenergetic and isoenthalpic conditions.  The fluctuation-dissipation
theorem defines gamma to be set equal to sigma\*sigma/(2 dpdTheta), where
dpdTheta is the average internal temperature for the pair. The particle
internal temperature is related to the particle internal energy through
a mesoparticle equation of state (see :doc:`fix eos <fix>`). The
differential internal conductive and mechanical energies are computed
within style *dpd/fdt/energy* as:

.. image:: Eqs/pair_dpd_energy.jpg
   :align: center

where

.. image:: Eqs/pair_dpd_energy_terms.jpg
   :align: center

Zeta\_ij\^q is a second Gaussian random number with zero mean and unit
variance that is used to compute the internal conductive energy. The
fluctuation-dissipation theorem defines alpha\*alpha to be set
equal to 2\*kB\*kappa, where kappa is the mesoparticle thermal
conductivity parameter.   The following coefficients must be defined for
each pair of atoms types via the :doc:`pair\_coeff <pair_coeff>`
command as in the examples above, or in the data file or restart files
read by the :doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands:

* A (force units)
* sigma (force\*time\^(1/2) units)
* kappa (energy\*temperature/time units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global DPD
cutoff is used.

The pairwise energy associated with styles *dpd/fdt* and
*dpd/fdt/energy* is only due to the conservative force term Fc, and is
shifted to be zero at the cutoff distance Rc.  The pairwise virial is
calculated using only the conservative term.

The forces computed through the *dpd/fdt* and *dpd/fdt/energy* styles
can be integrated with the velocity-Verlet integration scheme or the
Shardlow splitting integration scheme described by :ref:`(Lisal) <Lisal3>`.
In the cases when these pair styles are combined with the
:doc:`fix shardlow <fix_shardlow>`, these pair styles differ from the
other dpd styles in that the dissipative and random forces are split
from the force calculation and are not computed within the pair style.
Thus, only the conservative force is computed by the pair style,
while the stochastic integration of the dissipative and random forces
are handled through the Shardlow splitting algorithm approach.  The
Shardlow splitting algorithm is advantageous, especially when
performing DPD under isoenergetic conditions, as it allows
significantly larger timesteps to be taken.


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


These commands are part of the USER-DPD package.  They are only
enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Pair styles *dpd/fdt* and *dpd/fdt/energy* require use of the
:doc:`comm\_modify vel yes <comm_modify>` option so that velocities are
stored by ghost atoms.

Pair style *dpd/fdt/energy* requires :doc:`atom\_style dpd <atom_style>`
to be used in order to properly account for the particle internal
energies and temperatures.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`, :doc:`fix shardlow <fix_shardlow>`

**Default:** none


----------


.. _Lisal3:



**(Lisal)** M. Lisal, J.K. Brennan, J. Bonet Avalos, "Dissipative
particle dynamics at isothermal, isobaric, isoenergetic, and
isoenthalpic conditions using Shardlow-like splitting algorithms.",
J. Chem. Phys., 135, 204105 (2011).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
