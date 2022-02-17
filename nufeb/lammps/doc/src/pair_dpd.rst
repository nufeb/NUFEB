.. index:: pair\_style dpd

pair\_style dpd command
=======================

pair\_style dpd/gpu command
===========================

pair\_style dpd/intel command
=============================

pair\_style dpd/omp command
===========================

pair\_style dpd/tstat command
=============================

pair\_style dpd/tstat/gpu command
=================================

pair\_style dpd/tstat/omp command
=================================

Syntax
""""""


.. parsed-literal::

   pair_style dpd T cutoff seed
   pair_style dpd/tstat Tstart Tstop cutoff seed

* T = temperature (temperature units)
* Tstart,Tstop = desired temperature at start/end of run (temperature units)
* cutoff = global cutoff for DPD interactions (distance units)
* seed = random # seed (positive integer)

Examples
""""""""


.. parsed-literal::

   pair_style dpd 1.0 2.5 34387
   pair_coeff \* \* 3.0 1.0
   pair_coeff 1 1 3.0 1.0 1.0

   pair_style dpd/tstat 1.0 1.0 2.5 34387
   pair_coeff \* \* 1.0
   pair_coeff 1 1 1.0 1.0

Description
"""""""""""

Style *dpd* computes a force field for dissipative particle dynamics
(DPD) following the exposition in :ref:`(Groot) <Groot1>`.

Style *dpd/tstat* invokes a DPD thermostat on pairwise interactions,
which is equivalent to the non-conservative portion of the DPD force
field.  This pair-wise thermostat can be used in conjunction with any
:doc:`pair style <pair_style>`, and in leiu of per-particle thermostats
like :doc:`fix langevin <fix_langevin>` or ensemble thermostats like
Nose Hoover as implemented by :doc:`fix nvt <fix_nh>`.  To use
*dpd/tstat* as a thermostat for another pair style, use the :doc:`pair\_style hybrid/overlay <pair_hybrid>` command to compute both the desired
pair interaction and the thermostat for each pair of particles.

For style *dpd*\ , the force on atom I due to atom J is given as a sum
of 3 terms

.. image:: Eqs/pair_dpd.jpg
   :align: center

where Fc is a conservative force, Fd is a dissipative force, and Fr is
a random force.  Rij is a unit vector in the direction Ri - Rj, Vij is
the vector difference in velocities of the two atoms = Vi - Vj, alpha
is a Gaussian random number with zero mean and unit variance, dt is
the timestep size, and w(r) is a weighting factor that varies between
0 and 1.  Rc is the cutoff.  Sigma is set equal to sqrt(2 Kb T gamma),
where Kb is the Boltzmann constant and T is the temperature parameter
in the pair\_style command.

For style *dpd/tstat*\ , the force on atom I due to atom J is the same
as the above equation, except that the conservative Fc term is
dropped.  Also, during the run, T is set each timestep to a ramped
value from Tstart to Tstop.

For style *dpd*\ , the pairwise energy associated with style *dpd* is
only due to the conservative force term Fc, and is shifted to be zero
at the cutoff distance Rc.  The pairwise virial is calculated using
all 3 terms.  For style *dpd/tstat* there is no pairwise energy, but
the last two terms of the formula make a contribution to the virial.

For style *dpd*\ , the following coefficients must be defined for each
pair of atoms types via the :doc:`pair\_coeff <pair_coeff>` command as in
the examples above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands:

* A (force units)
* gamma (force/velocity units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global DPD
cutoff is used.  Note that sigma is set equal to sqrt(2 T gamma),
where T is the temperature set by the :doc:`pair\_style <pair_style>`
command so it does not need to be specified.

For style *dpd/tstat*\ , the coefficients defined for each pair of
atoms types via the :doc:`pair\_coeff <pair_coeff>` command is the same,
except that A is not included.

The GPU-accelerated versions of these styles are implemented based on
the work of :ref:`(Afshar) <Afshar>` and :ref:`(Phillips) <Phillips>`.

.. note::

   If you are modeling DPD polymer chains, you may want to use the
   :doc:`pair\_style srp <pair_srp>` command in conjunction with these pair
   styles.  It is a soft segmental repulsive potential (SRP) that can
   prevent DPD polymer chains from crossing each other.

.. note::

   The virial calculation for pressure when using this pair style
   includes all the components of force listed above, including the
   random force.


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

These pair styles do not support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

These pair styles do not support the :doc:`pair\_modify <pair_modify>`
shift option for the energy of the pair interaction.  Note that as
discussed above, the energy due to the conservative Fc term is already
shifted to be 0.0 at the cutoff distance Rc.

The :doc:`pair\_modify <pair_modify>` table option is not relevant
for these pair styles.

These pair style do not support the :doc:`pair\_modify <pair_modify>`
tail option for adding long-range tail corrections to energy and
pressure.

These pair styles writes their information to :doc:`binary restart files <restart>`, so pair\_style and pair\_coeff commands do not need
to be specified in an input script that reads a restart file.  Note
that the user-specified random number seed is stored in the restart
file, so when a simulation is restarted, each processor will
re-initialize its random number generator the same way it did
initially.  This means the random forces will be random, but will not
be the same as they would have been if the original simulation had
continued past the restart time.

These pair styles can only be used via the *pair* keyword of the
:doc:`run\_style respa <run_style>` command.  They do not support the
*inner*\ , *middle*\ , *outer* keywords.

The *dpd/tstat* style can ramp its target temperature over multiple
runs, using the *start* and *stop* keywords of the :doc:`run <run>`
command.  See the :doc:`run <run>` command for details of how to do
this.


----------


Restrictions
""""""""""""


The default frequency for rebuilding neighbor lists is every 10 steps
(see the :doc:`neigh\_modify <neigh_modify>` command). This may be too
infrequent for style *dpd* simulations since particles move rapidly
and can overlap by large amounts.  If this setting yields a non-zero
number of "dangerous" reneighborings (printed at the end of a
simulation), you should experiment with forcing reneighboring more
often and see if system energies/trajectories change.

These pair styles requires you to use the :doc:`comm\_modify vel yes <comm_modify>` command so that velocities are stored by ghost
atoms.

These pair styles will not restart exactly when using the
:doc:`read\_restart <read_restart>` command, though they should provide
statistically similar results.  This is because the forces they
compute depend on atom velocities.  See the
:doc:`read\_restart <read_restart>` command for more details.

Related commands
""""""""""""""""

:doc:`pair\_coeff <pair_coeff>`, :doc:`fix nvt <fix_nh>`, :doc:`fix langevin <fix_langevin>`, :doc:`pair\_style srp <pair_srp>`

**Default:** none


----------


.. _Groot1:



**(Groot)** Groot and Warren, J Chem Phys, 107, 4423-35 (1997).

.. _Afshar:



**(Afshar)** Afshar, F. Schmid, A. Pishevar, S. Worley, Comput Phys
Comm, 184, 1119-1128 (2013).

.. _Phillips:



**(Phillips)** C. L. Phillips, J. A. Anderson, S. C. Glotzer, Comput
Phys Comm, 230, 7191-7201 (2011).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
