.. index:: fix npt/body

fix npt/body command
====================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID npt/body keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* npt/body = style name of this fix command
* additional thermostat and barostat related keyword/value pairs from the :doc:`fix npt <fix_nh>` command can be appended

Examples
""""""""


.. parsed-literal::

   fix 1 all npt/body temp 300.0 300.0 100.0 iso 0.0 0.0 1000.0
   fix 2 all npt/body temp 300.0 300.0 100.0 x 5.0 5.0 1000.0
   fix 2 all npt/body temp 300.0 300.0 100.0 x 5.0 5.0 1000.0 drag 0.2
   fix 2 water npt/body temp 300.0 300.0 100.0 aniso 0.0 0.0 1000.0 dilate partial

Description
"""""""""""

Perform constant NPT integration to update position, velocity,
orientation, and angular velocity each timestep for body
particles in the group using a Nose/Hoover temperature
thermostat and Nose/Hoover pressure barostat.  P is pressure; T is
temperature.  This creates a system trajectory consistent with the
isothermal-isobaric ensemble.

This fix differs from the :doc:`fix npt <fix_nh>` command, which
assumes point particles and only updates their position and velocity.

The thermostat is applied to both the translational and rotational
degrees of freedom for the body particles, assuming a compute is
used which calculates a temperature that includes the rotational
degrees of freedom (see below).  The translational degrees of freedom
can also have a bias velocity removed from them before thermostatting
takes place; see the description below.

Additional parameters affecting the thermostat and barostat are
specified by keywords and values documented with the :doc:`fix npt <fix_nh>` command.  See, for example, discussion of the *temp*\ ,
*iso*\ , *aniso*\ , and *dilate* keywords.

The particles in the fix group are the only ones whose velocities and
positions are updated by the velocity/position update portion of the
NPT integration.

Regardless of what particles are in the fix group, a global pressure is
computed for all particles.  Similarly, when the size of the simulation
box is changed, all particles are re-scaled to new positions, unless the
keyword *dilate* is specified with a value of *partial*\ , in which case
only the particles in the fix group are re-scaled.  The latter can be
useful for leaving the coordinates of particles in a solid substrate
unchanged and controlling the pressure of a surrounding fluid.


----------


This fix computes a temperature and pressure each timestep.  To do
this, the fix creates its own computes of style "temp/body" and
"pressure", as if these commands had been issued:


.. parsed-literal::

   compute fix-ID_temp all temp/body
   compute fix-ID_press all pressure fix-ID_temp

See the :doc:`compute temp/body <compute_temp_body>` and :doc:`compute pressure <compute_pressure>` commands for details.  Note that the
IDs of the new computes are the fix-ID + underscore + "temp" or fix\_ID
+ underscore + "press", and the group for the new computes is "all"
since pressure is computed for the entire system.

Note that these are NOT the computes used by thermodynamic output (see
the :doc:`thermo\_style <thermo_style>` command) with ID = *thermo\_temp*
and *thermo\_press*.  This means you can change the attributes of this
fix's temperature or pressure via the
:doc:`compute\_modify <compute_modify>` command or print this temperature
or pressure during thermodynamic output via the :doc:`thermo\_style custom <thermo_style>` command using the appropriate compute-ID.
It also means that changing attributes of *thermo\_temp* or
*thermo\_press* will have no effect on this fix.

Like other fixes that perform thermostatting, this fix can be used
with :doc:`compute commands <compute>` that calculate a temperature
after removing a "bias" from the atom velocities.  E.g. removing the
center-of-mass velocity from a group of atoms or only calculating
temperature on the x-component of velocity or only calculating
temperature for atoms in a geometric region.  This is not done by
default, but only if the :doc:`fix\_modify <fix_modify>` command is used
to assign a temperature compute to this fix that includes such a bias
term.  See the doc pages for individual :doc:`compute commands <compute>` to determine which ones include a bias.  In
this case, the thermostat works in the following manner: the current
temperature is calculated taking the bias into account, bias is
removed from each atom, thermostatting is performed on the remaining
thermal degrees of freedom, and the bias is added back in.


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

**Restart, fix\_modify, output, run start/stop, minimize info:**

This fix writes the state of the Nose/Hoover thermostat and barostat
to :doc:`binary restart files <restart>`.  See the
:doc:`read\_restart <read_restart>` command for info on how to re-specify
a fix in an input script that reads a restart file, so that the
operation of the fix continues in an uninterrupted fashion.

The :doc:`fix\_modify <fix_modify>` *temp* and *press* options are
supported by this fix.  You can use them to assign a
:doc:`compute <compute>` you have defined to this fix which will be used
in its thermostatting or barostatting procedure.  If you do this, note
that the kinetic energy derived from the compute temperature should be
consistent with the virial term computed using all atoms for the
pressure.  LAMMPS will warn you if you choose to compute temperature
on a subset of atoms.

The :doc:`fix\_modify <fix_modify>` *energy* option is supported by this
fix to add the energy change induced by Nose/Hoover thermostatting and
barostatting to the system's potential energy as part of
:doc:`thermodynamic output <thermo_style>`.

This fix computes the same global scalar and global vector of
quantities as does the :doc:`fix npt <fix_nh>` command.

This fix can ramp its target temperature and pressure over multiple
runs, using the *start* and *stop* keywords of the :doc:`run <run>`
command.  See the :doc:`run <run>` command for details of how to do
this.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the BODY package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

This fix requires that atoms store torque and angular momentum and a
quaternion as defined by the :doc:`atom\_style body <atom_style>`
command.

Related commands
""""""""""""""""

:doc:`fix npt <fix_nh>`, :doc:`fix nve\_body <fix_nve_body>`, :doc:`fix nvt\_body <fix_nvt_body>`, :doc:`fix\_modify <fix_modify>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
