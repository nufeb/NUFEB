.. index:: fix wall/reflect

fix wall/reflect command
========================

fix wall/reflect/kk command
===========================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID wall/reflect face arg ... keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* wall/reflect = style name of this fix command
* one or more face/arg pairs may be appended
* face = *xlo* or *xhi* or *ylo* or *yhi* or *zlo* or *zhi*
  
  .. parsed-literal::
  
       *xlo*\ ,\ *ylo*\ ,\ *zlo* arg = EDGE or constant or variable
         EDGE = current lo edge of simulation box
         constant = number like 0.0 or -30.0 (distance units)
         variable = :doc:`equal-style variable <variable>` like v_x or v_wiggle
       *xhi*\ ,\ *yhi*\ ,\ *zhi* arg = EDGE or constant or variable
         EDGE = current hi edge of simulation box
         constant = number like 50.0 or 100.3 (distance units)
         variable = :doc:`equal-style variable <variable>` like v_x or v_wiggle

* zero or more keyword/value pairs may be appended
* keyword = *units*
  
  .. parsed-literal::
  
       *units* value = *lattice* or *box*
         *lattice* = the wall position is defined in lattice units
         *box* = the wall position is defined in simulation box units



Examples
""""""""


.. parsed-literal::

   fix xwalls all wall/reflect xlo EDGE xhi EDGE
   fix walls all wall/reflect xlo 0.0 ylo 10.0 units box
   fix top all wall/reflect zhi v_pressdown

Description
"""""""""""

Bound the simulation with one or more walls which reflect particles
in the specified group when they attempt to move through them.

Reflection means that if an atom moves outside the wall on a timestep
by a distance delta (e.g. due to :doc:`fix nve <fix_nve>`), then it is
put back inside the face by the same delta, and the sign of the
corresponding component of its velocity is flipped.

When used in conjunction with :doc:`fix nve <fix_nve>` and :doc:`run\_style verlet <run_style>`, the resultant time-integration algorithm is
equivalent to the primitive splitting algorithm (PSA) described by
:ref:`Bond <Bond1>`.  Because each reflection event divides
the corresponding timestep asymmetrically, energy conservation is only
satisfied to O(dt), rather than to O(dt\^2) as it would be for
velocity-Verlet integration without reflective walls.

Up to 6 walls or faces can be specified in a single command: *xlo*\ ,
*xhi*\ , *ylo*\ , *yhi*\ , *zlo*\ , *zhi*\ .  A *lo* face reflects particles
that move to a coordinate less than the wall position, back in the
*hi* direction.  A *hi* face reflects particles that move to a
coordinate higher than the wall position, back in the *lo* direction.

The position of each wall can be specified in one of 3 ways: as the
EDGE of the simulation box, as a constant value, or as a variable.  If
EDGE is used, then the corresponding boundary of the current
simulation box is used.  If a numeric constant is specified then the
wall is placed at that position in the appropriate dimension (x, y, or
z).  In both the EDGE and constant cases, the wall will never move.
If the wall position is a variable, it should be specified as v\_name,
where name is an :doc:`equal-style variable <variable>` name.  In this
case the variable is evaluated each timestep and the result becomes
the current position of the reflecting wall.  Equal-style variables
can specify formulas with various mathematical functions, and include
:doc:`thermo\_style <thermo_style>` command keywords for the simulation
box parameters and timestep and elapsed time.  Thus it is easy to
specify a time-dependent wall position.

The *units* keyword determines the meaning of the distance units used
to define a wall position, but only when a numeric constant or
variable is used.  It is not relevant when EDGE is used to specify a
face position.  In the variable case, the variable is assumed to
produce a value compatible with the *units* setting you specify.

A *box* value selects standard distance units as defined by the
:doc:`units <units>` command, e.g. Angstroms for units = real or metal.
A *lattice* value means the distance units are in lattice spacings.
The :doc:`lattice <lattice>` command must have been previously used to
define the lattice spacings.


----------


Here are examples of variable definitions that move the wall position
in a time-dependent fashion using equal-style
:doc:`variables <variable>`.


.. parsed-literal::

   variable ramp equal ramp(0,10)
   fix 1 all wall/reflect xlo v_ramp

   variable linear equal vdisplace(0,20)
   fix 1 all wall/reflect xlo v_linear

   variable wiggle equal swiggle(0.0,5.0,3.0)
   fix 1 all wall/reflect xlo v_wiggle

   variable wiggle equal cwiggle(0.0,5.0,3.0)
   fix 1 all wall/reflect xlo v_wiggle

The ramp(lo,hi) function adjusts the wall position linearly from lo to
hi over the course of a run.  The vdisplace(c0,velocity) function does
something similar using the equation position = c0 + velocity\*delta,
where delta is the elapsed time.

The swiggle(c0,A,period) function causes the wall position to
oscillate sinusoidally according to this equation, where omega = 2 PI
/ period:


.. parsed-literal::

   position = c0 + A sin(omega\*delta)

The cwiggle(c0,A,period) function causes the wall position to
oscillate sinusoidally according to this equation, which will have an
initial wall velocity of 0.0, and thus may impose a gentler
perturbation on the particles:


.. parsed-literal::

   position = c0 + A (1 - cos(omega\*delta))


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


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix\_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


Any dimension (xyz) that has a reflecting wall must be non-periodic.

A reflecting wall should not be used with rigid bodies such as those
defined by a "fix rigid" command.  This is because the wall/reflect
displaces atoms directly rather than exerts a force on them.  For
rigid bodies, use a soft wall instead, such as :doc:`fix wall/lj93 <fix_wall>`.  LAMMPS will flag the use of a rigid
fix with fix wall/reflect with a warning, but will not generate an
error.

Related commands
""""""""""""""""

:doc:`fix wall/lj93 <fix_wall>`, :doc:`fix oneway <fix_oneway>`

**Default:** none


----------


.. _Bond1:



**(Bond)** Bond and Leimkuhler, SIAM J Sci Comput, 30, p 134 (2007).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
