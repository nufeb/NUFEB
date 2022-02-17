.. index:: fix hyper/local

fix hyper/local command
=======================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID hyper/local cutbond qfactor Vmax Tequil Dcut alpha Btarget

* ID, group-ID are documented in :doc:`fix <fix>` command
* hyper/local = style name of this fix command
* cutbond = max distance at which a pair of atoms is considered bonded (distance units)
* qfactor = max strain at which bias potential goes to 0.0 (unitless)
* Vmax = estimated height of bias potential (energy units)
* Tequil = equilibration temperature (temperature units)
* Dcut = minimum distance between boosted bonds (distance units)
* alpha = boostostat relaxation time (time units)
* Btarget = desired time boost factor (unitless)
* zero or more keyword/value pairs may be appended
* keyword = *check/ghost* or *check/bias*
  
  .. parsed-literal::
  
       *check/ghost* values = none
       *check/bias* values = Nevery error/warn/ignore



Examples
""""""""


.. parsed-literal::

   fix 1 all hyper/local 1.0 0.3 0.8 300.0

Description
"""""""""""

This fix is meant to be used with the :doc:`hyper <hyper>` command to
perform a bond-boost local hyperdynamics (LHD) simulation.  The role
of this fix is to a select multiple pairs of atoms in the system at
each timestep to add a local bias potential to, which will alter the
dynamics of the system in a manner that effectively accelerates time.
This is in contrast to the :doc:`fix hyper/global <fix_hyper_global>`
command, which can be user to perform a global hyperdynamics (GHD)
simulation, by adding a global bias potential to a single pair of
atoms at each timestep.  GHD can time accelerate a small simulation
with up to a few 100 atoms.  For larger systems, LHD is needed to
achieve good time acceleration.

For a system that undergoes rare transition events, where one or more
atoms move over an energy barrier to a new potential energy basin, the
effect of the bias potential is to induce more rapid transitions.
This can lead to a dramatic speed-up in the rate at which events
occurs, without altering their relative frequencies, thus leading to
an overall increase in the elapsed real time of the simulation as
compared to running for the same number of timesteps with normal MD.
See the :doc:`hyper <hyper>` doc page for a more general discussion of
hyperdynamics and citations that explain both GHD and LHD.

The equations and logic used by this fix and described here to perform
LHD follow the description given in :ref:`(Voter2013) <Voter2013lhd>`.  The
bond-boost form of a bias potential for HD is due to Miron and
Fichthorn as described in :ref:`(Miron) <Mironlhd>`.

To understand this description, you should first read the description
of the GHD algorithm on the :doc:`fix hyper/global <fix_hyper_global>`
doc page.  This description of LHD builds on the GHD description.

The definition of bonds and Eij are the same for GHD and LHD.  The
formulas for Vij(max) and Fij(max) are also the same except for a
pre-factor Cij, explained below.

The bias energy Vij applied to a bond IJ with maximum strain is


.. parsed-literal::

   Vij(max) = Cij \* Vmax \* (1 - (Eij/q)\^2) for abs(Eij) < qfactor
            = 0 otherwise

The derivative of Vij(max) with respect to the position of each atom
in the IJ bond gives a bias force Fij(max) acting on the bond as


.. parsed-literal::

   Fij(max) = - dVij(max)/dEij = 2 Cij Vmax Eij / qfactor\^2   for abs(Eij) < qfactor
            = 0 otherwise

which can be decomposed into an equal and opposite force acting on
only the two I,J atoms in the IJ bond.

The key difference is that in GHD a bias energy and force is added (on
a particular timestep) to only one bond (pair of atoms) in the system,
which is the bond with maximum strain Emax.

In LHD, a bias energy and force can be added to multiple bonds
separated by the specified *Dcut* distance or more.  A bond IJ is
biased if it is the maximum strain bond within its local
"neighborhood", which is defined as the bond IJ plus any neighbor
bonds within a distance *Dcut* from IJ.  The "distance" between bond
IJ and bond KL is the minimum distance between any of the IK, IL, JK,
JL pairs of atoms.

For a large system, multiple bonds will typically meet this
requirement, and thus a bias potential Vij(max) will be applied to
many bonds on the same timestep.

In LHD, all bonds store a Cij prefactor which appears in the Vij(max)
and Fij(max) equations above.  Note that the Cij factor scales the
strength of the bias energy and forces whenever bond IJ is the maximum
strain bond in its neighborhood.

Cij is initialized to 1.0 when a bond between the I,J atoms is first
defined.  The specified *Btarget* factor is then used to adjust the
Cij prefactors for each bond every timestep in the following manner.

An instantaneous boost factor Bij is computed each timestep
for each bond, as


.. parsed-literal::

   Bij = exp(beta \* Vkl(max))

where Vkl(max) is the bias energy of the maxstrain bond KL within bond
IJ's neighborhood, beta = 1/kTequil, and *Tequil* is the temperature
of the system and an argument to this fix.

.. note::

   To run an LHD simulation, the input script must also use the
   :doc:`fix langevin <fix_langevin>` command to thermostat the atoms at
   the same *Tequil* as specified by this fix, so that the system is
   running constant-temperature (NVT) dynamics.  LAMMPS does not check
   that this is done.

Note that if IJ = KL, then bond IJ is a biased bond on that timestep,
otherwise it is not.  But regardless, the boost factor Bij can be
thought of an estimate of time boost currently being applied within a
local region centered on bond IJ.  For LHD, we want this to be the
specified *Btarget* value everywhere in the simulation domain.

To accomplish this, if Bij < Btarget, the Cij prefactor for bond IJ is
incremented on the current timestep by an amount proportional to the
inverse of the specified *alpha* and the difference (Bij - Btarget).
Conversely if Bij > Btarget, Cij is decremented by the same amount.
This procedure is termed "boostostatting" in
:ref:`(Voter2013) <Voter2013lhd>`.  It drives all of the individual Cij to
values such that when Vij\ *max* is applied as a bias to bond IJ, the
resulting boost factor Bij will be close to *Btarget* on average.
Thus the LHD time acceleration factor for the overall system is
effectively *Btarget*\ .

Note that in LHD, the boost factor *Btarget* is specified by the user.
This is in contrast to global hyperdynamics (GHD) where the boost
factor varies each timestep and is computed as a function of *Vmax*\ ,
Emax, and *Tequil*\ ; see the :doc:`fix hyper/global <fix_hyper_global>`
doc page for details.


----------


Here is additional information on the input parameters for LHD.

Note that the *cutbond*\ , *qfactor*\ , and *Tequil* arguments have the
same meaning as for GHD.  The *Vmax* argument is slightly different.
The *Dcut*\ , *alpha*\ , and *Btarget* parameters are unique to LHD.

The *cutbond* argument is the cutoff distance for defining bonds
between pairs of nearby atoms.  A pair of I,J atoms in their
equilibrium, minimum-energy configuration, which are separated by a
distance Rij < *cutbond*\ , are flagged as a bonded pair.  Setting
*cubond* to be ~25% larger than the nearest-neighbor distance in a
crystalline lattice is a typical choice for solids, so that bonds
exist only between nearest neighbor pairs.

The *qfactor* argument is the limiting strain at which the bias
potential goes to 0.0.  It is dimensionless, so a value of 0.3 means a
bond distance can be up to 30% larger or 30% smaller than the
equilibrium (quenched) R0ij distance and the two atoms in the bond
could still experience a non-zero bias force.

If *qfactor* is set too large, then transitions from one energy basin
to another are affected because the bias potential is non-zero at the
transition state (e.g. saddle point).  If *qfactor* is set too small
than little boost can be achieved because the Eij strain of some bond in
the system will (nearly) always exceed *qfactor*\ .  A value of 0.3 for
*qfactor* is typically a reasonable value.

The *Vmax* argument is a fixed prefactor on the bias potential.  There
is a also a dynamic prefactor Cij, driven by the choice of *Btarget*
as discussed above.  The product of these should be a value less than
the smallest barrier height for an event to occur.  Otherwise the
applied bias potential may be large enough (when added to the
interatomic potential) to produce a local energy basin with a maxima
in the center.  This can produce artificial energy minima in the same
basin that trap an atom.  Or if Cij\*\ *Vmax* is even larger, it may
induce an atom(s) to rapidly transition to another energy basin.  Both
cases are "bad dynamics" which violate the assumptions of LHD that
guarantee an accelerated time-accurate trajectory of the system.

.. note::

   It may seem that *Vmax* can be set to any value, and Cij will
   compensate to reduce the overall prefactor if necessary.  However the
   Cij are initialized to 1.0 and the boostostatting procedure typically
   operates slowly enough that there can be a time period of bad dynamics
   if *Vmax* is set too large.  A better strategy is to set *Vmax* to the
   smallest barrier height for an event (the same as for GHD), so that
   the Cij remain near unity.

The *Tequil* argument is the temperature at which the system is
simulated; see the comment above about the :doc:`fix langevin <fix_langevin>` thermostatting.  It is also part of the
beta term in the exponential factor that determines how much boost is
achieved as a function of the bias potential.  See the discussion of
the *Btarget* argument below.

As discussed above, the *Dcut* argument is the distance required
between two locally maxstrain bonds for them to both be selected as
biased bonds on the same timestep.  Computationally, the larger *Dcut*
is, the more work (computation and communication) must be done each
timestep within the LHD algorithm.  And the fewer bonds can be
simultaneously biased, which may mean the specified *Btarget* time
acceleration cannot be achieved.

Physically *Dcut* should be a long enough distance that biasing two
pairs of atoms that close together will not influence the dynamics of
each pair.  E.g. something like 2x the cutoff of the interatomic
potential.  In practice a *Dcut* value of ~10 Angstroms seems to work
well for many solid-state systems.

.. note::

   You should insure that ghost atom communication is performed for
   a distance of at least *Dcut* + *cutevent* = the distance one or more
   atoms move (between quenched states) to be considered an "event".  It
   is an argument to the "compute event/displace" command used to detect
   events.  By default the ghost communication distance is set by the
   pair\_style cutoff, which will typically be < *Dcut*\ .  The :doc:`comm\_modify cutoff <comm_modify>` command should be used to override the ghost
   cutoff explicitly, e.g.


.. parsed-literal::

   comm_modify cutoff 12.0

Note that this fix does not know the *cutevent* parameter, but uses
half the *cutbond* parameter as an estimate to warn if the ghost
cutoff is not long enough.

As described above the *alpha* argument is a pre-factor in the
boostostat update equation for each bond's Cij prefactor.  *Alpha* is
specified in time units, similar to other thermostat or barostat
damping parameters.  It is roughly the physical time it will take the
boostostat to adjust a Cij value from a too high (or too low) value to
a correct one.  An *alpha* setting of a few ps is typically good for
solid-state systems.  Note that the *alpha* argument here is the
inverse of the alpha parameter discussed in
:ref:`(Voter2013) <Voter2013lhd>`.

The *Btarget* argument is the desired time boost factor (a value > 1)
that all the atoms in the system will experience.  The elapsed time
t\_hyper for an LHD simulation running for *N* timesteps is simply


.. parsed-literal::

   t_hyper = Btarget \* N\*dt

where dt is the timestep size defined by the :doc:`timestep <timestep>`
command.  The effective time acceleration due to LHD is thus t\_hyper /
N\*dt = Btarget, where N\*dt is elapsed time for a normal MD run
of N timesteps.

You cannot choose an arbitrarily large setting for *Btarget*\ .  The
maximum value you should choose is


.. parsed-literal::

   Btarget = exp(beta \* Vsmall)

where Vsmall is the smallest event barrier height in your system, beta
= 1/kTequil, and *Tequil* is the specified temperature of the system
(both by this fix and the Langevin thermostat).

Note that if *Btarget* is set smaller than this, the LHD simulation
will run correctly.  There will just be fewer events because the hyper
time (t\_hyper equation above) will be shorter.

.. note::

   If you have no physical intuition as to the smallest barrier
   height in your system, a reasonable strategy to determine the largest
   *Btarget* you can use for an LHD model, is to run a sequence of
   simulations with smaller and smaller *Btarget* values, until the event
   rate does not change (as a function of hyper time).


----------


Here is additional information on the optional keywords for this fix.

The *check/ghost* keyword turns on extra computation each timestep to
compute statistics about ghost atoms used to determine which bonds to
bias.  The output of these stats are the vector values 14 and 15,
described below.  If this keyword is not enabled, the output
of the stats will be zero.

The *check/bias* keyword turns on extra computation and communication
to check if any biased bonds are closer than *Dcut* to each other,
which should not be the case if LHD is operating correctly.  Thus it
is a debugging check.  The *Nevery* setting determines how often the
check is made.  The *error*\ , *warn*\ , or *ignore* setting determines
what is done if the count of too-close bonds is not zero.  Either the
code will exit, or issue a warning, or silently tally the count.  The
count can be output as vector value 17, as described below.  If this
keyword is not enabled, the output of that statistic will be 0.

Note that both of these computations are costly, hence they are only
enabled by these keywords.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix\_modify <fix_modify>` *energy* option is supported by this
fix to add the energy of the bias potential to the system's
potential energy as part of :doc:`thermodynamic output <thermo_style>`.

This fix computes a global scalar and global vector of length 21,
which can be accessed by various :doc:`output commands <Howto_output>`.
The scalar is the magnitude of the bias potential (energy units)
applied on the current timestep, summed over all biased bonds.  The
vector stores the following quantities:

* 1 = # of biased bonds on this step
* 2 = max strain Eij of any bond on this step (absolute value, unitless)
* 3 = average bias coeff for all bonds on this step (unitless)
* 4 = average # of bonds/atom on this step
* 5 = average neighbor bonds/bond on this step within *Dcut*

* 6 = max bond length during this run (distance units)
* 7 = average # of biased bonds/step during this run
* 8 = fraction of biased bonds with no bias during this run
* 9 = fraction of biased bonds with negative strain during this run
* 10 = average bias coeff for all bonds during this run (unitless)
* 11 = min bias coeff for any bond during this run (unitless)
* 12 = max bias coeff for any bond during this run (unitless)

* 13 = max drift distance of any bond atom during this run (distance units)
* 14 = max distance from proc subbox of any ghost atom with maxstrain < qfactor during this run (distance units)
* 15 = max distance outside my box of any ghost atom with any maxstrain during this run (distance units)
* 16 = count of ghost atoms that could not be found on reneighbor steps during this run
* 17 = count of bias overlaps (< Dcut) found during this run

* 18 = cumulative hyper time since fix created (time units)
* 19 = cumulative count of event timesteps since fix created
* 20 = cumulative count of atoms in events since fix created
* 21 = cumulative # of new bonds formed since fix created

The first quantities (1-5) are for the current timestep.  Quantities
6-17 are for the current hyper run.  They are reset each time a new
hyper run is performed.  Quantities 18-21 are cumulative across
multiple runs (since the point in the input script the fix was
defined).

For value 8, the numerator is a count of all biased bonds on each
timestep whose bias energy = 0.0 due to Eij >= *qfactor*\ .  The
denominator is the count of all biased bonds on all timesteps.

For value 9, the numerator is a count of all biased bonds on each
timestep with negative strain.  The denominator is the count of all
biased bonds on all timesteps.

Values 13-17 are mostly useful for debugging and diagnostic purposes.

For value 13, drift is the distance an atom moves between two quenched
states when the second quench determines an event has occurred.  Atoms
involved in an event will typically move the greatest distance since
others typically remain near their original quenched position.

For values 14-16, neighbor atoms in the full neighbor list with cutoff
*Dcut* may be ghost atoms outside a processor's sub-box.  Before the
next event occurs they may move further than *Dcut* away from the
sub-box boundary.  Value 14 is the furthest (from the sub-box) any
ghost atom in the neighbor list with maxstrain < *qfactor* was
accessed during the run.  Value 15 is the same except that the ghost
atom's maxstrain may be >= *qfactor*\ , which may mean it is about to
participate in an event.  Value 16 is a count of how many ghost atoms
could not be found on reneighbor steps, presumably because they moved
too far away due to their participation in an event (which will likely
be detected at the next quench).

Typical values for 14 and 15 should be slightly larger than *Dcut*\ ,
which accounts for ghost atoms initially at a *Dcut* distance moving
thermally before the next event takes place.

Note that for values 14 and 15 to be computed, the optional keyword
*check/ghost* must be specified.  Otherwise these values will be zero.
This is because computing them incurs overhead, so the values are only
computed if requested.

Value 16 should be zero or small.  As explained above a small count
likely means some ghost atoms were participating in their own events
and moved a longer distance.  If the value is large, it likely means
the communication cutoff for ghosts is too close to *Dcut* leading to
many not-found ghost atoms before the next event.  This may lead to a
reduced number of bonds being selected for biasing, since the code
assumes those atoms are part of highly strained bonds.  As explained
above, the :doc:`comm\_modify cutoff <comm_modify>` command can be used
to set a longer cutoff.

For value 17, no two bonds should be biased if they are within a
*Dcut* distance of each other.  This value should be zero, indicating
that no pair of biased bonds are closer than *Dcut* from each other.

Note that for values 17 to be computed, the optional keyword
*check/bias* must be specified and it determines how often this check
is performed.  This is because performing the check incurs overhead,
so if only computed as often as requested.

The result at the end of the run is the cumulative total from every
timestep the check was made.  Note that the value is a count of atoms
in bonds which found other atoms in bonds too close, so it is almost
always an over-count of the number of too-close bonds.

Value 18 is simply the specified *boost* factor times the number of
timesteps times the timestep size.

For value 19, events are checked for by the :doc:`hyper <hyper>` command
once every *Nevent* timesteps.  This value is the count of those
timesteps on which one (or more) events was detected.  It is NOT the
number of distinct events, since more than one event may occur in the
same *Nevent* time window.

For value 20, each time the :doc:`hyper <hyper>` command checks for an
event, it invokes a compute to flag zero or more atoms as
participating in one or more events.  E.g. atoms that have displaced
more than some distance from the previous quench state.  Value 20 is
the cumulative count of the number of atoms participating in any of
the events that were found.

Value 21 tallies the number of new bonds created by the bond reset
operation.  Bonds between a specific I,J pair of atoms may persist for
the entire hyperdynamics simulation if neither I or J are involved in
an event.

The scalar and vector values calculated by this fix are all
"intensive".

This fix also computes a local vector of length the number of bonds
currently in the system.  The value for each bond is its Cij prefactor
(bias coefficient).  These values can be can be accessed by various
:doc:`output commands <Howto_output>`.  A particularly useful one is the
:doc:`fix ave/histo <fix_ave_histo>` command which can be used to
histogram the Cij values to see if they are distributed reasonably
close to 1.0, which indicates a good choice of *Vmax*\ .

The local values calculated by this fix are unitless.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the REPLICA package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`hyper <hyper>`, :doc:`fix hyper/global <fix_hyper_global>`

Default
"""""""

The check/ghost and check/bias keywords are not enabled by default.


----------


.. _Voter2013lhd:



**(Voter2013)** S. Y. Kim, D. Perez, A. F. Voter, J Chem Phys, 139,
144110 (2013).

.. _Mironlhd:



**(Miron)** R. A. Miron and K. A. Fichthorn, J Chem Phys, 119, 6210 (2003).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
