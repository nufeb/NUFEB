.. index:: pair\_style sdpd/taitwater/isothermal

pair\_style sdpd/taitwater/isothermal command
=============================================

Syntax
""""""


.. parsed-literal::

   pair_style sdpd/taitwater/isothermal temperature viscosity seed

* temperature = temperature of the fluid (temperature units)
* viscosity = dynamic viscosity of the fluid (mass\*distance/time units)
* seed = random number generator seed (positive integer, optional)

Examples
""""""""


.. parsed-literal::

   pair_style sdpd/taitwater/isothermal 300. 1. 28681
   pair_coeff \* \* 1000.0 1430.0 2.4

Description
"""""""""""

The sdpd/taitwater/isothermal style computes forces between mesoscopic
particles according to the Smoothed Dissipative Particle Dynamics model
described in this paper by :ref:`(Espanol and Revenga) <Espanol_Revenga>` under
the following assumptions:



#. The temperature is constant and uniform.
#. The shear viscosity is constant and uniform.
#. The volume viscosity is negligible before the shear viscosity.
#. The Boltzmann constant is negligible before the heat capacity of a
   single mesoscopic particle of fluid.

The third assumption is true for water in nearly incompressible flows.
The fourth holds true for water for any reasonable size one can
imagine for a mesoscopic particle.

The pressure forces between particles will be computed according to
Tait's equation of state:

.. image:: Eqs/pair_sph_tait.jpg
   :align: center

where gamma = 7 and B = c\_0\^2 rho\_0 / gamma, with rho\_0 being the
reference density and c\_0 the reference speed of sound.

The laminar viscosity and the random forces will be computed according
to formulas described in :ref:`(Espanol and Revenga) <Espanol_Revenga>`.

.. warning::

   Similar to :doc:`brownian <pair_brownian>` and
   :doc:`dpd <pair_dpd>` styles, the :doc:`newton <newton>` setting for
   pairwise interactions needs to be on when running LAMMPS in parallel
   if you want to ensure linear momentum conservation. Otherwise random
   forces generated for pairs straddling processor boundary will not be
   equal and opposite.

.. note::

   The actual random seed used will be a mix of what you specify
   and other parameters like the MPI ranks. This is to ensure that
   different MPI tasks have distinct seeds.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair\_coeff <pair_coeff>` command as in the examples
above.

* rho0 reference density (mass/volume units)
* c0 reference soundspeed (distance/time units)
* h kernel function cutoff (distance units)


----------


**Mixing, shift, table, tail correction, restart, rRESPA info**\ :

This style does not support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

This style does not support the :doc:`pair\_modify <pair_modify>`
shift, table, and tail options.

This style does not write information to :doc:`binary restart files <restart>`.  Thus, you need to re-specify the pair\_style and
pair\_coeff commands in an input script that reads a restart file.

This style can only be used via the *pair* keyword of the :doc:`run\_style respa <run_style>` command.  It does not support the *inner*\ ,
*middle*\ , *outer* keywords.

Restrictions
""""""""""""


This pair style is part of the USER-SDPD package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`pair coeff <pair_coeff>`, :doc:`pair sph/rhosum <pair_sph_rhosum>`

Default
"""""""

The default seed is 0 (before mixing).


----------


.. _Espanol\_Revenga:



**(Espanol and Revenga)** Espanol, Revenga, Physical Review E, 67, 026705 (2003).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
