.. index:: write\_data

write\_data command
===================

Syntax
""""""


.. parsed-literal::

   write_data file keyword value ...

* file = name of data file to write out
* zero or more keyword/value pairs may be appended
* keyword = *pair* or *nocoeff*
  
  .. parsed-literal::
  
       *nocoeff* = do not write out force field info
       *nofix* = do not write out extra sections read by fixes
       *pair* value = *ii* or *ij*
         *ii* = write one line of pair coefficient info per atom type
         *ij* = write one line of pair coefficient info per IJ atom type pair



Examples
""""""""


.. parsed-literal::

   write_data data.polymer
   write_data data.\*

Description
"""""""""""

Write a data file in text format of the current state of the
simulation.  Data files can be read by the :doc:`read data <read_data>`
command to begin a simulation.  The :doc:`read\_data <read_data>` command
also describes their format.

Similar to :doc:`dump <dump>` files, the data filename can contain a "\*"
wild-card character.  The "\*" is replaced with the current timestep
value.

.. note::

   The write-data command is not yet fully implemented in two
   respects.  First, most pair styles do not yet write their coefficient
   information into the data file.  This means you will need to specify
   that information in your input script that reads the data file, via
   the :doc:`pair\_coeff <pair_coeff>` command.  Second, a few of the :doc:`atom styles <atom_style>` (body, ellipsoid, line, tri) that store
   auxiliary "bonus" information about aspherical particles, do not yet
   write the bonus info into the data file.  Both these functionalities
   will be added to the write\_data command later.

Because a data file is in text format, if you use a data file written
out by this command to restart a simulation, the initial state of the
new run will be slightly different than the final state of the old run
(when the file was written) which was represented internally by LAMMPS
in binary format.  A new simulation which reads the data file will
thus typically diverge from a simulation that continued in the
original input script.

If you want to do more exact restarts, using binary files, see the
:doc:`restart <restart>`, :doc:`write\_restart <write_restart>`, and
:doc:`read\_restart <read_restart>` commands.  You can also convert
binary restart files to text data files, after a simulation has run,
using the :doc:`-r command-line switch <Run_options>`.

.. note::

   Only limited information about a simulation is stored in a data
   file.  For example, no information about atom :doc:`groups <group>` and
   :doc:`fixes <fix>` are stored.  :doc:`Binary restart files <read_restart>`
   store more information.

Bond interactions (angle, etc) that have been turned off by the :doc:`fix shake <fix_shake>` or :doc:`delete\_bonds <delete_bonds>` command will
be written to a data file as if they are turned on.  This means they
will need to be turned off again in a new run after the data file is
read.

Bonds that are broken (e.g. by a bond-breaking potential) are not
written to the data file.  Thus these bonds will not exist when the
data file is read.


----------


The *nocoeff* keyword requests that no force field parameters should
be written to the data file. This can be very helpful, if one wants
to make significant changes to the force field or if the parameters
are read in separately anyway, e.g. from an include file.

The *nofix* keyword requests that no extra sections read by fixes
should be written to the data file (see the *fix* option of the
:doc:`read\_data <read_data>` command for details). For example, this
option excludes sections for user-created per-atom properties
from :doc:`fix property/atom <fix_property_atom>`.

The *pair* keyword lets you specify in what format the pair
coefficient information is written into the data file.  If the value
is specified as *ii*\ , then one line per atom type is written, to
specify the coefficients for each of the I=J interactions.  This means
that no cross-interactions for I != J will be specified in the data
file and the pair style will apply its mixing rule, as documented on
individual :doc:`pair\_style <pair_style>` doc pages.  Of course this
behavior can be overridden in the input script after reading the data
file, by specifying additional :doc:`pair\_coeff <pair_coeff>` commands
for any desired I,J pairs.

If the value is specified as *ij*\ , then one line of coefficients is
written for all I,J pairs where I <= J.  These coefficients will
include any specific settings made in the input script up to that
point.  The presence of these I != J coefficients in the data file
will effectively turn off the default mixing rule for the pair style.
Again, the coefficient values in the data file can be overridden
in the input script after reading the data file, by specifying
additional :doc:`pair\_coeff <pair_coeff>` commands for any desired I,J
pairs.


----------


Restrictions
""""""""""""


This command requires inter-processor communication to migrate atoms
before the data file is written.  This means that your system must be
ready to perform a simulation before using this command (force fields
setup, atom masses initialized, etc).

Related commands
""""""""""""""""

:doc:`read\_data <read_data>`, :doc:`write\_restart <write_restart>`

Default
"""""""

The option defaults are pair = ii.


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
