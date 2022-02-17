.. index:: compute gyration/chunk

compute gyration/chunk command
==============================

Syntax
""""""


.. parsed-literal::

   compute ID group-ID gyration/chunk chunkID keyword value ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* gyration/chunk = style name of this compute command
* chunkID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command
* zero or more keyword/value pairs may be appended
* keyword = *tensor*
  
  .. parsed-literal::
  
       *tensor* value = none



Examples
""""""""


.. parsed-literal::

   compute 1 molecule gyration/chunk molchunk
   compute 2 molecule gyration/chunk molchunk tensor

Description
"""""""""""

Define a computation that calculates the radius of gyration Rg for
multiple chunks of atoms.

In LAMMPS, chunks are collections of atoms defined by a :doc:`compute chunk/atom <compute_chunk_atom>` command, which assigns each atom
to a single chunk (or no chunk).  The ID for this command is specified
as chunkID.  For example, a single chunk could be the atoms in a
molecule or atoms in a spatial bin.  See the :doc:`compute chunk/atom <compute_chunk_atom>` and :doc:`Howto chunk <Howto_chunk>`
doc pages for details of how chunks can be defined and examples of how
they can be used to measure properties of a system.

This compute calculates the radius of gyration Rg for each chunk,
which includes all effects due to atoms passing through periodic
boundaries.

Rg is a measure of the size of a chunk, and is computed by this
formula

.. image:: Eqs/compute_gyration.jpg
   :align: center

where M is the total mass of the chunk, Rcm is the center-of-mass
position of the chunk, and the sum is over all atoms in the
chunk.

Note that only atoms in the specified group contribute to the
calculation.  The :doc:`compute chunk/atom <compute_chunk_atom>` command
defines its own group; atoms will have a chunk ID = 0 if they are not
in that group, signifying they are not assigned to a chunk, and will
thus also not contribute to this calculation.  You can specify the
"all" group for this command if you simply want to include atoms with
non-zero chunk IDs.

If the *tensor* keyword is specified, then the scalar Rg value is not
calculated, but an Rg tensor is instead calculated for each chunk.
The formula for the components of the tensor is the same as the above
formula, except that (Ri - Rcm)\^2 is replaced by (Rix - Rcmx) \* (Riy -
Rcmy) for the xy component, etc.  The 6 components of the tensor are
ordered xx, yy, zz, xy, xz, yz.

.. note::

   The coordinates of an atom contribute to Rg in "unwrapped" form,
   by using the image flags associated with each atom.  See the :doc:`dump custom <dump>` command for a discussion of "unwrapped" coordinates.
   See the Atoms section of the :doc:`read\_data <read_data>` command for a
   discussion of image flags and how they are set for each atom.  You can
   reset the image flags (e.g. to 0) before invoking this compute by
   using the :doc:`set image <set>` command.

The simplest way to output the results of the compute gyration/chunk
calculation to a file is to use the :doc:`fix ave/time <fix_ave_time>`
command, for example:


.. parsed-literal::

   compute cc1 all chunk/atom molecule
   compute myChunk all gyration/chunk cc1
   fix 1 all ave/time 100 1 100 c_myChunk file tmp.out mode vector

**Output info:**

This compute calculates a global vector if the *tensor* keyword is not
specified and a global array if it is.  The length of the vector or
number of rows in the array = the number of chunks *Nchunk* as
calculated by the specified :doc:`compute chunk/atom <compute_chunk_atom>` command.  If the *tensor* keyword
is specified, the global array has 6 columns.  The vector or array can
be accessed by any command that uses global values from a compute as
input.  See the :doc:`Howto output <Howto_output>` doc page for an
overview of LAMMPS output options.

All the vector or array values calculated by this compute are
"intensive".  The vector or array values will be in distance
:doc:`units <units>`, since they are the square root of values
represented by the formula above.

Restrictions
""""""""""""
 none

**Related commands:** none

:doc:`compute gyration <compute_gyration>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
