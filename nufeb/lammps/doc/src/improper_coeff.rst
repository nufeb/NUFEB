.. index:: improper\_coeff

improper\_coeff command
=======================

Syntax
""""""


.. parsed-literal::

   improper_coeff N args

* N = improper type (see asterisk form below)
* args = coefficients for one or more improper types

Examples
""""""""


.. parsed-literal::

   improper_coeff 1 300.0 0.0
   improper_coeff \* 80.2 -1 2
   improper_coeff \*4 80.2 -1 2

Description
"""""""""""

Specify the improper force field coefficients for one or more improper
types.  The number and meaning of the coefficients depends on the
improper style.  Improper coefficients can also be set in the data
file read by the :doc:`read\_data <read_data>` command or in a restart
file.

N can be specified in one of two ways.  An explicit numeric value can
be used, as in the 1st example above.  Or a wild-card asterisk can be
used to set the coefficients for multiple improper types.  This takes
the form "\*" or "\*n" or "n\*" or "m\*n".  If N = the number of improper
types, then an asterisk with no numeric values means all types from 1
to N.  A leading asterisk means all types from 1 to n (inclusive).  A
trailing asterisk means all types from n to N (inclusive).  A middle
asterisk means all types from m to n (inclusive).

Note that using an improper\_coeff command can override a previous
setting for the same improper type.  For example, these commands set
the coeffs for all improper types, then overwrite the coeffs for just
improper type 2:


.. parsed-literal::

   improper_coeff \* 300.0 0.0
   improper_coeff 2 50.0 0.0

A line in a data file that specifies improper coefficients uses the
exact same format as the arguments of the improper\_coeff command in an
input script, except that wild-card asterisks should not be used since
coefficients for all N types must be listed in the file.  For example,
under the "Improper Coeffs" section of a data file, the line that
corresponds to the 1st example above would be listed as


.. parsed-literal::

   1 300.0 0.0

The :doc:`improper\_style class2 <improper_class2>` is an exception to
this rule, in that an additional argument is used in the input script
to allow specification of the cross-term coefficients.  See its doc
page for details.


----------


The list of all improper styles defined in LAMMPS is given on the
:doc:`improper\_style <improper_style>` doc page.  They are also listed
in more compact form on the :ref:`Commands improper <improper>` doc page.

On either of those pages, click on the style to display the formula it
computes and its coefficients as specified by the associated
improper\_coeff command.


----------


Restrictions
""""""""""""


This command must come after the simulation box is defined by a
:doc:`read\_data <read_data>`, :doc:`read\_restart <read_restart>`, or
:doc:`create\_box <create_box>` command.

An improper style must be defined before any improper coefficients are
set, either in the input script or in a data file.

Related commands
""""""""""""""""

:doc:`improper\_style <improper_style>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
