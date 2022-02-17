#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Get application directory
application=`getApplication`

./Allclean.sh

blockMesh > log.blockMesh

# ----------------------------------------------------------------- end-of-file
