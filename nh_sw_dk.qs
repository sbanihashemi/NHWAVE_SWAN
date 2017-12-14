#
# Template:  OpenMPI, Low-Bandwidth Ethernet Variant
# Revision:  $Id: openmpi-gige.qs,v 1.1 2012/01/11 21:33:26 frey Exp $
#
# Usage:
# 1. Modify "NPROC" in the -pe line to reflect the number
#    of processors desired
# 2. Modify the value of "MY_EXE" to be your MPI program and any
#    arguments to be passed to it.
# 3. Uncomment the SHOW_MPI_DEBUGGING line if you want very verbose
#    output written to the Grid Engine output file by OpenMPI.
# 
#$ -pe openmpi 16
#
# If you want an email message to be sent to you when your job ultimately
# finishes, edit the -M line to have your email address and change the
# next two lines to start with #$ instead of just #
# -m eas
# -M zdong@udel.edu
#

#
# Setup the environment; choose the OpenMPI version that's
# right for you:
#
source /opt/shared/valet/docs/valet.sh
vpkg_require openmpi/1.4.4-intel64

#
# The MPI program to execute and any arguments to it:
#
MY_EXE="./nhwave_swan"

#
# Uncomment to enable lots of additional information as OpenMPI executes
# your job:
#
#SHOW_MPI_DEBUGGING=YES

##
## You should NOT need to change anything after this comment.
##
OPENMPI_FLAGS="--display-map --mca btl ^openib --mca oob_tcp_if_exclude lo,ib0 --mca btl_tcp_if_exclude lo,ib0 --mca mtl ^psm
"
if [ "x$SHOW_MPI_DEBUGGING" = "xYES" ]; then
  OPENMPI_FLAGS="${OPENMPI_FLAGS} --debug-devel --debug-daemons --display-devel-map --display-devel-allocation --mca mca_verb
ose 1 --mca coll_base_verbose 1 --mca ras_base_verbose 1 --mca ras_gridengine_debug 1 --mca ras_gridengine_verbose 1 --mca bt
l_base_verbose 1 --mca mtl_base_verbose 1 --mca plm_base_verbose 1 --mca pls_rsh_debug 1"
fi

echo "GridEngine parameters:"
echo "  mpirun        = "`which mpirun`
echo "  nhosts        = $NHOSTS"
echo "  nproc         = $NSLOTS"
echo "  executable    = $MY_EXE"
echo "  MPI flags     = $OPENMPI_FLAGS"
echo "-- begin OPENMPI run --"
mpirun ${OPENMPI_FLAGS} $MY_EXE
echo "-- end OPENMPI run --"


