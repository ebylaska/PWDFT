#!/usr/bin/env bash -f
#
# HOW TO RUN in BASH: export NPROCS=3; ./dopwdft_tests.bash
#
if [ -z "${NPROCS}" ];  then
   NPROCS=2
fi

echo "Running QA tests: NPROCS="$NPROCS
echo " "

./runtest.bash -n $NPROCS C2 benzene
./runtest.bash -n $NPROCS ch3cl ccl4_born_test periodic_polarizability aperiodic_polarizability
