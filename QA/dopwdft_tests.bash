#!/usr/bin/env bash -f
#
#
echo "Running QA tests:"
echo
./runtest.bash -n 2 C2 benzene
./runtest.bash -n 2 ch3cl ccl4_born_test
