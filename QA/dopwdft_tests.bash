#!/usr/bin/env bash -f
#
# HOW TO RUN in BASH: export NPROCS=3; export PWDFT_BIN=../../build/pwdft;  ./dopwdft_tests.bash
# HOW TO RUN in BASH: export PWDFT_BIN=../../build/pwdft; export NPROCS=3;  ./dopwdft_tests.bash
# HOW TO RUN in BASH: export PWDFT_BIN=../../build/pwdft;  ./dopwdft_tests.bash
#
if [ -z "${NPROCS}" ];  then
   NPROCS=2
fi

if [ -z "${PWDFT_BIN}" ];  then
   TFILE1=../build/pwdft
   TFILE2=`which pwdft`
   if test -f "$TFILE1"; then
      export PWDFT_BIN=`readlink -f $TFILE1`
   elif test -f "$TFILE2"; then
      export PWDFT_BIN=`readlink -f $TFILE2`
   else
       echo "pwdft binary not found!! recompile?"
       exit 1
   fi
fi

echo "Running QA/tests: NPROCS="$NPROCS "PWDFT_BIN="$PWDFT_BIN
echo " "

echo "Short QA/tests:"
echo " "
./runtest.bash -n $NPROCS C2

echo " "
echo "Medium QA/tests:"
echo " "
./runtest.bash -n $NPROCS benzene ch3cl ccl4_born_test periodic_polarizability aperiodic_polarizability

echo " "
echo "Long QA/tests:"
echo " "
