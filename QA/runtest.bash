#!/bin/bash -f

#
#define the pwdft and python3 binaries
#
if [ -z "${PWDFT_BIN}" ];  then
   PWDFT_BIN=`which pwdft`
   if [ -z "${PWDFT_BIN}" ];  then PWDFT_BIN=../../build/pwdft; fi
fi

if [ -z "${PYTHON_BIN}" ];  then
   PYTHON_BIN=`which python3`
fi


if [[ $# -eq 0 ]] ; then
  echo "runtest.bash [-n nproc | -p nproc | -h] testdir1 [testdir2 ...]"
  echo " -procs nproc sets the number of processors to use"
  exit 0
fi

#NPROC=
while getopts 'n:p:h' opt; do
    case "$opt" in
     n)
      NPROC="$OPTARG"
      #echo "Processing option 'n' with '${OPTARG}' argument"
      shift; shift;
      ;;

     p)
      NPROC="$OPTARG"
      NOHEADER="printing"
      #echo "Processing option 'n' with '${OPTARG}' argument"
      shift; shift;
      ;;

    ?|h)
      echo "Usage: $(basename $0) [-h | -p nproc | -n nproc] testdir1 [testdir2 ...]" 
      exit 1
      ;;
  esac
done

if [ ! -z "${NOHEADER}" ];  then 
   echo "PWDFT_BIN="$PWDFT_BIN
   echo "PYTHON_BIN="$PYTHON_BIN
   echo  "NPROC="$NPROC
   if [ ! -z "${NPROC}" ];  then 
       echo "Shifted argument list!"; 
   fi
   NUMBER_OF_JOBS=$#
   echo "NUMBER_OF_JOBS="$NUMBER_OF_JOBS
   echo " "
fi

if [ -z "${NPROC}" ];  then 
    NPROC=1
fi

CWD=`pwd`
SCRATCHDIR=$CWD/scratchdir
TESTOUTPUTS=$CWD/testoutputs

for i in $@
do
   INPUT=$CWD/tests/$i/$i.nw
   echo "Running tests/$i/$i.nw"
   PWDFTPARSE=$CWD/pwparse.py
   RUNPWDFT="mpirun -np $NPROC $PWDFT_BIN "
   rm -rf $SCRATCHDIR
   mkdir -p $SCRATCHDIR $TESTOUTPUTS
   cd $SCRATCHDIR
   cp $INPUT .
   echo    "   Executing pwdft:"
   echo    "      Current directory:" `pwd`
   echo    "      Executable:" $RUNPWDFT  $i.nw ">" $i.out
   echo -n "      Running...."
   $RUNPWDFT $i.nw > $i.out
   echo "Finished"
   RUNSTATUS=$?
   if [ $RUNSTATUS -ne 0 ];then
      echo "   Execution failed!"
      exit 
   fi
   cp $i.out $TESTOUTPUTS
   cd $TESTOUTPUTS

   $PYTHON_BIN $PWDFTPARSE $CWD/tests/$i/$i.out $i.out > $i.check

   echo -n "   Verifying output: "
   if grep -Fxq "Failed" $1.check
   then
      echo "Failed"
   else
      echo "OK"
   fi
   echo " "
done




