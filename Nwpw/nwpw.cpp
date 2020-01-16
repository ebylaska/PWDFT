// A simple program that computes the square root of a number
#include <cmath>
//#include <cstdlib>
#include <iostream>
#include        <cstdio>
#include        <stdio.h>
#include <string>
#include "NwpwConfig.h"


extern int cpsd(int argc, char *argv[]);


int main(int argc, char* argv[])
{
  std::cout << argv[0] << " (NWChemEx) - Version " << Nwpw_VERSION_MAJOR << "."
            << Nwpw_VERSION_MINOR << std::endl;


  int ijk = cpsd(argc,argv);

  return 0;
}
