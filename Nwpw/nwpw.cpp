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
  if (argc < 2) {
    // report version
    std::cout << argv[0] << " Version " << Nwpw_VERSION_MAJOR << "."
              << Nwpw_VERSION_MINOR << std::endl;
    std::cout << "Usage: " << argv[0] << " number" << std::endl;
    std::cout << "Hello Eric"
            << std::endl;
    return 1;
  }

  std::cout << argv[0] << " (NWChemEx) - Version " << Nwpw_VERSION_MAJOR << "."
            << Nwpw_VERSION_MINOR << std::endl;


  // convert input to double
  //const double inputValue = atof(argv[1]);
  const double inputValue = std::stod(argv[1]);

  // calculate square root
  const double outputValue = sqrt(inputValue);
  std::cout << "The square root of " << inputValue << " is " << outputValue
            << std::endl;

  std::cout << "Hello Eric"
            << std::endl;

  int ijk = cpsd(argc,argv);

  return 0;
}
