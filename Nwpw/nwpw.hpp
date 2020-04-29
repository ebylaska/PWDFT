#ifndef _NWPW_HPP_
#define _NWPW_HPP_

#include <string>
#include "mpi.h"


extern int cpsd(MPI_Comm, std::string&);
extern int cpmd(MPI_Comm, std::string&);

#endif


