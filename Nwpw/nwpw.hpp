#ifndef _NWPW_HPP_
#define _NWPW_HPP_

#include <string>
#include "mpi.h"

namespace pwdft {


extern int cpsd(MPI_Comm, std::string&);
extern int cpmd(MPI_Comm, std::string&);
extern int pspw_minimizer(MPI_Comm, std::string&,std::ostream&);
extern int pspw_geovib(MPI_Comm, std::string&, std::ostream&);
extern int pspw_bomd(MPI_Comm, std::string&, std::ostream&);

extern int cpsd_debug(MPI_Comm, std::string&);
}

#endif


