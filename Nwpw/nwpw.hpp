#ifndef _NWPW_HPP_
#define _NWPW_HPP_

#pragma once

#include "mpi.h"
#include <string>

namespace pwdft {

extern int band_cpsd(MPI_Comm, std::string &);
extern int band_minimizer(MPI_Comm, std::string &, std::ostream &);

extern int cpsd(MPI_Comm, std::string &);
extern int cpmd(MPI_Comm, std::string &);
extern int pspw_minimizer(MPI_Comm, std::string &, std::ostream &);
extern int pspw_geovib(MPI_Comm, std::string &, std::ostream &);
extern int pspw_bomd(MPI_Comm, std::string &, std::ostream &);
extern int pspw_dplot(MPI_Comm, std::string &, std::ostream &);

extern int file_generate(std::string &);

extern int ctask_cpmd_start(MPI_Comm, std::string &, double *, double *,
                            double *, double *, double *, double *,
                            std::ostream &);
extern int ctask_cpmd_run(MPI_Comm, double *, double *, double *, double *,
                          double *, double *, std::ostream &);
extern int ctask_cpmd_stop(MPI_Comm, std::ostream &);

extern int cpsd_debug(MPI_Comm, std::string &);

} // namespace pwdft

#endif
