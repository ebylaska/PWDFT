#ifndef _VDW_HPP_
#define _VDW_HPP_

#include "Parallel.hpp" 

namespace pwdft {

extern double vdw_DF_kernel_phi_value(int, const double *, const double *, double *, 
                                      double *, const double *, double, double );

extern void vdw_DF_kernel_gen_phir(int, const double *, const double *, double *, double *, 
                                   const double *, double, double, int, double, double *);

extern void vdw_GaussLegendre(double, double, int, double *, double *);

extern double vdw_DF_Fsin(double, double, double, double, double, double, double);

extern void vdw_DF_kernel_gen_Wab(int, double *, double *, double *, double *, double *, double *);

extern void vdw_DF_kernel_gen_data(Parallel *, const std::string &, const std::string &);


} // namespace pwdft

#endif
