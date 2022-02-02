#ifndef _UTIL_LEGENDRE_HPP_
#define _UTIL_LEGENDRE_HPP_

namespace pwdft {
using namespace pwdft;

extern double util_legendre_lm(const int, const int, const double);
extern double util_rlegendre_lm(const int, const int, const double);
extern double util_legendre_lm_div(const int, const int, const double);
extern double util_dlegendre_lm(const int, const int, const double);
extern double util_rlegendre_lm_div(const int, const int, const double);

extern double util_theta_lm(const int, const int, const double);
extern double util_ytheta_lm(const int, const int, const double);
extern double util_rtheta_lm(const int, const int, const double);
extern double util_ytheta2_lm(const int, const int, const double);
extern double util_rtheta2_lm(const int, const int, const double);

extern double util_theta_lm_div(const int, const int, const double);
extern double util_ytheta_lm_div(const int, const int, const double);
extern double util_rtheta_lm_div(const int, const int, const double);

extern double util_drtheta_lm(const int, const int, const double);
extern double util_drtheta_cos_theta_lm(const int, const int, const double);
extern double util_dytheta_lm(const int, const int, const double);


}
#endif
