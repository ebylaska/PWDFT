#ifndef _UTIL_LINESEARCH_HPP_
#define _UTIL_LINESEARCH_HPP_

#pragma once

namespace pwdft {


extern void util_linesearch_init();
extern void util_linesearch_maxiter_set(const int);
extern int  util_linesearch_maxiter();
extern int  util_linesearch_counter();

extern double util_linesearch(const double, const double, const double, double,
                              double (*)(double), double (*)(double),
                              const double, 
                              double *, double *,
                              const int);
}
#endif
