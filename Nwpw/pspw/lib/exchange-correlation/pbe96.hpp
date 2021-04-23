#ifndef _PBE96_HPP_
#define _PBE96_HPP_

extern void gen_PBE96_BW_unrestricted(const int, 
                                      double *, double *,
                                      const double, const double,
                                      double *, double *, double *);

extern void gen_PBE96_BW_restricted(const int,
                                    double *, double *,
                                    const double, const double,
                                    double *, double *, double *);

extern void gen_PBE96_x_unrestricted(const double, const double,
                              double *, double *, double *);

extern void gen_PBE96_x_restricted(const double, const double,
                                   double *, double *, double *);
#endif
