#ifndef _PSEUDOPOTENTIAL_H_
#define _PSEUDOPOTENTIAL_H_

extern void pspsolve(int *, int *, int *, int *, double *, int *,
                     char *, int *, char *, int *, char *, int *, char *, int *);

extern void hgh_parse(int *, int *, int *, double *,
                     char *, int *, char *, int *, char *, int *, char *, int *);

extern void cpi_parse(int *, int *, int *, double *,
                      char *, int *, char *, int *, char *, int *, char *, int *, char *, int *);

extern void teter_parse(int *, int *, int *, double *,
                        char *, int *, char *, int *, char *, int *, char *, int *, char *, int *);

extern void paw_atom_driver(int *, int *, int *, double *,
                            char *, int *, char *, int *, char *, int *, char *, int *, char *, int *);

extern void qmmm_parse(int *, int *, int *, double *,
                       char *, int *, char *, int *, char *, int *, char *, int *, char *, int *);

extern void carter_parse(int *, int *, int *, double *,
                         char *, int *, char *, int *, char *, int *, char *, int *, char *, int *);


#endif
