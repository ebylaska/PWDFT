#ifndef _PSP_U_HPP_
#define _PSP_U_HPP_

#pragma once

//
//
// Pseudopotential.hpp - Header file for the Pseudopotential class.
// This class encapsulates operations and data related to pseudopotentials
// used in electronic structure calculations.
//

#include "Control2.hpp"
#include "Ion.hpp"
#include "Pneb.hpp"
#include "Strfac.hpp"
//#include "gdevice.hpp"

namespace pwdft {

/**
 * @class Pseudopotential
 * @brief Class for managing pseudopotential data and calculations.
 *
 * The Pseudopotential class encapsulates operations and data related to
 * pseudopotentials used in electronic structure calculations. It provides
 * methods for handling non-local and local pseudopotentials, semicore
 * corrections, and other properties used in electronic structure calculations.
 */
class Psp_U {

  // char **atomsym;
  Pneb *mypneb;
  Ion *myion;
  Strfac *mystrfac;

private:
  //void apply_pspspin_scaling(double *, int, int, int);

public:

  // uterm variables
  bool uterm = false;
  int nuterms = 0;

  bool   *uterm_ions;
  int    *uterm_l;
  double *uterm_uscale,*uterm_jscale;

  int    *uterm_vstart;

  /* Constructors */
  Psp_U(Ion *, Pneb *, Strfac *, Control2 &, std::ostream &);

    
  /* destructor */
  ~Psp_U() {
    if (uterm)
    {
       delete[] uterm_l;
       delete[] uterm_ions;
       delete[] uterm_uscale;
       delete[] uterm_jscale;
       delete[] uterm_vstart;
    }
  }

};

} // namespace pwdft
  
#endif



