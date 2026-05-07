#ifndef _EXCHANGE_CORRELATION_HPP_
#define _EXCHANGE_CORRELATION_HPP_

#include "Control2.hpp"
#include "Pneb.hpp"
#include "parsestring.hpp"
#include "vdw_DF.hpp"
#include <iostream>

namespace pwdft {

class XC_Operator {

  Pneb *mypneb;
  vdw_DF *myvdw = nullptr;

  double *xtmp;
  double *rho, *grx, *gry, *grz, *agr, *fn, *fdn;
  double *tau,*dfdtau;

  std::string xc_name,options_disp;
  int gga;
  bool use_lda, use_gga, use_mgga;

  bool has_disp = false;
  bool has_vdw  = false;
  bool is_grimme2,is_vdw2;

public:
  /* Constructors */
  XC_Operator(Pneb *, Control2 &);

  /* destructor */
  ~XC_Operator() 
  {
     // std::cout << "finishec XC_Operator" << std::endl;
     if (use_lda)
     {
        delete[] xtmp;
     }

     if (use_gga) 
     {
        delete[] rho;
    
        delete[] grx;
        delete[] gry;
        delete[] grz;
    
        delete[] agr;
        delete[] fn;
        delete[] fdn;
        delete[] fdn;
     }

     if (use_mgga) 
     {
        delete[] rho;
    
        delete[] grx;
        delete[] gry;
        delete[] grz;
    
        delete[] agr;
        delete[] fn;
        delete[] fdn;
        delete[] fdn;

        delete[] tau;
        delete[] dfdtau;
     }
    
     delete myvdw;
  }

  void v_exc_all(int, double *, double *, double *);

  friend std::ostream &operator<<(std::ostream &os, const XC_Operator &xc) 
  {
     os << "   exchange-correlation = ";
     if (xc.gga == 0)
        os << "LDA (Vosko et al) parameterization\n";
     if (xc.gga == 10)
        os << "PBE96 (White and Bird) parameterization\n";
     if (xc.gga == 11)
        os << "BLYP (White and Bird) parameterization\n";
     if (xc.gga == 12)
        os << "revPBE (White and Bird) parameterization\n";
     if (xc.gga == 13)
        os << "PBEsol (White and Bird) parameterization\n";
     // if (xc.gga==14) os << "HSER (White and Bird) remainder\n";
     // if (xc.gga==15) os << "B3LYPr (White and Bird) remainder\n";
     if (xc.gga == 16)
        os << "BEEF (White and Bird) parameterization\n";
     if (xc.gga == 17)
        os << "XBEEF-CPBE (White and Bird) parameterization\n";
    
     if (xc.gga == 110)
        os << "PBE0 (White and Bird) parameterization\n";
     if (xc.gga == 111)
        os << "BLYP0 (White and Bird) parameterization\n";
     if (xc.gga == 112)
        os << "revPBE0 (White and Bird) parameterization\n";
     if (xc.gga == 113)
        os << "BNL (White and Bird) parameterization\n";
     if (xc.gga == 114)
        os << "HSE (White and Bird) parameterization\n";
     if (xc.gga == 115)
        os << "B3LYP (White and Bird) parameterization\n";
    
     if (xc.gga == 200)
        os << "Hartree-Fock\n";

     if (xc.gga == 300)
        os << "VS98\n";
     if (xc.gga == 301)
        os << "TPSS03\n";
     if (xc.gga == 302)
        os << "SCAN\n";
     if (xc.gga == 303)
        os << "PKZB\n";
     if (xc.gga == 304)
        os << "M06-L\n";
     if (xc.gga == 305)
        os << "M06\n";
     if (xc.gga == 306)
        os << "M06-2X\n";
     if (xc.gga == 307)
        os << "r2SCAN\n";

     if (xc.has_vdw)
     {
        os << "   dispersion correction= ";
        if (xc.is_vdw2)
           os << "vdw2 Langreth functional\n";
        else
           os << "vdw Langreth functional\n";
     }
     if (xc.has_disp)
     {
        os << "   dispersion correction= ";
        if (mystring_contains(xc.options_disp,"-old"))  
           os << "Grimme2\n";
        else if (mystring_contains(xc.options_disp,"-zerom")) 
           os << "Grimme5\n";
        else if (mystring_contains(xc.options_disp,"-zero"))  
           os << "Grimme3\n";
        else if (mystring_contains(xc.options_disp,"-bjm"))  
           os << "Grimme6\n";
        else if (mystring_contains(xc.options_disp,"-bj"))  
           os << "Grimme4\n";
     }
    
     return os;
  }
};

} // namespace pwdft

#endif // _EXCHANGE_CORRELATION_HPP_
