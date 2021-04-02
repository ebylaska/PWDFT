#ifndef _EXCHANGE_CORRELATION_H_
#define _EXCHANGE_CORRELATION_H_


#include        <iostream>
#include        "Control2.hpp"
#include        "Pneb.hpp"

class   XC_Operator {

   Pneb   *mypneb;

   string xc_name;
   int  gga;
   bool use_lda,use_gga,use_mgga;

public:

   /* Constructors */
   XC_Operator(Pneb *, Control2&);

   /* destructor */
   ~XC_Operator() {
       //std::cout << "finishec XC_Operator" << std::endl;
    }

   friend ostream& operator<<(ostream& os, const XC_Operator& xc) {
      os << "   exchange-correlation = ";
      if (xc.gga==0)  os << "LDA (Vosko et al) parameterization\n";
      if (xc.gga==10) os << "PBE96 (White and Bird) parameterization\n";
      if (xc.gga==11) os << "BLYP (White and Bird) parameterization\n";
      if (xc.gga==12) os << "revPBE (White and Bird) parameterization\n";
      if (xc.gga==13) os << "PBEsol (White and Bird) parameterization\n";
      //if (xc.gga==14) os << "HSER (White and Bird) remainder\n";
      //if (xc.gga==15) os << "B3LYPr (White and Bird) remainder\n";
      if (xc.gga==16) os << "BEEF (White and Bird) parameterization\n";
      if (xc.gga==17) os << "XBEEF-CPBE (White and Bird) parameterization\n";

      if (xc.gga==110) os << "PBE0 (White and Bird) parameterization\n";
      if (xc.gga==111) os << "BLYP0 (White and Bird) parameterization\n";
      if (xc.gga==112) os << "revPBE0 (White and Bird) parameterization\n";
      if (xc.gga==113) os << "BNL (White and Bird) parameterization\n";
      if (xc.gga==114) os << "HSE (White and Bird) parameterization\n";
      if (xc.gga==115) os << "B3LYP (White and Bird) parameterization\n";

      if (xc.gga==200) os << "Hartree-Fock\n";

      return os;
   }

};

#endif
