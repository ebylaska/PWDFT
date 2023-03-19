#ifndef _EXCHANGE_CORRELATION_HPP_
#define _EXCHANGE_CORRELATION_HPP_


#include        <iostream>
#include        "Control2.hpp"
#include        "Pneb.hpp"

namespace pwdft {


class   XC_Operator {

   Pneb   *mypneb;

   double *xtmp;
   double *rho,*grx,*gry,*grz,*agr,*fn,*fdn;

   std::string xc_name;
   int  gga;
   bool use_lda,use_gga,use_mgga;

public:

   /* Constructors */
   XC_Operator(Pneb *, Control2&);

   /* destructor */
   ~XC_Operator() {
       //std::cout << "finishec XC_Operator" << std::endl;
       if (use_lda) delete [] xtmp;
       if (use_gga) {
          delete [] rho;

          delete [] grx;
          delete [] gry;
          delete [] grz;

          delete [] agr;
          delete [] fn;
          delete [] fdn;
       }
    }

   void v_exc_all(int, double *, double *, double *);


   friend std::ostream& operator<<(std::ostream& os, const XC_Operator& xc) {
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

} // namespace pwdft

#endif // _EXCHANGE_CORRELATION_HPP_
