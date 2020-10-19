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
       std::cout << "finishec XC_Operator" << std::endl;
    }

};

#endif
