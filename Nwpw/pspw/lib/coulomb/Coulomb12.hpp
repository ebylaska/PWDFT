#ifndef _COULOMB12_HPP_
#define _COULOMB12_HPP_

#pragma once

#include	"Control2.hpp"
#include	"Pneb.hpp"
#include	"Coulomb.hpp"
#include	"Coulomb2.hpp"

namespace pwdft {

class	Coulomb12_Operator {


public:

   bool has_coulomb1 = false;
   Coulomb_Operator  *mycoulomb1;

   bool has_coulomb2 = false;
   Coulomb2_Operator *mycoulomb2;

   /* Constructors */
   Coulomb12_Operator(Pneb *, Control2&);

   /* destructor */
   ~Coulomb12_Operator() {
       if (has_coulomb1) delete  mycoulomb1;
       if (has_coulomb2) delete  mycoulomb2;
    }

};

}

#endif
