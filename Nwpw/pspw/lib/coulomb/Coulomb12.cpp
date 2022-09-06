/* Coulomb12.cpp - 
   Author - Eric Bylaska
*/

#include	"Coulomb12.hpp"

namespace pwdft {

/*******************************************
 *                                         *
 *  Coulomb12_Operator::Coulomb12_Operator *
 *                                         *
 *******************************************/
Coulomb12_Operator::Coulomb12_Operator(Pneb *mygrid, Control2& control)
{
   if (control.version==3)
   {
      has_coulomb1 = true;
      mycoulomb1   = new (std::nothrow) Coulomb_Operator(mygrid);
   }
   
   if (control.version==4) 
   {
      has_coulomb2 = true;
      mycoulomb2   = new (std::nothrow) Coulomb2_Operator(mygrid);
   }
}




}
