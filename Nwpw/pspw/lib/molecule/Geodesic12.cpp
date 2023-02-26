/* Geodesic12.cpp - 
   Author - Eric Bylaska
*/

#include	"Geodesic12.hpp"

namespace pwdft {

/*******************************************
 *                                         *
 *        Geodesic12::Geodesic12           *
 *                                         *
 *******************************************/
Geodesic12::Geodesic12(int minimizer0, Molecule *mymolecule0, Control2& control)
{
   int minimizer = control.minimizer();
   if ((minimizer==1) || (minimizer==2))
   {
      has_geodesic1 = true;
      mygeodesic1   = new (std::nothrow) Geodesic(minimizer0,mymolecule0);
   }
   
   if ((minimizer==4) || (minimizer==7))
   {
      has_geodesic2 = true;
      mygeodesic2   = new (std::nothrow) Geodesic2(minimizer0,mymolecule0);
   }
}




}
