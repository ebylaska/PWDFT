

#include "Paw_gaunt.hpp"
#include "util_gaunt.hpp"


namespace pwdft {


static bool has_complex_gaunt = false;
static bool has_real_gaunt = false;
static Paw_gaunt  *complex_gaunt;
static Paw_gaunt  *real_gaunt;



/*******************************************
 *                                         *
 *            util_gaunt_init              *
 *                                         *
 *******************************************/
void util_gaunt_init(const bool iscmplx, const int lmax0)
{
   if (iscmplx)
   {
      has_complex_gaunt = true;
      complex_gaunt = new Paw_gaunt(true,lmax0);
   }
   else
   {
      has_real_gaunt = true;
      real_gaunt = new Paw_gaunt(false,lmax0);
   }
}

/*******************************************
 *                                         *
 *            util_gaunt_end               *
 *                                         *
 *******************************************/
void util_gaunt_end()
{
  if (has_complex_gaunt)
  {
     delete complex_gaunt;
     has_complex_gaunt = false;
  }

  if (has_real_gaunt)
  {
     delete real_gaunt;
     has_real_gaunt = false;
  }

}

/*******************************************
 *                                         *
 *            util_gaunt                   *
 *                                         *
 *******************************************/
double util_gaunt(const bool iscmplx,
                const int l,  const int m,
                const int l1, const int m1,
                const int l2, const int m2)
{
   if (iscmplx)
   {
      return complex_gaunt->gaunt(l,m,l1,m1,l2,m2);
   }
   else
   {
      return real_gaunt->gaunt(l,m,l1,m1,l2,m2);
   }
}

/*******************************************
 *                                         *
 *            util_gaunt2                  *
 *                                         *
 *******************************************/
double util_gaunt2(const bool iscmplx,
                const int l,  const int m,
                const int l1, const int m1,
                const int l2, const int m2)
{
   if (iscmplx)
      return complex_gaunt->gaunt2(l,m,l1,m1,l2,m2);
   else
      return real_gaunt->gaunt2(l,m,l1,m1,l2,m2);
}


/*******************************************
 *                                         *
 *            util_gaunt3                  *
 *                                         *
 *******************************************/
double util_gaunt3(const bool iscmplx,
                const int l,  const int m,
                const int l1, const int m1,
                const int l2, const int m2)
{
   if (iscmplx)
      return complex_gaunt->gaunt3(l,m,l1,m1,l2,m2);
   else
      return real_gaunt->gaunt3(l,m,l1,m1,l2,m2);
}



}


