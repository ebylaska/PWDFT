/* nwpw_Nose_Hoover.cpp -
   Author - Eric Bylaska
*/


#include	<iostream>
#include	<sstream>
#include	<cstring>
#include	<cmath>
//#include	"blas.h"

#include	"nwpw_Nose_Hoover.hpp"

#include	"json.hpp"
using json = nlohmann::json;

namespace pwdft {



/* Constructors */

/*******************************************
 *                                         *
 *   nwpw_Nose_Hoover::nwpw_Nose_Hoover    *
 *                                         *
 *******************************************/
nwpw_Nose_Hoover::nwpw_Nose_Hoover(Ion& myion, int nemax, double eke0, Control2& control)
{

   nose_on = control.Nose();

   if (nose_on)  
   {
      double kb = 3.16679e-6;

      nion = myion.nion;

      Te = control.Nose_Te();
      Tr = control.Nose_Tr();
      double Pe_tmp = control.Nose_Pe();
      double Pr_tmp = control.Nose_Pr();

      nosers  = control.Nose_restart(); 

      Mchain   = control.Nose_Mchain();
      Nchain   = control.Nose_Nchain();
      Ne_chain = control.Nose_Ne_chain();
      if (Ne_chain<=0.0) Ne_chain = 3.0*(nemax);

      eke0_init = control.Nose_eke0();
      if ((eke0_init<=0.0) || (!nosers))  eke0_init = eke0;

      fmass = control.fake_mass();


      Xem = new (std::nothrow) double [Mchain];
      Xe0 = new (std::nothrow) double [Mchain];
      Xe1 = new (std::nothrow) double [Mchain];
      Xe2 = new (std::nothrow) double [Mchain];
      Ee0 = new (std::nothrow) double [Mchain];
      Qe  = new (std::nothrow) double [Mchain];
      Pe  = new (std::nothrow) double [Mchain];
      for (auto m=0; m<Mchain; ++m)
      {
         Xem[m] = control.Nose_Xem(m);
         Xe0[m] = control.Nose_Xe0(m);
         Xe1[m] = control.Nose_Xe1(m);
         Xe2[m] = 0.0;

         Pe[m] = Pe_tmp;
         Qe[m] = control.Nose_Qe(m);
         if (Qe[m]<=0.0) nosers = false;
      }

      Xrm = new (std::nothrow) double [Nchain];
      Xr0 = new (std::nothrow) double [Nchain];
      Xr1 = new (std::nothrow) double [Nchain];
      Xr2 = new (std::nothrow) double [Nchain];
      Er0 = new (std::nothrow) double [Nchain];
      Qr  = new (std::nothrow) double [Nchain];
      Pr  = new (std::nothrow) double [Nchain];
      for (auto n=0; n<Nchain; ++n)
      {
         Xrm[n] = control.Nose_Xem(n);
         Xr0[n] = control.Nose_Xe0(n);
         Xr1[n] = control.Nose_Xe1(n);
         Xr2[n] = 0.0;

         Pr[n]  = Pr_tmp;
         Qr[n]  = control.Nose_Qr(n);
         if (Qr[n]<=0.0) nosers = false;
      }

      if (!nosers)
      {
         for (auto m=0; m<Mchain; ++m) Xem[m] = 0.0;
         for (auto m=0; m<Mchain; ++m) Xe0[m] = 0.0;
         for (auto m=0; m<Mchain; ++m) Xe1[m] = 0.0;

         for (auto n=0; n<Nchain; ++n) Xrm[n] = 0.0;
         for (auto n=0; n<Nchain; ++n) Xr0[n] = 0.0;
         for (auto n=0; n<Nchain; ++n) Xr1[n] = 0.0;
      }

      g_dof = 3.0*nion - 6.0;
      if (g_dof<1) g_dof=1.0;

      Er0[0] = 0.5*g_dof*kb*Tr;
      for (auto n=1; n<Nchain; ++n) Er0[n] = 0.5*kb*Tr;

      Ee0[0] = 4.0*kb*Te*fmass*(nion/mass_total)*eke0_init;;
      double betae = Ee0[0]/Ne_chain;
      for (auto m=1; m<Mchain; ++m) Ee0[m] = betae;

      double pi = 4.0*std::atan(1.0);
      if (!nosers)
      {
         for (auto m=0; m<Mchain; ++m) Qe[m] = Ee0[m]*std::pow((Pe[m]/pi),2);
         for (auto n=0; n<Nchain; ++n) Qr[n] = Er0[n]*std::pow((Pr[n]/pi),2);
      }
      else
      {
         for (auto m=0; m<Mchain; ++m) Pe[m] = pi*std::sqrt(Qe[m])/Ee0[m];
         for (auto n=0; n<Nchain; ++n) Pr[n] = pi*std::sqrt(Qr[n])/Er0[n];
      }

   }
}


/********************************************
 *                                          *
 *      nwpw_Nose_Hoover::writejsonstr      *
 *                                          *
 ********************************************/

void nwpw_Nose_Hoover::writejsonstr(std::string& rtdbstring)
{
   if (nose_on)  
   {
      /* put velocities in Xem and Xrm */
      for (auto m=0; m<Mchain; ++m) Xem[m] = (Xe2[m] - Xe0[m])/(2.0*time_step);
      for (auto n=0; n<Nchain; ++n) Xrm[n] = (Xr2[n] - Xr0[n])/(2.0*time_step);

      /* save restart information into json */
      auto rtdbjson = json::parse(rtdbstring);
      auto nosejson = rtdbjson["nwpw"]["car-parrinello"]["Nose-Hoover"];

      nosejson["eke0"]     = eke0_init;
      nosejson["Mchain"]   = Mchain;
      nosejson["Nchain"]   = Nchain;
      nosejson["Ne_chain"] = Ne_chain;
      nosejson["Te"]       = Te;
      nosejson["Tr"]       = Tr;
      nosejson["restart"]  = true;

      nosejson["Qe"]  = std::vector<double>(Qe, &Qe[Mchain]);
      nosejson["Xe1"] = std::vector<double>(Xe1,&Xe1[Mchain]);
      nosejson["Xe0"] = std::vector<double>(Xe0,&Xe0[Mchain]);
      nosejson["Xem"] = std::vector<double>(Xem,&Xem[Mchain]);

      nosejson["Qr"]  = std::vector<double>(Qr, &Qr[Nchain]);
      nosejson["Xr1"] = std::vector<double>(Xr1,&Xr1[Nchain]);
      nosejson["Xr0"] = std::vector<double>(Xr0,&Xr0[Nchain]);
      nosejson["Xrm"] = std::vector<double>(Xrm,&Xrm[Nchain]);

      rtdbstring = rtdbjson.dump();
   }
}




/*******************************************
 *                                         *
 *      nwpw_Nose_Hoover::shortprint       *
 *                                         *
 *******************************************/
std::string nwpw_Nose_Hoover::shortprint()
{
   std::string outstring = "";

   if (nose_on)  
   {
      std::stringstream stream;

      stream << std::endl; 
      stream << "APC Potential:" << std::endl;
      stream << std::endl
          << std::endl;

      outstring = stream.str();
   }
   return outstring;
}


}
