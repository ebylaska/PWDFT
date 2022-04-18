/* Ions.C - 
   Author - Eric Bylaska
*/


#include        <iostream>
#include        <cstdio>
#include        <cmath>
#include        <cstdlib>
#include	<cstring>
using namespace std;

#include "json.hpp"
using json = nlohmann::json;

#include	"Control2.hpp"
#include	"Ion.hpp"

namespace pwdft {
using namespace pwdft;

static void center_v_mass(int nion, double *mass, double *rion0, double *vx, double *vy, double *vz)
{
   double tmass = 0.0;
   double sx    = 0.0;
   double sy    = 0.0;
   double sz    = 0.0;
   for (auto ii=0; ii<nion; ++ii)
   {
      tmass += mass[ii];
      sx    += mass[ii]*rion0[3*ii];
      sy    += mass[ii]*rion0[3*ii+1];
      sz    += mass[ii]*rion0[3*ii+2];
   }
   *vx = sx/tmass;
   *vy = sy/tmass;
   *vz = sz/tmass;
}


/*******************************
 *                             *
 *        ion_find_nkatm       *
 *                             *
 *******************************/
static void ion_tmpsymbols(RTDB& myrtdb, const int nion, char *symb, char *symtmp)
{
   int j,ii,i,matype,nelem;
   char date[50];

   if (!myrtdb.get_info("geometry:geometry:tags",&matype,&nelem,date))
   { printf("rtdb.get_info error: geometry:geometry:tags not found\n"); exit(99);}

   else
   {
      if (!myrtdb.get("geometry:geometry:tags",matype,nelem,symtmp))
      { printf("rtdb_get error: geometry:geometry:tags not found\n"); exit(99);}

      j = 0;
      for (ii=0; ii<nion; ++ii)
      {
        i = 3*ii;
        while ((symtmp[j]!= ((char) 0)) && (symtmp[j]!= ((char) 10)))
           symb[i++] = symtmp[j++];
        symb[i] = '\0';
        if (symtmp[j]==((char) 10)) ++j;
      }
   }
}

static int ion_find_nkatm(RTDB& myrtdb, const int nion)
{
   int nkatm,i,ii,j,found;
   char *symtmp,*symb;
   symtmp = new char[3*nion];
   symb   = new char[3*nion];

   ion_tmpsymbols(myrtdb,nion,symb,symtmp);

   /* determine nkatm */
   nkatm = 0;
   for (ii=0; ii<nion; ++ii)
   {
      found = 0;
      for (i=0; i<nkatm; ++i)
         if (strcmp(&symb[3*ii],&symtmp[3*i])==0) found=1;
      if (!found)
      {
         strcpy(&symtmp[3*nkatm],&symb[3*ii]);
         ++nkatm;
      }
   }

   delete [] symb;
   delete [] symtmp;

   return nkatm;
}

/*******************************
 *                             *
 *    ion_find_atomkatmnatm    *
 *                             *
 *******************************/
static void ion_find_atomkatmnatm(RTDB& myrtdb, const int nion, const int nkatm,
                               char *atom, int *katm, int *natm)
{
   int i,ii,j,jj,found;
   char *symtmp,*symb;
   symtmp = new char[3*nion];
   symb   = new char[3*nion];
   ion_tmpsymbols(myrtdb,nion,symb,symtmp);

   /* determine atom */
   jj = 0;
   for (ii=0; ii<nion; ++ii)
   {
      found = 0;
      for (i=0; i<jj; ++i)
         if (strcmp(&symb[3*ii],&atom[3*i])==0) found=1;
      if (!found)
      {
         strcpy(&atom[3*jj],&symb[3*ii]);
         ++jj;
      }
   }

   /* determine katm */
   for (ii=0; ii<nion; ++ii)
   {
      j = 0;
      for (i=0; i<nkatm; ++i)
         if (strcmp(&symb[3*ii],&atom[3*i])==0) j=i;
      katm[ii] = j;
   }

   /* determine natm */
   for (i=0; i<nkatm; ++i)   natm[i] = 0;
   for (ii=0; ii<nion; ++ii) natm[katm[ii]] += 1;
}



/* Constructors */

/*********************************
 *                               *
 *          Ion::Ion             *
 *                               *
 *********************************/
Ion::Ion(RTDB& myrtdb, Control2& control) 
{

   int matype,nelem,ii,i,j,found;
   char *symtmp,*symb;
   char date[50];
   double amu_to_mass = 1822.89;

   /* get parallel mappings */
   if (!myrtdb.get("geometry:geometry:ncenter",rtdb_int,1,&nion)) nion = 1;

   nkatm = ion_find_nkatm(myrtdb,nion);

   time_step = control.time_step();

   charge    = new double[nion];
   mass      = new double[nion];
   dti       = new double[nion];
   rion0     = new double[3*nion];
   rion1     = new double[3*nion];
   rion2     = new double[3*nion];
   katm      = new int[nion];
   natm      = new int[nkatm];
   atomarray = new char[3*nkatm];
   zv_psp    = new double[nkatm];

   ion_find_atomkatmnatm(myrtdb,nion,nkatm,atomarray,katm,natm);

   /*  read in ion positions */
   if (!myrtdb.get("geometry:geometry:coords",rtdb_double,3*nion,rion1)) 
   { rion1[0] = 0.0; rion1[1] = 0.0; rion1[2] = 0.0; }
   if (!myrtdb.get("geometry:geometry:coords",rtdb_double,3*nion,rion2)) 
   { rion2[0] = 0.0; rion2[1] = 0.0; rion2[2] = 0.0; }
   if (!myrtdb.get("geometry:geometry:coords",rtdb_double,3*nion,rion0)) 
   { rion0[0] = 0.0; rion0[1] = 0.0; rion0[2] = 0.0; }

   /*  read in masses and charges */
   if (!myrtdb.get("geometry:geometry:masses",rtdb_double,nion,mass)) mass[0] = 1.0;
   if (!myrtdb.get("geometry:geometry:charges",rtdb_double,nion,charge)) charge[0] = 1.0;

   for (ii=0; ii<nion; ++ii)
   {
      mass[ii] *= amu_to_mass;
      dti[ii]  = (time_step*time_step)/mass[ii];
   }

}

Ion::Ion(string rtdbstring, Control2& control) 
{
   double amu_to_mass = 1822.89;

   auto rtdbjson =  json::parse(rtdbstring);


   string geomname = "geometry";
   if (rtdbjson["geometry"].is_string())
      geomname = rtdbjson["geometry"];


   json geomjson = rtdbjson["geometries"][geomname];

   nion = geomjson["nion"];
   auto symbols = geomjson["symbols"];

   vector<string> tmpsymbols;
   for (auto i=0; i<symbols.size(); ++i)
   {
      auto match = std::find( begin(tmpsymbols), end(tmpsymbols), symbols[i] );
      if (match == end(tmpsymbols)) tmpsymbols.push_back(symbols[i]);
   }
   nkatm = tmpsymbols.size();

   time_step = control.time_step();

   charge    = new double[nion];
   mass      = new double[nion];
   dti       = new double[nion];
   rion0     = new double[3*nion];
   rion1     = new double[3*nion];
   rion2     = new double[3*nion];
   katm      = new int[nion];
   natm      = new int[nkatm];
   atomarray = new char[3*nkatm];
   zv_psp    = new double[nkatm];

   for (auto ia=0; ia<nkatm; ++ia)
   {
       natm[ia] = 0.0;
       strcpy(&atomarray[3*ia],const_cast<char*>(tmpsymbols[ia].data()));
   }
   for (auto i=0; i<nion; ++i)
   {
      charge[i] = (double) geomjson["charges"][i];
      mass[i]   = ((double) geomjson["masses"][i])*amu_to_mass;
      dti[i]    = (time_step*time_step)/mass[i];

      rion0[3*i]  =(geomjson["velocities"][3*i].is_number_float())   ? (double) geomjson["velocities"][3*i]   : 0.0;
      rion0[3*i+1]=(geomjson["velocities"][3*i+1].is_number_float()) ? (double) geomjson["velocities"][3*i+1] : 0.0;
      rion0[3*i+2]=(geomjson["velocities"][3*i+2].is_number_float()) ? (double) geomjson["velocities"][3*i+2] : 0.0;

      rion1[3*i]   = (double) geomjson["coords"][3*i];
      rion1[3*i+1] = (double) geomjson["coords"][3*i+1];
      rion1[3*i+2] = (double) geomjson["coords"][3*i+2];
      rion2[3*i]   = (double) geomjson["coords"][3*i];
      rion2[3*i+1] = (double) geomjson["coords"][3*i+1];
      rion2[3*i+2] = (double) geomjson["coords"][3*i+2];

      auto match = std::find( begin(tmpsymbols), end(tmpsymbols), symbols[i] );
      if (match != end(tmpsymbols)) 
      {
         auto ia = std::distance( begin(tmpsymbols),match);
         katm[i] = ia;
         natm[ia] += 1;
      }
   }

   // generate random initial velocities
   double vgx,vgy,vgz,rr0,rr1,rr2,rr3,rr4,rr5;
   double twopi = 16.0*atan(1.0);
   int seed = -1;
   double Tf =  -1.0;
   double kb = 3.16679e-6;

   if (rtdbjson["nwpw"]["car-parrinello"]["initial_velocities"][0].is_number_float())
      Tf =  rtdbjson["nwpw"]["car-parrinello"]["initial_velocities"][0];
   else if (rtdbjson["nwpw"]["initial_velocities"][0].is_number_float())
      Tf =  rtdbjson["nwpw"]["initial_velocities"][0];

   if (rtdbjson["nwpw"]["car-parrinello"]["initial_velocities"][1].is_number_integer())
      seed =  rtdbjson["nwpw"]["car-parrinello"]["initial_velocities"][1];
   else if (rtdbjson["nwpw"]["initial_velocities"][1].is_number_integer())
      seed =  rtdbjson["nwpw"]["initial_velocities"][1];

   if ((Tf>=0.0) && (seed>0))
   {
      std::srand(seed);
      for (auto i=0; i<nion; ++i)
      {
         rr0 = ((double) std::rand())/((double) RAND_MAX);
         rr1 = ((double) std::rand())/((double) RAND_MAX);
         rr2 = ((double) std::rand())/((double) RAND_MAX);
         rr3 = ((double) std::rand())/((double) RAND_MAX);
         rr4 = ((double) std::rand())/((double) RAND_MAX);
         rr5 = ((double) std::rand())/((double) RAND_MAX);
         std::cout << "RANDS=" << rr0 << " " << rr1;
         std::cout <<      " " << rr2 << " " << rr3;
         std::cout <<      " " << rr4 << " " << rr5 << std::endl;
         std::cout <<      " seed=" << seed;
         std::cout <<      " Tf=" << Tf << std::endl;

         vgx = -(2.00*kb*Tf/mass[i])*log(rr0);
         vgy = cos(twopi*rr1);
         rion0[3*i] = sqrt(vgx)*vgy;

         vgx = -(2.00*kb*Tf/mass[i])*log(rr2);
         vgy = cos(twopi*rr3);
         rion0[3*i+1] = sqrt(vgx)*vgy;

         vgx = -(2.00*kb*Tf/mass[i])*log(rr4);
         vgy = cos(twopi*rr5);
         rion0[3*i+2] = sqrt(vgx)*vgy;
      }

      // rescale velocities
      center_v_mass(nion,mass,rion0,&vgx,&vgy,&vgz);
      for (auto i=0; i<nion; ++i)
      {
         rion0[3*i]   -= vgx;
         rion0[3*i+1] -= vgy;
         rion0[3*i+2] -= vgz;
      }
   }
 
      

/*  DEBUG CHECK
   std::cout << "NION=" << nion << std::endl;
   std::cout << "NKATM=" << nkatm << std::endl;
   for (auto i=0; i<nion; ++i)
   {
      std::cout << "I=" << i << std::endl;
      std::cout << "   KATM=" << katm[i] << std::endl;
      std::cout << "   MASS=" << mass[i] << std::endl;
      std::cout << "   CHARGE=" << charge[i] << std::endl;
   }
   std::cout << "NATM=" << natm[0] << std::endl;
   std::cout << "ATOMARRAY=" << atomarray << std::endl;
*/
}


/*******************************
 *                             *
 *      Ion::writejsonstr      *
 *                             *
 *******************************/

void Ion::writejsonstr(string& rtdbstring) 
{
   auto rtdbjson =  json::parse(rtdbstring);

   string geomname = "geometry";
   if (rtdbjson["geometry"].is_string())
      geomname = rtdbjson["geometry"];

   vector<double> coords;
   for (auto i=0; i<(3*nion); ++i)
      coords.push_back(rion1[i]);

   rtdbjson["geometries"][geomname]["coords"] = coords;

   rtdbstring = rtdbjson.dump();
}

}
