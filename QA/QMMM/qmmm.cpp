
#include        <iomanip>
#include        <cstring>
#include        <cmath>
#include        <iostream>
#include        <fstream>
#include        <cstdio>
#include        <string>
#include        <unistd.h>
#include        <sys/stat.h>
#include        <string>
#include        <vector>
#include        <set>
#include        <map>
#include        <algorithm>

#include "parsestring.hpp"

#include	"qmmm.hpp"


#define MASTER          0
#define ANGTOBOHR       1.88972687777

/*******************************************
 *                                         *
 *        QMMM_Operator::QMMM_Operator     *
 *                                         *
 *******************************************/
QMMM_Operator::QMMM_Operator(std::string nwinput)
{
   //set up qmmm fragments and sizes
   std::vector<std::string> fragments;
   if  (mystring_contains(mystring_lowercase(nwinput),"fragment"))
   {
      fragments = mystring_split(nwinput,"fragment");
      fragments.erase(fragments.begin());
      nkfrag = fragments.size();
      for (auto & frag: fragments)
      {
         int nfrag0=0;
         frag = mystring_rtrim(mystring_ltrim(mystring_split(frag,"end")[0]));

         if (mystring_contains(mystring_lowercase(frag),"index_start"))
            nfrag0 = mystring_split0(mystring_trim(mystring_split(mystring_split(frag,"index_start")[1],"\n")[0])).size();

         if (mystring_contains(mystring_lowercase(frag),"bond_spring"))
            nbond += (mystring_split(frag,"bond_spring").size() - 1);

         if (mystring_contains(mystring_lowercase(frag),"angle_spring"))
            nangle += (mystring_split(frag,"angle_spring").size() - 1);

         nfrag += nfrag0;
      }
   }

   // set QM and MM LJ parameters
   if  (mystring_contains(mystring_lowercase(nwinput),"lj_qm_parameters"))
   {
      std::vector<std::string> ss;
      std::vector<std::string> lj_fragments = mystring_split(nwinput,"lj_qm_parameters");
      lj_fragments.erase(lj_fragments.begin());
      for (auto & lj_frag: lj_fragments)
      {
         ss = mystring_split0(lj_frag);
         lj_qm_data[ss[0]].push_back(std::stod(ss[1]));
         lj_qm_data[ss[0]].push_back(std::stod(ss[2]));
      }
   }
   if  (mystring_contains(mystring_lowercase(nwinput),"lj_mm_parameters"))
   {
      std::vector<std::string> ss;
      std::vector<std::string> lj_fragments = mystring_split(nwinput,"lj_mm_parameters");
      lj_fragments.erase(lj_fragments.begin());
      for (auto & lj_frag: lj_fragments)
      {
         ss = mystring_split0(lj_frag);
         lj_mm_data[ss[0]].push_back(std::stod(ss[1]));
         lj_mm_data[ss[0]].push_back(std::stod(ss[2]));
      }
   }


   // switching parameters
   switch_Rin  = new (std::nothrow) double [nkfrag]; //  = { (2.0160/0.529177) };
   switch_Rout = new (std::nothrow) double [nkfrag]; //  = { (3.1287/0.529177) };

   // Allocate QMMM indexing arrays 
   size_frag        = new (std::nothrow) int [nkfrag];
   nfrag_frag       = new (std::nothrow) int [nkfrag];
   nbond_frag       = new (std::nothrow) int [nkfrag];
   nangle_frag      = new (std::nothrow) int [nkfrag];
   bond_start_frag  = new (std::nothrow) int [nkfrag];
   angle_start_frag = new (std::nothrow) int [nkfrag];

   indxfrag_start = new (std::nothrow) int [nfrag];
   kfrag          = new (std::nothrow) int [nfrag];

   // bond and angle spring parameters
   bond_indx  = new (std::nothrow) int [2*nbond],
   angle_indx = new (std::nothrow) int [3*nangle];
   Krb        = new (std::nothrow) double [2*nbond];
   Kra        = new (std::nothrow) double [2*nangle];

   double pi = 4.0*std::atan(1.0);

   int bond_start_sum  = 0;
   int angle_start_sum = 0;
   int ntf = 0;
   int ia  = 0;
   for (auto & frag: fragments)
   {
      frag = mystring_rtrim(mystring_ltrim(mystring_split(frag,"end")[0]));

      int tsize = std::stoi(mystring_split0(mystring_split(frag,"size")[1])[0]);

      size_frag[ia] = tsize;

      if (mystring_contains(mystring_lowercase(frag),"index_start"))
      {
         std::vector<std::string> ss = mystring_split0(mystring_trim(mystring_split(mystring_split(frag,"index_start")[1],"\n")[0]));
         nfrag_frag[ia] = ss.size();
         for (auto i=0; i<ss.size(); ++i)
         {
            indxfrag_start[ntf] = std::stoi(ss[i])-1;
            kfrag[ntf] = ia;
            ++ntf;
         }
      }

      if (mystring_contains(mystring_lowercase(frag),"bond_spring"))
      {
         std::string line;
         std::vector<std::string> ss;
         nbond_frag[ia]      = (mystring_split(frag,"bond_spring").size() - 1);
         bond_start_frag[ia] = bond_start_sum;
         bond_start_sum     += nbond_frag[ia];
         for (auto b=0; b<nbond_frag[ia]; ++b)
         {
            int b1 = b + bond_start_frag[ia];
            line = mystring_split(frag,"bond_spring")[b+1];
            ss   = mystring_split0(line);
            int i = std::stoi(ss[0]);
            int j = std::stoi(ss[1]);
            double k = std::stod(ss[2]);
            double r = std::stod(ss[3]);

            bond_indx[2*b1  ] = i-1;
            bond_indx[2*b1+1] = j-1;
            Krb[2*b1  ] = k/27.2116/23.06/ANGTOBOHR/ANGTOBOHR;
            Krb[2*b1+1] = r*ANGTOBOHR;
         }
      }

      if (mystring_contains(mystring_lowercase(frag),"angle_spring"))
      {
         std::string line;
         std::vector<std::string> ss;
         nangle_frag[ia]      = (mystring_split(frag,"angle_spring").size() - 1);
         angle_start_frag[ia] = angle_start_sum;
         angle_start_sum     += nangle_frag[ia];
         for (auto a=0; a<nangle_frag[ia]; ++a)
         {
            int a1 = a + angle_start_frag[ia];
            line = mystring_split(frag,"angle_spring")[a+1];
            ss   = mystring_split0(line);
            int i = std::stoi(ss[0]);
            int j = std::stoi(ss[1]);
            int k = std::stoi(ss[2]);
            double ks    = std::stod(ss[3]);
            double theta = std::stod(ss[4]);
            std::cout << "i=" << i << " j=" << j << " k=" << k
                      << " ks =" << ks << " theta=" << theta << std::endl;
            angle_indx[3*a1  ] = i-1;
            angle_indx[3*a1+1] = j-1;
            angle_indx[3*a1+2] = k-1;
            Kra[2*a1  ] = ks/27.2116/23.06;
            Kra[2*a1+1] = theta*pi/180.0;
         }
      }
      ++ia;
   }
}


