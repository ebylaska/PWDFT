#include        <iostream>
#include        <cstdio>
#include        <cmath>
#include        <cstdlib>
#include        <string>

#include        "Parallel.hpp"
#include        "Ion.hpp"
#include        "Pneb.hpp"
#include        "Molecule.hpp"

/******************************************
 *                                        *
 *            cgsd_noit_energy            *
 *                                        *
 ******************************************/
double cgsd_noit_energy(Molecule& mymolecule) 
{
   Parallel *parall = mymolecule.mygrid->d3db::parall;
   Pneb *mygrid = mymolecule.mygrid;
   Ion  *myion  = mymolecule.myion;

   int    ne[2],ispin,nion;
   double E[50],total_energy;

   nion  = mymolecule.myion->nion;
   ispin = mymolecule.mygrid->ispin;
   ne[0] = mymolecule.mygrid->ne[0];
   ne[1] = mymolecule.mygrid->ne[1];

   total_energy  = mymolecule.energy_eigenvalues();

   E[0] = total_energy;
   E[1] = mymolecule.eorbit();
   E[2] = mymolecule.ehartree();
   E[3] = mymolecule.exc();
   E[4] = mymolecule.eion();

   E[5] = mymolecule.eke();
   E[6] = mymolecule.vl_ave();
   E[7] = mymolecule.vnl_ave();
   E[8] = 2*E[2];
   E[9] = mymolecule.pxc();


   //                 |***************************|
   //****************** report summary of results **********************
   //                 |***************************|
   if (parall->is_master())
   {
      std::cout << std::endl << std::endl;
      std::cout << "          =============  summary of results  =================" << std::endl;;
      cout << "\n\n";
      printf(" total     energy    : %19.10le (%15.5le /ion)\n",      E[0],E[0]/nion);
      printf(" total orbital energy: %19.10le (%15.5le /electron)\n", E[1],E[1]/(ne[0]+ne[1]));
      printf(" hartree energy      : %19.10le (%15.5le /electron)\n", E[2],E[2]/(ne[0]+ne[1]));
      printf(" exc-corr energy     : %19.10le (%15.5le /electron)\n", E[3],E[3]/(ne[0]+ne[1]));
      printf(" ion-ion energy      : %19.10le (%15.5le /ion)\n",      E[4],E[4]/nion);
      printf("\n");
      printf(" K.S. kinetic energy : %19.10le (%15.5le /electron)\n",      E[5],E[5]/(ne[0]+ne[1]));
      printf(" K.S. V_l energy     : %19.10le (%15.5le /electron)\n",      E[6],E[6]/(ne[0]+ne[1]));
      printf(" K.S. V_nl energy    : %19.10le (%15.5le /electron)\n",      E[7],E[7]/(ne[0]+ne[1]));
      printf(" K.S. V_Hart energy  : %19.10le (%15.5le /electron)\n",      E[8],E[8]/(ne[0]+ne[1]));
      printf(" K.S. V_xc energy    : %19.10le (%15.5le /electron)\n",      E[9],E[9]/(ne[0]+ne[1]));
      double viral = (E[9]+E[8]+E[7]+E[6])/E[5];
      printf(" Viral Coefficient   : %19.10le\n",viral);

      printf("\n orbital energies:\n");
      int nn = ne[0] - ne[1];
      double ev = 27.2116;
      for (int i=0; i<nn; ++i)
      {
         printf("%18.7le",mymolecule.eig[i]); printf(" ("); printf("%8.3f",mymolecule.eig[i]*ev); printf("eV)\n");
      }
      for (int i=0; i<ne[1]; ++i)
      {
         printf("%18.7le",mymolecule.eig[i+nn]); printf(" ("); printf("%8.3lf",mymolecule.eig[i+nn]*ev); printf("eV) ");
         printf("%18.7le",mymolecule.eig[i+(ispin-1)*ne[0]]); printf(" ("); printf("%8.3lf",mymolecule.eig[i+(ispin-1)*ne[0]]*ev); printf("eV)\n");
      }
   }

   return total_energy;
}
