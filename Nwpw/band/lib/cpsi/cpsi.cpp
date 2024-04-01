

#include <cmath>
#include <cstring> //memset
#include <iostream>
//#include	"control.hpp"
// extern "C" {
//#include        "compressed_io.h"
//}
#include "compressed_io.hpp"

#include "Parallel.hpp"
#include "Cneb.hpp"
#include "cpsi.hpp"

namespace pwdft {

/*****************************************************
 *                                                   *
 *               wvfnc_expander_convert              *
 *                                                   *
 *****************************************************/

static void cwvfnc_expander_convert(int ngrid[], double *psi1, int dngrid[], double *psi2) 
{
   int indx, dindx, i2, j2, k2;
   int nfft3d = (ngrid[0]) * ngrid[1] * ngrid[2];
   int dnfft3d = (dngrid[0]) * dngrid[1] * dngrid[2];
   int n2ft3d = 2 * nfft3d;
   int dn2ft3d = 2 * dnfft3d;
   int inc2 = (ngrid[0]);
   int dinc2 = (dngrid[0]);
   int inc3 = (ngrid[0]) * ngrid[1];
   int dinc3 = (dngrid[0]) * dngrid[1];
 
   int n1 = ngrid[0];
   int n2 = ngrid[1];
   int n3 = ngrid[2];
   if (n1 > dngrid[0])
     n1 = dngrid[0];
   if (n2 > dngrid[1])
     n2 = dngrid[1];
   if (n3 > dngrid[2])
     n3 = dngrid[2];
 
   int idiff = dngrid[0] - ngrid[0];
   int jdiff = dngrid[1] - ngrid[1];
   int kdiff = dngrid[2] - ngrid[2];
   bool ireverse = (idiff < 0);
   bool jreverse = (jdiff < 0);
   bool kreverse = (kdiff < 0);
   if (jreverse)
     jdiff = -jdiff;
   if (kreverse)
     kdiff = -kdiff;
 
   std::memset(psi2, 0, dn2ft3d * sizeof(double));
   for (auto k = 0; k < n3; ++k)
     for (auto j = 0; j < n2; ++j)
       for (auto i = 0; i < (n1); ++i) {
         indx = i;
         dindx = i;
 
         if (k < (n3 / 2))
           k2 = k;
         else
           k2 = kdiff + k;
 
         if (j < (n2 / 2))
           j2 = j;
         else
           j2 = jdiff + j;

         if (i < (n1 / 2))
           i2 = i;
         else
           i2 = idiff + i;

         if (ireverse) {
           indx = indx + i2;
           dindx = dindx + i;
         } else {
           indx = indx + i;
           dindx = dindx + i2;
         }
 
         if (jreverse) {
           indx = indx + j2 * inc2;
           dindx = dindx + j * dinc2;
         } else {
           indx = indx + j * inc2;
           dindx = dindx + j2 * dinc2;
         }
 
         if (kreverse) {
           indx = indx + k2 * inc3;
           dindx = dindx + k * dinc3;
         } else {
           indx = indx + k * inc3;
           dindx = dindx + k2 * dinc3;
         }
 
         psi2[2 * dindx] = psi1[2 * indx];
         psi2[2 * dindx + 1] = psi1[2 * indx + 1];
       }
}

/*****************************************************
 *                                                   *
 *                dwvfnc_expander                    *
 *                                                   *
 *****************************************************/
static void cwvfnc_expander(Cneb *mycneb, char *filename, std::ostream &coutput)
{
   Parallel *myparall = mycneb->c3db::parall;
 
   bool lprint = myparall->base_stdio_print;
 
   int version, ispin, occupation, nfft[3], dnfft[3], ne[2], nbrillouin;
   double unita[9], dunita[9];
   char tmpfilename[256];
   strcpy(tmpfilename, filename);
   strcat(tmpfilename, ".wvfnc_expander");
 
   if (myparall->is_master()) {
     openfile(4, filename, "r");
     iread(4, &version, 1);
     iread(4, nfft, 3);
     dread(4, unita, 9);
     iread(4, &ispin, 1);
     iread(4, ne, 2);
     iread(4, &nbrillouin, 1);
     iread(4, &occupation, 1);
 
     dnfft[0] = mycneb->nx;
     dnfft[1] = mycneb->ny;
     dnfft[2] = mycneb->nz;
     dunita[0] = mycneb->lattice->unita1d(0);
     dunita[1] = mycneb->lattice->unita1d(1);
     dunita[2] = mycneb->lattice->unita1d(2);
     dunita[3] = mycneb->lattice->unita1d(3);
     dunita[4] = mycneb->lattice->unita1d(4);
     dunita[5] = mycneb->lattice->unita1d(5);
     dunita[6] = mycneb->lattice->unita1d(6);
     dunita[7] = mycneb->lattice->unita1d(7);
     dunita[8] = mycneb->lattice->unita1d(8);
 
     openfile(6, tmpfilename, "w");
     iwrite(6, &version, 1);
     iwrite(6, dnfft, 3);
     dwrite(6, dunita, 9);
     iwrite(6, &ispin, 1);
     iwrite(6, ne, 2);
     iwrite(6, &nbrillouin, 1);
     iwrite(6, &occupation, 1);
 
     int n2ft3d = (nfft[0] + 2) * nfft[1] * nfft[2];
     int dn2ft3d = (dnfft[0] + 2) * dnfft[1] * dnfft[2];
     double *psi1 = new double[n2ft3d];
     double *psi2 = new double[dn2ft3d];
     for (auto ms=0; ms<ispin; ++ms)
     for (auto n=0; n<ne[ms]; ++n) 
     {
        if (lprint) coutput << " converting .... psi:" << n + 1 << " spin:" << ms + 1 << std::endl;
        dread(4, psi1, n2ft3d);
        cwvfnc_expander_convert(nfft, psi1, dnfft, psi2);
        dwrite(6, psi2, dn2ft3d);
     }
     if (lprint) coutput << std::endl;
     if (occupation > 0) 
     {
        double *occ1 = new double[ne[0] + ne[1]];
        dread(4, occ1, (ne[0] + ne[1]));
        dwrite(6, occ1, (ne[0] + ne[1]));
        delete[] occ1;
     }
 
     delete[] psi1;
     delete[] psi2;
 
     closefile(4);
     closefile(6);
 
     /* copy and delete the tmpfilenane file */
     std::rename(tmpfilename, filename);
     // std::remove(tmpfilename);
   }
}

/*****************************************************
 *                                                   *
 *                cpsi_get_header                    *
 *                                                   *
 *****************************************************/
void cpsi_get_header(Parallel *myparall, int *version, int nfft[],
                     double unita[], int *ispin, int ne[], int *nbrillouin, char *filename) 
{
   if (myparall->is_master()) 
   {
      // char *fname = control_input_movecs_filename();
      openfile(4, filename, "r");
      iread(4, version, 1);
      iread(4, nfft, 3);
      dread(4, unita, 9);
      iread(4, ispin, 1);
      iread(4, ne, 2);
      iread(4, nbrillouin, 1);
      closefile(4);
   }
   myparall->Brdcst_iValue(0, 0, version);
   myparall->Brdcst_iValues(0, 0, 3, nfft);
   myparall->Brdcst_Values(0, 0, 9, unita);
   myparall->Brdcst_iValue(0, 0, ispin);
   myparall->Brdcst_iValues(0, 0, 2, ne);
   myparall->Brdcst_iValue(0, 0, nbrillouin);
}

/*****************************************************
 *                                                   *
 *                cpsi_check_convert                  *
 *                                                   *
 *****************************************************/
static bool cpsi_check_convert(Cneb *mycneb, char *filename, std::ostream &coutput) 
{
   Parallel *myparall = mycneb->c3db::parall;
   int version0, ispin0, nfft0[3], ne0[2], nbrillouin0;
   double unita0[9];
   bool converted = false;
 
   cpsi_get_header(myparall, &version0, nfft0, unita0, &ispin0, ne0, &nbrillouin0, filename);
   if ((nfft0[0] != mycneb->nx) || (nfft0[1] != mycneb->ny) || (nfft0[2] != mycneb->nz)) 
   {
      if (myparall->base_stdio_print)
         coutput << " psi grids are being converted: " << std::endl
                 << " -----------------------------: " << std::endl;
     
      cwvfnc_expander(mycneb, filename, coutput);
      converted = true;
   }
 
   return converted;
}

/*****************************************************
 *                                                   *
 *                cpsi_read0                         *
 *                                                   *
 *****************************************************/
/*
   Just reads psi and its header.

   Note - the file must exist

*/
void cpsi_read0(Cneb *mycneb, int *version, int nfft[], double unita[],
               int *ispin, int ne[], int *nbrillouin, double *psi, char *filename) 
{
   int occupation;
 
   Parallel *myparall = mycneb->c3db::parall;
 
   if (myparall->is_master()) 
   {
      openfile(4, filename, "r");
      iread(4, version, 1);
      iread(4, nfft, 3);
      dread(4, unita, 9);
      iread(4, ispin, 1);
      iread(4, ne, 2);
      iread(4, nbrillouin, 1);
      iread(4, &occupation, 1);
   }
   myparall->Brdcst_iValue(0, 0, version);
   myparall->Brdcst_iValues(0, 0, 3, nfft);
   myparall->Brdcst_Values(0, 0, 9, unita);
   myparall->Brdcst_iValue(0, 0, ispin);
   myparall->Brdcst_iValues(0, 0, 2, ne);
   myparall->Brdcst_iValue(0, 0, nbrillouin);

   myparall->Brdcst_iValue(0, 0, &occupation);
 
   /* reads in c format and automatically packs the result to g format */
   //mycneb->g_read(4,ispin,psi);
   mycneb->g_read_ne(4,ne,*nbrillouin,psi);
   
 
   if (myparall->is_master())
     closefile(4);
}

/*****************************************************
 *                                                   *
 *                cpsi_read                          *
 *                                                   *
 *****************************************************/
/*
   Reads psi and checks the header and check for orthonormalization.
   Entry - mycneb: Cneb grid structure
           filename: input filename
           wvfnc_initialize: force initialize wavefuntion
   Exit - psi2: complex wavefunction

   Uses - psi_filefind, psi_read0,
          mycneb->g_generate_random,
          mycneb->gg_traceall,
          mycneb->g_ortho
*/
bool cpsi_read(Cneb *mycneb, char *filename, bool wvfnc_initialize, double *psi2, std::ostream &coutput) 
{
   nwpw_timing_function ftimer(50);
   int version, ispin, nfft[3], ne[2],nbrillouin;
   double unita[9];
   Parallel *myparall = mycneb->c3db::parall;
   bool newpsi = true;
 
   /* read psi from file if psi_exist and not forcing wavefunction initialization */
   if (cpsi_filefind(mycneb,filename) && (!wvfnc_initialize)) 
   {
      newpsi = cpsi_check_convert(mycneb,filename,coutput);
 
      if (myparall->base_stdio_print)
         coutput << " input psi exists, reading from file: " << filename << std::endl;
 
      cpsi_read0(mycneb, &version, nfft, unita, &ispin, ne, &nbrillouin, psi2, filename);
   }
 
   /* generate new psi */
   else 
   {
      if (myparall->base_stdio_print) coutput << " generating random psi from scratch" << std::endl;
      mycneb->g_generate_random(psi2);
   }
   newpsi = newpsi || (ispin != mycneb->ispin)
                   || (nbrillouin != mycneb->nbrillouin)
                   || (ne[0] != mycneb->ne[0])
                   || (ne[1] != mycneb->ne[1])
                   || (std::abs(unita[0] - mycneb->lattice->unita1d(0)) > 1.0e-4) 
                   || (std::abs(unita[1] - mycneb->lattice->unita1d(1)) > 1.0e-4) 
                   || (std::abs(unita[2] - mycneb->lattice->unita1d(2)) > 1.0e-4) 
                   || (std::abs(unita[3] - mycneb->lattice->unita1d(3)) > 1.0e-4) 
                   || (std::abs(unita[4] - mycneb->lattice->unita1d(4)) > 1.0e-4) 
                   || (std::abs(unita[5] - mycneb->lattice->unita1d(5)) > 1.0e-4) 
                   || (std::abs(unita[6] - mycneb->lattice->unita1d(6)) > 1.0e-4) 
                   || (std::abs(unita[7] - mycneb->lattice->unita1d(7)) > 1.0e-4) 
                   || (std::abs(unita[8] - mycneb->lattice->unita1d(8)) > 1.0e-4);
 
   /* ortho check */
   double sum2 = mycneb->gg_traceall(psi2, psi2);
   double sum1 = mycneb->ne[0] + mycneb->ne[1];
 
   if ((mycneb->ispin) == 1)
      sum1 *= 2;
   if (std::fabs(sum2 - sum1) > 1.0e-10) 
   {
      if (myparall->base_stdio_print)
         coutput << " Warning - Gram-Schmidt being performed on psi2" << std::endl;
      mycneb->g_ortho(psi2);
 
      double sum3 = mycneb->gg_traceall(psi2, psi2);
      if (myparall->base_stdio_print)
        coutput << "         - exact norm = " << sum1 << " norm=" << sum2
                << " corrected norm=" << sum3
                << " (error=" << std::abs(sum2 - sum1) << ")" << std::endl;
   }
 
   return newpsi;
}

/*****************************************************
 *                                                   *
 *                cpsi_write                         *
 *                                                   *
 *****************************************************/
void cpsi_write(Cneb *mycneb, int *version, int nfft[], double unita[],
               int *ispin, int ne[], int *nbrillouin, double *psi, char *filename,
               std::ostream &coutput) 
{
   nwpw_timing_function ftimer(50);
   int occupation = -1;
 
   Parallel *myparall = mycneb->c3db::parall;
 
   if (myparall->base_stdio_print)
     coutput << " output psi to filename: " << filename << std::endl;
 
   if (myparall->is_master()) 
   {
      openfile(6, filename, "w");
      iwrite(6, version, 1);
      iwrite(6, nfft, 3);
      dwrite(6, unita, 9);
      iwrite(6, ispin, 1);
      iwrite(6, ne, 2);
      iwrite(6, nbrillouin, 1);
      iwrite(6, &occupation, 1);
   }
 
   mycneb->g_write(6, psi);

   if (occupation>0)
       std::cout << "Write the occupations here!" << std::endl;
 
   if (myparall->is_master())
      closefile(6);
}

/*****************************************************
 *                                                   *
 *                cpsi_filefind                      *
 *                                                   *
 *****************************************************/
bool cpsi_filefind(Cneb *mycneb, char *filename) 
{
   int ifound = 0;
   Parallel *myparall = mycneb->c3db::parall;
 
   if (myparall->is_master()) 
   {
      ifound = cfileexists(filename);
   }
 
   myparall->Brdcst_iValue(0, 0, &ifound);
 
   return (ifound > 0);
}

/*
void v_psi_read(Cneb *mycneb,int *version, int nfft[],
              double unita[], int *ispin, int ne[],
              double *psi)
{
   int occupation;

   Parallel *myparall = mycneb->c3db::parall;

   if (myparall->is_master())
   {
      openfile(4,control_input_v_movecs_filename(),"r");
      iread(4,version,1);
      iread(4,nfft,3);
      dread(4,unita,9);
      iread(4,ispin,1);
      iread(4,ne,2);
      iread(4,&occupation,1);
   }
   myparall->Brdcst_iValue(0,0,version);
   myparall->Brdcst_iValues(0,0,3,nfft);
   myparall->Brdcst_Values(0,0,9,unita);
   myparall->Brdcst_iValue(0,0,ispin);
   myparall->Brdcst_iValues(0,0,2,ne);

   mycneb->g_read(4,psi);

   if (myparall->is_master()) closefile(4);
}
*/

/*
void v_psi_write(Cneb *mycneb,int *version, int nfft[],
              double unita[], int *ispin, int ne[],
              double *psi)
{
   int occupation = -1;

   Parallel *myparall = mycneb->c3db::parall;

   if (myparall->is_master())
   {
      openfile(6,control_output_v_movecs_filename(),"w");
      iwrite(6,version,1);
      iwrite(6,nfft,3);
      dwrite(6,unita,9);
      iwrite(6,ispin,1);
      iwrite(6,ne,2);
      iwrite(6,&occupation,1);
   }

   mycneb->g_write(6,psi);

   if (myparall->is_master()) closefile(6);
}

*/

} // namespace pwdft
