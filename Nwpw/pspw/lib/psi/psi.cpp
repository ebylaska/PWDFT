

/*
#include        <cmath>
#include        <cstdio>
#include        <cstdlib>
#include        <stdio.h>

*/

#include <cmath>
#include <cstring> //memset
#include <iostream>
#include <string>
#include <limits>
//#include	"control.hpp"
// extern "C" {
//#include        "compressed_io.h"
//}
#include "compressed_io.hpp"

#include "Parallel.hpp"
#include "Pneb.hpp"
#include "psi.hpp"
#include "../../nwpwlib/utilities/util.hpp"

namespace pwdft {

/*****************************************************
 *                                                   *
 *               wvfnc_expander_convert              *
 *                                                   *
 *****************************************************/

static void wvfnc_expander_convert(int ngrid[], double *psi1, int dngrid[], double *psi2) 
{
   int indx, dindx, j2, k2;
   int nfft3d = (ngrid[0] / 2 + 1) * ngrid[1] * ngrid[2];
   int dnfft3d = (dngrid[0] / 2 + 1) * dngrid[1] * dngrid[2];
   int n2ft3d = 2 * nfft3d;
   int dn2ft3d = 2 * dnfft3d;
   int inc2 = (ngrid[0] / 2 + 1);
   int dinc2 = (dngrid[0] / 2 + 1);
   int inc3 = (ngrid[0] / 2 + 1) * ngrid[1];
   int dinc3 = (dngrid[0] / 2 + 1) * dngrid[1];
 
   int n1 = ngrid[0];
   int n2 = ngrid[1];
   int n3 = ngrid[2];
   if (n1 > dngrid[0])
      n1 = dngrid[0];
   if (n2 > dngrid[1])
      n2 = dngrid[1];
   if (n3 > dngrid[2])
      n3 = dngrid[2];
 
   int jdiff = dngrid[1] - ngrid[1];
   int kdiff = dngrid[2] - ngrid[2];
   bool jreverse = (jdiff < 0);
   bool kreverse = (kdiff < 0);
   if (jreverse)
      jdiff = -jdiff;
   if (kreverse)
      kdiff = -kdiff;

   std::memset(psi2, 0, dn2ft3d * sizeof(double));
   for (auto k = 0; k < n3; ++k)
   for (auto j = 0; j < n2; ++j)
   for (auto i = 0; i < (n1/2+1); ++i) 
   {
      indx = i;
      dindx = i;

      if (k < (n3/2))
         k2 = k;
      else
         k2 = kdiff + k;

      if (j < (n2/2))
         j2 = j;
      else
         j2 = jdiff + j;

      if (jreverse) 
      {
         indx = indx   + j2*inc2;
         dindx = dindx + j*dinc2;
      } 
      else 
      {
         indx = indx   + j*inc2;
         dindx = dindx + j2*dinc2;
      }

      if (kreverse) 
      {
         indx = indx   + k2*inc3;
         dindx = dindx + k*dinc3;
      } 
      else 
      {
         indx = indx   + k*inc3;
         dindx = dindx + k2*dinc3;
      }

      psi2[2*dindx]   = psi1[2*indx];
      psi2[2*dindx+1] = psi1[2*indx+1];
   }
}

/*****************************************************
 *                                                   *
 *                wvfnc_expander                     *
 *                                                   *
 *****************************************************/
static void wvfnc_expander(Pneb *mypneb, char *filename, std::ostream &coutput)
{
   Parallel *myparall = mypneb->d3db::parall;
 
   bool lprint = myparall->base_stdio_print;
 
   int version, ispin, occupation, nfft[3], dnfft[3], ne[2];
   double unita[9], dunita[9];
   char tmpfilename[256];
   strcpy(tmpfilename, filename);
   strcat(tmpfilename, ".wvfnc_expander");
 
   if (myparall->is_master())
   {
      openfile(4, filename, "r");
      iread(4, &version, 1);
      iread(4, nfft, 3);
      dread(4, unita, 9);
      iread(4, &ispin, 1);
      iread(4, ne, 2);
      iread(4, &occupation, 1);
     
      dnfft[0] = mypneb->nx;
      dnfft[1] = mypneb->ny;
      dnfft[2] = mypneb->nz;
      dunita[0] = mypneb->lattice->unita1d(0);
      dunita[1] = mypneb->lattice->unita1d(1);
      dunita[2] = mypneb->lattice->unita1d(2);
      dunita[3] = mypneb->lattice->unita1d(3);
      dunita[4] = mypneb->lattice->unita1d(4);
      dunita[5] = mypneb->lattice->unita1d(5);
      dunita[6] = mypneb->lattice->unita1d(6);
      dunita[7] = mypneb->lattice->unita1d(7);
      dunita[8] = mypneb->lattice->unita1d(8);
     
      openfile(6, tmpfilename, "w");
      iwrite(6, &version, 1);
      iwrite(6, dnfft, 3);
      dwrite(6, dunita, 9);
      iwrite(6, &ispin, 1);
      iwrite(6, ne, 2);
      iwrite(6, &occupation, 1);
     
      int n2ft3d = (nfft[0] + 2) * nfft[1] * nfft[2];
      int dn2ft3d = (dnfft[0] + 2) * dnfft[1] * dnfft[2];
      double *psi1 = new double[n2ft3d];
      double *psi2 = new double[dn2ft3d];
      for (auto ms = 0; ms < ispin; ++ms)
      for (auto n = 0;  n  < ne[ms]; ++n) 
      {
         if (lprint) coutput << " converting .... psi:" << n + 1 << " spin:" << ms + 1 << std::endl;
         dread(4, psi1, n2ft3d);
         wvfnc_expander_convert(nfft, psi1, dnfft, psi2);
         dwrite(6, psi2, dn2ft3d);
      }
      if (lprint) coutput << std::endl;
      if (occupation > 0) {
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

static bool isDescending(const double* arr, int size) {
    for (int i = 0; i < size - 1; ++i) {
        if (arr[i] < arr[i + 1]) {
            return false;
        }
    }
    return true;
}

/*****************************************************
 *                                                   *
 *                psi_get_header                     *
 *                                                   *
 *****************************************************/
void psi_get_header(Parallel *myparall, int *version, int nfft[],
                    double unita[], int *ispin, int ne[], char *filename) 
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
      closefile(4);
   }
   myparall->Brdcst_iValue(0, 0, version);
   myparall->Brdcst_iValues(0, 0, 3, nfft);
   myparall->Brdcst_Values(0, 0, 9, unita);
   myparall->Brdcst_iValue(0, 0, ispin);
   myparall->Brdcst_iValues(0, 0, 2, ne);
}

/*****************************************************
 *                                                   *
 *                psi_get_header_occupation          *
 *                                                   *
 *****************************************************/
/*   Reads header information from a file and broadcasts
 *   the data across all parallel processes. The master
 *   process performs the file reading, and the data
 *   includes version, grid dimensions, unit cell
 *   dimensions, spin information, electron counts,
 *   and occupation details.
 */
 void psi_get_header_occupation(Parallel *myparall, int *version, int nfft[],
                                double unita[], int *ispin, int ne[], int *occupation, 
                                char *filename)
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
      iread(4, occupation, 1);
      closefile(4);
   }
   myparall->Brdcst_iValue(0, 0, version);
   myparall->Brdcst_iValues(0, 0, 3, nfft);
   myparall->Brdcst_Values(0, 0, 9, unita);
   myparall->Brdcst_iValue(0, 0, ispin);
   myparall->Brdcst_iValues(0, 0, 2, ne);
   myparall->Brdcst_iValue(0, 0, occupation);
}

/*****************************************************
 *                                                   *
 *                psi_check_convert                  *
 *                                                   *
 *****************************************************/
static bool psi_check_convert(Pneb *mypneb, char *filename, std::ostream &coutput) 
{
   Parallel *myparall = mypneb->d3db::parall;
   int version0, ispin0, nfft0[3], ne0[2];
   double unita0[9];
   bool converted = false;
 
   psi_get_header(myparall, &version0, nfft0, unita0, &ispin0, ne0, filename);
   if ((nfft0[0] != mypneb->nx) || (nfft0[1] != mypneb->ny) || (nfft0[2] != mypneb->nz)) 
   {
      if (myparall->base_stdio_print)
         coutput << " psi grids are being converted: " << std::endl
                 << " -----------------------------: " << std::endl;
     
      wvfnc_expander(mypneb, filename, coutput);
      converted = true;
   }
 
   return converted;
}

// Helper: check for NaN/Inf in psi array
static bool psi_has_nan_or_inf(const double* psi, int ncoeff) {
    for (int i = 0; i < ncoeff; ++i) {
        if (std::isnan(psi[i]) || std::isinf(psi[i])) return true;
    }
    return false;
}

/*****************************************************
 *                                                   *
 *                psi_read0                          *
 *                                                   *
 *****************************************************/
/*
   Just reads psi and its header.

   Note - the file must exist

*/
void psi_read0(Pneb *mypneb, int *version, int nfft[], double unita[],
               int *ispin, int ne[], double *psi, 
               int *occupation, double occ[],  char *filename, bool reverse) 
{
   Parallel *myparall = mypneb->d3db::parall;
   int neout[2];
   bool header_mismatch = false;
 
   if (myparall->is_master()) 
   {
      openfile(4, filename, "r");
      iread(4, version, 1);
      iread(4, nfft, 3);
      dread(4, unita, 9);
      iread(4, ispin, 1);
      iread(4, neout, 2);
      iread(4, occupation, 1);
      // Header sanity check
      if (nfft[0] != mypneb->nx || nfft[1] != mypneb->ny || nfft[2] != mypneb->nz ||
          *ispin != mypneb->ispin || neout[0] != mypneb->ne[0] || neout[1] != mypneb->ne[1]) {
          header_mismatch = true;
      }
   }
   myparall->Brdcst_iValue(0, 0, version);
   myparall->Brdcst_iValues(0, 0, 3, nfft);
   myparall->Brdcst_Values(0, 0, 9, unita);
   myparall->Brdcst_iValue(0, 0, ispin);
   myparall->Brdcst_iValues(0, 0, 2, neout);
   myparall->Brdcst_iValue(0, 0, occupation);
   // Fix: Broadcast header_mismatch as int
   int header_mismatch_int = header_mismatch ? 1 : 0;
   myparall->Brdcst_iValue(0, 0, &header_mismatch_int);
   header_mismatch = (header_mismatch_int != 0);

   if (header_mismatch) {
      if (myparall->base_stdio_print)
         std::cout << "[PWDFT] Wavefunction file header does not match current run. Reinitializing psi." << std::endl;
      // Fill psi with nonsense so caller can detect and reinit
      for (int i = 0; i < 2 * (mypneb->ne[0] + mypneb->ne[1]) * mypneb->npack(1); ++i)
          psi[i] = std::numeric_limits<double>::quiet_NaN();
      if (myparall->is_master()) closefile(4);
      return;
   }

 
   /* reads in c format and automatically packs the result to g format */
   //mypneb->g_read(4,ispin,psi);
   if (reverse)
      mypneb->g_read_ne_reverse(4,ne,psi);
   else
      mypneb->g_read_ne(4,ne,psi);


   if (((*occupation) > 0) && (occ))
   {
      if (myparall->is_master()) 
      {
            dread(4,occ, ne[0]+ne[1]);
 
            if (reverse)
            {
               std::reverse(occ, occ + ne[0]); // Reverse the first section if descending
           
               // Check and reverse the second part
               if (*ispin>1)
                  std::reverse(occ + ne[0], occ + ne[0] + ne[1]); // Reverse the second section if descending and *ispin > 1
            }
      }
      myparall->Brdcst_Values(0, 0, ne[0]+ne[1], occ);
   }
   
   if (myparall->is_master())
      closefile(4);

   ne[0] = neout[0];
   ne[1] = neout[1];
}

void psi_read0(Pneb *mypneb, int *version, int nfft[], double unita[],
               int *ispin, int ne[], double *psi,
               char *filename)
{
   int occupation = 0;
   Parallel *myparall = mypneb->d3db::parall;

   if (myparall->is_master()) 
   {
      openfile(4, filename, "r");
      iread(4, version, 1);
      iread(4, nfft, 3);
      dread(4, unita, 9);
      iread(4, ispin, 1);
      iread(4, ne, 2);
      iread(4, &occupation, 1);
   }
   myparall->Brdcst_iValue(0, 0, version);
   myparall->Brdcst_iValues(0, 0, 3, nfft);
   myparall->Brdcst_Values(0, 0, 9, unita);
   myparall->Brdcst_iValue(0, 0, ispin);
   myparall->Brdcst_iValues(0, 0, 2, ne);
   //myparall->Brdcst_iValue(0, 0, occupation);

   /* reads in c format and automatically packs the result to g format */
   //mypneb->g_read(4,ispin,psi);
   mypneb->g_read_ne(4,ne,psi);

   if (myparall->is_master())
      closefile(4);
}



/*****************************************************
 *                                                   *
 *                psi_read                           *
 *                                                   *
 *****************************************************/
/*
   Reads psi and checks the header and check for orthonormalization.
   Entry - mypneb: Pneb grid structure
           filename: input filename
           wvfnc_initialize: force initialize wavefuntion
   Exit - psi2: complex wavefunction

   Uses - psi_filefind, psi_read0,
          mypneb->g_generate_random,
          mypneb->gg_traceall,
          mypneb->g_ortho
*/
bool psi_read(Pneb *mypneb, char *filename, bool wvfnc_initialize, double *psi2, 
              int *occupation, double occ2[], 
              std::ostream &coutput) 
{
   nwpw_timing_function ftimer(50);
   int version, ispin, nfft[3], ne[2];
   double unita[9];
   Parallel *myparall = mypneb->d3db::parall;
   bool newpsi = true;
   bool bad_restart = false;
 
   /* read psi from file if psi_exist and not forcing wavefunction initialization */
   if (psi_filefind(mypneb,filename) && (!wvfnc_initialize)) 
   {
      newpsi = psi_check_convert(mypneb,filename,coutput);
 
      if (myparall->base_stdio_print)
         coutput << " input psi exists, reading from file: " << filename;
 
      ne[0] = mypneb->ne[0];
      ne[1] = mypneb->ne[1];
      psi_read0(mypneb, &version, nfft, unita, &ispin, ne, psi2, occupation, occ2, filename,false);
      // Sanity check: header mismatch or NaN/Inf in psi2
      int ncoeff = 2 * (mypneb->ne[0] + mypneb->ne[1]) * mypneb->npack(1);
      if (psi_has_nan_or_inf(psi2, ncoeff)) {
         if (myparall->base_stdio_print)
            coutput << "\n[PWDFT] Wavefunction file contains NaN/Inf or header mismatch. Reinitializing psi." << std::endl;
         bad_restart = true;
      } else {
         // Norm check
         double norm = mypneb->gg_traceall(psi2, psi2);
         if (norm < 1e-8) {
            if (myparall->base_stdio_print)
               coutput << "\n[PWDFT] Wavefunction norm is zero or too small. Reinitializing psi." << std::endl;
            bad_restart = true;
         }
      }
      if (bad_restart) {
         std::string guess = get_initial_wavefunction_guess();
         if (guess == "atomic") {
            if (myparall->base_stdio_print) coutput << " generating atomic guess for psi" << std::endl;
            mypneb->g_generate_atomic_guess(psi2); // To be implemented
         } else {
            if (myparall->base_stdio_print) coutput << " generating random psi from scratch" << std::endl;
            mypneb->g_generate_random(psi2);
         }
         newpsi = true;
      }

      if ((*occupation>0) and (myparall->base_stdio_print))
         coutput << " ... reading occupations";

      if (myparall->base_stdio_print)
         coutput << std::endl;

      if ((occ2) && (*occupation > 0))
         if (isDescending(occ2, ne[0]))
         {
            if (myparall->base_stdio_print)
               coutput << " - reversing order of psi and occupation" << std::endl;
            psi_read0(mypneb, &version, nfft, unita, &ispin, ne, psi2, occupation, occ2, filename,true);
         }
   }
 
   /* generate new psi */
   else 
   {
      ispin = mypneb->ispin;
      ne[0] = mypneb->ne[0];
      ne[1] = mypneb->ne[1];
      std::string guess = get_initial_wavefunction_guess();
      if (guess == "atomic") {
         if (myparall->base_stdio_print) coutput << " generating atomic guess for psi" << std::endl;
         mypneb->g_generate_atomic_guess(psi2); // To be implemented
      } else {
         if (myparall->base_stdio_print) coutput << " generating random psi from scratch" << std::endl;
         mypneb->g_generate_random(psi2);
      }
   }

   newpsi = newpsi || (ispin != mypneb->ispin)
                   || (ne[0] != mypneb->ne[0])
                   || (ne[1] != mypneb->ne[1])
                   || (std::abs(unita[0] - mypneb->lattice->unita1d(0)) > 1.0e-4) 
                   || (std::abs(unita[1] - mypneb->lattice->unita1d(1)) > 1.0e-4) 
                   || (std::abs(unita[2] - mypneb->lattice->unita1d(2)) > 1.0e-4) 
                   || (std::abs(unita[3] - mypneb->lattice->unita1d(3)) > 1.0e-4) 
                   || (std::abs(unita[4] - mypneb->lattice->unita1d(4)) > 1.0e-4) 
                   || (std::abs(unita[5] - mypneb->lattice->unita1d(5)) > 1.0e-4) 
                   || (std::abs(unita[6] - mypneb->lattice->unita1d(6)) > 1.0e-4) 
                   || (std::abs(unita[7] - mypneb->lattice->unita1d(7)) > 1.0e-4) 
                   || (std::abs(unita[8] - mypneb->lattice->unita1d(8)) > 1.0e-4);
 
   // Add extra missing orbitals to psi2
   int nextra[2] = {(mypneb->ne[0]-ne[0]),(mypneb->ne[1]-ne[1])};
   for (auto ms=0; ms<ispin; ++ms)
   {
      if (nextra[ms]>0)
      {
         if (myparall->base_stdio_print) 
            coutput << " - adding " << nextra[ms] << " orbitals to ms=" << ms << " psi" << std::endl;
         int indxa = 2*mypneb->npack(1)*(ne[ms]+ms*mypneb->ne[0]);
         mypneb->g_generate_extra_random(nextra[ms],psi2+indxa);

         // project out filled psi
         mypneb->g_project_out_filled_extra(ms,nextra,psi2);
      }
   }

   /* ortho check */
   double sum2 = mypneb->gg_traceall(psi2, psi2);
   double sum1 = mypneb->ne[0] + mypneb->ne[1];
 
   if ((mypneb->ispin) == 1)
      sum1 *= 2;
   if (std::fabs(sum2 - sum1) > 1.0e-10) 
   {
      if (myparall->base_stdio_print)
         coutput << " Warning - Gram-Schmidt being performed on psi2" << std::endl;
      mypneb->g_ortho(-1,psi2);
 
      double sum3 = mypneb->gg_traceall(psi2, psi2);
      if (myparall->base_stdio_print)
        coutput << "         - exact norm = " << sum1 << " norm=" << sum2
                << " corrected norm=" << sum3
                << " (error=" << std::abs(sum2 - sum1) << ")" << std::endl;
   }
 
   return newpsi;
}

bool psi_read(Pneb *mypneb, char *filename, bool wvfnc_initialize, double *psi2, std::ostream &coutput)
{
   nwpw_timing_function ftimer(50);
   int version, ispin, nfft[3], ne[2];
   double unita[9];
   Parallel *myparall = mypneb->d3db::parall;
   bool newpsi = true;
      /* read psi from file if psi_exist and not forcing wavefunction initialization */
   if (psi_filefind(mypneb,filename) && (!wvfnc_initialize))
   {
      newpsi = psi_check_convert(mypneb,filename,coutput);

      if (myparall->base_stdio_print)
         coutput << " input psi exists, reading from file: " << filename << std::endl;

      psi_read0(mypneb, &version, nfft, unita, &ispin, ne, psi2, filename);
   }

   /* generate new psi */
   else
   {
      if (myparall->base_stdio_print) coutput << " generating random psi from scratch" << std::endl;
      mypneb->g_generate_random(psi2);
   }
   newpsi = newpsi || (ispin != mypneb->ispin)
                   || (ne[0] != mypneb->ne[0])
                   || (ne[1] != mypneb->ne[1])
                   || (std::abs(unita[0] - mypneb->lattice->unita1d(0)) > 1.0e-4)
                   || (std::abs(unita[1] - mypneb->lattice->unita1d(1)) > 1.0e-4)
                   || (std::abs(unita[2] - mypneb->lattice->unita1d(2)) > 1.0e-4)
                   || (std::abs(unita[3] - mypneb->lattice->unita1d(3)) > 1.0e-4)
                   || (std::abs(unita[4] - mypneb->lattice->unita1d(4)) > 1.0e-4)
                   || (std::abs(unita[5] - mypneb->lattice->unita1d(5)) > 1.0e-4)
                   || (std::abs(unita[6] - mypneb->lattice->unita1d(6)) > 1.0e-4)
                   || (std::abs(unita[7] - mypneb->lattice->unita1d(7)) > 1.0e-4)
                   || (std::abs(unita[8] - mypneb->lattice->unita1d(8)) > 1.0e-4);

   /* ortho check */
   double sum2 = mypneb->gg_traceall(psi2, psi2);
   double sum1 = mypneb->ne[0] + mypneb->ne[1];

   if ((mypneb->ispin) == 1)
      sum1 *= 2;
   if (std::fabs(sum2 - sum1) > 1.0e-10)
   {
      if (myparall->base_stdio_print)
         coutput << " Warning - Gram-Schmidt being performed on psi2" << std::endl;
      mypneb->g_ortho(-1,psi2);

      double sum3 = mypneb->gg_traceall(psi2, psi2);
      if (myparall->base_stdio_print)
        coutput << "         - exact norm = " << sum1 << " norm=" << sum2
                << " corrected norm=" << sum3
                << " (error=" << std::abs(sum2 - sum1) << ")" << std::endl;
   }

   return newpsi;
}

/*****************************************************
 *                                                   *
 *                psi_write                          *
 *                                                   *
 *****************************************************/
void psi_write(Pneb *mypneb, int *version, int nfft[], double unita[],
               int *ispin, int ne[], double *psi, 
               int *occupation, double occ[], 
               char *filename,
               std::ostream &coutput) 
{
   nwpw_timing_function ftimer(50);
 
   Parallel *myparall = mypneb->d3db::parall;
 
   if (myparall->base_stdio_print)
     coutput << " output psi to filename: " << filename;
 
   if (myparall->is_master()) {
     openfile(6, filename, "w");
     iwrite(6, version, 1);
     iwrite(6, nfft, 3);
     dwrite(6, unita, 9);
     iwrite(6, ispin, 1);
     iwrite(6, ne, 2);
     iwrite(6, occupation, 1);
   }
 
   // write out 3d blocs
   mypneb->g_write(6, psi);

   // write out occupations
   if (*occupation>0)
   {
      if (myparall->is_master())
      {
         coutput << " ... writing occupations";
         dwrite(6, occ, ne[0]+ne[1]);
      }
   }
 
   if (myparall->is_master())
     closefile(6);

   if (myparall->base_stdio_print) coutput << std::endl;
}

void psi_write(Pneb *mypneb, int *version, int nfft[], double unita[],
               int *ispin, int ne[], double *psi,
               char *filename,
               std::ostream &coutput)
{
   int occupation =0;
   nwpw_timing_function ftimer(50);

   Parallel *myparall = mypneb->d3db::parall;

   if (myparall->base_stdio_print)
     coutput << " output psi to filename: " << filename << std::endl;

   if (myparall->is_master()) {
     openfile(6, filename, "w");
     iwrite(6, version, 1);
     iwrite(6, nfft, 3);
     dwrite(6, unita, 9);
     iwrite(6, ispin, 1);
     iwrite(6, ne, 2);
     iwrite(6, &occupation, 1);
   }

   // write out 3d blocs
   mypneb->g_write(6, psi);


   if (myparall->is_master())
     closefile(6);
}

/*****************************************************
 *                                                   *
 *                epsi_read                          *
 *                                                   *
 *****************************************************/
/*
   Reads psi and checks the header and check for orthonormalization.
   Entry - mypneb: Pneb grid structure
           filename: input filename
           wvfnc_initialize: force initialize wavefuntion
   Exit - psi2: complex wavefunction

   Uses - psi_filefind, psi_read0,
          mypneb->g_generate_random,
          mypneb->gg_traceall,
          mypneb->g_ortho
*/
bool epsi_read(Pneb *mypneb, char *filename, bool wvfnc_initialize, const int nex[], double *psi2, std::ostream &coutput) 
{
   nwpw_timing_function ftimer(50);
   int version, ispin, nfft[3], nexin[2],occupation;
   double unita[9];
   Parallel *myparall = mypneb->d3db::parall;
   bool newpsi = true;
   bool generate_newpsi = true;

   /* read psi from file if psi_exist and not forcing wavefunction initialization */
   if (psi_filefind(mypneb,filename) && (!wvfnc_initialize))
   {
      newpsi = psi_check_convert(mypneb,filename,coutput);

      if (myparall->base_stdio_print)
         coutput << " reading from file " << std::endl << std::endl;

      if (myparall->is_master())
      {
         openfile(4, filename, "r");
         iread(4, &version, 1);
         iread(4, nfft, 3);
         dread(4, unita, 9);
         iread(4, &ispin, 1);
         iread(4, nexin, 2);
         iread(4, &occupation, 1);
      }
      myparall->Brdcst_iValue(0, 0, &version);
      myparall->Brdcst_iValues(0, 0, 3, nfft);
      myparall->Brdcst_Values(0, 0, 9, unita);
      myparall->Brdcst_iValue(0, 0, &ispin);
      myparall->Brdcst_iValues(0, 0, 2, nexin);
      myparall->Brdcst_iValue(0, 0, &occupation);

      mypneb->g_read_excited(4, nexin,psi2);
      if (myparall->is_master())
         closefile(4);

      //psi_read0(mypneb, &version, nfft, unita, &ispin, ne, psi2, filename);
      if ((nex[0] == nexin[0]) && (nex[1]==nexin[1]))
         generate_newpsi = false;
   }

   /* generate new psi */
   if (generate_newpsi)
   {
      if (myparall->base_stdio_print) coutput << " generating random psi from scratch" << std::endl << std::endl;
      mypneb->g_generate_excited_random(nex, psi2);
   }


   newpsi = newpsi || (ispin != mypneb->ispin)
                   || (nexin[0] != nex[0])
                   || (nexin[1] != nex[1])
                   || (std::abs(unita[0] - mypneb->lattice->unita1d(0)) > 1.0e-4)
                   || (std::abs(unita[1] - mypneb->lattice->unita1d(1)) > 1.0e-4)
                   || (std::abs(unita[2] - mypneb->lattice->unita1d(2)) > 1.0e-4)
                   || (std::abs(unita[3] - mypneb->lattice->unita1d(3)) > 1.0e-4)
                   || (std::abs(unita[4] - mypneb->lattice->unita1d(4)) > 1.0e-4)
                   || (std::abs(unita[5] - mypneb->lattice->unita1d(5)) > 1.0e-4)
                   || (std::abs(unita[6] - mypneb->lattice->unita1d(6)) > 1.0e-4)
                   || (std::abs(unita[7] - mypneb->lattice->unita1d(7)) > 1.0e-4)
                   || (std::abs(unita[8] - mypneb->lattice->unita1d(8)) > 1.0e-4);

   /* ortho check */
   double sum2 = mypneb->gg_traceall_excited(nex, psi2, psi2);
   double sum1 = nex[0] + nex[1];
   newpsi = newpsi && (std::abs(sum2-sum1)> 1.0e-4);

   return newpsi;
}



/*****************************************************
 *                                                   *
 *                epsi_write                         *
 *                                                   *
 *****************************************************/
void epsi_write(Pneb *mypneb, int *version, int nfft[], double unita[],
               int *ispin, int ne[], double *psi, char *filename,
               std::ostream &coutput)
{
   nwpw_timing_function ftimer(50);
   int occupation = -1;

   Parallel *myparall = mypneb->d3db::parall;

   if (myparall->base_stdio_print)
     coutput << " output epsi to filename: " << filename << std::endl;

   if (myparall->is_master()) {
     openfile(6, filename, "w");
     iwrite(6, version, 1);
     iwrite(6, nfft, 3);
     dwrite(6, unita, 9);
     iwrite(6, ispin, 1);
     iwrite(6, ne, 2);
     iwrite(6, &occupation, 1);
   }

   mypneb->g_write_excited(6, ne, psi);

   if (myparall->is_master())
     closefile(6);
}


/*****************************************************
 *                                                   *
 *                psi_filefind                       *
 *                                                   *
 *****************************************************/
bool psi_filefind(Pneb *mypneb, char *filename) 
{
   int ifound = 0;
   Parallel *myparall = mypneb->d3db::parall;
 
   if (myparall->is_master()) {
     ifound = cfileexists(filename);
   }
 
   myparall->Brdcst_iValue(0, 0, &ifound);
 
   return (ifound > 0);
}

/*
void v_psi_read(Pneb *mypneb,int *version, int nfft[],
              double unita[], int *ispin, int ne[],
              double *psi)
{
   int occupation;

   Parallel *myparall = mypneb->d3db::parall;

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

   mypneb->g_read(4,psi);

   if (myparall->is_master()) closefile(4);
}
*/

/*
void v_psi_write(Pneb *mypneb,int *version, int nfft[],
              double unita[], int *ispin, int ne[],
              double *psi)
{
   int occupation = -1;

   Parallel *myparall = mypneb->d3db::parall;

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

   mypneb->g_write(6,psi);

   if (myparall->is_master()) closefile(6);
}

*/

// Atomic guess implementation
void Pneb::g_generate_atomic_guess(double *psi) {
    // Generate atomic-like orbitals centered on each atom
    // This creates a more physically motivated initial guess
    
    int taskid = d1db::parall->taskid();
    util_random(taskid + 91);
    
    double tmp2[n2ft3d];
    int taskid_j = d1db::parall->taskid_j();
    
    // Get ion positions and lattice info (we need to access these from somewhere)
    // For now, we'll use a simplified approach with atomic-like functions
    
    for (auto ms = 0; ms < ispin; ++ms)
        for (auto n = 0; n < ne[ms]; ++n) 
        {
            int qj = msntoindex(ms, n);
            int pj = msntop(ms, n);
            if (pj == taskid_j) 
            {
                // Generate atomic-like orbital with exponential decay
                d3db::r_zero(tmp2);
                
                // Create atomic-like function: exp(-alpha*r) with some randomness
                double alpha = 2.0 + util_random(0) * 3.0; // Random alpha between 2-5
                double center_x = 0.5 * nx * (0.5 - util_random(0)); // Random center
                double center_y = 0.5 * ny * (0.5 - util_random(0));
                double center_z = 0.5 * nz * (0.5 - util_random(0));
                
                for (auto i = 0; i < n2ft3d; ++i) {
                    int ix = i % (nx + 2);
                    int iy = (i / (nx + 2)) % ny;
                    int iz = i / ((nx + 2) * ny);
                    
                    if (ix < nx && iy < ny && iz < nz) {
                        double dx = (ix - center_x) * lattice->unita(0,0) / nx;
                        double dy = (iy - center_y) * lattice->unita(1,1) / ny;
                        double dz = (iz - center_z) * lattice->unita(2,2) / nz;
                        double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                        
                        if (r > 1e-6) {
                            tmp2[i] = std::exp(-alpha * r) * (0.5 - util_random(0)) * 0.1;
                        } else {
                            tmp2[i] = (0.5 - util_random(0)) * 0.1;
                        }
                    }
                }
                
                d3db::r_zero_ends(tmp2);
                d3db::rc_fft3d(tmp2);
                
                PGrid::c_pack(1, tmp2);
                int indx = 2 * PGrid::npack(1) * qj;
                PGrid::cc_pack_copy(1, tmp2, psi + indx);
            }
        }
}

// Superposition of atomic orbitals
void Pneb::g_generate_superposition_guess(double *psi) {
    // Generate superposition of atomic orbitals
    // This creates a more physically motivated initial guess
    // by combining atomic orbitals with appropriate phases
    
    int taskid = d1db::parall->taskid();
    util_random(taskid + 92);
    
    double tmp2[n2ft3d];
    int taskid_j = d1db::parall->taskid_j();
    
    for (auto ms = 0; ms < ispin; ++ms)
        for (auto n = 0; n < ne[ms]; ++n) 
        {
            int qj = msntoindex(ms, n);
            int pj = msntop(ms, n);
            if (pj == taskid_j) 
            {
                d3db::r_zero(tmp2);
                
                // Create superposition of multiple atomic-like functions
                int n_centers = 2 + (n % 3); // 2-4 centers per orbital
                
                for (auto center = 0; center < n_centers; ++center) {
                    double alpha = 1.5 + util_random(0) * 2.0;
                    double center_x = nx * (0.3 + 0.4 * util_random(0));
                    double center_y = ny * (0.3 + 0.4 * util_random(0));
                    double center_z = nz * (0.3 + 0.4 * util_random(0));
                    double phase = util_random(0) * 2.0 * M_PI;
                    double amplitude = (0.5 - util_random(0)) * 0.2;
                    
                    for (auto i = 0; i < n2ft3d; ++i) {
                        int ix = i % (nx + 2);
                        int iy = (i / (nx + 2)) % ny;
                        int iz = i / ((nx + 2) * ny);
                        
                        if (ix < nx && iy < ny && iz < nz) {
                            double dx = (ix - center_x) * lattice->unita(0,0) / nx;
                            double dy = (iy - center_y) * lattice->unita(1,1) / ny;
                            double dz = (iz - center_z) * lattice->unita(2,2) / nz;
                            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                            
                            if (r > 1e-6) {
                                tmp2[i] += amplitude * std::exp(-alpha * r) * std::cos(phase);
                            } else {
                                tmp2[i] += amplitude * std::cos(phase);
                            }
                        }
                    }
                }
                
                d3db::r_zero_ends(tmp2);
                d3db::rc_fft3d(tmp2);
                
                PGrid::c_pack(1, tmp2);
                int indx = 2 * PGrid::npack(1) * qj;
                PGrid::cc_pack_copy(1, tmp2, psi + indx);
            }
        }
}

// Gaussian basis functions
void Pneb::g_generate_gaussian_guess(double *psi) {
    // Generate Gaussian basis functions as initial guess
    // This can be more stable than random for some systems
    
    int taskid = d1db::parall->taskid();
    util_random(taskid + 93);
    
    double tmp2[n2ft3d];
    int taskid_j = d1db::parall->taskid_j();
    
    for (auto ms = 0; ms < ispin; ++ms)
        for (auto n = 0; n < ne[ms]; ++n) 
        {
            int qj = msntoindex(ms, n);
            int pj = msntop(ms, n);
            if (pj == taskid_j) 
            {
                d3db::r_zero(tmp2);
                
                // Create Gaussian functions with different widths
                double sigma = 0.5 + util_random(0) * 2.0; // Random width
                double center_x = nx * (0.2 + 0.6 * util_random(0));
                double center_y = ny * (0.2 + 0.6 * util_random(0));
                double center_z = nz * (0.2 + 0.6 * util_random(0));
                double amplitude = (0.5 - util_random(0)) * 0.3;
                
                for (auto i = 0; i < n2ft3d; ++i) {
                    int ix = i % (nx + 2);
                    int iy = (i / (nx + 2)) % ny;
                    int iz = i / ((nx + 2) * ny);
                    
                    if (ix < nx && iy < ny && iz < nz) {
                        double dx = (ix - center_x) * lattice->unita(0,0) / nx;
                        double dy = (iy - center_y) * lattice->unita(1,1) / ny;
                        double dz = (iz - center_z) * lattice->unita(2,2) / nz;
                        double r2 = dx*dx + dy*dy + dz*dz;
                        
                        tmp2[i] = amplitude * std::exp(-r2 / (2.0 * sigma * sigma));
                    }
                }
                
                d3db::r_zero_ends(tmp2);
                d3db::rc_fft3d(tmp2);
                
                PGrid::c_pack(1, tmp2);
                int indx = 2 * PGrid::npack(1) * qj;
                PGrid::cc_pack_copy(1, tmp2, psi + indx);
            }
        }
}

// Mixed strategy: atomic-like for occupied, random for virtual
void Pneb::g_generate_mixed_guess(double *psi) {
    // Mixed strategy: atomic-like for occupied, random for virtual
    // This combines the best of both approaches
    
    int taskid = d1db::parall->taskid();
    util_random(taskid + 94);
    
    double tmp2[n2ft3d];
    int taskid_j = d1db::parall->taskid_j();
    
    for (auto ms = 0; ms < ispin; ++ms)
        for (auto n = 0; n < ne[ms]; ++n) 
        {
            int qj = msntoindex(ms, n);
            int pj = msntop(ms, n);
            if (pj == taskid_j) 
            {
                // Use atomic-like for first half of orbitals, random for second half
                if (n < ne[ms] / 2) {
                    // Atomic-like approach
                    d3db::r_zero(tmp2);
                    
                    double alpha = 2.0 + util_random(0) * 2.0;
                    double center_x = nx * (0.3 + 0.4 * util_random(0));
                    double center_y = ny * (0.3 + 0.4 * util_random(0));
                    double center_z = nz * (0.3 + 0.4 * util_random(0));
                    
                    for (auto i = 0; i < n2ft3d; ++i) {
                        int ix = i % (nx + 2);
                        int iy = (i / (nx + 2)) % ny;
                        int iz = i / ((nx + 2) * ny);
                        
                        if (ix < nx && iy < ny && iz < nz) {
                            double dx = (ix - center_x) * lattice->unita(0,0) / nx;
                            double dy = (iy - center_y) * lattice->unita(1,1) / ny;
                            double dz = (iz - center_z) * lattice->unita(2,2) / nz;
                            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                            
                            if (r > 1e-6) {
                                tmp2[i] = std::exp(-alpha * r) * (0.5 - util_random(0)) * 0.1;
                            } else {
                                tmp2[i] = (0.5 - util_random(0)) * 0.1;
                            }
                        }
                    }
                } else {
                    // Random approach
                    d3db::r_setrandom(tmp2);
                    d3db::r_zero_ends(tmp2);
                }
                
                d3db::rc_fft3d(tmp2);
                
                PGrid::c_pack(1, tmp2);
                int indx = 2 * PGrid::npack(1) * qj;
                PGrid::cc_pack_copy(1, tmp2, psi + indx);
            }
        }
}

} // namespace pwdft
