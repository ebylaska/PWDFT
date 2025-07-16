

#include <cmath>
#include <cstring> //memset
#include <iostream>
#include <string>
//#include	"control.hpp"
// extern "C" {
//#include        "compressed_io.h"
//}
#include "compressed_io.hpp"

#include "Parallel.hpp"
#include "Cneb.hpp"
#include "cpsi.hpp"
#include "../../nwpwlib/utilities/util.hpp"

namespace pwdft {

/*****************************************************
 *                                                   *
 *              cwvfnc_expander_convert              *
 *                                                   *
 *****************************************************/

static void cwvfnc_expander_convert(int ngrid[], double *psi1, int dngrid[], double *psi2) 
{
   int nfft3d  =  ngrid[0] *  ngrid[1] *  ngrid[2];
   int dnfft3d = dngrid[0] * dngrid[1] * dngrid[2];
   int n2ft3d  = 2*nfft3d;
   int dn2ft3d = 2*dnfft3d;
   int inc2  = ngrid[0];
   int dinc2 = dngrid[0];
   int inc3  = ngrid[0]*ngrid[1];
   int dinc3 = dngrid[0]*dngrid[1];
 
   int n1 = ngrid[0];
   int n2 = ngrid[1];
   int n3 = ngrid[2];
   if (n1 > dngrid[0]) n1 = dngrid[0];
   if (n2 > dngrid[1]) n2 = dngrid[1];
   if (n3 > dngrid[2]) n3 = dngrid[2];
 
   int idiff = dngrid[0] - ngrid[0];
   int jdiff = dngrid[1] - ngrid[1];
   int kdiff = dngrid[2] - ngrid[2];
   bool ireverse = (idiff < 0);
   bool jreverse = (jdiff < 0);
   bool kreverse = (kdiff < 0);
   if (ireverse) idiff = -idiff;
   if (jreverse) jdiff = -jdiff;
   if (kreverse) kdiff = -kdiff;
 
   std::memset(psi2, 0, dn2ft3d*sizeof(double));

   for (auto k=0; k<n3; ++k)
   for (auto j=0; j<n2; ++j)
   for (auto i=0; i<n1; ++i) 
   {
      int i2 = (i < n1 / 2) ? i : i + idiff;
      int j2 = (j < n2 / 2) ? j : j + jdiff;
      int k2 = (k < n3 / 2) ? k : k + kdiff;

      int indx = 0, dindx = 0;

      indx  += ireverse ? i2 : i;
      dindx += ireverse ? i  : i2;

      indx  += (jreverse ? j2 : j)  * inc2;
      dindx += (jreverse ? j  : j2) * dinc2;

      indx  += (kreverse ? k2 : k)  * inc3;
      dindx += (kreverse ? k  : k2) * dinc3;

      // SAFETY: Bounds check before accessing memory
      if (indx < 0 || indx >= nfft3d || dindx < 0 || dindx >= dnfft3d)
         continue;

      psi2[2 * dindx]     = psi1[2 * indx];
      psi2[2 * dindx + 1] = psi1[2 * indx + 1];
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
 *                cwvfnc_expander                    *
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
 
     int n2ft3d  = 2* nfft[0] *  nfft[1] *  nfft[2];
     int dn2ft3d = 2*dnfft[0] * dnfft[1] * dnfft[2];
     double *psi1 = new double[n2ft3d];
     double *psi2 = new double[dn2ft3d];
     for (auto nb=0; nb<nbrillouin; ++nb)
     for (auto ms=0; ms<ispin; ++ms)
     for (auto n=0; n<ne[ms]; ++n) 
     {
        if (lprint) coutput << " converting .... cpsi: nb=" << nb+1 << " n=" << n+1 << " spin=" << ms + 1 << std::endl;
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
 *                cpsi_check_convert                 *
 *                                                   *
 *****************************************************/
/**
 * @brief Checks and converts the psi grid data if necessary based on input file parameters.
 *
 * This static function validates the compatibility of psi grid data dimensions from an input file
 * (`filename`) with the current Cneb object (`_mycneb`). It compares grid dimensions and, if 
 * mismatched, initiates a conversion process using `cwvfnc_expander`. The function generates debug 
 * output and supports both standard output and custom stream output through `coutput`.
 *
 * Operational Steps:
 * 1. Retrieves parallel and file header information, extracting version number, grid dimensions,
 *    unit cell vectors, spin details, electron counts, and Brillouin zone data from the file.
 * 2. Compares file-retrieved grid dimensions (`nfft0`) with those in `_mycneb` (`nx`, `ny`, `nz`).
 * 3. If dimensions do not match, engages `cwvfnc_expander` to convert the psi data, logging
 *    conversion activity through the provided `coutput`.
 *
 * @param  mycneb Cneb object containing current simulation parameters.
 * @param filename Character array representing the path to the input file with psi data.
 * @param coutput Output stream for logging conversion messages and activities.
 * @return Returns `true` if conversion is necessary and has been executed, otherwise `false`.
 *
 * Note:
 * - Ensure proper MPI and parallel setup in `mycneb` to accurately invoke file header operations
 *   with `cpsi_get_header`.
 * - This function relies on version-specific checks for compatibility, which should be validated
 *   according to the latest implementation or file structures.
 */
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

   // also convert if ne and nbrillouin are wrong
 
   return converted;
}


/*****************************************************
 *                                                   *
 *                cpsi_read0                         *
 *                                                   *
 *****************************************************/
/**
 * @brief Reads psi data and its header information from a file.
 *
 * This function is responsible for reading psi grid data and associated header metadata from a 
 * specified file (`filename`). It efficiently manages file input operations, distributing the 
 * read data across MPI tasks within the `Cneb` context. The function facilitates both standard 
 * and reverse reading modes, leveraging grid-specific data handling operations and optional 
 * occupation data adjustments.
 *
 * Operational Steps:
 * 1. Master node opens the specified file for reading and retrieves header details, including:
 *    - Version, grid dimensions (`nfft`), unit cell vectors (`unita`), spin configuration, electron counts (`ne`),
 *      and Brillouin zones, as well as initial occupation data.
 * 2. Utilizes MPI broadcast (`Brdcst_`) to distribute header information across all tasks.
 * 3. Conditionally reads psi data:
 *    - Uses reverse packing logic (`g_read_reverse`) if `reverse` is true.
 *    - Otherwise applies standard reading logic (`g_read`).
 * 4. If occupation data is present (`occupation > 0`), retrieves occupation numbers using `g_read_occ`,
 *    optionally reversing the data ordering based on the `reverse` flag settings and spin configurations.
 * 5. Closes the file on the master node after completing all read operations.
 *
 * @param mycneb Cneb object integrating simulation parameters and framework.
 * @param version Integer representing version information as file metadata.
 * @param nfft Array denoting grid dimensions (x, y, z) for the psi data.
 * @param unita Array representing unit cell vectors for contextualizing data layout.
 * @param ispin Integer specifying spin details for the psi dataset.
 * @param ne Array containing electron counts per spin configuration.
 * @param nbrillouin Integer specifying the number of Brillouin zones.
 * @param psi Pointer to the psi data array for input.
 * @param occupation Integer holding occupation state details.
 * @param occ Array for potential occupation number reading and adjustments.
 * @param filename Character array with the name of the file to read from.
 * @param reverse Boolean flag indicating whether to apply reverse logic during reads.
 *
 * Note:
 * - Ensure file existence and accessibility prior to invocation, as the function assumes a pre-existing file.
 * - Verify MPI initialization and task role consistency (via `myparall`) for accurate distribution and data handling.
 */
void cpsi_read0(Cneb *mycneb, int *version, int nfft[], double unita[],
                int *ispin, int ne[], int *nbrillouin, double *psi, 
                int *occupation, double occ[], char *filename, bool reverse)
{
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
      iread(4, occupation, 1);
   }
   myparall->Brdcst_iValue(0, 0, version);
   myparall->Brdcst_iValues(0, 0, 3, nfft);
   myparall->Brdcst_Values(0, 0, 9, unita);
   myparall->Brdcst_iValue(0, 0, ispin);
   myparall->Brdcst_iValues(0, 0, 2, ne);
   myparall->Brdcst_iValue(0, 0, nbrillouin);
   myparall->Brdcst_iValue(0, 0, occupation);

   /* reads in c format and automatically packs the result to g format */
   // need to make sure psi is the proper size before reading
   if (reverse)
      mycneb->g_read_ne_reverse(4,ne,*nbrillouin,psi);
   else
      mycneb->g_read_ne(4,ne,*nbrillouin,psi);

   if ((*occupation) > 0)
   {
      mycneb->g_read_occ(4,occ);
      if (reverse)
      {
         std::reverse(occ, occ + ne[0]); // Reverse the first section if descending

         // Check and reverse the second part
         if (*ispin>1)
            std::reverse(occ + ne[0], occ + ne[0] + ne[1]); // Reverse the second section if descending and *ispin > 1
      }
   }

   if (myparall->is_master())
     closefile(4);
}

/*****************************************************
 *                                                   *
 *                cpsi_read0   (overloaded)          *
 *                                                   *
 *****************************************************/
/**
 * @brief Overloaded function that reads psi data and its header information from a file.
 *
 * This overloaded variant of `cpsi_read0` focuses on reading psi grid data and header metadata 
 * from a specified file (`filename`). It effectively manages file input operations and distributes 
 * the read data across MPI tasks within the `Cneb` context. Simplified compared to its counterpart, 
 * this version does not handle reverse logic nor occupation data adjustments via an external `occ` array.
 *
 * Operational Steps:
 * 1. Master node opens the specified file for reading and retrieves header details, including:
 *    - Version, grid dimensions (`nfft`), unit cell vectors (`unita`), spin configuration, electron counts (`ne`),
 *      and Brillouin zones, along with occupation state data.
 * 2. Utilizes MPI broadcast (`Brdcst_`) to distribute header information across all tasks.
 * 3. Reads psi data using the g-format packing logic of `g_read`.
 * 4. Closes the file on the master node after completing all read operations.
 *
 * @param mycneb Cneb object integrating simulation parameters and framework.
 * @param version Integer representing version information as file metadata.
 * @param nfft Array denoting grid dimensions (x, y, z) for the psi data.
 * @param unita Array representing unit cell vectors for contextualizing data layout.
 * @param ispin Integer specifying spin details for the psi dataset.
 * @param ne Array containing electron counts per spin configuration.
 * @param nbrillouin Integer specifying the number of Brillouin zones.
 * @param psi Pointer to the psi data array for input.
 * @param filename Character array with the name of the file to read from.
 *
 * Note:
 * - Ensure file existence and accessibility prior to function invocation, as it assumes a pre-existing file.
 * - Verify MPI setup and consistency across tasks (via `myparall`) to accurately distribute header and data information.
 * - This function does not apply occupation number reversal nor additional occupation array adjustments.
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
   //mycneb->g_read(4,psi);
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
bool cpsi_read(Cneb *mycneb, char *filename, bool wvfnc_initialize, double *psi2,
              int *occupation, double occ2[],
              std::ostream &coutput)
{  
   nwpw_timing_function ftimer(50);
   int version, ispin, nfft[3], ne[2], nbrillouin;
   double unita[9];
   Parallel *myparall = mycneb->c3db::parall;
   bool newpsi = true;
                   
   /* read psi from file if psi_exist and not forcing wavefunction initialization */
   if (cpsi_filefind(mycneb,filename) && (!wvfnc_initialize)) 
   {               
      newpsi = cpsi_check_convert(mycneb,filename,coutput); // also convert if ne and nbrillouin are wrong
                   
      if (myparall->base_stdio_print)
         coutput << " input psi exists, reading from file: " << filename;
   
      ne[0] = mycneb->ne[0]; 
      ne[1] = mycneb->ne[1];
      cpsi_read0(mycneb, &version, nfft, unita, &ispin, ne, &nbrillouin, psi2, occupation, occ2, filename,false);
      
      if ((*occupation>0) and (myparall->base_stdio_print))
         coutput << " ... reading occupations";

      if (myparall->base_stdio_print)
         coutput << std::endl;


      if ((occ2) && (*occupation > 0))
      {
         if (isDescending(occ2, ne[0]))
         {
            if (myparall->base_stdio_print)
               coutput << " - reversing order of psi and occupation" << std::endl;
            cpsi_read0(mycneb, &version, nfft, unita, &ispin, ne, &nbrillouin, psi2, occupation, occ2, filename,true);
         }
      }
   }

   /* generate new psi */
   else
   {
      ispin = mycneb->ispin;
      ne[0] = mycneb->ne[0];
      ne[1] = mycneb->ne[1];
      std::string guess = get_initial_wavefunction_guess();
      if (guess == "atomic") {
         if (myparall->base_stdio_print) coutput << " generating atomic guess for cpsi" << std::endl;
         mycneb->g_generate_atomic_guess(psi2); // To be implemented
      } else {
         if (myparall->base_stdio_print) coutput << " generating random cpsi from scratch" << std::endl;
         mycneb->g_generate_random(psi2);
      }
   }

   newpsi = newpsi || (ispin != mycneb->ispin)
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

   // Add extra missing orbitals to psi2
   //   Current orbital structure is:
   //      psi[0] ... psi[ne-1]      ← old orbitals
   //      psi[ne] ... psi[ne+nextra-1] ← not set
   //   Orbitals should be restructured to be:
   //      psi[0] ... psi[nextra-1]  ← new orbitals set to be (high energy, low occ)
   //      psi[nextra] ... psi[ne-1] ← old orbitals
   int nextra[2] = {(mycneb->ne[0]-ne[0]),(mycneb->ne[1]-ne[1])};
   int orb_size = 2 * mycneb->npack1_max();  // size of a single orbital
   for (auto ms=0; ms<ispin; ++ms)
   {
      if (nextra[ms]>0)
      {
         if (myparall->base_stdio_print)
            coutput << " - adding " << nextra[ms] << " to ms=" << ms << " psi" << std::endl;

         int nold = ne[ms];
         int ntotal = mycneb->ne[ms];

         // Shift existing psi upward to make room for new orbitals at the bottom
         for (int n = nold - 1; n >= 0; --n)
         {
             int src = orb_size * (n + ms * mycneb->ne[0]);
             int dst = orb_size * (n + nextra[ms] + ms * mycneb->ne[0]);
             std::memmove(psi2 + dst, psi2 + src, sizeof(double) * orb_size);
         }

         // Generate new random orbitals at the start of this spin channel
         int start = orb_size * (ms * mycneb->ne[0]);
         mycneb->g_generate_extra_random(nextra[ms], psi2 + start);

         // need to orthopsi and  belowo
         // project out filled psi
         mycneb->g_project_out_filled_extra(ms,nextra,psi2);

         // remove nwpw":{"virtual":[4,0]}
      }
   }
      
   /* ortho check */
   double sum2 = mycneb->gg_traceall(psi2,psi2);
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
         coutput << " input cpsi exists, reading from file: " << filename << std::endl;
 
      cpsi_read0(mycneb, &version, nfft, unita, &ispin, ne, &nbrillouin, psi2, filename);
   }
 
   /* generate new psi */
   else 
   {
      if (myparall->base_stdio_print) coutput << " generating random cpsi from scratch" << std::endl;
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
/**
 * @brief Writes psi data and associated computational parameters to a file.
 *
 * This function handles the writing of psi grid data, occupation information, and other relevant 
 * simulation parameters to an output file, specified by `filename`. It ensures data is recorded 
 * correctly using collaborative MPI operations within the `Cneb` context. The master node initiates 
 * the file-writing sequence, followed by distributed grid and occupation data writes from `mycneb`.
 *
 * Operational Steps:
 * 1. Initializes function timing with `nwpw_timing_function` for performance tracking.
 * 2. Logs the output operation, detailing the filename, to the given output stream `coutput`.
 * 3. If current node is the master:
 *    - Opens the specified file (`filename`).
 *    - Writes header information including version, grid dimensions (`nfft`), unit cell vectors 
 *      (`unita`), spin data, electron counts (`ne`), Brillouin zone data, and occupation data.
 * 4. Delegates the writing of psi data to `mycneb->g_write`.
 * 5. If occupation data is present, executes `mycneb->g_write_occ` to record occupation numbers.
 * 6. Closes the file on the master node once all data writing operations are complete.
 *
 * @param mycneb  Cneb object containing simulation parameters and grid data.
 * @param version Integer representing the version of the data or simulation specifications.
 * @param nfft Array indicating grid dimensions (x, y, z) for the psi data.
 * @param unita Array representing unit cell vectors necessary for positional context.
 * @param ispin Integer specifying spin data configuration.
 * @param ne Array holding electron counts per spin configuration.
 * @param nbrillouin Integer specifying the number of Brillouin zones.
 * @param psi Pointer to psi data array to be written.
 * @param occupation Pointer to integer representing occupation state or configuration.
 * @param occ Array containing occupation numbers if applicable.
 * @param filename Character array holding the name of the file to write to.
 * @param coutput Output stream for logging actions and file operations.
 *
 * Note:
 * - Ensure consistent MPI setup and master node recognition to manage file operations and streamlining data writes.
 * - Validate that `mycneb` contains all necessary methods (`g_write`, `g_write_occ`) for efficient distributed data handling.
 */
void cpsi_write(Cneb *mycneb, int *version, int nfft[], double unita[],
               int *ispin, int ne[], int *nbrillouin, double *psi, 
               int *occupation, double occ[],
               char *filename, std::ostream &coutput)
{
   nwpw_timing_function ftimer(50);

   Parallel *myparall = mycneb->c3db::parall;

   if (myparall->base_stdio_print)
     coutput << " output cpsi to filename: " << filename << std::endl;

   if (myparall->is_master())
   {
      openfile(6, filename, "w");
      iwrite(6, version, 1);
      iwrite(6, nfft, 3);
      dwrite(6, unita, 9);
      iwrite(6, ispin, 1);
      iwrite(6, ne, 2);
      iwrite(6, nbrillouin, 1);
      iwrite(6, occupation, 1);
   }

   mycneb->g_write(6, psi);

   if (*occupation>0)
      mycneb->g_write_occ(6, occ);

   if (myparall->is_master())
      closefile(6);
}

/*****************************************************
 *                                                   *
 *                cpsi_write  (overloaded)           *
 *                                                   *
 *****************************************************/
/**
 * @brief Overloaded function that writes psi data and associated simulation parameters to a file 
 * without occupation data.
 *
 * As an overloaded version of `cpsi_write`, this function records psi grid data and associated metadata 
 * to an output file specified by `filename`. The master node handles file operations while distributing 
 * data writes via MPI processes within the `Cneb` context. Notably, although an occupation variable is 
 * defined, no actual occupation data is written to the file, and its default value is set to -1.
 *
 * Operational Steps:
 * 1. Initializes function timing with `nwpw_timing_function` for performance assessment.
 * 2. Logs the intended output operation, including filename details, if applicable, to the `coutput` stream.
 * 3. If the current node is the master node:
 *    - Opens the specified file (`filename`) for writing.
 *    - Writes header information including version, grid dimensions (`nfft`), unit cell vectors (`unita`),
 *      spin configuration, electron counts (`ne`), and Brillouin zones.
 *    - Writes occupation data with a default placeholder value (-1), without involving actual occupation data logic.
 * 4. Delegates the writing of psi data to `mycneb->g_write`.
 * 5. Evaluates if occupation numbers are greater than zero but does not perform the write.
 * 6. Closes the file on the master node once all data writing operations are complete.
 *
 * @param mycneb Cneb object containing simulation parameters and grid data.
 * @param version Integer representing the version of the data or simulation specifications.
 * @param nfft Array indicating grid dimensions (x, y, z) for the psi data.
 * @param unita Array representing unit cell vectors necessary for structural context.
 * @param ispin Integer specifying spin configuration for the data.
 * @param ne Array holding electron counts per spin configuration.
 * @param nbrillouin Integer denoting the number of Brillouin zones.
 * @param psi Pointer to psi data array to be written.
 * @param filename Character array holding the name of the file to write to.
 * @param coutput Output stream for logging actions and file operations.
 *
 * Note:
 * - Confirm consistent MPI setup and master node recognition to effectively manage file operations and data writes.
 * - The function, as an overloaded version, does not perform occupation data writes, utilizing a placeholder (-1) instead.
 * - Ensure `mycneb` provides necessary methods (`g_write`) for distributed data handling.
 */
void cpsi_write(Cneb *mycneb, int *version, int nfft[], double unita[],
               int *ispin, int ne[], int *nbrillouin, double *psi, char *filename,
               std::ostream &coutput) 
{
   nwpw_timing_function ftimer(50);
   int occupation = -1;
 
   Parallel *myparall = mycneb->c3db::parall;
 
   if (myparall->base_stdio_print)
 
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
       coutput << "Write the occupations here!" << std::endl;
 
   if (myparall->is_master())
      closefile(6);
}

/*****************************************************
 *                                                   *
 *                ecpsi_read                         *
 *                                                   *
 *****************************************************/
bool ecpsi_read(Cneb *mycneb, char *filename, bool wvfnc_initialize, const int nex[], double *psi2, std::ostream &coutput)
{
   nwpw_timing_function ftimer(50);
   int version, ispin, nfft[3], nexin[2],occupation,nbrillouin;
   double unita[9];
   Parallel *myparall = mycneb->c3db::parall;
   bool newpsi = true;
   bool generate_newpsi = true;

   /* read psi from file if psi_exist and not forcing wavefunction initialization */
   if (cpsi_filefind(mycneb,filename) && (!wvfnc_initialize))
   {
      newpsi = cpsi_check_convert(mycneb,filename,coutput);

      if (myparall->base_stdio_print)
         coutput << " reading from file " << std::endl << std::endl;

      //cpsi_read0(mycneb, &version, nfft, unita, &ispin, ne, &nbrillouin, psi2, filename);
      if (myparall->is_master())
      {
         openfile(4, filename, "r");
         iread(4, &version, 1);
         iread(4, nfft, 3);
         dread(4, unita, 9);
         iread(4, &ispin, 1);
         iread(4, nexin, 2);
         iread(4, &nbrillouin, 1);
         iread(4, &occupation, 1);
      }
      myparall->Brdcst_iValue(0, 0, &version);
      myparall->Brdcst_iValues(0, 0, 3, nfft);
      myparall->Brdcst_Values(0, 0, 9, unita);
      myparall->Brdcst_iValue(0, 0, &ispin);
      myparall->Brdcst_iValues(0, 0, 2, nexin);
      myparall->Brdcst_iValue(0, 0, &nbrillouin);
      myparall->Brdcst_iValue(0, 0, &occupation);

      mycneb->g_read_excited(4, nexin, nbrillouin,psi2);
      if (myparall->is_master())
         closefile(4);

      if ((nex[0] == nexin[0]) && (nex[1]==nexin[1]) & (nbrillouin == mycneb->nbrillouin))
         generate_newpsi = false;
   }

   /* generate new psi */
   if (generate_newpsi)
   {
      if (myparall->base_stdio_print) coutput << " generating random cpsi from scratch" << std::endl << std::endl;
      mycneb->g_generate_excited_random(nex, psi2);
   }

   newpsi = newpsi || (ispin != mycneb->ispin)
                   || (nbrillouin != mycneb->nbrillouin)
                   || (nexin[0] != nex[0])
                   || (nexin[1] != nex[1])
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
   double sum2 = mycneb->gg_traceall_excited(nex, psi2, psi2);
   double sum1 = nex[0] + nex[1];
   newpsi = newpsi && (std::abs(sum2-sum1)> 1.0e-4);

   return newpsi;
}

/*****************************************************
 *                                                   *
 *               ecpsi_write                         *
 *                                                   *
 *****************************************************/
void ecpsi_write(Cneb *mycneb, int *version, int nfft[], double unita[],
               int *ispin, int ne[], int *nbrillouin, double *psi, char *filename,
               std::ostream &coutput)
{
   nwpw_timing_function ftimer(50);
   int occupation = -1;

   Parallel *myparall = mycneb->c3db::parall;

   if (myparall->base_stdio_print)
     coutput << " output ecpsi to filename: " << filename << std::endl;

   if (myparall->is_master()) {
     openfile(6, filename, "w");
     iwrite(6, version, 1);
     iwrite(6, nfft, 3);
     dwrite(6, unita, 9);
     iwrite(6, ispin, 1);
     iwrite(6, ne, 2);
     iwrite(6, nbrillouin, 1);
     iwrite(6, &occupation, 1);
   }

   mycneb->g_write_excited(6, ne, *nbrillouin, psi);

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

// Stub for atomic guess (to be implemented)
void Cneb::g_generate_atomic_guess(double *psi) {
    // TODO: Implement atomic guess projection
    g_generate_random(psi); // fallback to random for now
}

} // namespace pwdft
