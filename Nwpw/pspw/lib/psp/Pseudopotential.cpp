/* Pseudopotential.cpp - 
   Author - Eric Bylaska
*/

/**
 * @class Pseudopotential
 * @brief Class for managing pseudopotential data and calculations.
 *
 * The Pseudopotential class encapsulates operations and data related to
 * pseudopotentials used in electronic structure calculations. It provides
 * methods for handling non-local and local pseudopotentials, semicore
 * corrections, and other properties used in electronic structure calculations.
 */

#include <cmath>
#include <cstring>
#include <iostream>

#include "Psp1d_Hamann.hpp"
#include "Psp1d_pawppv1.hpp"
//#include "gdevice.hpp"
#include "nwpw_timing.hpp"
#include "util.hpp"

#include "blas.h"

#include "Pseudopotential.hpp"
#include "compressed_io.hpp"

namespace pwdft {

/*******************************************
 *                                         *
 *              vpp_read_header            *
 *                                         *
 *******************************************/
/**
 * @brief Reads header information from a VPP (Variational Pseudopotential) file.
 *
 * This function reads header information from a VPP file, including comments, PSP type, version, FFT dimensions,
 * unit cell vectors, atom symbol, atomic mass, and valence charge. The information is stored in the provided
 * variables.
 *
 * @param fname The filename of the VPP file to read.
 * @param comment A character array to store the comment from the VPP file.
 * @param psp_type Pointer to an integer variable to store the PSP (Pseudopotential) type.
 * @param version Pointer to an integer variable to store the version of the VPP file.
 * @param nfft Pointer to an integer array to store FFT dimensions (3 elements).
 * @param unita Pointer to a double array to store unit cell vectors (9 elements).
 * @param atom A character array to store the atom symbol.
 * @param amass Pointer to a double variable to store the atomic mass.
 * @param zv Pointer to a double variable to store the valence charge.
 * @return `true` if the file was found and read successfully, `false` otherwise.
 */
static bool vpp_read_header(char *fname, char *comment, int *psp_type,
                            int *version, int nfft[], double unita[],
                            char *atom, double *amass, double *zv) {
  int i, ifound;

  ifound = cfileexists(fname);

  if (ifound > 0) {
    openfile(5, fname, "r");
    cread(5, comment, 80);
    comment[79] = '\0';
    i = 78;
    while (comment[i] == ' ')
      comment[i--] = '\0';

    iread(5, psp_type, 1);
    iread(5, version, 1);
    iread(5, nfft, 3);
    dread(5, unita, 9);
    cread(5, atom, 2);
    dread(5, amass, 1);
    dread(5, zv, 1);
    closefile(5);
  }

  return (ifound > 0);
}

/*****************************************************
 *                                                   *
 *                vpp_formatter_check                *
 *                                                   *
 *****************************************************/
/**
 * @brief Checks if a VPP (Variational Pseudopotential) file requires reformatting for compatibility.
 *
 * This function checks if a VPP file needs reformatting to be compatible with the current code. It reads the header
 * information from the VPP file and compares it to the parameters of the provided PGrid, including lattice vectors,
 * FFT dimensions, and PSP version. If any of these parameters differ or the file doesn't exist, it signals
 * that reformatting is needed.
 *
 * @param mygrid Pointer to the PGrid object associated with the grid.
 * @param fname The filename of the VPP file to check.
 * @param psp_version The expected PSP (Pseudopotential) version for compatibility.
 * @return `true` if reformatting is needed, `false` otherwise.
 */
static bool vpp_formatter_check(PGrid *mygrid, char *fname,
                                const int psp_version) {
  char comment[80], atom[2];
  int psp_type, version, nfft[3];
  double unita[9], amass, zv;
  double tol = 1.0e-9;
  int ireformat;

  ireformat = 1;
  if (mygrid->parall->is_master()) {
    if (vpp_read_header(fname, comment, &psp_type, &version, nfft, unita, atom,
                        &amass, &zv)) {
      bool reformat = false;
      for (auto i = 0; i < 9; ++i)
        reformat = reformat ||
                   (std::fabs(mygrid->lattice->unita1d(i) - unita[i]) > tol);
      reformat = reformat || (mygrid->nx != nfft[0]);
      reformat = reformat || (mygrid->ny != nfft[1]);
      reformat = reformat || (mygrid->nz != nfft[2]);
      reformat = reformat || (psp_version != version);
      if (!reformat)
        ireformat = 0;
    }
  }
  mygrid->parall->Brdcst_iValue(0, 0, &ireformat);

  return (ireformat == 1);
}


/*******************************************
 *                                         *
 *           Multiply_Gijl_sw1             *
 *                                         *
 *******************************************/
/**
 * @brief Multiply Gijl matrix with sw1 to obtain sw2.
 *
 * This function multiplies the Gijl matrix with sw1 to obtain sw2. It loops over the
 * projectors and calculates the product for matching l and m values. The result is stored
 * in sw2.
 *
 * @param nn Number of grid points.
 * @param nprj Number of projectors.
 * @param nmax Maximum value for n.
 * @param lmax Maximum value for l.
 * @param n_prj Array of n values for projectors.
 * @param l_prj Array of l values for projectors.
 * @param m_prj Array of m values for projectors.
 * @param G Matrix containing Gijl elements.
 * @param sw1 Array containing sw1 values.
 * @param sw2 Array to store the result (sw2 = Gijl * sw1).
 */
static void Multiply_Gijl_sw1(int nn, const int nprj, const int nmax,
                              const int lmax, int *n_prj, int *l_prj,
                              int *m_prj, double *G, double *sw1, double *sw2) {
  int a, b, na, nb;
  int one = 1;
  int zero = 0;
  double rzero = 0.0;
  int nmax2 = nmax * nmax;
  int nnn = nn * nprj;
  int nna = nn;

  // DCOPY_PWDFT(nnn,&rzero,zero,sw2,one);
  std::memset(sw2, 0, nnn * sizeof(double));

  for (b = 0; b < nprj; ++b)
    for (a = 0; a < nprj; ++a)
      if ((l_prj[a] == l_prj[b]) && (m_prj[a] == m_prj[b])) {
        na = n_prj[a] - 1;
        nb = n_prj[b] - 1;
        // na = n_prj[a];
        // nb = n_prj[b];

        DAXPY_PWDFT(nna, (G[nb + na * nmax + nmax2 * l_prj[a]]), &(sw1[a * nn]),
                    one, &(sw2[b * nn]), one);
      }
}

/*******************************************
 *                                         *
 *                vpp_read                 *
 *                                         *
 *******************************************/
/**
 * @brief Read data from a formatted PSP (Pseudopotential) file.
 *
 * This function reads various data from a formatted PSP file, including lattice information,
 * atomic properties, wavefunction-related data, and potential-related data. The data is
 * distributed among different variables and arrays passed as function parameters.
 *
 * @param mygrid A pointer to the grid structure.
 * @param fname The filename of the PSP file to read.
 * @param comment A character array to store comments from the file.
 * @param psp_type An integer to store the PSP type.
 * @param version An integer to store the PSP version.
 * @param nfft An integer array to store FFT grid dimensions.
 * @param unita A double array to store lattice vectors.
 * @param atom A character array to store the atomic symbol.
 * @param amass A double to store atomic mass.
 * @param zv A double to store the valence charge.
 * @param lmmax An integer to store (lmax + 1)^2 - (2 * locp + 1).
 * @param lmax An integer to store the maximum l value.
 * @param locp An integer to store the local pseudopotential flag.
 * @param nmax An integer to store the maximum n value.
 * @param rc A double array to store rc values.
 * @param nprj An integer to store the number of projectors.
 * @param n_projector An integer array to store n values for projectors.
 * @param l_projector An integer array to store l values for projectors.
 * @param m_projector An integer array to store m values for projectors.
 * @param b_projector An integer array to store b values for projectors.
 * @param Gijl A double array to store Gijl elements.
 * @param rlocal A double to store the local radius.
 * @param semicore A boolean flag indicating whether the PSP is semicore.
 * @param rcore A double to store the core radius.
 * @param ncore A double array to store core data.
 * @param vl A double array to store the vl data.
 * @param vnl A double array to store the vnl data.
 * @param log_amesh A double to store the log of the amesh value.
 * @param r1 A double to store the r1 value.
 * @param rmax A double to store the rmax value.
 * @param sigma A double to store the sigma value.
 * @param zion A double to store the zion value.
 * @param n1dgrid An integer to store the 1D grid dimension.
 * @param n1dbasis An integer to store the 1D basis dimension.
 * @param nae An integer array to store nae values.
 * @param nps An integer array to store nps values.
 * @param lps An integer array to store lps values.
 * @param icut An integer to store the icut value.
 * @param eig A double array to store eig values.
 * @param phi_ae A double array to store phi_ae values.
 * @param dphi_ae A double array to store dphi_ae values.
 * @param phi_ps A double array to store phi_ps values.
 * @param dphi_ps A double array to store dphi_ps values.
 * @param core_ae A double array to store core_ae values.
 * @param core_ps A double array to store core_ps values.
 * @param core_ae_prime A double array to store core_ae_prime values.
 * @param core_ps_prime A double array to store core_ps_prime values.
 * @param rgrid A double array to store rgrid values.
 * @param core_kin_energy A double to store the core kinetic energy.
 * @param core_ion_energy A double to store the core ionization energy.
 * @param hartree_matrix A double array to store the Hartree matrix.
 * @param comp_charge_matrix A double array to store the compensation charge matrix.
 * @param comp_pot_matrix A double array to store the compensation potential matrix.
 * @param coutput The output stream for printing.
 */
static void vpp_read(PGrid *mygrid, char *fname, char *comment, int *psp_type, int *version,
                     int *nfft, double *unita, char *atom, double *amass, double *zv,
                     int *lmmax, int *lmax, int *locp, int *nmax, double **rc, int *nprj,
                     int **n_projector, int **l_projector, int **m_projector,
                     int **b_projector, double **Gijl, double *rlocal, bool *semicore,
                     double *rcore, double **ncore, double *vl, double **vnl,
                     double *log_amesh, double *r1, double *rmax, double *sigma,
                     double *zion, int *n1dgrid, int *n1dbasis, int **nae, int **nps,
                     int **lps, int *icut, double **eig, double **phi_ae, double **dphi_ae,
                     double **phi_ps, double **dphi_ps, double **core_ae, double **core_ps,
                     double **core_ae_prime, double **core_ps_prime, double **rgrid,
                     double *core_kin_energy, double *core_ion_energy,
                     double **hartree_matrix, double **comp_charge_matrix,
                     double **comp_pot_matrix, std::ostream &coutput) 
{
   nwpw_timing_function ftimer(50);
   int i, nn;
   double *tmp2, *prj;
   Parallel *parall = mygrid->parall;
 
   *rlocal = 0.0;
 
   if (parall->base_stdio_print)
     coutput << std::endl
             << " reading formatted psp filename: " << fname << std::endl;
 
   if (parall->is_master()) {
     openfile(5, fname, "r");
     cread(5, comment, 80);
     comment[79] = '\0';
     i = 78;
     while (comment[i] == ' ')
       comment[i--] = '\0';
 
     iread(5, psp_type, 1);
     iread(5, version, 1);
     iread(5, nfft, 3);
     dread(5, unita, 9);
     cread(5, atom, 2);
     dread(5, amass, 1);
     dread(5, zv, 1);
     iread(5, lmax, 1);
     iread(5, locp, 1);
     iread(5, nmax, 1);
   }
   parall->Brdcst_cValues(0, 0, 80, comment);
   parall->Brdcst_iValue(0, 0, psp_type);
   parall->Brdcst_iValue(0, 0, version);
   parall->Brdcst_iValues(0, 0, 3, nfft);
   parall->Brdcst_Values(0, 0, 9, unita);
   parall->Brdcst_cValues(0, 0, 2, atom);
   parall->Brdcst_Values(0, 0, 1, amass);
   parall->Brdcst_Values(0, 0, 1, zv);
   parall->Brdcst_iValue(0, 0, lmax);
   parall->Brdcst_iValue(0, 0, locp);
   parall->Brdcst_iValue(0, 0, nmax);
   *lmmax = ((*lmax) + 1) * ((*lmax) + 1) - (2 * (*locp) + 1);
 
   if (*psp_type == 4)
     *n1dbasis = *locp;
   else
     *n1dbasis = *lmax + 1;
 
   *rc = new (std::nothrow) double[*lmax + 1]();
   if (parall->is_master()) {
     dread(5, *rc, *lmax + 1);
     iread(5, nprj, 1);
   }
   parall->Brdcst_Values(0, 0, *lmax + 1, *rc);
   parall->Brdcst_iValue(0, 0, nprj);
   if (*nprj > 0) {
     *n_projector = new int[*nprj]();
     *l_projector = new int[*nprj]();
     *m_projector = new int[*nprj]();
     *b_projector = new int[*nprj]();
     if (parall->is_master()) {
       iread(5, *n_projector, *nprj);
       iread(5, *l_projector, *nprj);
       iread(5, *m_projector, *nprj);
       iread(5, *b_projector, *nprj);
     }
     parall->Brdcst_iValues(0, 0, *nprj, *n_projector);
     parall->Brdcst_iValues(0, 0, *nprj, *l_projector);
     parall->Brdcst_iValues(0, 0, *nprj, *m_projector);
     parall->Brdcst_iValues(0, 0, *nprj, *b_projector);
 
     nn = (*nmax) * (*nmax) * (*lmax + 1);
     if (*psp_type == 4)
       nn *= 5;
     *Gijl = new (std::nothrow) double[nn]();
     if (parall->is_master()) {
       dread(5, *Gijl, nn);
     }
     parall->Brdcst_Values(0, 0, nn, *Gijl);
   }
 
   if (parall->is_master()) {
     if (*version == 4)
       dread(5, rlocal, 1);
     dread(5, rcore, 1);
   }
   parall->Brdcst_Values(0, 0, 1, rlocal);
   parall->Brdcst_Values(0, 0, 1, rcore);
   if (*rcore > 0.0)
     *semicore = true;
   else
     *semicore = false;
 
   /* read in miscellaneous paw energies and 1d wavefunctions */
   if (*psp_type == 4) {
     nn = (*n1dbasis) * (*n1dbasis) * (*n1dbasis) * (*n1dbasis) *
          (2 * (*lmax) + 1);
     *hartree_matrix = new (std::nothrow) double[nn]();
     if (parall->is_master())
       dread(5, *hartree_matrix, nn);
     parall->Brdcst_Values(0, 0, nn, *hartree_matrix);
 
     nn = (*n1dbasis) * (*n1dbasis) * (2 * (*lmax) + 1);
     *comp_charge_matrix = new (std::nothrow) double[nn]();
     if (parall->is_master())
       dread(5, *comp_charge_matrix, nn);
     parall->Brdcst_Values(0, 0, nn, *comp_charge_matrix);
 
     *comp_pot_matrix = new (std::nothrow) double[nn]();
     if (parall->is_master())
       dread(5, *comp_pot_matrix, nn);
     parall->Brdcst_Values(0, 0, nn, *comp_pot_matrix);
 
     if (parall->is_master()) {
       dread(5, core_kin_energy, 1);
       dread(5, core_ion_energy, 1);
     }
     parall->Brdcst_Values(0, 0, 1, core_kin_energy);
     parall->Brdcst_Values(0, 0, 1, core_ion_energy);
 
     if (parall->is_master()) {
       iread(5, n1dgrid, 1);
       iread(5, icut, 1);
       dread(5, log_amesh, 1);
       dread(5, r1, 1);
       dread(5, rmax, 1);
       dread(5, sigma, 1);
       dread(5, zion, 1);
     }
     parall->Brdcst_iValue(0, 0, n1dgrid);
     parall->Brdcst_iValue(0, 0, icut);
     parall->Brdcst_Values(0, 0, 1, log_amesh);
     parall->Brdcst_Values(0, 0, 1, r1);
     parall->Brdcst_Values(0, 0, 1, rmax);
     parall->Brdcst_Values(0, 0, 1, sigma);
     parall->Brdcst_Values(0, 0, 1, zion);
 
     if ((*n1dbasis) > 0) {
       nn = (*n1dbasis);
       *nae = new (std::nothrow) int[nn]();
       *nps = new (std::nothrow) int[nn]();
       *lps = new (std::nothrow) int[nn]();
 
       *eig = new (std::nothrow) double[nn]();
 
       nn = (*n1dbasis) * (*n1dgrid);
       *phi_ae = new (std::nothrow) double[nn]();
       *dphi_ae = new (std::nothrow) double[nn]();
       *phi_ps = new (std::nothrow) double[nn]();
       *dphi_ps = new (std::nothrow) double[nn]();
     }
     nn = (*n1dgrid);
     *core_ae = new (std::nothrow) double[nn]();
     *core_ps = new (std::nothrow) double[nn]();
     *core_ae_prime = new (std::nothrow) double[nn]();
     *core_ps_prime = new (std::nothrow) double[nn]();
     *rgrid = new (std::nothrow) double[nn]();
 
     nn = (*n1dbasis);
     if (parall->is_master())
       dread(5, *eig, nn);
     parall->Brdcst_Values(0, 0, nn, *eig);
 
     if (parall->is_master())
       iread(5, *nae, nn);
     parall->Brdcst_iValues(0, 0, nn, *nae);
 
     if (parall->is_master())
       iread(5, *nps, nn);
     parall->Brdcst_iValues(0, 0, nn, *nps);
 
     if (parall->is_master())
       iread(5, *lps, nn);
     parall->Brdcst_iValues(0, 0, nn, *lps);
 
     nn = (*n1dgrid);
     if (parall->is_master())
       dread(5, *rgrid, nn);
     parall->Brdcst_Values(0, 0, nn, *rgrid);
 
     nn = (*n1dbasis) * (*n1dgrid);
     if (parall->is_master())
       dread(5, *phi_ae, nn);
     parall->Brdcst_Values(0, 0, nn, *phi_ae);
 
     if (parall->is_master())
       dread(5, *dphi_ae, nn);
     parall->Brdcst_Values(0, 0, nn, *dphi_ae);
 
     if (parall->is_master())
       dread(5, *phi_ps, nn);
     parall->Brdcst_Values(0, 0, nn, *phi_ps);
 
     if (parall->is_master())
       dread(5, *dphi_ps, nn);
     parall->Brdcst_Values(0, 0, nn, *dphi_ps);
 
     nn = (*n1dbasis) * (*n1dgrid);
     if (parall->is_master())
       dread(5, *core_ae, nn);
     parall->Brdcst_Values(0, 0, nn, *core_ae);
 
     if (parall->is_master())
       dread(5, *core_ps, nn);
     parall->Brdcst_Values(0, 0, nn, *core_ps);
 
     if (parall->is_master())
       dread(5, *core_ae_prime, nn);
     parall->Brdcst_Values(0, 0, nn, *core_ae_prime);
 
     if (parall->is_master())
       dread(5, *core_ps_prime, nn);
     parall->Brdcst_Values(0, 0, nn, *core_ps_prime);
   }
 
   /* readin vl 3d block */
   tmp2 = new (std::nothrow) double[mygrid->nfft3d]();
   mygrid->t_read(5, tmp2, -1);
   mygrid->t_pack(0, tmp2);
   mygrid->tt_pack_copy(0, tmp2, vl);
 
   /* reading vnl 3d block */
   if (*nprj > 0) 
   {
      *vnl = new (std::nothrow) double[(*nprj) * (mygrid->npack(1))]();
      prj = *vnl;
      for (i = 0; i < (*nprj); ++i) {
         mygrid->t_read(5, tmp2, -1);
         mygrid->t_pack(1, tmp2);
         mygrid->tt_pack_copy(1, tmp2, &prj[i * mygrid->npack(1)]);
      }
   }
   if (*semicore) {
     nn = 5 * mygrid->npack(0);
     *ncore = new (std::nothrow) double[nn]();
     prj = *ncore;
 
     mygrid->t_read(5, tmp2, -1);
     mygrid->t_pack(0, tmp2);
     mygrid->tt_pack_copy(0, tmp2, prj);
 
     mygrid->t_read(5, tmp2, -1);
     mygrid->t_pack(0, tmp2);
     mygrid->tt_pack_copy(0, tmp2, &prj[2 * mygrid->npack(0)]);
 
     mygrid->t_read(5, tmp2, -1);
     mygrid->t_pack(0, tmp2);
     mygrid->tt_pack_copy(0, tmp2, &prj[3 * mygrid->npack(0)]);
 
     mygrid->t_read(5, tmp2, -1);
     mygrid->t_pack(0, tmp2);
     mygrid->tt_pack_copy(0, tmp2, &prj[4 * mygrid->npack(0)]);
   }
 
   delete[] tmp2;
 
   if (parall->is_master())
     closefile(5);
}


/*******************************************
 *                                         *
 *                vpp_write                *
 *                                         *
 *******************************************/
/**
 * @brief Write data to a formatted PSP (Pseudopotential) file.
 *
 * This function writes various data to a formatted PSP file, including lattice information,
 * atomic properties, wavefunction-related data, and potential-related data. The data is
 * taken from the provided input parameters and saved to the specified file.
 *
 * @param mygrid A pointer to the grid structure.
 * @param fname The filename of the PSP file to write.
 * @param comment A character array containing comments to include in the file.
 * @param psp_type An integer specifying the PSP type.
 * @param version An integer specifying the PSP version.
 * @param nfft An integer array containing FFT grid dimensions.
 * @param unita A double array containing lattice vectors.
 * @param atom A character array specifying the atomic symbol.
 * @param amass A double specifying the atomic mass.
 * @param zv A double specifying the valence charge.
 * @param lmmax An integer specifying (lmax + 1)^2 - (2 * locp + 1).
 * @param lmax An integer specifying the maximum l value.
 * @param locp An integer specifying the local pseudopotential flag.
 * @param nmax An integer specifying the maximum n value.
 * @param rc A double array containing rc values.
 * @param nprj An integer specifying the number of projectors.
 * @param n_projector An integer array containing n values for projectors.
 * @param l_projector An integer array containing l values for projectors.
 * @param m_projector An integer array containing m values for projectors.
 * @param b_projector An integer array containing b values for projectors.
 * @param Gijl A double array containing Gijl elements.
 * @param rlocal A double specifying the local radius.
 * @param semicore A boolean flag indicating whether the PSP is semicore.
 * @param rcore A double specifying the core radius.
 * @param ncore A double array containing core data.
 * @param vl A double array containing the vl data.
 * @param vnl A double array containing the vnl data.
 * @param log_amesh A double specifying the log of the amesh value.
 * @param r1 A double specifying the r1 value.
 * @param rmax A double specifying the rmax value.
 * @param sigma A double specifying the sigma value.
 * @param zion A double specifying the zion value.
 * @param n1dgrid An integer specifying the 1D grid dimension.
 * @param n1dbasis An integer specifying the 1D basis dimension.
 * @param nae An integer array containing nae values.
 * @param nps An integer array containing nps values.
 * @param lps An integer array containing lps values.
 * @param icut An integer specifying the icut value.
 * @param eig A double array containing eig values.
 * @param phi_ae A double array containing phi_ae values.
 * @param dphi_ae A double array containing dphi_ae values.
 * @param phi_ps A double array containing phi_ps values.
 * @param dphi_ps A double array containing dphi_ps values.
 * @param core_ae A double array containing core_ae values.
 * @param core_ps A double array containing core_ps values.
 * @param core_ae_prime A double array containing core_ae_prime values.
 * @param core_ps_prime A double array containing core_ps_prime values.
 * @param rgrid A double array containing rgrid values.
 * @param core_kin_energy A double specifying the core kinetic energy.
 * @param core_ion_energy A double specifying the core ionization energy.
 * @param hartree_matrix A double array containing the Hartree matrix.
 * @param comp_charge_matrix A double array containing the compensation charge matrix.
 * @param comp_pot_matrix A double array containing the compensation potential matrix.
 * @param coutput The output stream for printing.
 */
static void vpp_write(PGrid *mygrid, char *fname, char *comment, int psp_type, int version,
                      int *nfft, double *unita, char *atom, double amass, double zv, int lmmax,
                      int lmax, int locp, int nmax, double *rc, int nprj, int *n_projector,
                      int *l_projector, int *m_projector, int *b_projector, double *Gijl,
                      double rlocal, bool semicore, double rcore, double *ncore, double *vl,
                      double *vnl, double log_amesh, double r1, double rmax, double sigma,
                      double zion, int n1dgrid, int n1dbasis, int *nae, int *nps, int *lps,
                      int icut, double *eig, double *phi_ae, double *dphi_ae, double *phi_ps,
                      double *dphi_ps, double *core_ae, double *core_ps, double *core_ae_prime,
                      double *core_ps_prime, double *rgrid, double core_kin_energy,
                      double core_ion_energy, double *hartree_matrix, double *comp_charge_matrix,
                      double *comp_pot_matrix, std::ostream &coutput) 
{
   nwpw_timing_function ftimer(50);
   int i, nn;
   double *prj;
   Parallel *parall = mygrid->parall;
 
   // double tmp2[mygrid->nfft3d];
   double *tmp2 = new (std::nothrow) double[mygrid->nfft3d]();
 
   if (parall->base_stdio_print)
      coutput << std::endl
              << " writing formatted psp filename: " << fname << std::endl;
 
   if (parall->is_master()) 
   {
      openfile(6, fname, "w");
      comment[79] = '\0';
      cwrite(6, comment, 80);
     
      iwrite(6, &psp_type, 1);
      iwrite(6, &version, 1);
      iwrite(6, nfft, 3);
      dwrite(6, unita, 9);
      cwrite(6, atom, 2);
      dwrite(6, &amass, 1);
      dwrite(6, &zv, 1);
      iwrite(6, &lmax, 1);
      if (psp_type == 4)
         iwrite(6, &n1dbasis, 1);
      else
         iwrite(6, &locp, 1);
      iwrite(6, &nmax, 1);
     
      dwrite(6, rc, lmax + 1);
      iwrite(6, &nprj, 1);
     
      if (nprj > 0) 
      {
         iwrite(6, n_projector, nprj);
         iwrite(6, l_projector, nprj);
         iwrite(6, m_projector, nprj);
         iwrite(6, b_projector, nprj);
      }
      nn = (nmax) * (nmax) * (lmax + 1);
      if (psp_type == 4)
         nn *= 5;
      dwrite(6, Gijl, nn);
      if (version == 4)
         dwrite(6, &rlocal, 1);
      dwrite(6, &rcore, 1);
      // double x = mygrid->parall->SumAll(0,1.0); // Probably not needed!!
     
      // ***** Miscellaneous paw energies and 1d wavefunctions ****
      if (psp_type == 4) 
      {
         nn = n1dbasis * n1dbasis * n1dbasis * n1dbasis * (2 * lmax + 1);
         dwrite(6, hartree_matrix, nn);
        
         nn = n1dbasis * n1dbasis * (2 * lmax + 1);
         dwrite(6, comp_charge_matrix, nn);
         dwrite(6, comp_pot_matrix, nn);
        
         dwrite(6, &core_kin_energy, 1);
         dwrite(6, &core_ion_energy, 1);
        
         /* write 1d-wavefunctions */
         iwrite(6, &n1dgrid, 1);
         iwrite(6, &icut, 1);
         dwrite(6, &log_amesh, 1);
         dwrite(6, &r1, 1);
         dwrite(6, &rmax, 1);
         dwrite(6, &sigma, 1);
         dwrite(6, &zion, 1);
         dwrite(6, eig, n1dbasis);
         iwrite(6, nae, n1dbasis);
         iwrite(6, nps, n1dbasis);
         iwrite(6, lps, n1dbasis);
        
         dwrite(6, rgrid, n1dgrid);
         dwrite(6, phi_ae, n1dgrid * n1dbasis);
         dwrite(6, dphi_ae, n1dgrid * n1dbasis);
         dwrite(6, phi_ps, n1dgrid * n1dbasis);
         dwrite(6, dphi_ps, n1dgrid * n1dbasis);
         dwrite(6, core_ae, n1dgrid);
         dwrite(6, core_ps, n1dgrid);
         dwrite(6, core_ae_prime, n1dgrid);
         dwrite(6, core_ps_prime, n1dgrid);
      }
   }
 
   /* readin vl 3d block */
   mygrid->tt_pack_copy(0, vl, tmp2);
   mygrid->t_unpack(0, tmp2);
   mygrid->t_write_buffer(6, tmp2, 0);
 
   /* reading vnl 3d block */
   if (nprj > 0) 
   {
      prj = vnl;
      for (i = 0; i < (nprj); ++i) 
      {
         mygrid->tt_pack_copy(1, &prj[i * mygrid->npack(1)], tmp2);
         mygrid->t_unpack(1, tmp2);
         mygrid->t_write_buffer(6, tmp2, 0);
      }
   }
 
   if (semicore) 
   {
      prj = ncore;
     
      mygrid->tt_pack_copy(0, prj, tmp2);
      mygrid->t_unpack(0, tmp2);
      mygrid->t_write_buffer(6, tmp2, 0);
     
      mygrid->tt_pack_copy(0, &prj[2 * mygrid->npack(0)], tmp2);
      mygrid->t_unpack(0, tmp2);
      mygrid->t_write_buffer(6, tmp2, 0);
     
      mygrid->tt_pack_copy(0, &prj[3 * mygrid->npack(0)], tmp2);
      mygrid->t_unpack(0, tmp2);
      mygrid->t_write_buffer(6, tmp2, 0);
     
      mygrid->tt_pack_copy(0, &prj[4 * mygrid->npack(0)], tmp2);
      mygrid->t_unpack(0, tmp2);
      mygrid->t_write_buffer(6, tmp2, 0);
   }
 
   delete[] tmp2;
 
   if (parall->is_master())
      closefile(6);
}


/*******************************************
 *                                         *
 *            semicore_check               *
 *                                         *
 *******************************************/
/**
 * @brief Calculate a quantity related to the semicore pseudopotential.
 *
 * This function calculates a quantity related to the semicore pseudopotential
 * based on the provided grid and pseudopotential data. The result depends on
 * whether the pseudopotential is marked as semicore and the parameters such as
 * the core radius and core density.
 *
 * @param mygrid A pointer to the grid structure.
 * @param semicore A boolean flag indicating whether the pseudopotential is semicore.
 * @param rcore A double specifying the core radius.
 * @param ncore A double array containing core density data.
 * @return The calculated quantity related to the semicore pseudopotential.
 */
static double semicore_check(PGrid *mygrid, bool semicore, double rcore,
                             double *ncore) {
  double sum = 0.0;
  if (semicore) {
    double omega = mygrid->lattice->omega();
    double scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
    // double scal2 = 1.0/lattice_omega();
    // double dv    = lattice_omega()*scal1;
    double scal2 = 1.0 / omega;
    double dv = omega * scal1;
    double *tmp = mygrid->r_alloc();

    /* put sqrt(core-density) at atom position */
    mygrid->tc_pack_copy(0, ncore, tmp);
    mygrid->c_pack_SMul(0, scal2, tmp);

    /* Put put tmp into real space */
    mygrid->c_unpack(0, tmp);
    mygrid->cr_fft3d(tmp);
    mygrid->r_zero_ends(tmp);

    /*  square it  */
    mygrid->r_sqr(tmp);

    /* integrate it */
    sum = mygrid->r_dsum(tmp) * dv;

    mygrid->r_dealloc(tmp);
  }
  return sum;
}


/*******************************************
 *                                         *
 *          convert_psp_type               *
 *                                         *
 *******************************************/
/**
 * @brief Convert a character to an integer pseudopotential type.
 *
 * This function takes a character representing a pseudopotential type and
 * converts it to an integer value. The character '0' corresponds to an integer
 * value of 0, '1' corresponds to 1, '2' to 2, and so on up to '9' corresponding
 * to 9. Any other character input will result in a return value of 0.
 *
 * @param test A character representing the pseudopotential type.
 * @return The integer value corresponding to the pseudopotential type.
 */
static int convert_psp_type(char *test) {
  int psp_type = 0;
  if (test[0] == '0')
    psp_type = 0;
  if (test[0] == '1')
    psp_type = 1;
  if (test[0] == '2')
    psp_type = 2;
  if (test[0] == '3')
    psp_type = 3;
  if (test[0] == '4')
    psp_type = 4;
  if (test[0] == '5')
    psp_type = 5;
  if (test[0] == '6')
    psp_type = 6;
  if (test[0] == '7')
    psp_type = 7;
  if (test[0] == '8')
    psp_type = 8;
  if (test[0] == '9')
    psp_type = 9;

  return psp_type;
}


/*******************************************
 *                                         *
 *            vpp_get_psp_type             *
 *                                         *
 *******************************************/
/**
 * @brief Get the pseudopotential type from a pseudopotential file.
 *
 * This function reads the pseudopotential file specified by the `pspname` parameter
 * and extracts the pseudopotential type. The pseudopotential type is determined by
 * reading the first character in the file, which represents the type.
 *
 * @param myparall A pointer to the parallel processing context.
 * @param pspname  The name of the pseudopotential file to read.
 * @return The pseudopotential type as an integer value.
 */
static int vpp_get_psp_type(Parallel *myparall, char *pspname) {
   int psp_type;
   char atom[3];
 
   if (myparall->is_master()) 
   {
      FILE *fp = std::fopen(pspname, "r");
      std::fscanf(fp, "%s", atom);
      psp_type = convert_psp_type(atom);
      // if (psp_type>0) std::fscanf(fp,"%s",atom);
      fclose(fp);
   }
   myparall->Brdcst_iValue(0, 0, &psp_type);
 
   return psp_type;
}


/*******************************************
 *                                         *
 *                vpp_generate             *
 *                                         *
 *******************************************/
/**
 * @brief Generate a pseudopotential in various formats based on the specified type.
 *
 * This function generates a pseudopotential in various formats based on the provided
 * pseudopotential type. The generated pseudopotential data is stored in different
 * output arrays and variables, depending on the type of pseudopotential.
 *
 * @param mygrid            A pointer to the grid object.
 * @param pspname           The name of the pseudopotential file.
 * @param fname             The name of the output file.
 * @param comment           A comment to be included in the output file.
 * @param psp_type          A pointer to the pseudopotential type.
 * @param psp_version_in    The input pseudopotential version.
 * @param version           A pointer to the output pseudopotential version.
 * @param nfft              A pointer to the output array for FFT grid dimensions.
 * @param unita             A pointer to the output array for lattice unit vectors.
 * @param atom              A pointer to the output array for the atomic symbol.
 * @param amass             A pointer to the output atomic mass.
 * @param zv                A pointer to the output atomic number.
 * @param lmmax             A pointer to the output lmmax value.
 * @param lmax              A pointer to the output lmax value.
 * @param locp              A pointer to the output locp value.
 * @param nmax              A pointer to the output nmax value.
 * @param rc                A pointer to the output array for radial cutoff values.
 * @param nprj              A pointer to the output number of projectors.
 * @param n_projector       A pointer to the output array for n-projector values.
 * @param l_projector       A pointer to the output array for l-projector values.
 * @param m_projector       A pointer to the output array for m-projector values.
 * @param b_projector       A pointer to the output array for b-projector values.
 * @param Gijl              A pointer to the output array for Gijl values.
 * @param rlocal            A pointer to the output rlocal value.
 * @param semicore          A pointer to the output semicore flag.
 * @param rcore             A pointer to the output rcore value.
 * @param ncore             A pointer to the output array for core density values.
 * @param vl                A pointer to the output array for local potential values.
 * @param vlpaw             A pointer to the output array for local PAW potential values.
 * @param vnl               A pointer to the output array for non-local potential values.
 * @param log_amesh         A pointer to the output log of the amesh value.
 * @param r1                A pointer to the output r1 value.
 * @param rmax              A pointer to the output rmax value.
 * @param sigma             A pointer to the output sigma value.
 * @param zion              A pointer to the output zion value.
 * @param n1dgrid           A pointer to the output 1D grid size.
 * @param n1dbasis          A pointer to the output 1D basis size.
 * @param nae               A pointer to the output array for AE occupancies.
 * @param nps               A pointer to the output array for PS occupancies.
 * @param lps               A pointer to the output array for PS angular momenta.
 * @param icut              A pointer to the output icut value.
 * @param eig               A pointer to the output array for eigenvalues.
 * @param phi_ae            A pointer to the output array for AE wavefunctions.
 * @param dphi_ae           A pointer to the output array for AE wavefunction derivatives.
 * @param phi_ps            A pointer to the output array for PS wavefunctions.
 * @param dphi_ps           A pointer to the output array for PS wavefunction derivatives.
 * @param core_ae           A pointer to the output array for AE core densities.
 * @param core_ps           A pointer to the output array for PS core densities.
 * @param core_ae_prime     A pointer to the output array for AE core density derivatives.
 * @param core_ps_prime     A pointer to the output array for PS core density derivatives.
 * @param rgrid             A pointer to the output array for radial grid values.
 * @param core_kin_energy   A pointer to the output AE core kinetic energy.
 * @param core_ion_energy   A pointer to the output AE core ionization energy.
 * @param hartree_matrix    A pointer to the output array for Hartree matrix elements.
 * @param comp_charge_matrix A pointer to the output array for compensated charge matrix elements.
 * @param comp_pot_matrix   A pointer to the output array for compensated potential matrix elements.
 * @param coutput           The output stream for diagnostic messages.
 */
static void vpp_generate(PGrid *mygrid, char *pspname, char *fname, char *comment, int *psp_type,
                         int psp_version_in, int *version, int *nfft, double *unita, char *atom,
                         double *amass, double *zv, int *lmmax, int *lmax, int *locp, int *nmax,
                         double **rc, int *nprj, int **n_projector, int **l_projector,
                         int **m_projector, int **b_projector, double **Gijl, double *rlocal,
                         bool *semicore, double *rcore, double **ncore, double *vl, double *vlpaw,
                         double **vnl, double *log_amesh, double *r1, double *rmax, double *sigma,
                         double *zion, int *n1dgrid, int *n1dbasis, int **nae, int **nps, int **lps,
                         int *icut, double **eig, double **phi_ae, double **dphi_ae, double **phi_ps,
                         double **dphi_ps, double **core_ae, double **core_ps,
                         double **core_ae_prime, double **core_ps_prime, double **rgrid,
                         double *core_kin_energy, double *core_ion_energy, double **hartree_matrix,
                         double **comp_charge_matrix, double **comp_pot_matrix,
                         std::ostream &coutput) 
{
   int i, nn;
   double *tmp2, *prj;
   Parallel *myparall = mygrid->parall;
 
   *psp_type = vpp_get_psp_type(myparall,pspname);
   //std::cout << "psp_type=" << *psp_type << std::endl;
   //if (*psp_type==2) *psp_type= 0;
 
   if ((*psp_type == 0) || (*psp_type == 9)) 
   {
      int nray = mygrid->n_ray();
      Psp1d_Hamann psp1d(myparall,pspname,psp_version_in);
     
      nfft[0] = mygrid->nx;
      nfft[1] = mygrid->ny;
      nfft[2] = mygrid->nz;
      for (auto i = 0; i < 9; ++i)
         unita[i] = mygrid->lattice->unita1d(i);
     
      atom[0] = psp1d.atom[0];
      atom[1] = psp1d.atom[1];
      *amass = psp1d.amass;
      *zv = psp1d.zv;
      for (auto i = 0; i < 80; ++i)
         comment[i] = psp1d.comment[i];
     
      *psp_type = psp1d.psp_type;
      *version = psp1d.version;
      *lmax = psp1d.lmax;
      *locp = psp1d.locp;
      *nmax = psp1d.nmax;
      *lmmax = ((*lmax) + 1) * ((*lmax) + 1) - (2 * (*locp) + 1);
     
      *nprj = psp1d.nprj;
      *semicore = psp1d.semicore;
      *rcore = psp1d.rcore;
      *rlocal = psp1d.rlocal;
     
      *rc = new (std::nothrow) double[*lmax + 1]();
      for (auto l = 0; l <= (*lmax); ++l)
         (*rc)[l] = psp1d.rc[l];
     
      /* allocate Gijl and copy from psp1d */
      int nn = (psp1d.nmax) * (psp1d.nmax) * (psp1d.lmax + 1);
      *Gijl = new (std::nothrow) double[nn]();
      for (auto l = 0; l < nn; ++l) 
      {
         (*Gijl)[l] = psp1d.vnlnrm[l];
      }
     
      /* allocate n_projector, l_projector, m_projector, and b_projector and copy
       * from psp1d */
      if (psp1d.nprj > 0) 
      {
         *n_projector = new (std::nothrow) int[psp1d.nprj]();
         *l_projector = new (std::nothrow) int[psp1d.nprj]();
         *m_projector = new (std::nothrow) int[psp1d.nprj]();
         *b_projector = new (std::nothrow) int[psp1d.nprj]();
        
         for (auto l = 0; l < psp1d.nprj; ++l) 
         {
            (*n_projector)[l] = psp1d.n_prj[l];
            (*l_projector)[l] = psp1d.l_prj[l];
            (*m_projector)[l] = psp1d.m_prj[l];
            (*b_projector)[l] = psp1d.b_prj[l];
         }
      }
     
      /*  allocate and generate ray formatted grids */
      double *G_ray = mygrid->generate_G_ray();
      double *vl_ray = new (std::nothrow) double[nray]();
      double *vnl_ray = new (std::nothrow) double[(psp1d.lmax + 1 + psp1d.n_extra) * nray]();
      double *rho_sc_k_ray = new (std::nothrow) double[2 * nray]();
      psp1d.vpp_generate_ray(myparall, nray, G_ray, vl_ray, vnl_ray, rho_sc_k_ray);
     
      /* filter the ray formatted grids */
      double ecut = mygrid->lattice->ecut();
      double wcut = mygrid->lattice->wcut();
      util_filter(nray, G_ray, ecut, vl_ray);
      for (auto l = 0; l < (psp1d.lmax + 1 + psp1d.n_extra); ++l)
         util_filter(nray, G_ray, wcut, &(vnl_ray[l * nray]));
      if (*semicore) 
      {
         util_filter(nray, G_ray, ecut, rho_sc_k_ray);
         util_filter(nray, G_ray, ecut, &(rho_sc_k_ray[nray]));
      }
     
      /* allocate vnl and ncore  generate formated grids */
      *vnl = new (std::nothrow) double[(psp1d.nprj) * (mygrid->npack(1))]();
      if (*semicore)
         *ncore = new (std::nothrow) double[5 * mygrid->npack(0)]();
     
      /*  generate formatted grids using splines */
      psp1d.vpp_generate_spline(mygrid, nray, G_ray, vl_ray, vnl_ray,
                                rho_sc_k_ray, vl, *vnl, *ncore);
     
      /* deallocate ray formatted grids */
      delete[] rho_sc_k_ray;
      delete[] vnl_ray;
      delete[] vl_ray;
      delete[] G_ray;
 
   } 
   else if (*psp_type == 1) 
   {
      if (myparall->base_stdio_print)
         coutput << "in vpp_generate Not finished, hghppv1 psp_type = " << *psp_type << std::endl;
   } 
   else if (*psp_type == 2) 
   {
      std::cout << "im here kbppv3e psp_version_in=" << psp_version_in <<  std::endl;
      int nray = mygrid->n_ray();
      Psp1d_Hamann psp1d(myparall,pspname,psp_version_in);

      nfft[0] = mygrid->nx;
      nfft[1] = mygrid->ny;
      nfft[2] = mygrid->nz;
      for (auto i = 0; i < 9; ++i)
         unita[i] = mygrid->lattice->unita1d(i);

      atom[0] = psp1d.atom[0];
      atom[1] = psp1d.atom[1];
      *amass = psp1d.amass;
      *zv = psp1d.zv;
      for (auto i = 0; i < 80; ++i)
         comment[i] = psp1d.comment[i];
      
      *psp_type = psp1d.psp_type;
      *version = psp1d.version;
      *lmax = psp1d.lmax;
      *locp = psp1d.locp;
      *nmax = psp1d.nmax;
      *lmmax = ((*lmax) + 1) * ((*lmax) + 1) - (2 * (*locp) + 1);
      
      *nprj = psp1d.nprj;
      *semicore = psp1d.semicore;
      *rcore = psp1d.rcore;
      *rlocal = psp1d.rlocal;
      
      *rc = new (std::nothrow) double[*lmax + 1]();
      for (auto l = 0; l <= (*lmax); ++l)
         (*rc)[l] = psp1d.rc[l];

      /* allocate Gijl and copy from psp1d */
      int nn = (psp1d.nmax) * (psp1d.nmax) * (psp1d.lmax + 1);
      *Gijl = new (std::nothrow) double[nn]();
      for (auto l = 0; l < nn; ++l)  
      {
         (*Gijl)[l] = psp1d.vnlnrm[l];
      }

      /* allocate n_projector, l_projector, m_projector, and b_projector and copy
       * from psp1d */
      if (psp1d.nprj > 0) 
      {
         *n_projector = new (std::nothrow) int[psp1d.nprj]();
         *l_projector = new (std::nothrow) int[psp1d.nprj]();
         *m_projector = new (std::nothrow) int[psp1d.nprj]();
         *b_projector = new (std::nothrow) int[psp1d.nprj]();
        
         for (auto l = 0; l < psp1d.nprj; ++l)
         {
            (*n_projector)[l] = psp1d.n_prj[l];
            (*l_projector)[l] = psp1d.l_prj[l];
            (*m_projector)[l] = psp1d.m_prj[l];
            (*b_projector)[l] = psp1d.b_prj[l];
         }
      }
     
      /*  allocate and generate ray formatted grids */
      double *G_ray = mygrid->generate_G_ray();
      double *vl_ray = new (std::nothrow) double[nray]();
      double *vnl_ray = new (std::nothrow) double[(psp1d.lmax + 1 + psp1d.n_extra) * nray]();
      double *rho_sc_k_ray = new (std::nothrow) double[2 * nray]();
      psp1d.vpp_generate_ray(myparall, nray, G_ray, vl_ray, vnl_ray, rho_sc_k_ray);
            
      /* filter the ray formatted grids */
      double ecut = mygrid->lattice->ecut();
      double wcut = mygrid->lattice->wcut();
      util_filter(nray, G_ray, ecut, vl_ray);
      for (auto l = 0; l < (psp1d.lmax + 1 + psp1d.n_extra); ++l)
         util_filter(nray, G_ray, wcut, &(vnl_ray[l * nray]));
      if (*semicore)
      {
         util_filter(nray, G_ray, ecut, rho_sc_k_ray);
         util_filter(nray, G_ray, ecut, &(rho_sc_k_ray[nray]));
      }
      
      /* allocate vnl and ncore  generate formated grids */
      *vnl = new (std::nothrow) double[(psp1d.nprj) * (mygrid->npack(1))]();
      if (*semicore)
         *ncore = new (std::nothrow) double[5 * mygrid->npack(0)]();

      /*  generate formatted grids using splines */
      psp1d.vpp_generate_spline(mygrid, nray, G_ray, vl_ray, vnl_ray,
                                rho_sc_k_ray, vl, *vnl, *ncore);

      /* deallocate ray formatted grids */
      delete[] rho_sc_k_ray;
      delete[] vnl_ray;
      delete[] vl_ray;
      delete[] G_ray;


      if (myparall->base_stdio_print)
         coutput << "in vpp_generate Not finished, kbppv3e psp_type = " << *psp_type << std::endl;
   } 
   else if (*psp_type == 4) 
   {
      if (myparall->base_stdio_print)
         coutput << "in vpp_generate Not finished, pawppv1 psp_type = " << *psp_type << std::endl;
     
      int nray = mygrid->n_ray();
      Psp1d_pawppv1 paw1d(myparall,pspname);
     
      nfft[0] = mygrid->nx;
      nfft[1] = mygrid->ny;
      nfft[2] = mygrid->nz;
      for (auto i=0; i<9; ++i)
         unita[i] = mygrid->lattice->unita1d(i);
     
      atom[0] = paw1d.atom[0];
      atom[1] = paw1d.atom[1];
      *amass = paw1d.amass;
      *zv = paw1d.zv;
      for (auto i = 0; i < 80; ++i)
         comment[i] = paw1d.comment[i];
     
      *psp_type = paw1d.psp_type;
      *version = paw1d.version;
      *lmax = paw1d.lmax;
      *locp = paw1d.locp;
      *nmax = paw1d.nmax;
      //*lmmax=((*lmax)+1)*((*lmax)+1) - (2*(*locp)+1);
     
      *nprj = paw1d.nprj;
      *semicore = false;
      *rcore = 0.0;
      *rlocal = paw1d.rlocal;
      *rc = new (std::nothrow) double[*lmax + 1]();
      for (auto l = 0; l <= (*lmax); ++l)
         (*rc)[l] = paw1d.rc[l];
     
      /* allocate n_projector, l_projector, m_projector, and b_projector and copy
       * from psp1d */
      if (paw1d.nprj > 0) 
      {
         *n_projector = new int[paw1d.nprj]();
         *l_projector = new int[paw1d.nprj]();
         *m_projector = new int[paw1d.nprj]();
         *b_projector = new int[paw1d.nprj]();
        
         for (auto l = 0; l < paw1d.nprj; ++l) 
         {
            (*n_projector)[l] = paw1d.n_prj[l];
            (*l_projector)[l] = paw1d.l_prj[l];
            (*m_projector)[l] = paw1d.m_prj[l];
            (*b_projector)[l] = paw1d.b_prj[l];
         }
      }
     
      *log_amesh = paw1d.log_amesh;
      *r1 = paw1d.r1;
      *rmax = paw1d.rmax;
      *sigma = paw1d.sigma;
      *zion = paw1d.zion;
      *n1dgrid = paw1d.n1dgrid;
      *n1dbasis = paw1d.nbasis;
      *icut = paw1d.icut;
     
      *core_kin_energy = paw1d.core_kin_energy;
      *core_ion_energy = paw1d.core_ion_energy;
     
      if (paw1d.nbasis > 0) 
      {
         *nae = new (std::nothrow) int[paw1d.nbasis]();
         *nps = new (std::nothrow) int[paw1d.nbasis]();
         *lps = new (std::nothrow) int[paw1d.nbasis]();
        
         *eig = new double[paw1d.nbasis]();
        
         *phi_ae  = new (std::nothrow) double[paw1d.nbasis * paw1d.n1dgrid]();
         *dphi_ae = new (std::nothrow) double[paw1d.nbasis * paw1d.n1dgrid]();
         *phi_ps  = new (std::nothrow) double[paw1d.nbasis * paw1d.n1dgrid]();
         *dphi_ps = new (std::nothrow) double[paw1d.nbasis * paw1d.n1dgrid]();
      }
      *core_ae       = new (std::nothrow) double[paw1d.n1dgrid]();
      *core_ps       = new (std::nothrow) double[paw1d.n1dgrid]();
      *core_ae_prime = new (std::nothrow) double[paw1d.n1dgrid]();
      *core_ps_prime = new (std::nothrow) double[paw1d.n1dgrid]();
      *rgrid         = new (std::nothrow) double[paw1d.n1dgrid]();
     
      for (auto l=0; l<paw1d.nbasis; ++l) 
      {
         (*nae)[l] = paw1d.nae[l];
         (*nps)[l] = paw1d.nps[l];
         (*lps)[l] = paw1d.lps[l];
        
         (*eig)[l] = paw1d.eig[l];
      }
      for (auto il = 0; il < (paw1d.nbasis * paw1d.n1dgrid); ++il) 
      {
         (*phi_ae)[il] = paw1d.phi_ae[il];
         (*dphi_ae)[il] = paw1d.dphi_ae[il];
         (*phi_ps)[il] = paw1d.phi_ps[il];
         (*dphi_ps)[il] = paw1d.dphi_ps[il];
      }
      for (auto i = 0; i < (paw1d.n1dgrid); ++i) 
      {
         (*core_ae)[i] = paw1d.core_ae[i];
         (*core_ps)[i] = paw1d.core_ps[i];
         (*core_ae_prime)[i] = paw1d.core_ae_prime[i];
         (*core_ps_prime)[i] = paw1d.core_ps_prime[i];
         (*rgrid)[i] = paw1d.rgrid[i];
      }
     
      /*  allocate paw matrices */
      *Gijl               = new (std::nothrow) double[(paw1d.nmax)*(paw1d.nmax)*(paw1d.lmax+1)*5]();
      *comp_charge_matrix = new (std::nothrow) double[(paw1d.nbasis)*(paw1d.nbasis)*(2*paw1d.lmax+1)]();
      *comp_pot_matrix    = new (std::nothrow) double[(paw1d.nbasis)*(paw1d.nbasis)*(2*paw1d.lmax+1)]();
      *hartree_matrix     = new (std::nothrow) double[(paw1d.nbasis)*(paw1d.nbasis)*(paw1d.nbasis) 
                                                      *(paw1d.nbasis)*(2*paw1d.lmax+1)]();
     
      paw1d.vpp_generate_paw_matrices(myparall,*Gijl,*comp_charge_matrix,*comp_pot_matrix,*hartree_matrix);
     
      /*  allocate and generate ray formatted grids */
      double *G_ray = mygrid->generate_G_ray();
      double *vl_ray = new (std::nothrow) double[nray]();
      double *vlpaw_ray = new (std::nothrow) double[nray]();
      double *vnl_ray = new (std::nothrow) double[(paw1d.nbasis) * nray]();
     
      paw1d.vpp_generate_ray(myparall, nray, G_ray, vl_ray, vlpaw_ray, vnl_ray);
     
      /* filter the ray formatted grids */
      double ecut = mygrid->lattice->ecut();
      double wcut = mygrid->lattice->wcut();
      util_filter(nray,G_ray,ecut,vl_ray);
      util_filter(nray,G_ray,ecut,vlpaw_ray);
      for (auto l=0; l<(paw1d.nbasis); ++l)
         util_filter(nray,G_ray,wcut,&(vnl_ray[l*nray]));
     
      /* allocate vnl and other paw data then generate formated grids */
      *vnl = new (std::nothrow) double[(paw1d.nprj) * (mygrid->npack(1))]();
     
      /*  generate formatted grids using splines */
      paw1d.vpp_generate_spline(mygrid, nray, G_ray, vl_ray, vlpaw_ray, vnl_ray,
                                vl, vlpaw, *vnl);
     
      /* deallocate ray formatted grids */
      delete[] vnl_ray;
      delete[] vlpaw_ray;
      delete[] vl_ray;
      delete[] G_ray;
   } 
   else 
   {
      if (myparall->base_stdio_print)
         coutput << "in vpp_generate Not finished, psp_type = " << *psp_type << std::endl;
   }
   std::cout << "vpp_generate done" << std::endl;
}


/* Constructors */

/*******************************************
 *                                         *
 *     Pseudopotential::Pseudopotential    *
 *                                         *
 *******************************************/
/**
 * @brief Constructor for the Pseudopotential class.
 *
 * This constructor initializes and populates various data structures and objects related to pseudopotentials
 * for electronic structure calculations.
 *
 * @param myionin Pointer to the Ion object.
 * @param mypnebin Pointer to the Pneb object.
 * @param mystrfacin Pointer to the Strfac object.
 * @param control Reference to the Control2 object.
 * @param coutput Reference to the output stream.
 *
 * @note This constructor performs initialization, data reading, and setup for pseudopotentials and PAW pseudopotentials, if applicable.
 */
Pseudopotential::Pseudopotential(Ion *myionin, Pneb *mypnebin,
                                 Strfac *mystrfacin, Control2 &control,
                                 std::ostream &coutput) 
{
   int ia,version,nfft[3];
   int *n_ptr,*l_ptr,*m_ptr,*b_ptr;

   double *rc_ptr,*G_ptr,*vnl_ptr,*ncore_ptr;
   int *nae_ptr,*nps_ptr,*lps_ptr;

   double *eig_ptr,*phi_ae_ptr,*dphi_ae_ptr,*phi_ps_ptr,*dphi_ps_ptr;
   double *core_ae_ptr, *core_ps_ptr, *core_ae_prime_ptr, *core_ps_prime_ptr,*rgrid_ptr;
   double *hartree_matrix_ptr,*comp_charge_matrix_ptr,*comp_pot_matrix_ptr;
 
   double unita[9];
   char fname[256], pspname[256], aname[2];
   char fname2[256];
 
   myion = myionin;
   mypneb = mypnebin;
   mystrfac = mystrfacin;
 
   myefield = new nwpw_efield(myion,mypneb,mystrfac,control,coutput);
   myapc    = new nwpw_apc(myion,mypneb,mystrfac,control,coutput);
   mydipole = new nwpw_dipole(myion,mypneb,mystrfac,control);
 
   psp_version = control.version;

   // pspspin
   pspspin = control.pspspin();
   if (pspspin) 
   {
      int nion = myion->nion;
      pspspin_upions   = new bool[nion];
      pspspin_downions = new bool[nion];
      pspspin_upl = new int[nion];
      pspspin_upm = new int[nion];
      pspspin_downl = new int[nion];
      pspspin_downm = new int[nion];
      pspspin_upscale   = new double[nion];
      pspspin_downscale = new double[nion];
      for (auto ii=0; ii<nion; ++ii)
      { 
         pspspin_upions[ii] = false;
         pspspin_downions[ii] = false;
         pspspin_upl[ii] = -999;
         pspspin_upm[ii] = -999;
         pspspin_downl[ii] = -999;
         pspspin_downm[ii] = -999;
         pspspin_upscale[ii] = 1.0;
         pspspin_downscale[ii] = 1.0;
      }
      std::string pspspin_output = control.set_pspspin(nion,
                                                       pspspin_upscale,pspspin_downscale,
                                                       pspspin_upl, pspspin_downl,
                                                       pspspin_upm, pspspin_downm,
                                                       pspspin_upions, pspspin_downions);
      std::cout << "pspspin=" << pspspin << std::endl;
      std::cout << "UPIONS =" << pspspin_upions[0] << std::endl;
      std::cout << "DNIONS =" << pspspin_downions[0] << std::endl;
      if (mypneb->d3db::parall->is_master())
         coutput << pspspin_output << std::endl;
   }
 
   npsp = myion->nkatm;
   nprj_max = 0;
 
   psp_type = new int[npsp]();
   lmax     = new int[npsp]();
   lmmax    = new int[npsp]();
   locp     = new int[npsp]();
   nmax     = new int[npsp]();
   nprj     = new int[npsp]();
   semicore = new bool[npsp + 1]();
 
   n_projector = new int *[npsp]();
   l_projector = new int *[npsp]();
   m_projector = new int *[npsp]();
   b_projector = new int *[npsp]();
 
   zv        = new (std::nothrow) double[npsp]();
   amass     = new (std::nothrow) double[npsp]();
   rlocal    = new (std::nothrow) double[npsp]();
   rcore     = new (std::nothrow) double[npsp]();
   ncore_sum = new (std::nothrow) double[npsp]();
   rc        = new (std::nothrow) double *[npsp]();
 
   vl = new (std::nothrow) double *[npsp]();
   for (ia=0; ia<npsp; ++ia)
      vl[ia] = new (std::nothrow) double[mypneb->npack(0)]();
 
   Gijl       = new (std::nothrow) double *[npsp]();
   vnl        = new (std::nothrow) double *[npsp]();
   ncore_atom = new (std::nothrow) double *[npsp]();
   semicore_density = mypneb->r_alloc();
 
   comment = new char *[npsp]();
   for (ia=0; ia<npsp; ++ia)
      comment[ia] = new char[80]();
 
   semicore[npsp] = false;
 
   // *** paw data  ***
   pawexist = false;
   vlpaw = new (std::nothrow) double *[npsp]();
   for (ia=0; ia<npsp; ++ia)
      vlpaw[ia] = new (std::nothrow) double[mypneb->npack(0)]();
 
   hartree_matrix     = new (std::nothrow) double *[npsp]();
   comp_charge_matrix = new (std::nothrow) double *[npsp]();
   comp_pot_matrix    = new (std::nothrow) double *[npsp]();
 
   rgrid   = new (std::nothrow) double *[npsp]();
   eig     = new (std::nothrow) double *[npsp]();
   phi_ae  = new (std::nothrow) double *[npsp]();
   dphi_ae = new (std::nothrow) double *[npsp]();
   phi_ps  = new (std::nothrow) double *[npsp]();
   dphi_ps = new (std::nothrow) double *[npsp]();
   core_ae = new (std::nothrow) double *[npsp]();
   core_ps = new (std::nothrow) double *[npsp]();
   core_ae_prime = new (std::nothrow) double *[npsp]();
   core_ps_prime = new (std::nothrow) double *[npsp]();
 
   log_amesh = new (std::nothrow) double[npsp]();
   r1        = new (std::nothrow) double[npsp]();
   rmax      = new (std::nothrow) double[npsp]();
   sigma     = new (std::nothrow) double[npsp]();
   zion      = new (std::nothrow) double[npsp]();
   core_kin  = new (std::nothrow) double[npsp]();
   core_ion  = new (std::nothrow) double[npsp]();
 
   n1dgrid  = new (std::nothrow) int[npsp]();
   n1dbasis = new (std::nothrow) int[npsp]();
   icut     = new (std::nothrow) int[npsp]();
   nae      = new (std::nothrow) int *[npsp]();
   nps      = new (std::nothrow) int *[npsp]();
   lps      = new (std::nothrow) int *[npsp]();
 
   for (ia = 0; ia < npsp; ++ia) 
   {
      strcpy(fname, myion->atom(ia));
      strcat(fname, ".vpp");
      control.add_permanent_dir(fname);
     
      if (vpp_formatter_check(mypneb, fname, psp_version)) 
      {
         strcpy(pspname, myion->atom(ia));
         strcat(pspname, ".psp");
         control.add_permanent_dir(pspname);
         vpp_generate(mypneb, pspname, fname, comment[ia], &psp_type[ia], psp_version,
                      &version, nfft, unita, aname, &amass[ia], &zv[ia], &lmmax[ia],
                      &lmax[ia], &locp[ia], &nmax[ia], &rc_ptr, &nprj[ia], &n_ptr, &l_ptr,
                      &m_ptr, &b_ptr, &G_ptr, &rlocal[ia], &semicore[ia], &rcore[ia],
                      &ncore_ptr, vl[ia], vlpaw[ia], &vnl_ptr, &log_amesh[ia], &r1[ia],
                      &rmax[ia], &sigma[ia], &zion[ia], &n1dgrid[ia], &n1dbasis[ia],
                      &nae_ptr, &nps_ptr, &lps_ptr, &icut[ia], &eig_ptr, &phi_ae_ptr,
                      &dphi_ae_ptr, &phi_ps_ptr, &dphi_ps_ptr, &core_ae_ptr, &core_ps_ptr,
                      &core_ae_prime_ptr, &core_ps_prime_ptr, &rgrid_ptr, &core_kin[ia],
                      &core_ion[ia], &hartree_matrix_ptr, &comp_charge_matrix_ptr,
                      &comp_pot_matrix_ptr, coutput);
        
         // writing .vpp file to fname
         vpp_write(mypneb, fname, comment[ia], psp_type[ia], version, nfft, unita,
                   aname, amass[ia], zv[ia], lmmax[ia], lmax[ia], locp[ia],
                   nmax[ia], rc_ptr, nprj[ia], n_ptr, l_ptr, m_ptr, b_ptr, G_ptr,
                   rlocal[ia], semicore[ia], rcore[ia], ncore_ptr, vl[ia], vnl_ptr,
                   log_amesh[ia], r1[ia], rmax[ia], sigma[ia], zion[ia],
                   n1dgrid[ia], n1dbasis[ia], nae_ptr, nps_ptr, lps_ptr, icut[ia],
                   eig_ptr, phi_ae_ptr, dphi_ae_ptr, phi_ps_ptr, dphi_ps_ptr,
                   core_ae_ptr, core_ps_ptr, core_ae_prime_ptr, core_ps_prime_ptr,
                   rgrid_ptr, core_kin[ia], core_ion[ia], hartree_matrix_ptr,
                   comp_charge_matrix_ptr, comp_pot_matrix_ptr, coutput);
      } 
      else 
      {
         vpp_read(mypneb, fname, comment[ia], &psp_type[ia], &version, nfft, unita,
                  aname, &amass[ia], &zv[ia], &lmmax[ia], &lmax[ia], &locp[ia],
                  &nmax[ia], &rc_ptr, &nprj[ia], &n_ptr, &l_ptr, &m_ptr, &b_ptr,
                  &G_ptr, &rlocal[ia], &semicore[ia], &rcore[ia], &ncore_ptr,
                  vl[ia], &vnl_ptr, &log_amesh[ia], &r1[ia], &rmax[ia], &sigma[ia],
                  &zion[ia], &n1dgrid[ia], &n1dbasis[ia], &nae_ptr, &nps_ptr,
                  &lps_ptr, &icut[ia], &eig_ptr, &phi_ae_ptr, &dphi_ae_ptr,
                  &phi_ps_ptr, &dphi_ps_ptr, &core_ae_ptr, &core_ps_ptr,
                  &core_ae_prime_ptr, &core_ps_prime_ptr, &rgrid_ptr,
                  &core_kin[ia], &core_ion[ia], &hartree_matrix_ptr,
                  &comp_charge_matrix_ptr, &comp_pot_matrix_ptr, coutput);
      }
     
      rc[ia]          = rc_ptr;
      n_projector[ia] = n_ptr;
      l_projector[ia] = l_ptr;
      m_projector[ia] = m_ptr;
      b_projector[ia] = b_ptr;
      Gijl[ia]        = G_ptr;
      vnl[ia]         = vnl_ptr;
      if (nprj[ia] > nprj_max)
         nprj_max = nprj[ia];
     
      if (semicore[ia]) 
      {
         ncore_atom[ia] = ncore_ptr;
         ncore_sum[ia]  = semicore_check(mypneb,semicore[ia],rcore[ia],ncore_atom[ia]);
         semicore[npsp] = true;
      }
     
      if (psp_type[ia] == 4) 
      {
         pawexist = true;
         nae[ia] = nae_ptr;
         nps[ia] = nps_ptr;
         lps[ia] = lps_ptr;
         eig[ia] = eig_ptr;
         phi_ae[ia]        = phi_ae_ptr;
         dphi_ae[ia]       = dphi_ae_ptr;
         phi_ps[ia]        = phi_ps_ptr;
         dphi_ps[ia]       = dphi_ps_ptr, core_ae[ia] = core_ae_ptr;
         core_ps[ia]       = core_ps_ptr;
         core_ae_prime[ia] = core_ae_prime_ptr;
         core_ps_prime[ia] = core_ps_prime_ptr;
         rgrid[ia]              = rgrid_ptr;
         hartree_matrix[ia]     = hartree_matrix_ptr;
         comp_charge_matrix[ia] = comp_charge_matrix_ptr;
         comp_pot_matrix[ia]    = comp_pot_matrix_ptr;
      }
   }
 
   /* define the maximum number of projectors  */
   nprj_max *= 10;
   // nprj_max = 0;
   // for (auto ii=0; ii < myion->nion; ++ii)
   //     nprj_max += nprj[myion->katm[ii]];
 
   if (pawexist) 
   {
      // call psp_paw_init()
      mypaw_compcharge = new Paw_compcharge(myion,mypneb,control,nprj,n1dbasis,psp_type,
                                           lmax,sigma,nprj_max,l_projector,m_projector,
                                            b_projector,comp_charge_matrix,hartree_matrix);
      mypaw_xc = new Paw_xc(myion,mypneb,control,nprj,n1dbasis,n1dgrid,psp_type,
                            lmax,l_projector,m_projector,b_projector);
     
      // Paw_compcharge_init(ion_nion(),npsp,
      //                           int_mb(nprj(1)),
      //                           int_mb(n1dbasis(1)),
      //                           int_mb(psp_type(1)),
      //                           int_mb(lmax(1)),
      //                           dbl_mb(sigma(1)),
      //                           nmax_max*lmmax_max,
      //                           int_mb(l_projector(1)),
      //                           int_mb(m_projector(1)),
      //                           int_mb(b_projector(1)),
      //                           int_mb(comp_charge_matrix(1)),
      //                           int_mb(hartree_matrix(1)))
     
      // nwpw_xc_init(ion_nion(),npsp,
      //                 int_mb(nprj(1)),
      //                 int_mb(n1dbasis(1)),
      //                 int_mb(n1dgrid(1)),
      //                 int_mb(psp_type(1)),
      //                 int_mb(lmax(1)),
      //                 nmax_max*lmmax_max,
      //                 int_mb(l_projector(1)),
      //                 int_mb(m_projector(1)),
      //                 int_mb(b_projector(1)))
   }
}

 

/*******************************************
 *                                         *
 *     Pseudopotential::v_nonlocal_orb     *
 *                                         *
 *******************************************/
void Pseudopotential::v_nonlocal_orb(double *orb, double *Horb)
{
   nwpw_timing_function ftimer(6);
   bool done;
   int ii, ia, l, nshift0, sd_function, i;
   int jj, ll, jstart, jend, nprjall;
   double *exi;
   double *prjtmp, *sw1, *sw2, *prj, *vnlprj;
   Parallel *parall;
   double omega = mypneb->lattice->omega();
   // double scal = 1.0/lattice_omega();
   double scal = 1.0 / omega;
   int one = 1;
   int ntmp, nshift, nn;
   double rone = 1.0;
   double rmone = -1.0;

   //nn = mypneb->neq[0] + mypneb->neq[1];
   nn = 1;
   nshift0 = mypneb->npack(1);
   nshift = 2 * mypneb->npack(1);
   exi = new (std::nothrow) double[nshift]();
   prjtmp = new (std::nothrow) double[nprj_max * nshift]();
   sw1 = new (std::nothrow) double[nn * nprj_max]();
   sw2 = new (std::nothrow) double[nn * nprj_max]();

   parall = mypneb->d3db::parall;

   ii = 0;
   while (ii < (myion->nion))
   {
      ia = myion->katm[ii];
      nprjall = 0;
      jstart = ii;
      done = false;
      while (!done)
      {
         // generate projectors
         if (nprj[ia] > 0)
         {
            mystrfac->strfac_pack(1, ii, exi);
            for (l = 0; l < nprj[ia]; ++l)
            {
               sd_function = !(l_projector[ia][l] & 1);
               prj = prjtmp + ((l+nprjall)*nshift);
               vnlprj = vnl[ia] + (l*nshift0);
               if (sd_function)
                  mypneb->tcc_pack_Mul(1, vnlprj, exi, prj);
               else
                  mypneb->tcc_pack_iMul(1, vnlprj, exi, prj);
            }
            nprjall += nprj[ia];
         }
         ++ii;
         if (ii < (myion->nion))
         {
            ia = myion->katm[ii];
            done = ((nprjall + nprj[ia]) > nprj_max);
         }
         else
         {
            done = true;
         }
      }
      jend = ii;
      mypneb->cc_pack_inprjdot(1, nn, nprjall, orb, prjtmp, sw1);

   
      /* sw2 = Gijl*sw1 */
      ll = 0;
      for (jj = jstart; jj < jend; ++jj) {
        ia = myion->katm[jj];
        if (nprj[ia] > 0) {
          Multiply_Gijl_sw1(nn, nprj[ia], nmax[ia], lmax[ia], n_projector[ia],
                            l_projector[ia], m_projector[ia], Gijl[ia],
                            sw1+(ll*nn), sw2+(ll*nn));
          ll += nprj[ia];
        }
      }
   
      ntmp = nn * nprjall;
      DSCAL_PWDFT(ntmp, scal, sw2, one);
   
      mypneb->d3db::mygdevice.NT_dgemm(nshift, nn, nprjall, rmone, prjtmp, sw2, rone, Horb);
   }

   delete[] sw2;
   delete[] sw1;
   delete[] prjtmp;
   delete[] exi;
}




/*******************************************
 *                                         *
 *     Pseudopotential::v_nonlocal         *
 *                                         *
 *******************************************/
/**
 * @brief Calculate the nonlocal contribution to the Hamiltonian for the Pseudopotential class.
 *
 * This function computes the nonlocal contribution to the Hamiltonian matrix for electronic structure calculations.
 *
 * @param psi Pointer to the wavefunction.
 * @param Hpsi Pointer to the result Hamiltonian*wavefunction.
 *
 * @note This function performs matrix multiplications and calculations related to nonlocal contributions to the Hamiltonian.
 */
void Pseudopotential::v_nonlocal(double *psi, double *Hpsi) 
{
   nwpw_timing_function ftimer(6);
   bool done;
   int ii, ia, l, nshift0, sd_function, i;
   int jj, ll, jstart, jend, nprjall;
   double *exi;
   double *prjtmp, *sw1, *sw2, *prj, *vnlprj;
   Parallel *parall;
   double omega = mypneb->lattice->omega();
   // double scal = 1.0/lattice_omega();
   double scal = 1.0 / omega;
   int one = 1;
   int ntmp, nshift, nn;
   double rone = 1.0;
   double rmone = -1.0;
 
   nn = mypneb->neq[0] + mypneb->neq[1];
   nshift0 = mypneb->npack(1);
   nshift = 2 * mypneb->npack(1);
   exi = new (std::nothrow) double[nshift]();
   prjtmp = new (std::nothrow) double[nprj_max * nshift]();
   sw1 = new (std::nothrow) double[nn * nprj_max]();
   sw2 = new (std::nothrow) double[nn * nprj_max]();
 
   parall = mypneb->d3db::parall;

#if 1

   // Copy psi to device
   mypneb->d3db::mygdevice.psi_copy_host2gpu(nshift0, nn, psi);
   mypneb->d3db::mygdevice.hpsi_copy_host2gpu(nshift0, nn, Hpsi);
 
   ii = 0;
   while (ii < (myion->nion)) 
   {
      ia = myion->katm[ii];
      nprjall = 0;
      jstart = ii;
      done = false;
      while (!done) 
      {
         // generate projectors
         if (nprj[ia] > 0) 
         {
            mystrfac->strfac_pack(1, ii, exi);
            for (l = 0; l < nprj[ia]; ++l) 
            {
               sd_function = !(l_projector[ia][l] & 1);
               prj = prjtmp + ((l+nprjall)*nshift);
               vnlprj = vnl[ia] + (l*nshift0);
               if (sd_function)
                  mypneb->tcc_pack_Mul(1, vnlprj, exi, prj);
               else
                  mypneb->tcc_pack_iMul(1, vnlprj, exi, prj);
            }
            nprjall += nprj[ia];
         }
         ++ii;
         if (ii < (myion->nion)) 
         {
            ia = myion->katm[ii];
            done = ((nprjall + nprj[ia]) > nprj_max);
         } 
         else 
         {
            done = true;
         }
      }
      jend = ii;
      mypneb->cc_pack_inprjdot(1, nn, nprjall, psi, prjtmp, sw1);
      //parall->Vector_SumAll(1, nn*nprjall, sw1);
     
      /* sw2 = Gijl*sw1 */
      ll = 0;
      for (jj = jstart; jj < jend; ++jj) {
        ia = myion->katm[jj];
        if (nprj[ia] > 0) {
          Multiply_Gijl_sw1(nn, nprj[ia], nmax[ia], lmax[ia], n_projector[ia],
                            l_projector[ia], m_projector[ia], Gijl[ia],
                            sw1+(ll*nn), sw2+(ll*nn));
          ll += nprj[ia];
        }
      }
     
      ntmp = nn * nprjall;
      DSCAL_PWDFT(ntmp, scal, sw2, one);

      //scale (sw2) psp by factor - used for generating antiferromagnetic structures 
      //pwdft input: pspspin up/down scale l ion_numbers
      if (pspspin) apply_pspspin_scaling(sw2,nn,jstart,jend,ii);


      // DGEMM_PWDFT((char*) "N",(char*) "T",nshift,nn,nprjall,
      //                rmone,
      //                prjtmp,nshift,
      //                sw2,   nn,
      //                rone,
      //                Hpsi,nshift);
      mypneb->d3db::mygdevice.NT_dgemm(nshift, nn, nprjall, rmone, prjtmp, sw2, rone, Hpsi);
   }
   mypneb->d3db::mygdevice.hpsi_copy_gpu2host(nshift0, nn, Hpsi);
#else

  for (ii = 0; ii < (myion->nion); ++ii) {
    ia = myion->katm[ii];
    if (nprj[ia] > 0) {
      /* structure factor */
      mystrfac->strfac_pack(1, ii, exi);

      /* generate sw1's and projectors */
      for (l = 0; l < nprj[ia]; ++l) {
        sd_function = !(l_projector[ia][l] & 1);
        prj = &(prjtmp[l * nshift]);
        vnlprj = &(vnl[ia][l * nshift0]);
        if (sd_function)
          mypneb->tcc_pack_Mul(1, vnlprj, exi, prj);
        else
          mypneb->tcc_pack_iMul(1, vnlprj, exi, prj);
        // mypneb->cc_pack_indot(1,nn,psi,prj,&(sw1[l*nn]));
      }
      ntmp = nprj[ia];
      mypneb->cc_pack_inprjdot(1, nn, ntmp, psi, prjtmp, sw1);
      parall->Vector_SumAll(1, nn * nprj[ia], sw1);

      /* sw2 = Gijl*sw1 */
      Multiply_Gijl_sw1(nn, nprj[ia], nmax[ia], lmax[ia], n_projector[ia],
                        l_projector[ia], m_projector[ia], Gijl[ia], sw1, sw2);

      /* do Kleinman-Bylander Multiplication */
      ntmp = nn * nprj[ia];
      DSCAL_PWDFT(ntmp, scal, sw2, one);

      ntmp = nprj[ia];

      // DGEMM_PWDFT((char*) "N",(char*) "T",nshift,nn,ntmp,
      //       rmone,
      //       prjtmp,nshift,
      //       sw2,   nn,
      //       rone,
      //       Hpsi,nshift);
      mypneb->d3db::mygdevice.NT_dgemm(nshift, nn, ntmp, rmone, prjtmp, sw2, rone, Hpsi);

    } /*if nprj>0*/
  }   /*ii*/
#endif

   delete[] sw2;
   delete[] sw1;
   delete[] prjtmp;
   delete[] exi;
}


/*******************************************
 *                                         *
 *     Pseudopotential::v_nonlocal_fion    *
 *                                         *
 *******************************************/
/**
 * @brief Calculate the nonlocal contribution to the Hamiltonian for the Pseudopotential class with additional features.
 *
 * This function computes the nonlocal contribution to the Hamiltonian matrix for electronic structure calculations.
 * It also supports additional features such as handling movement (move) and calculating forces (fion).
 *
 * @param psi Pointer to the wavefunction.
 * @param Hpsi Pointer to the result Hamiltonian*wavefunction.
 * @param move Flag indicating whether to perform additional movement-related calculations.
 * @param fion Pointer to the array to store forces (if move is true).
 *
 * @note This function performs matrix multiplications and calculations related to nonlocal contributions to the Hamiltonian.
 * It optionally computes forces and supports movement-related calculations based on the 'move' flag.
 */
void Pseudopotential::v_nonlocal_fion(double *psi, double *Hpsi,
                                      const bool move, double *fion, double *occ) 
{
   nwpw_timing_function ftimer(6);
   bool done;
   int ii, ia, l, nshift0, sd_function, i, n;
   int jj, ll, jstart, jend, nprjall;
 
   double *exi;
   double *prjtmp, *sw1, *sw2, *prj, *vnlprj;
   double *Gx, *Gy, *Gz, *xtmp, *sum;
   double ff[3];
   Parallel *parall;
   double omega = mypneb->lattice->omega();
   // double scal = 1.0/lattice_omega();
   double scal = 1.0 / omega;
   int one = 1;
   int three = 3;
   int ntmp, nshift, nn, ispin;
 
   double rone = 1.0;
   double rmone = -1.0;
 
   nn = mypneb->neq[0] + mypneb->neq[1];
   ispin = mypneb->ispin;
   nshift0 = mypneb->npack(1);
   nshift = 2 * mypneb->npack(1);
   exi = new (std::nothrow) double[nshift]();
   prjtmp = new (std::nothrow) double[nprj_max*nshift]();
   sw1 = new (std::nothrow) double[nn*nprj_max]();
   sw2 = new (std::nothrow) double[nn*nprj_max]();
 
   // Copy psi to device
   mypneb->d3db::mygdevice.psi_copy_host2gpu(nshift0, nn, psi);
   mypneb->d3db::mygdevice.hpsi_copy_host2gpu(nshift0, nn, Hpsi);
 
   if (move) 
   {
      xtmp = new (std::nothrow) double[nshift0]();
      sum = new (std::nothrow) double[3*nn*nprj_max]();
     
      Gx = mypneb->Gpackxyz(1, 0);
      Gy = mypneb->Gpackxyz(1, 1);
      Gz = mypneb->Gpackxyz(1, 2);
   }
 
   parall = mypneb->d3db::parall;

   ii = 0;
   while (ii < (myion->nion)) 
   {
      ia = myion->katm[ii];
      nprjall = 0;
      jstart = ii;
      done = false;
      while (!done) 
      {
         // generate projectors
         if (nprj[ia] > 0) 
         {
            mystrfac->strfac_pack(1, ii, exi);
            for (l=0; l<nprj[ia]; ++l) 
            {
               sd_function = !(l_projector[ia][l] & 1);
               prj = &(prjtmp[(l + nprjall) * nshift]);
               vnlprj = &(vnl[ia][l * nshift0]);
               if (sd_function)
                  mypneb->tcc_pack_Mul(1, vnlprj, exi, prj);
               else
                  mypneb->tcc_pack_iMul(1, vnlprj, exi, prj);
            }
            nprjall += nprj[ia];
         }
         ++ii;
         if (ii < (myion->nion)) 
         {
            ia = myion->katm[ii];
            done = ((nprjall + nprj[ia]) > nprj_max);
         } 
         else 
         {
            done = true;
         }
      }
      jend = ii;
     
    
      mypneb->cc_pack_inprjdot(1, nn, nprjall, psi, prjtmp, sw1);

      if (move)
      {
         nwpw_timing_function f2timer(61);
         mypneb->n2ccttt_pack_i3ndot(1,nn,nprjall,psi,prjtmp,Gx,Gy,Gz,sum);
         if (occ)
         {
            int count3 = 0;
            for (auto l=0; l<nprjall; ++l)
            for (auto ms=0; ms<ispin; ++ms)
            for (auto q=0; q<mypneb->neq[ms]; ++q)
            {
                double wf = occ[mypneb->msntoindex(ms,q)];
                sum[count3]   *= wf;
                sum[count3+1] *= wf;
                sum[count3+2] *= wf;
                count3 += 3;
            }

         }
      }
     
      /* sw2 = Gijl*sw1 */
      ll = 0;
      for (jj=jstart; jj<jend; ++jj) 
      {
         ia = myion->katm[jj];
         if (nprj[ia] > 0) 
         {
            Multiply_Gijl_sw1(nn, nprj[ia], nmax[ia], lmax[ia], n_projector[ia],
                              l_projector[ia], m_projector[ia], Gijl[ia],
                              sw1 + ll*nn, sw2 + ll*nn);
            ll += nprj[ia];
         }
      }

     
      ntmp = nn*nprjall;
      DSCAL_PWDFT(ntmp, scal, sw2, one);

      //scale (sw2) psp by factor - used for generating antiferromagnetic structures 
      //pwdft input: pspspin up/down scale l ion_numbers
      if (pspspin) apply_pspspin_scaling(sw2,nn,jstart,jend,ii);

     
      mypneb->d3db::mygdevice.NT_dgemm(nshift, nn, nprjall, rmone, prjtmp, sw2, rone, Hpsi);
     
      if (move) 
      {
         nwpw_timing_function f2timer(61);
         // for (ll=0; ll<nprjall; ++ll)
         ll = 0;
         for (jj = jstart; jj < jend; ++jj) 
         {
            ia = myion->katm[jj];
            for (l=0; l<nprj[ia]; ++l) 
            {
               ff[0] = 2.0*DDOT_PWDFT(nn, sw2 + ll*nn, one, sum + 3*nn*ll,     three);
               ff[1] = 2.0*DDOT_PWDFT(nn, sw2 + ll*nn, one, sum + 3*nn*ll + 1, three);
               ff[2] = 2.0*DDOT_PWDFT(nn, sw2 + ll*nn, one, sum + 3*nn*ll + 2, three);
               parall->Vector_SumAll(2,3,ff);

               fion[3*jj]   += (3-ispin)*ff[0];
               fion[3*jj+1] += (3-ispin)*ff[1];
               fion[3*jj+2] += (3-ispin)*ff[2];
               ++ll;
            }
         }
      }
   }
 
   mypneb->d3db::mygdevice.hpsi_copy_gpu2host(nshift0, nn, Hpsi);


   if (move) 
   {
      delete[] xtmp;
      delete[] sum;
   }
   delete[] sw2;
   delete[] sw1;
   delete[] prjtmp;
   delete[] exi;
}






/*******************************************
 *                                         *
 *     Pseudopotential::f_nonlocal_fion    *
 *                                         *
 *******************************************/
/**
 * @brief Calculate the nonlocal forces for the Pseudopotential class with additional features.
 *
 * This function computes the nonlocal forces for electronic structure calculations with additional features.
 * It is used to calculate forces on ions based on the wavefunction.
 *
 * @param psi Pointer to the wavefunction.
 * @param fion Pointer to the array to store forces.
 *
 * @note This function calculates nonlocal forces by performing matrix multiplications and vector operations.
 * The resulting forces are stored in the 'fion' array.
 */
void Pseudopotential::f_nonlocal_fion(double *psi, double *fion, double *occ) 
{
   nwpw_timing_function ftimer(6);
   bool done;
   int ii, ia, l, nshift0, sd_function, i, n;
   int jj, ll, jstart, jend, nprjall;
 
   double *exi;
   double *prjtmp, *sw1, *sw2, *prj, *vnlprj;
   double *Gx, *Gy, *Gz, *xtmp, *sum;
   double ff[3];
   Parallel *parall;
   double omega = mypneb->lattice->omega();
   // double scal = 1.0/lattice_omega();
   double scal = 1.0 / omega;
   int one = 1;
   int three = 3;
   int ntmp, nshift, nn, ispin;
 
   double rone = 1.0;
   double rmone = -1.0;
 
   nn = mypneb->neq[0] + mypneb->neq[1];
   ispin = mypneb->ispin;
   nshift0 = mypneb->npack(1);
   nshift = 2 * mypneb->npack(1);
   exi = new (std::nothrow) double[nshift]();
   prjtmp = new (std::nothrow) double[nprj_max * nshift]();
   sw1 = new (std::nothrow) double[nn * nprj_max]();
   sw2 = new (std::nothrow) double[nn * nprj_max]();
 
   // Copy psi to device
   mypneb->d3db::mygdevice.psi_copy_host2gpu(nshift0, nn, psi);
   // mypneb->d3db::mygdevice.hpsi_copy_host2gpu(nshift0,nn,Hpsi);
 
   xtmp = new (std::nothrow) double[nshift0]();
   sum  = new (std::nothrow) double[3*nn*nprj_max]();
   // Gx = new (std::nothrow) double [mypneb->nfft3d]();
   // Gy = new (std::nothrow) double [mypneb->nfft3d]();
   // Gz = new (std::nothrow) double [mypneb->nfft3d]();
   // mypneb->tt_copy(mypneb->Gxyz(0),Gx);
   // mypneb->tt_copy(mypneb->Gxyz(1),Gy);
   // mypneb->tt_copy(mypneb->Gxyz(2),Gz);
   // mypneb->t_pack(1,Gx);
   // mypneb->t_pack(1,Gy);
   // mypneb->t_pack(1,Gz);
   Gx = mypneb->Gpackxyz(1, 0);
   Gy = mypneb->Gpackxyz(1, 1);
   Gz = mypneb->Gpackxyz(1, 2);
 
   parall = mypneb->d3db::parall;
 
   ii = 0;
   while (ii < (myion->nion)) 
   {
      ia = myion->katm[ii];
      nprjall = 0;
      jstart = ii;
      done = false;
      while (!done) 
      {
         // generate projectors
         if (nprj[ia] > 0) 
         {
            mystrfac->strfac_pack(1, ii, exi);
            for (l = 0; l < nprj[ia]; ++l) 
            {
               sd_function = !(l_projector[ia][l] & 1);
               prj = prjtmp + ((l+nprjall)*nshift);
               vnlprj = vnl[ia] + (l*nshift0);
               if (sd_function)
                  mypneb->tcc_pack_Mul(1, vnlprj, exi, prj);
               else
                  mypneb->tcc_pack_iMul(1, vnlprj, exi, prj);
              
               //mypneb->ch3t_pack_i3ndot(1,prj,nn,psi,Gx,Gy,Gz,sum+3*nn(l+nprjall));
               /*
               for (n = 0; n < nn; ++n) 
               {
                  mypneb->cct_pack_iconjgMul(1, prj, psi + n*nshift, xtmp);
                  sum[3*n +     3*nn*(l+nprjall)] = mypneb->tt_pack_idot(1, Gx, xtmp);
                  sum[3*n + 1 + 3*nn*(l+nprjall)] = mypneb->tt_pack_idot(1, Gy, xtmp);
                  sum[3*n + 2 + 3*nn*(l+nprjall)] = mypneb->tt_pack_idot(1, Gz, xtmp);
               }
               */
            }
            nprjall += nprj[ia];
         }
         ++ii;
         if (ii < (myion->nion)) 
         {
            ia = myion->katm[ii];
            done = ((nprjall + nprj[ia]) > nprj_max);
         } 
         else 
         {
            done = true;
         }
      }
      jend = ii;
      mypneb->cc_pack_inprjdot(1, nn, nprjall, psi, prjtmp, sw1);
      //parall->Vector_SumAll(1, nn*nprjall, sw1);

      mypneb->n2ccttt_pack_i3ndot(1,nn,nprjall,psi,prjtmp,Gx,Gy,Gz,sum);
      //parall->Vector_SumAll(1, 3*nn*nprjall, sum);
      if (occ)
      {
         int count3 = 0;
         for (auto l=0; l<nprjall; ++l)
         for (auto ms=0; ms<ispin; ++ms)
         for (auto q=0; q<mypneb->neq[ms]; ++q)
         {
             double wf = occ[mypneb->msntoindex(ms,q)];
             sum[count3]   *= wf;
             sum[count3+1] *= wf;
             sum[count3+2] *= wf;
             count3 += 3;
         }
      }
     
      /* sw2 = Gijl*sw1 */
      ll = 0;
      for (jj = jstart; jj < jend; ++jj) 
      {
         ia = myion->katm[jj];
         if (nprj[ia] > 0) 
         {
            Multiply_Gijl_sw1(nn, nprj[ia], nmax[ia], lmax[ia], n_projector[ia],
                              l_projector[ia], m_projector[ia], Gijl[ia],
                              sw1+(ll*nn),sw2+(ll*nn));
            ll += nprj[ia];
         }
      }
     
      ntmp = nn * nprjall;
      DSCAL_PWDFT(ntmp, scal, sw2, one);
     
      mypneb->d3db::mygdevice.T_free();


      //scale (sw2) psp by factor - used for generating antiferromagnetic structures 
      //pwdft input: pspspin up/down scale l ion_numbers
      if (pspspin) apply_pspspin_scaling(sw2,nn,jstart,jend,ii);


      // for (ll=0; ll<nprjall; ++ll)
      ll = 0;
      for (jj = jstart; jj < jend; ++jj) 
      {
         ia = myion->katm[jj];
         for (l = 0; l < nprj[ia]; ++l) 
         {
            //fion[3*jj]   += (3-ispin)*2.0*DDOT_PWDFT(nn, sw2+ll*nn, one, sum+(3*nn*ll),   three);
            //fion[3*jj+1] += (3-ispin)*2.0*DDOT_PWDFT(nn, sw2+ll*nn, one, sum+(3*nn*ll+1), three);
            //fion[3*jj+2] += (3-ispin)*2.0*DDOT_PWDFT(nn, sw2+ll*nn, one, sum+(3*nn*ll+2), three);
            ff[0] = 2.0*DDOT_PWDFT(nn, sw2 + ll*nn, one, sum + 3*nn*ll,     three);
            ff[1] = 2.0*DDOT_PWDFT(nn, sw2 + ll*nn, one, sum + 3*nn*ll + 1, three);
            ff[2] = 2.0*DDOT_PWDFT(nn, sw2 + ll*nn, one, sum + 3*nn*ll + 2, three);
            parall->Vector_SumAll(2,3,ff);

            fion[3*jj]   += (3-ispin)*ff[0];
            fion[3*jj+1] += (3-ispin)*ff[1];
            fion[3*jj+2] += (3-ispin)*ff[2];

            ++ll;
         }
      }
   }
   // mypneb->d3db::mygdevice.hpsi_copy_gpu2host(nshift0,nn,Hpsi);
 
   delete[] xtmp;
   delete[] sum;
   // delete [] Gx;
   // delete [] Gy;
   // delete [] Gz;
 
   delete[] sw2;
   delete[] sw1;
   delete[] prjtmp;
   delete[] exi;
}


/*******************************************
 *                                         *
 *     Pseudopotential::e_nonlocal         *
 *                                         *
 *******************************************/
/**
 * @brief Calculate the nonlocal energy contribution for the Pseudopotential class.
 *
 * This function computes the nonlocal energy contribution in the Pseudopotential class.
 * It is used to calculate the energy associated with nonlocal potentials based on the wavefunction.
 *
 * @param psi Pointer to the wavefunction.
 *
 * @return The nonlocal energy contribution.
 *
 * @note This function calculates the energy contribution by performing matrix multiplications and vector operations.
 */
double Pseudopotential::e_nonlocal(double *psi, double *occ) 
{
   nwpw_timing_function ftimer(6);

   Parallel *parall = mypneb->d3db::parall;
   double omega = mypneb->lattice->omega();
   double scal = 1.0 / omega;
   int one = 1;
   double rone = 1.0;
   double rmone = -1.0;
   double esum = 0.0;
 
   int nn = mypneb->neq[0] + mypneb->neq[1];
   int nshift0 = mypneb->npack(1);
   int nshift = 2 * mypneb->npack(1);
   double *exi    = new (std::nothrow) double[nshift]();
   double *prjtmp = new (std::nothrow) double[nprj_max * nshift]();
   double *sw1    = new (std::nothrow) double[nn * nprj_max]();
   double *sw2    = new (std::nothrow) double[nn * nprj_max]();
 
   // Copy psi to device
   mypneb->d3db::mygdevice.psi_copy_host2gpu(nshift0, nn, psi);
 
   auto ii = 0;
   while (ii<(myion->nion)) 
   {
      auto ia = myion->katm[ii];
      auto nprjall = 0;
      auto jstart = ii;
      auto done = false;
      while (!done) 
      {
         // generate projectors
         if (nprj[ia] > 0) 
         {
            mystrfac->strfac_pack(1, ii, exi);
            for (auto l=0; l<nprj[ia]; ++l) 
            {
               auto sd_function = !(l_projector[ia][l] & 1);
               auto prj    = prjtmp+(l+nprjall)*nshift;
               auto vnlprj = vnl[ia]+l*nshift0;
               if (sd_function)
                  mypneb->tcc_pack_Mul(1, vnlprj, exi, prj);
               else
                  mypneb->tcc_pack_iMul(1, vnlprj, exi, prj);
            }
            nprjall += nprj[ia];
         }
         ++ii;
         if (ii<(myion->nion)) 
         {
            ia = myion->katm[ii];
            done = ((nprjall + nprj[ia]) > nprj_max);
         } 
         else 
         {
            done = true;
         }
      }
      auto jend = ii;
      mypneb->cc_pack_inprjdot(1, nn, nprjall, psi, prjtmp, sw1);
      //parall->Vector_SumAll(1, nn * nprjall, sw1);

      /* sw2 = Gijl*sw1 */
      auto ll = 0;
      for (auto jj=jstart; jj<jend; ++jj) 
      {
         ia = myion->katm[jj];
         if (nprj[ia] > 0) 
         {
            Multiply_Gijl_sw1(nn, nprj[ia], nmax[ia], lmax[ia], n_projector[ia],
                              l_projector[ia], m_projector[ia], Gijl[ia],
                              sw1+ll*nn, sw2+ll*nn);
            ll += nprj[ia];
         }
      }
     
      auto ntmp = nn*nprjall;
      DSCAL_PWDFT(ntmp, scal, sw2, one);

      if (occ)
      {
         int count = 0;
         for (auto p=0; p<nprjall; ++p)
         {
            for (auto ms=0; ms<mypneb->ispin; ++ms)
            for (auto q=0; q<mypneb->neq[ms]; ++q)
            {
               sw1[count] *= occ[mypneb->msntoindex(ms,q)];
               ++count;
            }
         }
      }
      esum += DDOT_PWDFT(ntmp, sw1, one, sw2, one);
      mypneb->d3db::mygdevice.T_free();
   }


   esum = parall->SumAll(2,esum);

   if (mypneb->ispin == 1)
      esum *= 2.0;
 
   delete[] exi;
   delete[] prjtmp;
   delete[] sw1;
   delete[] sw2;
 
   return esum;
}


/*******************************************
 *                                         *
 *  Pseudopotential::apply_pspspin_scaling *
 *                                         *
 *******************************************/
/**
 * @brief Apply spin-dependent scaling factors to pseudopotential projectors.
 *
 * This function modifies the `sw2` array by applying scaling factors to
 * the pseudopotential projectors associated with specific angular momentum
 * channels (l, m) and atomic ion indices. It supports separate scaling for 
 * spin-up and spin-down components, as defined by the `pspspin_upscale` and 
 * `pspspin_downscale` arrays, respectively.
 *
 * Matching is performed based on the orbital angular momentum and magnetic 
 * quantum numbers (`l_projector`, `m_projector`) and ion index masks 
 * (`pspspin_upions`, `pspspin_downions`). The optional wildcard value `999` 
 * for `m` enables applying scaling across all m values for a given l.
 *
 * The function is typically used to construct antiferromagnetic or 
 * ferrimagnetic configurations by applying user-defined spin asymmetries
 * to pseudopotential projectors during initialization.
 *
 * @param sw2    Pointer to the projector array to be scaled (e.g., real- or reciprocal-space).
 *               This array is expected to be indexed as sw2[n + (l + ll) * nn].
 * @param nn     Leading dimension or stride between (l, m) channels (typically grid points or plane-waves).
 * @param jstart Start index (inclusive) of local ions assigned to this process.
 * @param jend   End index (exclusive) of local ions assigned to this process.
 * @param ii     Index into the user-defined `pspspin_*` arrays specifying which scaling entry to apply.
 */
void Pseudopotential::apply_pspspin_scaling(double* sw2, int nn, int jstart, int jend, int ii)
{
   int nup = mypneb->neq[0];
   int ndn = mypneb->neq[1];
   int ll = 0;

   for (int jj=jstart; jj<jend; ++jj)
   {
      int ia = myion->katm[jj];
      if (nprj[ia] > 0)
      {
         for (int l=0; l<nprj[ia]; ++l)
         {
            int l_prj = l_projector[ia][l];
            int m_prj = m_projector[ia][l];

            if (pspspin_upions[ii] &&
                (l_prj == pspspin_upl[ii]) &&
                (m_prj == pspspin_upm[ii] || pspspin_upm[ii] > 999))
            {
               for (int n=0; n<nup; ++n)
                  sw2[n + (l+ll)*nn] *= pspspin_upscale[ii];
            }

            if (pspspin_downions[ii] &&
                (l_prj == pspspin_downl[ii]) &&
                (m_prj == pspspin_downm[ii] || pspspin_downm[ii] > 999))
            {
               for (int n=0; n<ndn; ++n)
                  sw2[n + nup + (l+ll)*nn] *= pspspin_downscale[ii];
            }
         }
         ll += nprj[ia];
      }
   }
}



/*******************************************
 *                                         *
 *       Pseudopotential::v_local          *
 *                                         *
 *******************************************/
/**
 * @brief Calculate the local potential energy and forces for the Pseudopotential class.
 *
 * This function computes the local potential energy and, if `move` is set to true, also calculates
 * the forces on the ions. It is used to calculate the energy and forces associated with the local
 * potential based on the electron density gradient and ion positions.
 *
 * @param vout Pointer to the array to store the local potential energy.
 * @param move A boolean flag indicating whether to calculate forces on ions (true) or only energy (false).
 * @param dng Pointer to the electron density gradient.
 * @param fion Pointer to the array to store the forces on ions (if `move` is true).
 *
 * @note This function performs matrix multiplications and vector operations to calculate the energy
 * and forces.
 */
void Pseudopotential::v_local(double *vout, const bool move, double *dng, double *fion) 
{
   nwpw_timing_function ftimer(5);
   int ii, ia, nshift, npack0;
   double *exi, *vtmp, *xtmp, *Gx, *Gy, *Gz;
 
   npack0 = mypneb->npack(0);
   xtmp = new (std::nothrow) double[npack0]();
   if (move) 
   {
      // xtmp = new (std::nothrow) double[npack0]();
     
      // Gx = new (std::nothrow) double [mypneb->nfft3d]();
      // Gy = new (std::nothrow) double [mypneb->nfft3d]();
      // Gz = new (std::nothrow) double [mypneb->nfft3d]();
      // mypneb->tt_copy(mypneb->Gxyz(0),Gx);
      // mypneb->tt_copy(mypneb->Gxyz(1),Gy);
      // mypneb->tt_copy(mypneb->Gxyz(2),Gz);
      // mypneb->t_pack(0,Gx);
      // mypneb->t_pack(0,Gy);
      // mypneb->t_pack(0,Gz);
     
      Gx = mypneb->Gpackxyz(0, 0);
      Gy = mypneb->Gpackxyz(0, 1);
      Gz = mypneb->Gpackxyz(0, 2);
   }
 
   mypneb->c_pack_zero(0, vout);
   nshift = 2 * npack0;
   exi = new (std::nothrow) double[nshift]();
   vtmp = new (std::nothrow) double[nshift]();
   for (ii = 0; ii < (myion->nion); ++ii) 
   {
      ia = myion->katm[ii];
      mystrfac->strfac_pack(0, ii, exi);
      // mypneb->tcc_pack_MulSum2(0,vl[ia],exi,vout);
      mypneb->tcc_pack_Mul(0, vl[ia], exi, vtmp);
      mypneb->cc_pack_Sum2(0, vtmp, vout);
     
      if (move) {
        // double xx =  mypneb->cc_pack_dot(0,dng,dng);
        // double yy =  mypneb->cc_pack_dot(0,vtmp,vtmp);
     
        mypneb->cct_pack_iconjgMulb(0, dng, vtmp, xtmp);
        // double zz =  mypneb->tt_pack_dot(0,xtmp,xtmp);
     
        fion[3 * ii] = mypneb->tt_pack_dot(0, Gx, xtmp);
        fion[3 * ii + 1] = mypneb->tt_pack_dot(0, Gy, xtmp);
        fion[3 * ii + 2] = mypneb->tt_pack_dot(0, Gz, xtmp);
      }
   }
   delete[] exi;
   delete[] vtmp;
   delete[] xtmp;
   // if (move)
   //{
   // delete [] xtmp;
   //}
}


/*******************************************
 *                                         *
 *     Pseudopotential::f_local            *
 *                                         *
 *******************************************/
/**
 * @brief Calculate the forces on ions due to the local potential gradient.
 *
 * This function computes the forces experienced by ions based on the gradient of the local
 * potential energy. It is used to calculate the forces on ions in the context of density
 * functional theory calculations.
 *
 * @param dng Pointer to the electron density gradient.
 * @param fion Pointer to the array to store the forces on ions.
 *
 * @note This function performs matrix multiplications and vector operations to calculate the forces.
 */
void Pseudopotential::f_local(double *dng, double *fion) 
{
   nwpw_timing_function ftimer(5);
   int ii, ia, nshift, npack0;
   double *exi, *vtmp, *xtmp, *Gx, *Gy, *Gz;
 
   npack0 = mypneb->npack(0);
 
   xtmp = new (std::nothrow) double[npack0]();
   // Gx = new (std::nothrow) double [mypneb->nfft3d]();
   // Gy = new (std::nothrow) double [mypneb->nfft3d]();
   // Gz = new (std::nothrow) double [mypneb->nfft3d]();
   // mypneb->tt_copy(mypneb->Gxyz(0),Gx);
   // mypneb->tt_copy(mypneb->Gxyz(1),Gy);
   // mypneb->tt_copy(mypneb->Gxyz(2),Gz);
   // mypneb->t_pack(0,Gx);
   // mypneb->t_pack(0,Gy);
   // mypneb->t_pack(0,Gz);
 
   Gx = mypneb->Gpackxyz(0, 0);
   Gy = mypneb->Gpackxyz(0, 1);
   Gz = mypneb->Gpackxyz(0, 2);
 
   nshift = 2 * npack0;
   exi = new (std::nothrow) double[nshift]();
   vtmp = new (std::nothrow) double[nshift]();
   for (ii = 0; ii < (myion->nion); ++ii) 
   {
      ia = myion->katm[ii];
      mystrfac->strfac_pack(0, ii, exi);
      // mypneb->tcc_pack_MulSum2(0,vl[ia],exi,vout);
      mypneb->tcc_pack_Mul(0, vl[ia], exi, vtmp);
     
      // double xx =  mypneb->cc_pack_dot(0,dng,dng);
      // double yy =  mypneb->cc_pack_dot(0,vtmp,vtmp);
     
      mypneb->cct_pack_iconjgMulb(0, dng, vtmp, xtmp);
      // double zz =  mypneb->tt_pack_dot(0,xtmp,xtmp);
     
      fion[3 * ii] = mypneb->tt_pack_dot(0, Gx, xtmp);
      fion[3 * ii + 1] = mypneb->tt_pack_dot(0, Gy, xtmp);
      fion[3 * ii + 2] = mypneb->tt_pack_dot(0, Gz, xtmp);
   }
   delete[] exi;
   delete[] vtmp;
 
   delete[] xtmp;
   // delete [] Gx;
   // delete [] Gy;
   // delete [] Gz;
}


/*******************************************
 *                                         *
 *     Pseudopotential::v_lr_local         *
 *                                         *
 *******************************************/
/*
      This routine calculates the long-range part of the
      local pseudopotential (used by version4)
*/
/**
 * @brief Calculate the long-range part of the local pseudopotential.
 *
 * This function calculates the long-range part of the local pseudopotential used in a specific
 * version of the pseudopotential calculation (version4). The long-range part is computed based on
 * the positions of ions and their charges.
 *
 * @param vlr_out Pointer to the array to store the long-range part of the local pseudopotential.
 *
 * @note This function uses error function calculations for efficiency when available.
 */
void Pseudopotential::v_lr_local(double *vlr_out) 
{
   nwpw_timing_function ftimer(5);
 
   int n2ft3d = mypneb->n2ft3d;
   int nion = myion->nion;
   double sqrt_pi = std::sqrt(4.0 * std::atan(1.0));
 
   std::memset(vlr_out, 0, n2ft3d * sizeof(double));
 
   if (mypneb->lattice->fast_erf()) 
   {
      // Error function parameters
      double xerf, yerf;
      double c1 = 0.07052307840;
      double c2 = 0.04228201230;
      double c3 = 0.00927052720;
      double c4 = 0.00015201430;
      double c5 = 0.00027656720;
      double c6 = 0.00004306380;
     
      for (auto ii = 0; ii < nion; ++ii) {
        int ia = myion->katm[ii];
        auto c = 1.0 / rlocal[ia];
        auto q = -zv[ia];
        auto xii = myion->rion(0, ii);
        auto yii = myion->rion(1, ii);
        auto zii = myion->rion(2, ii);
     
        for (auto i = 0; i < n2ft3d; ++i) {
          auto dx = mypneb->PGrid::r_grid[3 * i] - xii;
          auto dy = mypneb->PGrid::r_grid[3 * i + 1] - yii;
          auto dz = mypneb->PGrid::r_grid[3 * i + 2] - zii;
          auto r = std::sqrt(dx * dx + dy * dy + dz * dz);
          if (r > 1.0 - 15) {
            xerf = r * c;
            auto yerf1 =
                (1.0 +
                 xerf *
                     (c1 +
                      xerf *
                          (c2 +
                           xerf * (c3 + xerf * (c4 + xerf * (c5 + xerf * c6))))));
            auto yerf4 = yerf1 * yerf1 * yerf1 * yerf1;
            auto yerf16 = yerf4 * yerf4 * yerf4 * yerf4;
            yerf = (1.0 - 1.0 / yerf16);
            vlr_out[i] += (q / r) * yerf;
          } else
            vlr_out[i] += (2.0 * q * c / sqrt_pi);
        }
      }
   } 
   else 
   {
      for (auto ii = 0; ii < nion; ++ii) 
      {
         int ia = myion->katm[ii];
         auto c = 1.0 / rlocal[ia];
         auto q = -zv[ia];
         auto xii = myion->rion(0, ii);
         auto yii = myion->rion(1, ii);
         auto zii = myion->rion(2, ii);
        
         for (auto i = 0; i < n2ft3d; ++i) 
         {
            auto dx = mypneb->PGrid::r_grid[3 * i] - xii;
            auto dy = mypneb->PGrid::r_grid[3 * i + 1] - yii;
            auto dz = mypneb->PGrid::r_grid[3 * i + 2] - zii;
            auto r = std::sqrt(dx * dx + dy * dy + dz * dz);
            vlr_out[i] += (r > 1.0e-15) ? ((q / r) * std::erf(r * c))
                                        : (2.0 * q * c / sqrt_pi);
         }
      }
   }
}


/*******************************************
 *                                         *
 *    Pseudopotential::grad_v_lr_local     *
 *                                         *
 *******************************************/
/*
    This routine calculates the gradient of the long-range part of the
    local pseudopotential (used by version 4)
*/
/**
 * @brief Calculate the gradient of the long-range part of the local pseudopotential.
 *
 * This function computes the gradient of the long-range part of the local pseudopotential used
 * in a specific version of the pseudopotential calculation (version4). The gradient is calculated
 * based on the electron density `rho` and ion positions and charges.
 *
 * @param rho    Pointer to the electron density array.
 * @param fion   Pointer to the array to store the gradient of the long-range part of the
 *               pseudopotential.
 *
 * @note This function uses error function calculations for efficiency when available.
 */
void Pseudopotential::grad_v_lr_local(const double *rho, double *fion) 
{
   nwpw_timing_function ftimer(5);
 
   int n2ft3d = mypneb->n2ft3d;
   int nion = myion->nion;
   int nion3 = 3 * nion;
   int one = 1;
   double rone = 1.0;
   double sqrt_pi = std::sqrt(4.0 * std::atan(1.0));
   double omega = mypneb->lattice->omega();
   double scal1 = 1.0 / ((double)((mypneb->nx) * (mypneb->ny) * (mypneb->nz)));
   double dv = omega * scal1;
 
   double ftmp[nion3];
 
   if (mypneb->lattice->fast_erf()) {
     // Error function parameters
     double xerf, yerf;
     double c1 = 0.07052307840;
     double c2 = 0.04228201230;
     double c3 = 0.00927052720;
     double c4 = 0.00015201430;
     double c5 = 0.00027656720;
     double c6 = 0.00004306380;
 
     for (auto ii = 0; ii < nion; ++ii) {
       int ia = myion->katm[ii];
       auto c = 1.0 / rlocal[ia];
       auto q = -zv[ia];
       auto xii = myion->rion(0, ii);
       auto yii = myion->rion(1, ii);
       auto zii = myion->rion(2, ii);
 
       auto fx = 0.0;
       auto fy = 0.0;
       auto fz = 0.0;
       for (auto i = 0; i < n2ft3d; ++i) {
         auto rx = xii - mypneb->PGrid::r_grid[3 * i];
         auto ry = yii - mypneb->PGrid::r_grid[3 * i + 1];
         auto rz = zii - mypneb->PGrid::r_grid[3 * i + 2];
         auto r = std::sqrt(rx * rx + ry * ry + rz * rz);
         double v;
         if (r > 1.0e-8) {
           xerf = r * c;
           auto yerf1 =
               (1.0 +
                xerf *
                    (c1 +
                     xerf *
                         (c2 +
                          xerf * (c3 + xerf * (c4 + xerf * (c5 + xerf * c6))))));
           auto yerf4 = yerf1 * yerf1 * yerf1 * yerf1;
           auto yerf16 = yerf4 * yerf4 * yerf4 * yerf4;
           yerf = (1.0 - 1.0 / yerf16);
           v = (q *
                ((2.0 / sqrt_pi) * (r * c) * std::exp(-(r * c) * (r * c)) -
                 yerf) /
                (r * r * r));
         } else
           v = 0.0;
 
         fx += rho[i] * rx * v;
         fy += rho[i] * ry * v;
         fz += rho[i] * rz * v;
       }
       ftmp[3 * ii] = -fx * dv;
       ftmp[3 * ii + 1] = -fy * dv;
       ftmp[3 * ii + 2] = -fz * dv;
     }
   } 
   else 
   {
     for (auto ii = 0; ii < nion; ++ii) {
       int ia = myion->katm[ii];
       auto c = 1.0 / rlocal[ia];
       auto q = -zv[ia];
       auto xii = myion->rion(0, ii);
       auto yii = myion->rion(1, ii);
       auto zii = myion->rion(2, ii);
 
       auto fx = 0.0;
       auto fy = 0.0;
       auto fz = 0.0;
       for (auto i = 0; i < n2ft3d; ++i) {
         auto rx = xii - mypneb->PGrid::r_grid[3 * i];
         auto ry = yii - mypneb->PGrid::r_grid[3 * i + 1];
         auto rz = zii - mypneb->PGrid::r_grid[3 * i + 2];
         auto r = std::sqrt(rx * rx + ry * ry + rz * rz);
         auto v =
             (r > 1.0e-8)
                 ? (q *
                    ((2.0 / sqrt_pi) * (r * c) * std::exp(-(r * c) * (r * c)) -
                     std::erf(r * c)) /
                    (r * r * r))
                 : 0.0;
 
         fx += rho[i] * rx * v;
         fy += rho[i] * ry * v;
         fz += rho[i] * rz * v;
       }
       ftmp[3 * ii] = -fx * dv;
       ftmp[3 * ii + 1] = -fy * dv;
       ftmp[3 * ii + 2] = -fz * dv;
     }
   }
   mypneb->PGrid::parall->Vector_SumAll(1, 3 * nion, ftmp);
   DAXPY_PWDFT(nion3, rone, ftmp, one, fion, one);
}


/********************************************
 *                                          *
 * Pseudopotential::semicore_density_update *
 *                                          *
 ********************************************/
/**
 * @brief Update the semicore electron density.
 *
 * This function updates the electron density associated with semicore electrons
 * based on the current atomic positions and pseudopotential information. It calculates
 * the semicore density in real space, squares it, and scales it by a factor.
 *
 * @note This function assumes the availability of various class-level variables and functions
 *       for the pseudopotential and lattice calculations.
 */
void Pseudopotential::semicore_density_update() 
{
   int ii, ia;
   double omega = mypneb->lattice->omega();
   double scal2 = 1.0 / omega;
   // double scal1 = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
   // double dv    = omega*scal1;
   double *exi = mypneb->c_pack_allocate(0);
   double *tmp = mypneb->r_alloc();
 
   mypneb->r_zero(semicore_density);
   for (ii = 0; ii < (myion->nion); ++ii) {
     ia = myion->katm[ii];
     if (semicore[ia]) {
       mystrfac->strfac_pack(0, ii, exi);
       mypneb->tcc_pack_Mul(0, ncore_atom[ia], exi, tmp);
 
       /* Put put tmp into real space */
       mypneb->c_unpack(0, tmp);
       mypneb->cr_fft3d(tmp);
       mypneb->r_zero_ends(tmp);
 
       /*  square it  */
       mypneb->r_sqr(tmp);
       mypneb->rr_Sum(tmp, semicore_density);
     }
   }
   mypneb->r_SMul(scal2 * scal2, semicore_density);
 
   mypneb->r_dealloc(tmp);
   mypneb->c_pack_deallocate(exi);
}


/********************************************
 *                                          *
 *    Pseudopotential::semicore_xc_fion     *
 *                                          *
 ********************************************/
/* this routine need to be check!! */
/**
 * @brief Calculate the exchange-correlation contribution to the forces on semicore electrons.
 *
 * This function calculates the exchange-correlation contribution to the forces on semicore electrons
 * based on the electron density, the exchange-correlation potential, and atomic positions.
 *
 * @param vxc       Exchange-correlation potential.
 * @param fion      Forces on atoms, including semicore electrons.
 *
 * @note This function assumes the availability of various class-level variables and functions
 *       for pseudopotential and lattice calculations.
 * @note This routine needs to be checked for correctness.
 */
void Pseudopotential::semicore_xc_fion(double *vxc, double *fion) 
{
   int ia;
   int ispin = mypneb->ispin;
   int n2ft3d = mypneb->n2ft3d;
   double omega = mypneb->lattice->omega();
   double scal2 = 1.0 / omega;
   double scal1 = 1.0 / ((double)((mypneb->nx) * (mypneb->ny) * (mypneb->nz)));
 
   double *exi = mypneb->c_pack_allocate(0);
   double *tmp = mypneb->r_alloc();
   double *tmpx = mypneb->r_alloc();
   double *tmpy = mypneb->r_alloc();
   double *tmpz = mypneb->r_alloc();
   double *vxcG = mypneb->r_alloc();
 
   // double *Gx = new (std::nothrow) double [mypneb->nfft3d]();
   // double *Gy = new (std::nothrow) double [mypneb->nfft3d]();
   // double *Gz = new (std::nothrow) double [mypneb->nfft3d]();
   // mypneb->tt_copy(mypneb->Gxyz(0),Gx);
   // mypneb->tt_copy(mypneb->Gxyz(1),Gy);
   // mypneb->tt_copy(mypneb->Gxyz(2),Gz);
   // mypneb->t_pack(1,Gx);
   // mypneb->t_pack(1,Gy);
   // mypneb->t_pack(1,Gz);
   double *Gx = mypneb->Gpackxyz(0, 0);
   double *Gy = mypneb->Gpackxyz(0, 1);
   double *Gz = mypneb->Gpackxyz(0, 2);
 
   mypneb->rrr_Sum(vxc, vxc + (ispin-1)*n2ft3d, vxcG);
 
   for (int ii = 0; ii < (myion->nion); ++ii) {
     ia = myion->katm[ii];
     if (semicore[ia]) {
       /* put sqrt(core-density) at atom position */
       mystrfac->strfac_pack(0, ii, exi);
       mypneb->tcc_pack_Mul(0, ncore_atom[ia], exi, tmp);
       mypneb->tcc_pack_iMul(0, Gx, tmp, tmpx);
       mypneb->tcc_pack_iMul(0, Gy, tmp, tmpy);
       mypneb->tcc_pack_iMul(0, Gz, tmp, tmpz);
 
       /* Put put tmp,tmpx,tmpy,tmpz into real space */
       mypneb->c_unpack(0, tmp);
       mypneb->c_unpack(0, tmpx);
       mypneb->c_unpack(0, tmpy);
       mypneb->c_unpack(0, tmpz);
       mypneb->cr_fft3d(tmp);
       mypneb->cr_fft3d(tmpx);
       mypneb->cr_fft3d(tmpy);
       mypneb->cr_fft3d(tmpz);
       mypneb->r_zero_ends(tmpx);
       mypneb->r_zero_ends(tmpy);
       mypneb->r_zero_ends(tmpz);
 
       mypneb->rr_Mul(tmp, tmpx);
       mypneb->rr_Mul(tmp, tmpy);
       mypneb->rr_Mul(tmp, tmpz);
 
       fion[3 * ii] += scal1 * scal2 * mypneb->rr_dot(tmpx, vxcG);
       fion[3 * ii + 1] += scal1 * scal2 * mypneb->rr_dot(tmpy, vxcG);
       fion[3 * ii + 2] += scal1 * scal2 * mypneb->rr_dot(tmpz, vxcG);
     }
   }
 
   mypneb->r_dealloc(tmp);
   mypneb->r_dealloc(tmpx);
   mypneb->r_dealloc(tmpy);
   mypneb->r_dealloc(tmpz);
   mypneb->r_dealloc(vxcG);
   mypneb->c_pack_deallocate(exi);
   // delete [] Gx;
   // delete [] Gy;
   // delete [] Gz;
}


/*******************************************
 *                                         *
 *      Pseudopotential::print_pspall      *
 *                                         *
 *******************************************/
/**
 * @brief Generate a string representation of the pseudopotential details for all elements in the cluster.
 *
 * This function generates a string representation of the pseudopotential details for all elements
 * involved in the cluster, including atomic positions, pseudopotential type, loggrid parameters,
 * augmentation sphere radius, compensation sigma, total number of projectors, and basis set details.
 *
 * @return A string containing the pseudopotential details for all elements in the cluster.
 */
std::string Pseudopotential::print_pspall() 
{
   std::stringstream stream;
 
   std::ios init(NULL);
   init.copyfmt(stream);
 
   stream << std::endl << " elements involved in the cluster:" << std::endl;
   for (int ia = 0; ia < (myion->nkatm); ++ia) 
   {
      if (psp_type[ia] == 4) 
      {
         stream << std::setw(7) << (ia + 1) << ": " << std::left << std::setw(4)
                << myion->atom(ia) << "valence charge =" << std::right
                << std::fixed << std::setprecision(1) << std::setw(5) << zv[ia]
                << "  core charge =" << std::fixed << std::setprecision(1)
                << std::setw(6) << zion[ia] - zv[ia] << std::endl;
         stream.copyfmt(init);
         stream << "             comment = " << comment[ia] << std::endl
                << "             pseudopotential type           :" << std::setw(10)
                << psp_type[ia] << std::endl
                << "             loggrid parameter r0           :"
                << std::scientific << std::setw(10) << std::setprecision(3)
                << r1[ia] << std::endl
                << "             loggrid parameter rmax         :"
                << std::scientific << std::setw(10) << std::setprecision(3)
                << rmax[ia] << std::endl
                << "             loggrid parameter npoints      :" << std::setw(10)
                << n1dgrid[ia] << std::endl
                << "             augmentation sphere radius     :" << std::fixed
                << std::setw(10) << std::setprecision(3) << sphere_radius(ia)
                << " ( " << icut[ia] << " npoints " << icut[ia] << " per task)"
                << std::endl
                << "             compensation sigma             :" << std::fixed
                << std::setw(10) << std::setprecision(3) << sigma[ia] << std::endl
                << "             total number of projectors     :" << std::setw(10)
                << nprj[ia] << std::endl;
         if (psp_version == 4)
           stream << "             aperiodic cutoff radius         : "
                  << std::fixed << std::setprecision(3) << std::setw(6)
                  << rlocal[ia] << std::endl;
         stream << "             n_ps (n) l          eig    #projector"
                << std::endl;
         for (auto i = 0; i < n1dbasis[ia]; ++i)
           stream << std::setw(17) << nps[ia][i] << " (" << std::setw(1)
                  << nae[ia][i] << ") " << std::setw(1) << spdf_name(lps[ia][i])
                  << std::fixed << std::setw(13) << std::setprecision(6)
                  << eig[ia][i] << std::setw(14) << 2 * lps[ia][i] + 1
                  << std::endl;
      } 
      else 
      {
         stream << std::setw(7) << (ia + 1) << ": " << std::left << std::setw(4)
                << myion->atom(ia) << "valence charge =" << std::right
                << std::fixed << std::setprecision(1) << std::setw(5) << zv[ia]
                << "  lmax =" << std::setw(1) << lmax[ia] << std::endl;
         stream.copyfmt(init);
         stream << "             comment = " << comment[ia] << std::endl
                << "             pseudopotential type            =" << std::setw(3)
                << psp_type[ia] << std::endl
                << "             highest angular component       =" << std::setw(3)
                << lmax[ia] << std::endl
                << "             local potential used            =" << std::setw(3)
                << locp[ia] << std::endl
                << "             number of non-local projections =" << std::setw(3)
                << nprj[ia] << std::endl;
         if (psp_version == 4)
            stream << "             aperiodic cutoff radius         = "
                   << std::fixed << std::setprecision(3) << std::setw(6)
                   << rlocal[ia] << std::endl;
         if (semicore[ia])
            stream << "             semicore corrections inlcuded   = "
                   << std::fixed << std::setprecision(3) << std::setw(6)
                   << rcore[ia] << " (radius) " << ncore(ia) << " (charge)"
                   << std::endl;
         stream << "             cutoff = ";
         for (auto l = 0; l <= lmax[ia]; ++l)
            stream << std::fixed << std::setprecision(3) << std::setw(8) << rc[ia][l];
         stream << std::endl;
      }
   }
 
   return stream.str();
}

} // namespace pwdft
