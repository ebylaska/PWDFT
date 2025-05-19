/* Ions.C -
   Author - Eric Bylaska
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <sstream>


#include "json.hpp"
using json = nlohmann::json;

#include "iofmt.hpp"
#include "Control2.hpp"
#include "ion_bondings.hpp"
#include "symmetry_elements.hpp"
#include "Ion.hpp"


//extern "C" void nwpwxc_vdw3_dftd3_(char *, int *, int *, double *, double *, double *, double *, double *);

#define xyzstream(S, X, Y, Z, VX, VY, VZ)                                      \
  std::left << std::setw(3) << (S) << E124 << (X) << E124 << (Y) << E124       \
            << (Z) << E124 << (VX) << E124 << (VY) << E124 << (VZ)


namespace pwdft {

/*******************************
 *                             *
 *        incell1              *
 *                             *
 *******************************/
/**
 * @brief Wrap atomic positions into the unit cell using lattice vectors.
 *
 * This function adjusts Cartesian atomic positions to lie within the unit cell,
 * by converting to fractional coordinates, wrapping them into the [-0.5, 0.5)
 * range, and converting them back to Cartesian coordinates.
 *
 * @param nion Number of ions.
 * @param a    3x3 lattice matrix in row-major order (Cartesian lattice vectors).
 * @param r1   Array of atomic positions in Cartesian coordinates (modified in-place).
 */
static void incell1(const int nion, const double* unita, double* r1) 
{
   // Compute reciprocal lattice matrix b from lattice matrix a
   double ub[9];
   ub[0] = unita[4]*unita[8] - unita[5]*unita[7];
   ub[1] = unita[5]*unita[6] - unita[3]*unita[8];
   ub[2] = unita[3]*unita[7] - unita[4]*unita[6];
   ub[3] = unita[7]*unita[2] - unita[8]*unita[1];
   ub[4] = unita[8]*unita[0] - unita[6]*unita[2];
   ub[5] = unita[6]*unita[1] - unita[7]*unita[0];
   ub[6] = unita[1]*unita[5] - unita[2]*unita[4];
   ub[7] = unita[2]*unita[3] - unita[0]*unita[5];
   ub[8] = unita[0]*unita[4] - unita[1]*unita[3];
   double volume = unita[0]*ub[0] + unita[1]*ub[1] + unita[2]*ub[2];
   for (auto i=0; i<9; ++i)
      ub[i] /= volume;


   std::vector<double> sion1(3 * nion, 0.0);

   // Convert to fractional coordinates
   for (auto i=0; i<nion; ++i)
   {
      for (auto j=0; j<3; ++j)
      {
         for (auto k=0; k<3; ++k)
            sion1[3*i+j] += r1[3*i+k]*ub[3*k+j];
      }
   }

   // Wrap fractional coordinates to [-0.5, 0.5)
   for (auto i=0; i<nion; ++i) 
   {
      for (auto j=0; j<3; ++j) 
      {
         while (sion1[3*i+j] > 0.5)  sion1[3*i+j] -= 1.0;
         while (sion1[3*i+j] < -0.5) sion1[3*i+j] += 1.0;
      }
   }

   // Convert back to Cartesian coordinates
   for (auto i=0; i<nion; ++i)
   {
      for (auto j=0; j<3; ++j)
      {
         r1[3*i+j] = 0.0;
         for (auto k=0; k<3; ++k)
            r1[3*i+j] += unita[3*k+j]*sion1[3*i+k];
      }
   }
}


/*******************************
 *                             *
 *        solve_3by3           *
 *                             *
 *******************************/
/**
 * @brief Solves a 3x3 linear system of equations using a custom approach.
 *
 * Solves the system Im * omega = L, where:
 * - Im is a symmetric 3x3 moment of inertia tensor, given as a flat array.
 * - L is a 3D vector representing angular momentum.
 * - omega is the resulting angular velocity vector.
 *
 * This implementation assumes that Im is symmetric and uses explicit algebraic
 * expressions to compute the components of omega, avoiding full matrix inversion.
 *
 * @param[in] Im   The 3x3 symmetric matrix (flattened to length 9).
 * @param[in] L    The right-hand side vector (length 3).
 * @param[out] omega The solution vector (length 3), containing angular velocity.
 */
static void solve_3by3(const double Im[], const double L[], double omega[]) {
  double a = Im[0];
  double b = Im[4];
  double c = Im[8];
  double d = Im[3];
  double e = Im[6];
  double f = Im[7];
  double o = L[0];
  double p = L[1];
  double q = L[2];

  double af_de = a * f - d * e;
  double aq_eo = a * q - e * o;
  double ab_dd = a * b - d * d;
  double ac_ee = a * c - e * e;

  double z = (af_de * (a * p - d * o) - ab_dd * aq_eo) /
             (af_de * af_de - ab_dd * ac_ee);
  double y = (aq_eo - z * ac_ee) / af_de;
  double x = (o - d * y - e * z) / a;

  omega[0] = x;
  omega[1] = y;
  omega[2] = z;
}

/*******************************
 *                             *
 *        center_v_mass        *
 *                             *
 *******************************/
/**
 * @brief Computes the center of mass position for a system of ions.
 *
 * Calculates the weighted average of atomic positions based on their masses,
 * resulting in the center of mass in the x, y, and z directions.
 *
 * @param[in]  nion   Number of ions (atoms).
 * @param[in]  mass   Array of atomic masses (length nion).
 * @param[in]  rion0  Cartesian coordinates of atoms (length 3*nion).
 * @param[out] vx     X-coordinate of the center of mass.
 * @param[out] vy     Y-coordinate of the center of mass.
 * @param[out] vz     Z-coordinate of the center of mass.
 */
static void center_v_mass(int nion, double *mass, double *rion0, double *vx,
                          double *vy, double *vz) {
  double tmass = 0.0;
  double sx = 0.0;
  double sy = 0.0;
  double sz = 0.0;
  for (auto ii = 0; ii < nion; ++ii) {
    tmass += mass[ii];
    sx += mass[ii] * rion0[3 * ii];
    sy += mass[ii] * rion0[3 * ii + 1];
    sz += mass[ii] * rion0[3 * ii + 2];
  }
  *vx = sx / tmass;
  *vy = sy / tmass;
  *vz = sz / tmass;
}

/*******************************
 *                             *
 *        ion_find_nkatm       *
 *                             *
 *******************************/
// static void ion_tmpsymbols(RTDB& myrtdb, const int nion, char *symb, char
// *symtmp)
//{
//    int j,ii,i,matype,nelem;
//    char date[50];
//
//    if (!myrtdb.get_info("geometry:geometry:tags",&matype,&nelem,date))
//    { printf("rtdb.get_info error: geometry:geometry:tags not found\n");
//    exit(99);}
//
//    else
//    {
//       if (!myrtdb.get("geometry:geometry:tags",matype,nelem,symtmp))
//       { printf("rtdb_get error: geometry:geometry:tags not found\n");
//       exit(99);}
//
//       j = 0;
//       for (ii=0; ii<nion; ++ii)
//       {
//         i = 3*ii;
//         while ((symtmp[j]!= ((char) 0)) && (symtmp[j]!= ((char) 10)))
//            symb[i++] = symtmp[j++];
//         symb[i] = '\0';
//         if (symtmp[j]==((char) 10)) ++j;
//       }
//    }
// }

// static int ion_find_nkatm(RTDB& myrtdb, const int nion)
//{
//    int nkatm,i,ii,j,found;
//    char *symtmp,*symb;
//    symtmp = new char[3*nion];
//    symb   = new char[3*nion];
//
//    ion_tmpsymbols(myrtdb,nion,symb,symtmp);
//
//    /* determine nkatm */
//    nkatm = 0;
//    for (ii=0; ii<nion; ++ii)
//    {
//       found = 0;
//       for (i=0; i<nkatm; ++i)
//          if (strcmp(&symb[3*ii],&symtmp[3*i])==0) found=1;
//       if (!found)
//       {
//          strcpy(&symtmp[3*nkatm],&symb[3*ii]);
//          ++nkatm;
//       }
//    }
//
//    delete [] symb;
//    delete [] symtmp;
//
//    return nkatm;
// }

/*******************************
 *                             *
 *    ion_find_atomkatmnatm    *
 *                             *
 *******************************/
// static void ion_find_atomkatmnatm(RTDB& myrtdb, const int nion, const int
// nkatm,
//                                char *atom, int *katm, int *natm)
//{
//    int i,ii,j,jj,found;
//    char *symtmp,*symb;
//    symtmp = new char[3*nion];
//    symb   = new char[3*nion];
//    ion_tmpsymbols(myrtdb,nion,symb,symtmp);
//
//    /* determine atom */
//    jj = 0;
//    for (ii=0; ii<nion; ++ii)
//    {
//       found = 0;
//       for (i=0; i<jj; ++i)
//          if (strcmp(&symb[3*ii],&atom[3*i])==0) found=1;
//       if (!found)
//       {
//          strcpy(&atom[3*jj],&symb[3*ii]);
//          ++jj;
//       }
//    }
//
//    /* determine katm */
//    for (ii=0; ii<nion; ++ii)
//    {
//       j = 0;
//       for (i=0; i<nkatm; ++i)
//          if (strcmp(&symb[3*ii],&atom[3*i])==0) j=i;
//       katm[ii] = j;
//    }
//
//    /* determine natm */
//    for (i=0; i<nkatm; ++i)   natm[i] = 0;
//    for (ii=0; ii<nion; ++ii) natm[katm[ii]] += 1;
// }

/* Constructors */

/*********************************
 *                               *
 *          Ion::Ion             *
 *                               *
 *********************************/
// Ion::Ion(RTDB& myrtdb, Control2& control)
//{
//
//    int matype,nelem,ii,i,j,found;
//    char *symtmp,*symb;
//    char date[50];
//    double amu_to_mass = 1822.89;
//
//    /* get parallel mappings */
//    if (!myrtdb.get("geometry:geometry:ncenter",rtdb_int,1,&nion)) nion = 1;
//
//    nkatm = ion_find_nkatm(myrtdb,nion);
//
//    time_step = control.time_step();
//
//    charge    = new double[nion];
//    mass      = new double[nion];
//    dti       = new double[nion];
//    rion0     = new double[3*nion];
//    rion1     = new double[3*nion];
//    rion2     = new double[3*nion];
//    vionhalf  = new double[3*nion];
//    fion1     = new double[3*nion];
//    katm      = new int[nion];
//    natm      = new int[nkatm];
//    atomarray = new char[3*nkatm];
//    zv_psp    = new double[nkatm];
//
//    ion_find_atomkatmnatm(myrtdb,nion,nkatm,atomarray,katm,natm);
//
//    /*  read in ion positions */
//    if (!myrtdb.get("geometry:geometry:coords",rtdb_double,3*nion,rion1))
//    { rion1[0] = 0.0; rion1[1] = 0.0; rion1[2] = 0.0; }
//    if (!myrtdb.get("geometry:geometry:coords",rtdb_double,3*nion,rion2))
//    { rion2[0] = 0.0; rion2[1] = 0.0; rion2[2] = 0.0; }
//    if (!myrtdb.get("geometry:geometry:coords",rtdb_double,3*nion,rion0))
//    { rion0[0] = 0.0; rion0[1] = 0.0; rion0[2] = 0.0; }
//
//    /*  read in masses and charges */
//    if (!myrtdb.get("geometry:geometry:masses",rtdb_double,nion,mass)) mass[0]
//    = 1.0; if
//    (!myrtdb.get("geometry:geometry:charges",rtdb_double,nion,charge))
//    charge[0] = 1.0;
//
//    for (ii=0; ii<nion; ++ii)
//    {
//       mass[ii] *= amu_to_mass;
//       dti[ii]  = (time_step*time_step)/mass[ii];
//    }
//
// }

Ion::Ion(std::string rtdbstring, Control2 &control) 
{
   double amu_to_mass = 1822.89;
 
   auto rtdbjson = json::parse(rtdbstring);
 
   std::string geomname = "geometry";
   if (rtdbjson["geometry"].is_string())
     geomname = rtdbjson["geometry"];
 
   json geomjson = rtdbjson["geometries"][geomname];
 
   nion = geomjson["nion"];
 
   auto symbols = geomjson["symbols"];
 
   std::vector<std::string> tmpsymbols;
   for (auto i = 0; i < symbols.size(); ++i) {
     auto match = std::find(begin(tmpsymbols), end(tmpsymbols), symbols[i]);
     if (match == end(tmpsymbols))
       tmpsymbols.push_back(symbols[i]);
   }
   nkatm = tmpsymbols.size();
 
   time_step = control.time_step();
 
   charge = new double[nion];
   mass = new double[nion];
   dti = new double[nion];
   rion0 = new double[3 * nion];
   rion1 = new double[3 * nion];
   rion2 = new double[3 * nion];
   rion_incell0 = new double[3 * nion];
   vionhalf = new double[3 * nion];
   fion1 = new double[3 * nion];
   katm = new int[nion];
   natm = new int[nkatm];
   atomarray = new char[3 * nkatm];
   zv_psp = new double[nkatm];
 
   for (auto ia = 0; ia < nkatm; ++ia) 
   {
      natm[ia] = 0.0;
      strcpy(&atomarray[3*ia], const_cast<char *>(tmpsymbols[ia].data()));
   }
   for (auto i=0; i<nion; ++i) 
   {
      charge[i] = (double)geomjson["charges"][i];
      mass[i] = ((double)geomjson["masses"][i]) * amu_to_mass;
      dti[i] = (time_step * time_step) / mass[i];
    
      rion0[3*i] = (geomjson["velocities"][3 * i].is_number_float())
                         ? (double)geomjson["velocities"][3 * i]
                         : 0.0;
      rion0[3*i+1] = (geomjson["velocities"][3 * i + 1].is_number_float())
                             ? (double)geomjson["velocities"][3 * i + 1]
                             : 0.0;
      rion0[3*i+2] = (geomjson["velocities"][3 * i + 2].is_number_float())
                             ? (double)geomjson["velocities"][3 * i + 2]
                             : 0.0;
    
      rion1[3*i]   = (double)geomjson["coords"][3*i];
      rion1[3*i+1] = (double)geomjson["coords"][3*i + 1];
      rion1[3*i+2] = (double)geomjson["coords"][3*i + 2];
    
      rion2[3*i]   = (double)geomjson["coords"][3*i];
      rion2[3*i+1] = (double)geomjson["coords"][3*i + 1];
      rion2[3*i+2] = (double)geomjson["coords"][3*i + 2];
    
      fion1[3*i] = 0.0;
      fion1[3*i+1] = 0.0;
      fion1[3*i+2] = 0.0;
    
      auto match = std::find(begin(tmpsymbols), end(tmpsymbols), symbols[i]);
      if (match != end(tmpsymbols)) {
        auto ia = std::distance(begin(tmpsymbols), match);
        katm[i] = ia;
        natm[ia] += 1;
      }
   }

   //Check for Point group here????
   double *rion_sym = fion1; //temporary
   sym_tolerance = (geomjson["symmetry_tolerance"].is_number_float())
                           ? (double) geomjson["symmetry_tolerance"] 
                           : 0.001;

   if (control.is_crystal())
   {
      is_crystal = true;
      rotation_type = "crystal";
      group_name = "P1";
      group_rank = 1;
   }
   else
   {
      determine_point_group(rion1,mass,nion,
                            sym_tolerance,
                            group_name,group_rank,rotation_type,
                            inertia_tensor,inertia_moments,inertia_axes,
                            rion_sym);
   }

    std::fill(rion_sym,rion_sym+3*nion,0.0); //re-zero the array
 
   // generate random initial velocities  (temperature, seed) - only set with
   // random velocities if seed > 0
   seed = -1;
   Tf = -1.0;
   double vgx, vgy, vgz, rr0, rr1, rr2, rr3, rr4, rr5;
   double twopi = 16.0 * atan(1.0);
 
   if (rtdbjson["nwpw"]["car-parrinello"]["initial_velocities"][0]
           .is_number_float())
     Tf = rtdbjson["nwpw"]["car-parrinello"]["initial_velocities"][0];
   else if (rtdbjson["nwpw"]["initial_velocities"][0].is_number_float())
     Tf = rtdbjson["nwpw"]["initial_velocities"][0];
 
   if (rtdbjson["nwpw"]["car-parrinello"]["initial_velocities"][1]
           .is_number_integer())
     seed = rtdbjson["nwpw"]["car-parrinello"]["initial_velocities"][1];
   else if (rtdbjson["nwpw"]["initial_velocities"][1].is_number_integer())
     seed = rtdbjson["nwpw"]["initial_velocities"][1];
 
   if ((Tf >= 0.0) && (seed > 0)) {
     std::srand(seed);
     for (auto i = 0; i < nion; ++i) {
       rr0 = ((double)std::rand()) / ((double)RAND_MAX);
       rr1 = ((double)std::rand()) / ((double)RAND_MAX);
       rr2 = ((double)std::rand()) / ((double)RAND_MAX);
       rr3 = ((double)std::rand()) / ((double)RAND_MAX);
       rr4 = ((double)std::rand()) / ((double)RAND_MAX);
       rr5 = ((double)std::rand()) / ((double)RAND_MAX);
       // std::cout << "RANDS=" << rr0 << " " << rr1;
       // std::cout <<      " " << rr2 << " " << rr3;
       // std::cout <<      " " << rr4 << " " << rr5 << std::endl;
       // std::cout <<      " seed=" << seed;
       // std::cout <<      " Tf=" << Tf << std::endl;
 
       vgx = -(2.00 * kb * Tf / mass[i]) * log(rr0);
       vgy = cos(twopi * rr1);
       rion0[3 * i] = sqrt(vgx) * vgy;
 
       vgx = -(2.00 * kb * Tf / mass[i]) * log(rr2);
       vgy = cos(twopi * rr3);
       rion0[3 * i + 1] = sqrt(vgx) * vgy;
 
       vgx = -(2.00 * kb * Tf / mass[i]) * log(rr4);
       vgy = cos(twopi * rr5);
       rion0[3 * i + 2] = sqrt(vgx) * vgy;
     }
 
     // rescale velocities
     center_v_mass(nion, mass, rion0, &vgx, &vgy, &vgz);
     for (auto i = 0; i < nion; ++i) {
       rion0[3 * i] -= vgx;
       rion0[3 * i + 1] -= vgy;
       rion0[3 * i + 2] -= vgz;
     }
     eki0 = ke();
     double Tscale;
     if (nion > 2)
       Tscale = 2.0 * eki0 / (3.0 * nion - 6.0) / kb;
     else
       Tscale = 2.0 * eki0 / kb;
     Tscale = sqrt(Tf / Tscale);
     for (auto i = 0; i < (3 * nion); ++i)
       rion0[i] *= Tscale;
   }
 
   // generate initial kinetic energies
   ekg = ke_com();
   eki0 = ke();
 
   // shift by velocity COM
   bool do_com_shift = true;
   if (rtdbjson["nwpw"]["car-parrinello"]["com_shift"].is_boolean())
     do_com_shift = rtdbjson["nwpw"]["car-parrinello"]["com_shift"];
 
   if (do_com_shift) {
     center_v_mass(nion, mass, rion0, &vgx, &vgy, &vgz);
     for (auto i = 0; i < nion; ++i) {
       rion0[3 * i] -= vgx;
       rion0[3 * i + 1] -= vgy;
       rion0[3 * i + 2] -= vgz;
     }
   }
 
   // scale velocities then find kinetic energy
   // rr0 = control.scaling(1);
   rr0 = control.ion_scaling();
   for (auto i = 0; i < (3 * nion); ++i)
     rion0[i] *= rr0;
   eki1 = ke();
 
   // set intitial vionhalf
   std::memcpy(vionhalf, rion0, 3 * nion * sizeof(double));
 
   // set ke_count, ke_total,and kg_total and g_dof
   ke_count = 0;
   ke_total = 0.0;
   kg_total = 0.0;
   if (rtdbjson["nwpw"]["ke_count"].is_number_integer()) ke_count = rtdbjson["nwpw"]["ke_count"];
   if (rtdbjson["nwpw"]["ke_total"].is_number_float())   ke_total = rtdbjson["nwpw"]["ke_total"];
   if (rtdbjson["nwpw"]["kg_total"].is_number_float())   kg_total = rtdbjson["nwpw"]["kg_total"];
 
   if (rtdbjson["nwpw"]["fix_translation"].is_boolean()) fix_translation = rtdbjson["nwpw"]["fix_translation"];
   if (rtdbjson["nwpw"]["fix_rotation"].is_boolean())    fix_rotation    = rtdbjson["nwpw"]["fix_rotation"];
 
   dof_translation = !fix_translation;
 
   if (rtdbjson["nwpw"]["dof_translation"].is_boolean()) dof_translation = rtdbjson["nwpw"]["dof_translation"];
   if (rtdbjson["nwpw"]["dof_rotation"].is_boolean())    dof_rotation    = rtdbjson["nwpw"]["dof_rotation"];
 
   g_dof = 3.0*nion - 6.0;
   if (dof_translation) g_dof += 3;
   if (dof_rotation)    g_dof += 3;
   if (g_dof < 1)       g_dof = 1.0;
 
   mybond     = new (std::nothrow) ion_bond(rion1,control);     has_bond_constraints = mybond->has_bond();
   mybondings = new (std::nothrow) ion_bondings(rion1,control); has_bondings_constraints = mybondings->has_bondings();
   has_constraints = has_bond_constraints || has_bondings_constraints;

   // dispersion stuff
   disp_on = control.has_disp();
   if (disp_on)
   {
      disp_options = control.options_disp();
      if (control.version==3)
         disp_options += " -pbc";
      ua_disp = control.unita_ptr();
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
 *      Ion::set_rion_incell   *
 *                             *
 *******************************/
void Ion::set_rion_incell(const int n, double *unita)
{
   if (n==0)      std::memcpy(rion_incell0,rion0,3*nion*sizeof(double));
   else if (n==2) std::memcpy(rion_incell0,rion2,3*nion*sizeof(double));
   else           std::memcpy(rion_incell0,rion1,3*nion*sizeof(double));
   incell1(nion,unita,rion_incell0);
}


/*******************************
 *                             *
 *      Ion::writefilename     *
 *                             *
 *******************************/
void Ion::writefilename(std::string &filename) 
{
   // open xyz file
   /*
   std::ofstream xyz(filename.data(),std::ios::app);
         
   double AACONV = 0.529177;

   xyz << nion << std::endl << std::endl;
   for (auto ii=0; ii<nion; ++ii)
      xyz << xyzstream(symbol(ii), 
                       rion(0,ii)*AACONV,rion(1,ii)*AACONV,rion(2,ii)*AACONV,
                       vion(0,ii)*AACONV,vion(1,ii)*AACONV,vion(2,ii)*AACONV) << std::endl;
                       */
}


/*******************************
 *                             *
 *      Ion::writejsonstr      *
 *                             *
 *******************************/

void Ion::writejsonstr(std::string &rtdbstring) {
  auto rtdbjson = json::parse(rtdbstring);

  std::string geomname = "geometry";
  if (rtdbjson["geometry"].is_string())
    geomname = rtdbjson["geometry"];

  // write coordinates
  std::vector<double> coords;
  for (auto i = 0; i < (3 * nion); ++i)
    coords.push_back(rion1[i]);
  rtdbjson["geometries"][geomname]["coords"] = coords;

  // write velocities and ke running averages
  if (ke_count > 0) {
    std::vector<double> vcoords;
    for (auto i = 0; i < (3 * nion); ++i)
      vcoords.push_back(rion0[i]);
    rtdbjson["geometries"][geomname]["velocities"] = vcoords;

    rtdbjson["nwpw"]["ke_count"] = ke_count;
    rtdbjson["nwpw"]["ke_total"] = ke_total;
    rtdbjson["nwpw"]["kg_total"] = kg_total;
  }

  rtdbstring = rtdbjson.dump();
}

/*******************************
 *                             *
 * Ion::remove_com_translation *
 *                             *
 *******************************/
void Ion::remove_com_translation() {
  double am0 = 0.0;
  double hx = 0.0;
  double gx = 0.0;
  double hy = 0.0;
  double gy = 0.0;
  double hz = 0.0;
  double gz = 0.0;
  for (auto ii = 0; ii < nion; ++ii) {
    am0 += mass[ii];
    gx += mass[ii] * rion1[3 * ii];
    gy += mass[ii] * rion1[3 * ii + 1];
    gz += mass[ii] * rion1[3 * ii + 2];
    hx += mass[ii] * rion2[3 * ii];
    hy += mass[ii] * rion2[3 * ii + 1];
    hz += mass[ii] * rion2[3 * ii + 2];
  }
  hx /= am0;
  hy /= am0;
  hz /= am0;
  gx /= am0;
  gy /= am0;
  gz /= am0;

  for (auto ii = 0; ii < nion; ++ii) {
    rion2[3 * ii] += gx - hx;
    rion2[3 * ii + 1] += gy - hy;
    rion2[3 * ii + 2] += gz - hz;
  }
}

/*******************************
 *                             *
 *    Ion::remove_rotation     *
 *                             *
 *******************************/
void Ion::remove_rotation() {
  double h = 1.0 / (2.0 * time_step);

  // center of mass
  double tmass = 0.0;
  double cm[3] = {0.0, 0.0, 0.0};
  for (auto ii = 0; ii < nion; ++ii) {
    tmass += mass[ii];
    cm[0] += mass[ii] * rion1[3 * ii];
    cm[1] += mass[ii] * rion1[3 * ii + 1];
    cm[2] += mass[ii] * rion1[3 * ii + 2];
  }
  cm[0] /= tmass;
  cm[1] /= tmass;
  cm[2] /= tmass;

  // total angular momentum and inertia
  double L[3] = {0.0, 0.0, 0.0};
  double Im[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (auto ii = 0; ii < nion; ++ii) {
    double temp[3] = {rion1[3 * ii] - cm[0], rion1[3 * ii + 1] - cm[1],
                      rion1[3 * ii + 2] - cm[2]};
    double v[3] = {h * (rion2[3 * ii] - rion0[3 * ii]),
                   h * (rion2[3 * ii + 1] - rion0[3 * ii + 1]),
                   h * (rion2[3 * ii + 2] - rion0[3 * ii + 2])};
    L[0] += mass[ii] * (temp[1] * v[2] - temp[2] * v[1]);
    L[1] += mass[ii] * (temp[2] * v[0] - temp[0] * v[2]);
    L[2] += mass[ii] * (temp[0] * v[1] - temp[1] * v[0]);
    for (auto j = 0; j < 3; ++j)
      for (auto i = 0; i < 3; ++i)
        Im[i + 3 * j] -= mass[ii] * temp[i] * temp[j];
  }

  tmass = Im[0] + Im[4] + Im[8];
  Im[0] -= tmass;
  Im[4] -= tmass;
  Im[8] -= tmass;
  double L2 = L[0] * L[0] + L[1] * L[1] + L[2] * L[2];

  if (L2 > 1.0e-12) {
    double omega[9];
    solve_3by3(Im, L, omega);
    double hinv = 1.0 / h;

    for (auto ii = 0; ii < nion; ++ii) {
      double temp[3] = {rion1[3 * ii] - cm[0], rion1[3 * ii + 1] - cm[1],
                        rion1[3 * ii + 2] - cm[2]};
      double v[3] = {(omega[1] * temp[2] - omega[2] * temp[1]),
                     (omega[2] * temp[0] - omega[0] * temp[2]),
                     (omega[0] * temp[1] - omega[1] * temp[0])};
      rion2[3 * ii] -= v[0] * hinv;
      rion2[3 * ii + 1] -= v[1] * hinv;
      rion2[3 * ii + 2] -= v[2] * hinv;
    }
  }
}


/*******************************************
 *                                         *
 *         Ion::print_constraints          *
 *                                         *
 *******************************************/
std::string Ion::print_constraints(const int opt) 
{
    std::string tmp = "";
    if (has_constraints)
    {
       if (opt==0) tmp = " ion constraints:\n";
       if (opt==1) tmp = " == Ion Constraints ==\n\n";
       if (has_bond_constraints)     tmp += mybond->print_all(opt);
       if (has_bondings_constraints) tmp += mybondings->print_all(opt);
    }
    return tmp;
}



/*******************************************
 *                                         *
 *        Ion::print_symmetry_group        *
 *                                         *
 *******************************************/
std::string Ion::print_symmetry_group()
{
   std::stringstream stream;

   std::ios init(NULL);
   init.copyfmt(stream);

   stream << " symmetry information: (symmetry_tolerance = " << Efmt(8,2) << this->sym_tolerance << ")" <<  std::endl;
   stream << "      group name   : " << this->group_name
                                     << "  (group rank = "  << this->group_rank
                                     << " rotation type : " << this->rotation_type <<")" <<  std::endl;
   if (!is_crystal)
   {
      stream << "      inertia axes : e1 = <" << Ffmt(8,3) << this->inertia_axes[0] << " "
                                              << Ffmt(8,3) << this->inertia_axes[1] << " "
                                              << Ffmt(8,3) << this->inertia_axes[2] << " > - "
                                              << "moment =" << Efmt(14,7) << this->inertia_moments[0] << std::endl;
      stream << "                     e2 = <" << Ffmt(8,3) << this->inertia_axes[3] << " "
                                              << Ffmt(8,3) << this->inertia_axes[4] << " "
                                              << Ffmt(8,3) << this->inertia_axes[5] << " > - "
                                              << "moment =" << Efmt(14,7) << this->inertia_moments[1] << std::endl;
      stream << "                     e3 = <" << Ffmt(8,3) << this->inertia_axes[6] << " "
                                              << Ffmt(8,3) << this->inertia_axes[7] << " "
                                              << Ffmt(8,3) << this->inertia_axes[8] << " > - "
                                              << "moment =" << Efmt(14,7) << this->inertia_moments[2] << std::endl;
   }
   return stream.str();
}




/*******************************************
 *                                         *
 *           Ion::disp_energy              *
 *                                         *
 *******************************************/
 /*
double Ion::disp_energy() 
{
   double edisp = 0.0; 
   if (disp_on)
   {
      double g[3*nion];
      double g_lat[9];
      int icharge[nion];
      for (auto ii=0; ii<nion; ++ii)
          icharge[ii] = static_cast<int>(charge[ii]);
      nwpwxc_vdw3_dftd3_( (char *)disp_options.c_str(), &nion,icharge,rion1,ua_disp,&edisp,g,g_lat);
   }
   return edisp;
}
*/
   


} // namespace pwdft
