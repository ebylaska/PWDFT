/* Ions.C -
   Author - Eric Bylaska
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "json.hpp"
using json = nlohmann::json;

#include "Control2.hpp"
#include "ion_bondings.hpp"
#include "Ion.hpp"

namespace pwdft {

/*******************************
 *                             *
 *        solve_3by3           *
 *                             *
 *******************************/
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
   vionhalf = new double[3 * nion];
   fion1 = new double[3 * nion];
   katm = new int[nion];
   natm = new int[nkatm];
   atomarray = new char[3 * nkatm];
   zv_psp = new double[nkatm];
 
   for (auto ia = 0; ia < nkatm; ++ia) {
     natm[ia] = 0.0;
     strcpy(&atomarray[3 * ia], const_cast<char *>(tmpsymbols[ia].data()));
   }
   for (auto i = 0; i < nion; ++i) {
     charge[i] = (double)geomjson["charges"][i];
     mass[i] = ((double)geomjson["masses"][i]) * amu_to_mass;
     dti[i] = (time_step * time_step) / mass[i];
 
     rion0[3 * i] = (geomjson["velocities"][3 * i].is_number_float())
                        ? (double)geomjson["velocities"][3 * i]
                        : 0.0;
     rion0[3 * i + 1] = (geomjson["velocities"][3 * i + 1].is_number_float())
                            ? (double)geomjson["velocities"][3 * i + 1]
                            : 0.0;
     rion0[3 * i + 2] = (geomjson["velocities"][3 * i + 2].is_number_float())
                            ? (double)geomjson["velocities"][3 * i + 2]
                            : 0.0;
 
     rion1[3 * i] = (double)geomjson["coords"][3 * i];
     rion1[3 * i + 1] = (double)geomjson["coords"][3 * i + 1];
     rion1[3 * i + 2] = (double)geomjson["coords"][3 * i + 2];
 
     rion2[3 * i] = (double)geomjson["coords"][3 * i];
     rion2[3 * i + 1] = (double)geomjson["coords"][3 * i + 1];
     rion2[3 * i + 2] = (double)geomjson["coords"][3 * i + 2];
 
     fion1[3 * i] = 0.0;
     fion1[3 * i + 1] = 0.0;
     fion1[3 * i + 2] = 0.0;
 
     auto match = std::find(begin(tmpsymbols), end(tmpsymbols), symbols[i]);
     if (match != end(tmpsymbols)) {
       auto ia = std::distance(begin(tmpsymbols), match);
       katm[i] = ia;
       natm[ia] += 1;
     }
   }
 
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
 
   mybond     = new (std::nothrow) ion_bond(rion1,control);
   mybondings = new (std::nothrow) ion_bondings(rion1,control);
   has_constraints = mybond->has_bond() || mybondings->has_bondings();

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
std::string Ion::print_constraints() 
{
    std::string tmp = "";
    if (has_constraints)
    {
       tmp = " ion constraints:\n";
       tmp += mybond->print_all();
       tmp += mybondings->print_all();
    }
    return tmp;
}



} // namespace pwdft
