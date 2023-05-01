/* nwpw_aimd_running_data.cpp -
   Author - Eric Bylaska
*/

#include "nwpw_aimd_running_data.hpp"
#include "blas.h"
#include "gdevice.hpp"
#include "nwpw_timing.hpp"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>

// Option for C++17
//#include	<filesystem>
// namespace fs = std::filesystem;

// Option for C++ before C++17
namespace fs {
inline bool exists(const std::string &filename) {
   struct stat buffer;
   return (stat(filename.c_str(), &buffer) == 0);
}
inline void remove(const std::string &filename) {
   int result = std::remove(filename.c_str());
}
inline void copy(const std::string &source_filename,
                 const std::string &dest_filename) {
   std::ifstream src(source_filename, std::ios::binary);
   std::ofstream dst(dest_filename, std::ios::binary);
   dst << src.rdbuf();
}
} // namespace fs

#define E124                                                                   \
  std::right << std::setw(12) << std::setprecision(4) << std::scientific
#define E156                                                                   \
  std::right << std::setw(15) << std::setprecision(6) << std::scientific
#define E1910                                                                  \
  std::right << std::setw(19) << std::setprecision(10) << std::scientific
#define F206 std::right << std::setw(20) << std::setprecision(6) << std::fixed
#define F105 std::right << std::setw(10) << std::setprecision(5) << std::fixed

#define hxyzstream(a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z)                \
  E124 << (a1x) << E124 << (a1y) << E124 << (a1z) << E124 << (a2x) << E124     \
       << (a2y) << E124 << (a2z) << E124 << (a3x) << E124 << (a3y) << E124     \
       << (a3z)

#define xyzstream(S, X, Y, Z, VX, VY, VZ)                                      \
  std::left << std::setw(3) << (S) << E124 << (X) << E124 << (Y) << E124       \
            << (Z) << E124 << (VX) << E124 << (VY) << E124 << (VZ)

#define hionstream(t, N, V, a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z)       \
  E156 << (t) << std::setw(6) << (N) << E156 << (V) << E156 << (a1x) << E156   \
       << (a1y) << E156 << (a1z) << E156 << (a2x) << E156 << (a2y) << E156     \
       << (a2z) << E156 << (a3x) << E156 << (a3y) << E156 << (a3z)

#define ionstream(II, S1, S2, X, Y, Z, VX, VY, VZ)                             \
  std::setw(6) << (II) << std::setw(3) << (S1) << std::setw(5) << (S2) << E156 \
               << (X) << E156 << (Y) << E156 << (Z) << E156 << (VX) << E156    \
               << (VY) << E156 << (VZ)

namespace pwdft {

/* Constructors */

/****************************************************
 *                                                  *
 *  nwpw_aimd_running_data::nwpw_aimd_running_data  *
 *                                                  *
 ****************************************************/
nwpw_aimd_running_data::nwpw_aimd_running_data(Control2 &control, Parallel *inparall, Pneb *inpneb, Ion *inion, 
                                               double *Ein, double *hmlin, double *psiin, double *dnin) 
{
   myparall = inparall;
   myion = inion;
   mypneb = inpneb;
   E = Ein;
   hml = hmlin;
   psi = psiin;
   dn = dnin;
 
   calculate_fei      = control.Fei_on();
   calculate_cif      = control.CIF_on();
   calculate_mulliken = control.Mulliken_on();
   calculate_dipole   = control.dipole_on();
   use_nose_output    = control.Nose();
 
   bool includedatafiles = true;
   if (control.task()==5) 
      includedatafiles = control.geometry_optimize();

   if (includedatafiles) 
   {
      if (control.task()==7) 
      {
         dt = control.bo_time_step();
         dt_inner = dt*control.bo_steps(0);
      } 
      else 
      {
         dt = control.time_step();
         dt_inner = dt*control.loop(0);
      }
  
      // set running data filenames
      xyz_filename           = control.permanent_dir_str + "/" + control.xyz_filename;
      ion_motion_filename    = control.permanent_dir_str + "/" + control.ion_motion_filename;
      emotion_filename       = control.permanent_dir_str + "/" + control.emotion_filename;
      fei_filename           = control.permanent_dir_str + "/" + control.fei_filename;
      cif_filename           = control.permanent_dir_str + "/" + control.cif_filename;
      dipole_motion_filename = control.permanent_dir_str + "/" + control.dipole_motion_filename;
      eigmotion_filename     = control.permanent_dir_str + "/" + control.eigmotion_filename;
      omotion_filename       = control.permanent_dir_str + "/" + control.omotion_filename;
      hmotion_filename       = control.permanent_dir_str + "/" + control.hmotion_filename;
  
      // set running data backup filenames
      xyz_bakfile       = control.permanent_dir_str + "/" + "XYZ99-bak";
      motion_bakfile    = control.permanent_dir_str + "/" + "MOTION99-bak";
      emotion_bakfile   = control.permanent_dir_str + "/" + "EMOTION99-bak";
      fei_bakfile       = control.permanent_dir_str + "/" + "FEI99-bak";
      cif_bakfile       = control.permanent_dir_str + "/" + "CIF99-bak";
      dipole_bakfile    = control.permanent_dir_str + "/" + "DIPOLE99-bak";
      eigmotion_bakfile = control.permanent_dir_str + "/" + "EIGMOTION999-bak";
      omotion_bakfile   = control.permanent_dir_str + "/" + "OMOTION99-bak";
      hmotion_bakfile   = control.permanent_dir_str + "/" + "HMOTION99-bak";
  
      // standard datafiles
      if ((myparall->is_master()) && includedatafiles) 
      {
         // check for backup files
         if (fs::exists(xyz_bakfile))     fs::copy(xyz_bakfile, xyz_filename);
         if (fs::exists(motion_bakfile))  fs::copy(motion_bakfile, ion_motion_filename);
         if (fs::exists(emotion_bakfile)) fs::copy(emotion_bakfile, emotion_filename);
        
         // make new backup files
         if (fs::exists(xyz_filename))        fs::copy(xyz_filename, xyz_bakfile);
         if (fs::exists(ion_motion_filename)) fs::copy(ion_motion_filename, motion_bakfile);
         if (fs::exists(emotion_filename))    fs::copy(emotion_filename, emotion_bakfile);
        
         // open xyz file
         xyz_open = true;
         xyz = new (std::nothrow) std::ofstream;
         xyz->open(const_cast<char *>(xyz_filename.data()), std::ios::app);
        
         // ion_motion file
         ion_motion_open = true;
         ion_motion = new (std::nothrow) std::ofstream;
         ion_motion->open(const_cast<char *>(ion_motion_filename.data()), std::ios::app);
        
         // read emotion data
         emotion_ishift = 0;
         emotion_time_shift = 0.0;
         if (fs::exists(emotion_filename))
         {
            std::ifstream emotiontmp;
            emotiontmp.open(const_cast<char *>(emotion_filename.data()));
           
            std::string line;
            while (std::getline(emotiontmp, line)) 
            {
               double w, a, b, c, d;
               std::istringstream iss(line);
              
               if (!(iss >> emotion_time_shift >> w >> a >> b >> c >> d)) 
               {
                  break;
               } // error
              
               // take care of running sums
               E[24] += a;
               E[25] += a * a;
               E[26] += (a + b + c);
               E[27] += (a + b + c) * (a + b + c);
               E[22] += d;
               E[23] += d * d;
               ++emotion_ishift;
           }
           emotiontmp.close();
         }
        
         // emotion file
         emotion_open = true;
         emotion = new (std::nothrow) std::ofstream;
         emotion->open(const_cast<char *>(emotion_filename.data()), std::ios::app);
        
         ion_motion_ishift = emotion_ishift;
      }
  
      // fei file
      if (calculate_fei && (myparall->is_master())) 
      {
         // check for a backup file
         if (fs::exists(fei_bakfile)) fs::copy(fei_bakfile, fei_filename);
        
         // make a new backup file
         if (fs::exists(fei_filename)) fs::copy(fei_filename, fei_bakfile);
        
         // fei file
         fei_open = true;
         fei = new (std::nothrow) std::ofstream;
         fei->open(const_cast<char *>(fei_filename.data()), std::ios::app);
      }
  
      // cif file
      if (calculate_cif && (myparall->is_master())) 
      {
         // check for a backup file
         if (fs::exists(cif_bakfile)) fs::copy(cif_bakfile, cif_filename);
     
         // make a new backup file
         if (fs::exists(cif_filename)) fs::copy(cif_filename, cif_bakfile);
     
         // cif file
         cif_open = true;
         cif = new (std::nothrow) std::ofstream;
         cif->open(const_cast<char *>(cif_filename.data()), std::ios::app);
      }
  
      // dipole file
      if (calculate_dipole && (myparall->is_master())) 
      {
         // check for a backup file
         if (fs::exists(dipole_bakfile)) fs::copy(dipole_bakfile, dipole_motion_filename);
        
         // make a new backup file
         if (fs::exists(dipole_motion_filename)) fs::copy(dipole_motion_filename, dipole_bakfile);
        
         // dipole_motion file
         dipole_motion_open = true;
         dipole_motion = new (std::nothrow) std::ofstream;
         dipole_motion->open(const_cast<char *>(dipole_motion_filename.data()),std::ios::app);
      }
  
      // Mulliken data
      if (calculate_mulliken && (myparall->is_master())) 
      {
         // check for backup files
         if (fs::exists(eigmotion_bakfile)) fs::copy(eigmotion_bakfile, eigmotion_filename);
         if (fs::exists(omotion_bakfile))   fs::copy(omotion_bakfile, omotion_filename);
         if (fs::exists(hmotion_bakfile))   fs::copy(hmotion_bakfile, hmotion_filename);
        
         // make new backup files
         if (fs::exists(eigmotion_filename)) fs::copy(eigmotion_filename, eigmotion_bakfile);
         if (fs::exists(omotion_filename))   fs::copy(omotion_filename, omotion_bakfile);
         if (fs::exists(hmotion_filename))   fs::copy(hmotion_filename, hmotion_bakfile);
        
         // eigmotion file
         eigmotion_open = true;
         eigmotion = new (std::nothrow) std::ofstream;
         eigmotion->open(const_cast<char *>(eigmotion_filename.data()),
                         std::ios::app);
        
         // hmotion file
         hmotion_open = true;
         hmotion = new (std::nothrow) std::ofstream;
         hmotion->open(const_cast<char *>(hmotion_filename.data()), std::ios::app);
        
         // omotion file
         omotion_open = true;
         omotion = new (std::nothrow) std::ofstream;
         omotion->open(const_cast<char *>(omotion_filename.data()), std::ios::app);
      }
   }
}

/********************************************
 *                                          *
 * nwpw_aimd_running_data::update_iteration *
 *                                          *
 ********************************************/
void nwpw_aimd_running_data::update_iteration(const int icount) 
{
   // Running xyz
   if (xyz_open) 
   {
      double AACONV = 0.529177;
      
      double a1x = mypneb->lattice->unita(0,0)*AACONV;
      double a1y = mypneb->lattice->unita(1,0)*AACONV;
      double a1z = mypneb->lattice->unita(2,0)*AACONV;
      double a2x = mypneb->lattice->unita(0,1)*AACONV;
      double a2y = mypneb->lattice->unita(1,1)*AACONV;
      double a2z = mypneb->lattice->unita(2,1)*AACONV;
      double a3x = mypneb->lattice->unita(0,2)*AACONV;
      double a3y = mypneb->lattice->unita(1,2)*AACONV;
      double a3z = mypneb->lattice->unita(2,2)*AACONV;
      
      *xyz << myion->nion << std::endl
           << hxyzstream(a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z)
           << std::endl;
      
      for (auto ii = 0; ii < myion->nion; ++ii)
        *xyz << xyzstream(
                    myion->symbol(ii), 
                    myion->rion(0,ii)*AACONV,
                    myion->rion(1,ii)*AACONV, 
                    myion->rion(2,ii)*AACONV,
                    myion->vion(0,ii)*AACONV, 
                    myion->vion(1,ii)*AACONV,
                    myion->vion(2,ii)*AACONV) << std::endl;
   }
 
   // Running ion_motion
   if (ion_motion_open) 
   {
      double current_time = dt_inner * (icount + ion_motion_ishift);
     
      double omega = mypneb->lattice->omega();
      double a1x = mypneb->lattice->unita(0,0);
      double a1y = mypneb->lattice->unita(1,0);
      double a1z = mypneb->lattice->unita(2,0);
      double a2x = mypneb->lattice->unita(0,1);
      double a2y = mypneb->lattice->unita(1,1);
      double a2z = mypneb->lattice->unita(2,1);
      double a3x = mypneb->lattice->unita(0,2);
      double a3y = mypneb->lattice->unita(1,2);
      double a3z = mypneb->lattice->unita(2,2);
     
      *ion_motion << hionstream(current_time,myion->nion,omega,
                                a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z)
                  << std::endl;
     
      for (auto ii = 0; ii < myion->nion; ++ii)
        *ion_motion << ionstream(ii,myion->symbol(ii),myion->symbol(ii),
                                 myion->rion(0,ii), myion->rion(1,ii), myion->rion(2,ii), 
                                 myion->vion(0,ii), myion->vion(1,ii), myion->vion(2,ii))
                    << std::endl;
   }
 
   // Running emotion
   if (emotion_open) 
   {
      double current_time = dt_inner*(icount+ emotion_ishift);
      double fac = ((double)(icount+emotion_ishift));
      double pressure = 0.0;
     
      eave = E[24]/fac;
      evar = E[25]/fac;
      evar -= eave*eave;
      have = E[26]/fac;
      hvar = E[27]/fac;
      hvar -= have*have;
     
      qave = E[22]/fac;
      qvar = E[23]/fac;
      qvar -= qave*qave;
     
      if (use_nose_output)
         *emotion << E1910 << current_time << E1910 << E[0] << E1910 << E[1]
                  << E1910 << E[2] << E1910 << E[3] << E1910 << E[13] << E1910
                  << E[4] << E1910 << E[5] << E1910 << E[6] << E1910 << E[7]
                  << E1910 << E[8] << E1910 << E[9] << E1910 << eave << E1910
                  << evar << E1910 << have << E1910 << hvar << E1910 << qave
                  << E1910 << qvar << E1910 << myion->Temperature() << E1910
                  << pressure << std::endl;
      else
         *emotion << E1910 << current_time << E1910 << E[0] << E1910 << E[1]
                  << E1910 << E[2] << E1910 << E[3] << E1910 << E[13] << E1910
                  << E[4] << E1910 << E[5] << E1910 << E[6] << E1910 << E[7]
                  << E1910 << eave << E1910 << evar << E1910 << have << E1910
                  << hvar << E1910 << qave << E1910 << qvar << E1910
                  << myion->Temperature() << E1910 << pressure << std::endl;
   }
 
   // Running fei
   if (fei_open) 
   {
      double a1x = mypneb->lattice->unita(0,0);
      double a1y = mypneb->lattice->unita(1,0);
      double a1z = mypneb->lattice->unita(2,0);
      double a2x = mypneb->lattice->unita(0,1);
      double a2y = mypneb->lattice->unita(1,1);
      double a2z = mypneb->lattice->unita(2,1);
      double a3x = mypneb->lattice->unita(0,2);
      double a3y = mypneb->lattice->unita(1,2);
      double a3z = mypneb->lattice->unita(2,2);
      *fei << myion->nion << std::endl
           << F206 << E[1] << std::endl
           << F105 << a1x << F105 << a1y << F105 << a1z << std::endl
           << F105 << a2x << F105 << a2y << F105 << a2z << std::endl
           << F105 << a3x << F105 << a3y << F105 << a3z << std::endl;
      for (auto ii = 0; ii < myion->nion; ++ii) 
      {
         *fei << std::setw(5) << (ii + 1) << std::setw(3) << myion->symbol(ii)
              << std::setw(5) << myion->symbol(ii) << F105 << myion->rion(0, ii)
              << F105 << myion->rion(1, ii) << F105 << myion->rion(2, ii) << F105
              << myion->fion(0, ii) << F105 << myion->fion(1, ii) << F105
              << myion->fion(2, ii) << std::endl;
      }
   }
 
   // Running cif
   if (cif_open) { }
 
   // Running dipole
   if (calculate_dipole) { }
 
   // Running Mulliken
   if (calculate_mulliken) { }
}

/*******************************************
 *                                         *
 * nwpw_aimd_running_data::remove_bakfiles *
 *                                         *
 *******************************************/
void nwpw_aimd_running_data::remove_bakfiles() 
{
   if (xyz_open)        if (fs::exists(xyz_bakfile))     fs::remove(xyz_bakfile);
   if (ion_motion_open) if (fs::exists(motion_bakfile))  fs::remove(motion_bakfile);
   if (emotion_open)    if (fs::exists(emotion_bakfile)) fs::remove(emotion_bakfile);
 
   if (fei_open) if (fs::exists(fei_bakfile)) fs::remove(fei_bakfile);
   if (cif_open) if (fs::exists(cif_bakfile)) fs::remove(cif_bakfile);
 
   if (dipole_motion_open) if (fs::exists(dipole_bakfile)) fs::remove(dipole_bakfile);
 
   if (eigmotion_open) if (fs::exists(eigmotion_bakfile)) fs::remove(eigmotion_bakfile);
   if (omotion_open)   if (fs::exists(omotion_bakfile))   fs::remove(omotion_bakfile);
   if (hmotion_open)   if (fs::exists(hmotion_bakfile))   fs::remove(hmotion_bakfile);
}

} // namespace pwdft
