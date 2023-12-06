#ifndef _NWPW_AIMD_RUNNING_DATA_HPP_
#define _NWPW_AIMD_RUNNING_DATA_HPP_

// ********************************************************************
// *                                                                  *
// *       nwpw_aimd_running_data :      ldjfdlsjf                    *
// *                                                                  *
// *                                                                  *
// ********************************************************************
#include "blas.h"
#include <fstream>
#include <iostream>

#include "Control2.hpp"
#include "Ion.hpp"
#include "Parallel.hpp"
#include "Lattice.hpp"
//#include "gdevice.hpp"

namespace pwdft {

class nwpw_aimd_running_data {

  Parallel *myparall;
  Lattice *mylattice;

  Ion *myion;
  double *E, *hml, *psi, *dn;

  bool ismaster, use_nose_output, use_bo_output;

  int emotion_ishift = 0;
  int ion_motion_ishift = 0;
  bool xyz_open = false;
  std::ofstream *xyz;
  std::string xyz_filename, xyz_bakfile;
  bool ion_motion_open = false;
  std::ofstream *ion_motion;
  std::string ion_motion_filename, motion_bakfile;
  bool emotion_open = false;
  std::ofstream *emotion;
  std::string emotion_filename, emotion_bakfile;

  bool calculate_fei = false;
  bool fei_open = false;
  std::ofstream *fei;
  std::string fei_filename, fei_bakfile;

  bool calculate_cif = false;
  bool cif_open = false;
  std::ofstream *cif;
  std::string cif_filename, cif_bakfile;
  bool cif_shift_cell = true;

  bool calculate_dipole = false;
  bool dipole_motion_open = false;
  std::ofstream *dipole_motion;
  std::string dipole_motion_filename, dipole_bakfile;

  bool calculate_mulliken = false;
  bool eigmotion_open = false;
  std::ofstream *eigmotion;
  std::string eigmotion_filename, eigmotion_bakfile;
  bool omotion_open = false;
  std::ofstream *omotion;
  std::string omotion_filename, omotion_bakfile;
  bool hmotion_open = false;
  std::ofstream *hmotion;
  std::string hmotion_filename, hmotion_bakfile;

public:
  double dt, dt_inner, eave, evar, have, hvar, qave, qvar, emotion_time_shift;

  /* constructor */
  nwpw_aimd_running_data(Control2 &, Parallel *, Lattice *, Ion *, double *, double *, double *, double *);

  /* destructor */
  ~nwpw_aimd_running_data() 
  {
     // remove bakfiles
     this->remove_bakfiles();
    
     // close open files
     if (xyz_open) 
     {
        xyz->close();
        delete xyz;
     }
     if (ion_motion_open) 
     {
        ion_motion->close();
        delete ion_motion;
     }
     if (emotion_open) 
     {
        emotion->close();
        delete emotion;
     }
    
     if (fei_open) 
     {
        fei->close();
        delete fei;
     }
     if (cif_open) 
     {
        cif->close();
        delete cif;
     }
    
     if (dipole_motion_open) 
     {
        dipole_motion->close();
        delete dipole_motion;
     }
    
     if (hmotion_open) 
     {
        hmotion->close();
        delete hmotion;
     }
     if (omotion_open) 
     {
        omotion->close();
        delete omotion;
     }
     if (eigmotion_open) 
     {
        eigmotion->close();
        delete eigmotion;
     }
  }

  void update_iteration(const int);
  void remove_bakfiles();
};

} // namespace pwdft

#endif
