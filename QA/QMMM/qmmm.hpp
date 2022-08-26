#ifndef	_QMMM_HPP_
#define _QMMM_HPP_

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



class	QMMM_Operator {

   double pi = 4.0*std::atan(1.0);

   double qmmm_lmbda = 1.0;
   int    nkfrag = 0;
   int    nfrag  = 0;
   int    nbond  = 0;
   int    nangle = 0;

   //set up qmmm lj stuff
   std::map<std::string,std::vector<double>> lj_qm_data,lj_mm_data;

   // switching parameters
   double *switch_Rin;
   double *switch_Rout;

   int *size_frag;
   int *nfrag_frag;
   int *nbond_frag;
   int *nangle_frag;
   int *bond_start_frag;
   int *angle_start_frag;

   int *indxfrag_start;
   int *kfrag;

   bool *self_interaction;

   // bond and angle spring parameters
   int    *bond_indx,*angle_indx;
   double *Krb;
   double *Kra;


public:
   int    nion,nion_qm,nion_mm,nkatm;
   int    *katm;
   double *epsilon;
   double *sigma;
   double *rion;
   double *qion;
   double *mass;
   std::string *aname;
   std::string *symbol;

   /* Constructors */
   QMMM_Operator(std::string nwinput);

   /* destructor */
   ~QMMM_Operator() {
         delete [] katm;
         delete [] epsilon;
         delete [] sigma;
         delete [] aname;
         delete [] symbol;
         delete [] rion;
         delete [] qion;
         delete [] mass;
         delete [] epsilon;
         delete [] sigma;

         delete [] switch_Rin;
         delete [] switch_Rout;

         delete [] size_frag;
         delete [] nfrag_frag;
         delete [] nbond_frag;
         delete [] nangle_frag;
         delete [] bond_start_frag;
         delete [] angle_start_frag;

         delete [] indxfrag_start;
         delete [] kfrag;

         delete [] self_interaction;

         delete [] bond_indx;
         delete [] angle_indx;
         delete [] Krb;
         delete [] Kra;
    }

    double spring_Energy(const double []);
    void spring_Force(const double [], double []);

    double Coulomb_Energy(const double [], const double []);
    void   Coulomb_Force(const double [], const double [], double []);

    double Coulomb_Energy_qmmm(const double [], const double []);
    void   Coulomb_Force_qmmm(const double [], const double [], double []);

    double MMMM_electrostatic_energy(const double [], const double []);

    double QMMM_electrostatic_energy(const double [], const double []);
    void QMMM_electrostatic_potential(const double [], const double [], double []);
    void QMMM_electrostatic_force(const double [], const double [], double []);

    double QMQM_electrostatic_energy(const double [], const double []);
    void QMQM_electrostatic_force(const double [], const double [], double []);



    double QMMM_LJ_energy(const double []);
    void QMMM_LJ_force(const double [], double []);

    double MMMM_LJ_energy(const double []);
    void MMMM_LJ_force(const double [], double []);



    //void   ke(double *, double *);
    //double ke_ave(double *);
};



#endif

