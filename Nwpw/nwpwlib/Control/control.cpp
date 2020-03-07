#include        <iostream>
#include        <cmath>
#include        <cstdlib>
#include	"Int64.h"
using namespace std;

static double unita[9],tolerances[3],scaling[2];
static double time_step,fake_mass,ks_alpha,ecut,wcut,rcut;
static double bo_time_step;
static Int64 bo_steps[2],bo_algorithm;
static Int64 loop[2],ngrid[3],npsp,ncut,mapping,mapping1d;
static Int64 np_dimensions[3],ewald_grid[3];
static Int64 code,gga;
static Int64 ispin,multiplicity;
static Int64 move,frac_coord,SA,fei,fei_quench,gram_schmidt;
static Int64 rotation,translation,balance,spin_orbit;
static Int64 maxit_orb,maxit_orbs,scf_algorithm,ks_algorithm;
static Int64 symm_number;
static Int64 qsize;

static char input_movecs_filename[80];




#include	"json.hpp"
#include	"rtdb.hpp"
#include	"control.hpp"

using json = nlohmann::json;

void   control_check_unita_for_default(RTDB& myrtdb)
{

}

/***********************************
 *                                 *
 *         control_read_json       *
 *                                 *
 ***********************************/
void control_read_json(string jsonstr)
{
   json j = json::parse(jsonstr);
}

/***********************************
 *                                 *
 *         control_read            *
 *                                 *
 ***********************************/
void control_read(RTDB &myrtdb)
{
   //RTDB myrtdb(myparallel, "eric.db", "old");

   Int64 idum[10];
   int matype,nelem;
   char date[99];

   myrtdb.parallel_mode = 1;
   int np = myrtdb.parall->np();

   /* get parallel mappings */
   if (!myrtdb.get("nwpw:mapping",rtdb_int,1,&mapping)) mapping = 1;

   /* set mapping1d */
   if (!myrtdb.get("nwpw:mapping1d",rtdb_int,1,&mapping1d)) mapping1d = 1;

   /* qsize */
   if (!myrtdb.get("nwpw:pfft3_qsize",rtdb_int,1,&qsize)) qsize = 4;

   /* get np_dimensions */
   if (!myrtdb.get("nwpw:np_dimensions",rtdb_int,3,np_dimensions))
   {
      np_dimensions[0] = np;
      np_dimensions[1] = 1;
      np_dimensions[2] = 1;
   }
   else
   {
      np_dimensions[2] = 1;
      if (np_dimensions[1]<1) np_dimensions[1] = 1;
      if (np_dimensions[0]<1) np_dimensions[0] = 1;

      /* reset np_dimensions[1] if it is not a  multiple of np2 */
      np = np/np_dimensions[2];
      while ( ((np%np_dimensions[1])!=0) && (np_dimensions[1]>1) )
            np_dimensions[1] = np_dimensions[1] - 1;
      np_dimensions[0] = np/np_dimensions[1];
    }

   if (!myrtdb.get("nwpw:balance",rtdb_log,1,&balance)) balance = 1;

   if (!myrtdb.get_info("cpsd:input_wavefunction_filename",&matype,&nelem,date))
      strcpy(input_movecs_filename,"eric.movecs");
   else
      if (!myrtdb.get("cpsd:input_wavefunction_filename",matype,nelem,input_movecs_filename))
         strcpy(input_movecs_filename,"eric.movecs");

   if (!myrtdb.get("cpsd:fake_mass",rtdb_double,1,&fake_mass)) fake_mass = 400000.0;
   if (!myrtdb.get("cpsd:time_step",rtdb_double,1,&time_step)) time_step = 5.8;
   if (!myrtdb.get("cpsd:loop",rtdb_int,2,loop)) 
      { loop[0]=10 ; loop[1]=100; }
   if (!myrtdb.get("cpsd:tolerances",rtdb_double,3,tolerances))
      { tolerances[0] = 1.0e-9; tolerances[1] = 1.0e-9; tolerances[2] = 1.0e-4; }

   if (!myrtdb.get("cpsd:ecut",rtdb_double,1,&ecut))      ecut=9000.0;
   if (!myrtdb.get("cpsd:wcut",rtdb_double,1,&wcut))      wcut=ecut;
   if (!myrtdb.get("cpsd:rcut",rtdb_double,1,&rcut))      rcut = 0.0;
   if (!myrtdb.get("cpsd:ncut",rtdb_int,1,&ncut))         ncut = 1;
   if (!myrtdb.get("cpsd:mult",rtdb_int,1,&multiplicity)) multiplicity = 1;
   if (!myrtdb.get("cpsd:ispin",rtdb_int,1,&ispin))       ispin=1;

   if (!myrtdb.get("cell_default:unita",rtdb_double,9,unita))
   {
       unita[0] = 20.0; unita[1] =  0.0; unita[2] =  0.0;
       unita[3] =  0.0; unita[4] = 20.0; unita[5] =  0.0;
       unita[6] =  0.0; unita[7] =  0.0; unita[8] = 20.0;
   }
   if (!myrtdb.get("cell_default:ngrid",rtdb_int,3,ngrid))
   {
      ngrid[0] = 32; ngrid[1] = 32; ngrid[2] = 32;
   }
   if (!myrtdb.get("nwpw:ewald_ngrid",rtdb_int,3,ewald_grid))
   {
      ewald_grid[0] = ngrid[0]; ewald_grid[1] = ngrid[1]; ewald_grid[2] = ngrid[2];
   }

}
      

int control_balance()          { return (int) balance; }
int control_mapping()          { return (int) mapping; }
int control_mapping1d()        { return (int) mapping1d; }
int control_np_orbital()       { return (int) np_dimensions[1]; }
int control_ngrid(const int i) { return (int) ngrid[i]; }
int control_ewald_ngrid(const int i) { return (int) ewald_grid[i]; }
int control_ewald_ncut() { return (int) ncut; }
int control_loop(const int i) { return (int) loop[i]; }
int control_pfft3_qsize() { return (int) qsize; }

double control_ewald_rcut() { return (double) rcut; }
double control_unita(const int i,const int j) { return unita[i+j*3]; }
double control_ecut() { return ecut; }
double control_wcut() { return wcut; }
double control_time_step() { return time_step; }
double control_fake_mass() { return fake_mass; }
double control_tolerances(const int i) { return tolerances[i];}

char   *control_input_movecs_filename() { return input_movecs_filename; }
