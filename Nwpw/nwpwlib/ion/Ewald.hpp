#ifndef	_EWALD_H_
#define _EWALD_H_


#include        "Parallel.hpp"
#include        "Lattice.hpp"
#include	"Ion.hpp"
#include	"Control2.hpp"
//#include	"Pseudopotential.hpp"

class	Ewald {
   int    encut,enx,eny,enz,enshl3d,enpack,enpack_all,enida;
   int    *i_indx,*j_indx,*k_indx;
   float *vg,*rcell,*eG,*vcx,*zv,*ss,*exi,*tmp3,*ftmp;
   float *ewx1,*ewy1,*ewz1;
   float unita[9],unitg[9],ercut,cewald,alpha;
   float eecut;

public:
   Parallel  *ewaldparall;
   Ion	     *ewaldion;
   Lattice  *ewaldlattice;

   /* Constructors */
   //Ewald(Parallel *, Ion *, Pseudopotential *);
   Ewald(Parallel *, Ion *, Lattice *, Control2&, float *);

   /* destructor */
   ~Ewald() {
            delete [] i_indx;
            delete [] j_indx;
            delete [] k_indx;
            delete [] vg;
            delete [] vcx;
            delete [] ss;
            delete [] exi;
            delete [] tmp3;
            delete [] ftmp;
            delete [] rcell;
            delete [] eG;
            delete [] zv;
            delete [] ewx1;
            delete [] ewy1;
            delete [] ewz1;
         }

    void phafac();
    int ncut()  {return encut;}
    int nida()  {return enida;}
    int npack() {return enpack;}
    int npack_all() {return enpack_all;}
    int nshl3d() {return enshl3d;}
    int nx() {return enx;}
    int ny() {return eny;}
    int nz() {return enz;}
    float ecut() {return eecut;}
    float rcut() {return ercut;}
    float mandelung() {return alpha;}
    float energy();
    void   force(float *);

};

#endif
