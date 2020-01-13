#ifndef _MAPPING1d_H_
#define _MAPPING1d_H_
/* mapping.h
   Author - Eric Bylaska

	this class is used defining 1d parallel maps
*/
using namespace std;

class Mapping1 {

        int np,taskid;

        int *qmap[2],*pmap[2],*kmap[2],*nqarray[2];

public:
        int ispin,ne[2];
        int maptype1,neq[2];

	/* Constructors */
	Mapping1();
	Mapping1(const int, const int, const int, const int, int *);

        /* destructor */
	~Mapping1() {
           for (int ms=0; ms<ispin; ++ms)
           {
              delete [] qmap[ms];
              delete [] pmap[ms];
              delete [] kmap[ms];
              delete [] nqarray[ms];
           }
        }

        int msntoindex(const int ms, const int n) { return qmap[ms][n]; }
        int msntop(const int ms, const int n)     { return pmap[ms][n]; }
};

#endif
