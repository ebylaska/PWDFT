#ifndef _MAPPING_H_
#define _MAPPING_H_
/* mapping.h
   Author - Eric Bylaska

	this class is used defining 3d parallel maps
*/

class Mapping3 {

        int *qmap[3],*pmap[3],*kmap;

public:
        int np,taskid;
        int maptype,n2ft3d,nfft3d;
        int nx,ny,nz;
        int nq,nq1,nq2,nq3;

	/* Constructors */
	Mapping3();
	Mapping3(const int, const int, const int, const int, const int, const int);

        /* copy constructor */
        //Mapping3(const Mapping3&);

        /* destructor */
	~Mapping3();

        int ijktoindex(const int i, const int j, const int k) {
           if (maptype==1)
              return (i + j*(nx/2+1) + qmap[0][k]*(nx/2+1)*ny);
           else
              return (k + qmap[2][i+j*(nx/2+1)]*nz);
        }
        int ijktoindex1(const int i, const int j, const int k) {
           if (maptype==1)
              return (i + k*(nx/2+1) + qmap[0][j]*(nx/2+1)*nz);
           else
              return (j + qmap[1][k+i*nz]*ny);
        }
        int ijktoindex2(const int i, const int j, const int k) {
           if (maptype==1)
              return (i + k*(nx+2) + qmap[0][j]*(nx+2)*ny);
           else
              return (i + qmap[0][j+k*ny]*(nx+2));
        }
        int ijktoindex2t(const int i, const int j, const int k) {
           if (maptype==1)
              return (i + k*(nx+2) + qmap[0][j]*(nx/2+1)*ny);
           else
              return (i + qmap[0][j+k*ny]*(nx/2+1));
        }

        int ijktop(const int i, const int j, const int k) {
           if (maptype==1)
              return (pmap[0][k]);
           else
              return (pmap[2][i+j*(nx/2+1)]);
        }
        int ijktop1(const int i, const int j, const int k) {
           if (maptype==1)
              return (pmap[0][j]);
           else
              return (pmap[1][k+i*nz]);
        }
        int ijktop2(const int i, const int j, const int k) {
           if (maptype==1)
              return (pmap[0][j]);
           else
              return (pmap[0][j+k*ny]);
        }

        int ijktoq(const int i, const int j, const int k) {
           if (maptype==1)
              return (qmap[0][k]);
           else
              return (qmap[2][i+j*(nx/2+1)]);
        }

        int ijktoq1(const int i, const int j, const int k) {
           if (maptype==1)
              return (qmap[0][j]);
           else
              return (qmap[1][k+i*nz]);
        }

        int ijktoq2(const int i, const int j, const int k) {
           if (maptype==1)
              return (qmap[0][j]);
           else
              return (qmap[0][j+k*ny]);
        }

};

#endif
