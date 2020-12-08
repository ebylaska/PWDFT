#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void Generate_projectors(const int ii, const int ng0, const int nprj, const int nprjall, 
				  const int nx, const int ny, const int nz,
                                  const __global int    *indxi,
                                  const __global int    *indxj,
                                  const __global int    *indxk,
                                  const __global float *phfacx,
                                  const __global float *phfacy,
                                  const __global float *phfacz,
                                  const __global int    *sdfunction, 
                                  const __global float *vnl,
                                  __global float *prj) {
   int shftx = 2*ii*nx;
   int shfty = 2*ii*ny;
   int shftz = 2*ii*ny;
   int ng    = 2*ng0;

   int i = get_global_id(0); //ng
   //int l = get_global_id(1); //nprj
   
   float ai = phfacx[shftx+2*indxi[i]]; float bi = phfacx[shftx+2*indxi[i]+1];
   float aj = phfacy[shfty+2*indxj[i]]; float bj = phfacy[shfty+2*indxj[i]+1];
   float ak = phfacz[shftz+2*indxk[i]]; float bk = phfacz[shftz+2*indxk[i]+1];
   float c  = aj*ak - bj*bk;
   float d  = aj*bk + ak*bj;
   float rexi = (ai*c - bi*d);
   float iexi = (ai*d + bi*c);

   for (int l=0; l<nprj; ++l)
   {
      if (sdfunction[l])
      {
         prj[2*i   + (l+nprjall)*ng] = rexi * vnl[i+l*ng0];
         prj[2*i+1 + (l+nprjall)*ng] = iexi * vnl[i+l*ng0];
      } else { 
         prj[2*i   + (l+nprjall)*ng] = -iexi * vnl[i+l*ng0];
         prj[2*i+1 + (l+nprjall)*ng] =  rexi * vnl[i+l*ng0];
      }
   }
}

