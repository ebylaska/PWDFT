
#include <cmath>
#include <cstdlib>

//#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


namespace pwdft {

static void density_color(double v, double *r, double *g, double *b)
{
    *r = std::fmin(1.0,2.0*v);
    *g = 2.0*(1.0 - std::fabs(v-0.5));
    *b = std::fmin(1.0,2.0*(1.0-v));
}

void nwpw_density_volume_rgb(const char *fname,
                             const double *rho,
                             int nx, int ny, int nz,
                             int imgw, int imgh)
{
    size_t n = (size_t) nx*ny*nz;
    if(n == 0) return;

    double tilt = 0.3; // viewing tilt
    double cx = nx*0.5;
    double cy = ny*0.5;
    double scale = imgw*0.7;


    unsigned char *img = (unsigned char*) calloc(3*imgw*imgh,1);
    if(!img) return;

    double rmax = std::fabs(rho[0]);

    for(size_t i=0; i<n; ++i)
        if(std::fabs(rho[i]) > rmax)
            rmax = std::fabs(rho[i]);

    double denom = std::log(1.0 + rmax);
    if(denom==0.0) denom = 1.0;

    for(int k=1; k<nz; ++k)
    for(int j=0; j<ny; ++j)
    {

        size_t plane = (size_t)nx*ny*k;
        size_t row   = (size_t)nx*j;

        for(int i=0; i<nx; ++i)
        {
            size_t idx   = i + row + plane;

            //double r = log(1.0 + rho[idx]) / denom;
            double r = std::log(1.0 + std::fabs(rho[idx])) / denom;
            if(r<1e-4) continue;


            double kz = (double)k;
            double x = (i-cx)/kz;
            double y = (j-cy)/kz + tilt*(kz/(double)nz);

            int xp = (int)(imgw*0.5 + scale*x);
            int yp = (int)(imgh*0.5 + scale*y);

            if(xp<0 || xp>=imgw || yp<0 || yp>=imgh) continue;

            int p = xp + imgw*yp;

            double alpha = r*(1.0-(double)k/nz)*0.6;

            double cr,cg,cb;
            density_color(r,&cr,&cg,&cb);

            //img[3*p+0] = fmin(255,img[3*p+0] + 255*cr*alpha);
            //img[3*p+1] = fmin(255,img[3*p+1] + 255*cg*alpha);
            //img[3*p+2] = fmin(255,img[3*p+2] + 255*cb*alpha);

            img[3*p+0] = (unsigned char) std::fmin(255.0, img[3*p+0] + 255.0*cr*alpha);
            img[3*p+1] = (unsigned char) std::fmin(255.0, img[3*p+1] + 255.0*cg*alpha);
            img[3*p+2] = (unsigned char) std::fmin(255.0, img[3*p+2] + 255.0*cb*alpha);
        }
    }

    stbi_write_png(fname,imgw,imgh,3,img,3*imgw);

    free(img);
}
  
}
