/*
 $Id$
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include "typesf2c.h"
#include "get_word.h"

extern float cpi_Splint();
extern void   cpi_Spline();


void cpi_parse(debug_ptr,lmax_ptr,locp_ptr,rlocal_ptr,
               sdir_name,n9,dir_name,n0,in_filename,n1,out_filename,n2,atom,n3)
int	*debug_ptr;
int	*lmax_ptr;
int	*locp_ptr;
float 	*rlocal_ptr;
char	sdir_name[];
int	*n9;
char	dir_name[];
int	*n0;
char	in_filename[];
int	*n1;
char	out_filename[];
int	*n2;
char	atom[];
int	*n3;
{

    int      debug;
    int      lmax_out,locp_out;
    float   rlocal_out;

    int      lmax;


    float   Zion;      /* local psp parameters          */
    float over_fourpi;

    int      *nl;
    int      i,j,k,p,p1;
    int      Ngrid,nrl;
    float   *rgrid,*psi,*psp;
    float       *rl, *tmp, *tmp2, *sc_rho, *sc_rhol, *sc_drho, *sc_drhol,
    **psil,
    **pspl;
    float   r,ul,vl,amesh,al,drl,r_semicore,rmax;

    int      lmaxp;
    float   dum1,dum2,dum3,dum4,r0;
    int idum;


    char   *w,*tc;
    FILE   *fp;

    char     comment[255];
    int      argc,value;



    int		m9 = ((int) (*n9));
    int		m0 = ((int) (*n0));
    int		m1 = ((int) (*n1));
    int		m2 = ((int) (*n2));
    int		m3 = ((int) (*n3));
    char *infile  = (char *) malloc(m9+m1+5);
    char *outfile = (char *) malloc(m0+m2+5);
    char *atom_out = (char *) malloc(m3+5);

    char *full_filename = (char *) malloc(m9+25+5);


    debug = *debug_ptr;
    lmax_out   = *lmax_ptr;
    locp_out   = *locp_ptr;
    rlocal_out = *rlocal_ptr;

    (void) strncpy(infile, sdir_name, m9);
    infile[m9] = '\0';
    strcat(infile,"/");
    infile[m9+1] = '\0';
    strncat(infile,in_filename,m1);
    infile[m9+m1+1] = '\0';

    (void) strncpy(outfile, dir_name, m0);
    outfile[m0] = '\0';
    (void) strcat(outfile,"/");
    outfile[m0+1] = '\0';
    (void) strncat(outfile,out_filename,m2);
    outfile[m0+m2+1] = '\0';

    (void) strncpy(atom_out, atom, m3);
    atom_out[m3] = '\0';

    over_fourpi = 1.0/(16.0*atan(1.0));


    /* find the comment */
    strcpy(comment,"CPI formatted  pseudopotential");
    fp = fopen(infile,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<comment>",w)!=0))
        w = get_word(fp);

    if (w!=NIL)
    {
        w = get_word(fp);
        p  = 0;
        tc = comment;
        while ((w!=NIL)&&(strcmp("<end>",w) != 0))
        {
            p = (strlen(w));
            strcpy(tc, w);
            for (p1=0;p1<p; ++p1) ++tc;
            strcpy(tc, " ");
            ++tc;

            w = get_word(fp);
        }
    }
    fclose(fp);


    /* define linear grid */
    nrl  = 2001;
    rmax = 40.0;
    drl  = rmax/((float) (nrl-1));

    fp = fopen(infile,"r+");
    w = get_word(fp);
    while ((w != ((char *) EOF)) && (strcmp("<linear>",w) != 0))
        w = get_word(fp);
    if (w!=((char *) EOF))
    {
        fscanf(fp,"%d %f",&nrl,&drl);
        rmax = ((float) (nrl-1))*drl;
    }
    fclose(fp);




    /* Read CPI psp */
    fp = fopen(infile,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<CPI>",w)!=0))
        w = get_word(fp);

    /* Error occured */
    if (w==NIL)
    {
        printf("Error: <CPI> section not found\n");
        fclose(fp);
        exit(99);
    }

    argc = to_eoln(fp);


    fscanf(fp,"%f %d",&Zion,&lmaxp);
    lmax = lmaxp-1;

    fscanf(fp,"%f %f %f %f",&dum1,&dum2,&dum3,&dum4);
    fscanf(fp,"%f %f %f ",&dum1,&dum2,&dum3);
    fscanf(fp,"%f %f %f ",&dum1,&dum2,&dum3);
    fscanf(fp,"%f %f %f ",&dum1,&dum2,&dum3);
    fscanf(fp,"%f %f %f ",&dum1,&dum2,&dum3);
    fscanf(fp,"%f %f %f ",&dum1,&dum2,&dum3);
    fscanf(fp,"%f %f %f ",&dum1,&dum2,&dum3);
    fscanf(fp,"%f %f %f ",&dum1,&dum2,&dum3);
    fscanf(fp,"%f %f %f ",&dum1,&dum2,&dum3);
    fscanf(fp,"%f %f %f ",&dum1,&dum2,&dum3);

    fscanf(fp,"%d %f",&Ngrid,&amesh);
    al = log(amesh);

    psi     = (float *) malloc(Ngrid*sizeof(float));
    psp     = (float *) malloc(Ngrid*sizeof(float));
    rgrid   = (float *) malloc(Ngrid*sizeof(float));
    tmp     = (float *) malloc(Ngrid*sizeof(float));
    tmp2    = (float *) malloc(Ngrid*sizeof(float));
    sc_rho  = (float *) malloc(Ngrid*sizeof(float));
    sc_drho = (float *) malloc(Ngrid*sizeof(float));

    for (i=0; i<Ngrid; ++i)
    {
        fscanf(fp,"%d %f %f %f",&j, &r,&ul,&vl);
        rgrid[i]  = r;
        psi[i] = ul;
        psp[i] = vl;
    }


    /* check linear grid and redefine if necessary */
    if (rmax > rgrid[Ngrid-5])
    {
        rmax = rgrid[Ngrid-5];
        drl = rmax/((float) (nrl-1));
    }



    /* generate linear meshes */
    rl       = (float *) malloc(nrl*sizeof(float));
    nl       = (int *)    malloc(nrl*sizeof(int));
    psil     = (float **) malloc(lmaxp*sizeof(float*));
    pspl     = (float **) malloc(lmaxp*sizeof(float*));
    sc_rhol  = (float *) malloc(nrl*sizeof(float));
    sc_drhol = (float *) malloc(nrl*sizeof(float));

    r0    = rgrid[0];
    rl[0] = rgrid[0];
    for (i=1; i<nrl; ++i)
    {
        rl[i] = drl*((float) i);
        nl[i] = rint(log(rl[i]/r0)/al -0.5);
    }


    psil[0] = (float *) malloc(nrl*sizeof(float));
    pspl[0] = (float *) malloc(nrl*sizeof(float));

    cpi_Spline(rgrid,psp,Ngrid-4,0.0,0.0,tmp,tmp2);
    pspl[0][0] = psp[0];
    for (i=1; i<nrl; ++i)
    {
        pspl[0][i] = cpi_Splint(rgrid,psp,tmp,Ngrid-4,nl[i],rl[i]);
    }

    cpi_Spline(rgrid,psi,Ngrid-4,0.0,0.0,tmp,tmp2);
    psil[0][0] = psi[0];
    for (i=1; i<nrl; ++i)
    {
        psil[0][i] = cpi_Splint(rgrid,psi,tmp,Ngrid-4,nl[i],rl[i]);
    }

    for (p=1; p<lmaxp; ++p)
    {
        fscanf(fp,"%d %f",&idum,&dum1);

        for (i=0; i<Ngrid; ++i)
        {
            fscanf(fp,"%d %f %f %f",&j, &r,&ul,&vl);
            rgrid[i]  = r;
            psi[i] = ul;
            psp[i] = vl;
        }

        psil[p] = (float *) malloc(nrl*sizeof(float));
        pspl[p] = (float *) malloc(nrl*sizeof(float));

        cpi_Spline(rgrid,psp,Ngrid-4,0.0,0.0,tmp,tmp2);
        pspl[p][0] = psp[0];
        for (i=1; i<nrl; ++i)
        {
            pspl[p][i] = cpi_Splint(rgrid,psp,tmp,Ngrid-4,nl[i],rl[i]);
        }

        cpi_Spline(rgrid,psi,Ngrid-4,0.0,0.0,tmp,tmp2);
        psil[p][0] = psi[0];
        for (i=1; i<nrl; ++i)
        {
            psil[p][i] = cpi_Splint(rgrid,psi,tmp,Ngrid-4,nl[i],rl[i]);
        }


    }

    /* read semi-core */
    r_semicore = 0.0;
    value = fscanf(fp,"%f   %f %f %f", &r,&ul,&vl,&dum1);
    if (value!=EOF)
    {
        r_semicore =  99.99; /* not known?? */
        rgrid[0]  = r; sc_rho[0] = ul; sc_drho[0] = vl;
        for (i=1; i<Ngrid; ++i)
        {
            fscanf(fp,"%f   %f %f %f", &r,&ul,&vl,&dum1);
            rgrid[i]   = r;
            sc_rho[i]  = ul;
            sc_drho[i] = vl;
        }

        cpi_Spline(rgrid,sc_rho,Ngrid-4,0.0,0.0,tmp,tmp2);
        sc_rhol[0] = sc_rho[0];
        for (i=1; i<nrl; ++i)
        {
            sc_rhol[i] = cpi_Splint(rgrid,sc_rho,tmp,Ngrid-4,nl[i],rl[i]);
        }

        cpi_Spline(rgrid,sc_drho,Ngrid-4,0.0,0.0,tmp,tmp2);
        sc_drhol[0] = sc_drho[0];
        for (i=1; i<nrl; ++i)
        {
            sc_drhol[i] = cpi_Splint(rgrid,sc_drho,tmp,Ngrid-4,nl[i],rl[i]);
        }


    }
    free(nl);
    free(rgrid);
    free(psi);
    free(psp);
    free(tmp);
    free(tmp2);
    free(sc_rho);
    free(sc_drho);



    fclose(fp);


    /* write outfile */
    fp = fopen(outfile,"w+");
    fprintf(fp,"%s\n",atom_out);
    fprintf(fp,"%f %f %d   %d %d %f\n",Zion,0.0,lmax,lmax_out,locp_out,rlocal_out);
    for (p=0; p<=lmax; ++p)
        fprintf(fp,"%f ", -1.0);
    fprintf(fp,"\n");
    fprintf(fp,"%d %f\n",nrl,drl);
    fprintf(fp,"%s\n",comment);

    /* appending pseudopotentials */
    for (k=0; k<nrl; ++k)
    {
        fprintf(fp,"%12.8lf", rl[k]);
        for (p=0; p<=lmax; ++p)
            fprintf(fp," %12.8lf", pspl[p][k]);
        fprintf(fp,"\n");
    }
    for (p=0; p<=lmax; ++p) free(pspl[p]);
    free(pspl);

    /* appending pseudowavefunctions */
    for (k=0; k<nrl; ++k)
    {
        fprintf(fp,"%12.8lf", rl[k]);
        for (p=0; p<=lmax; ++p)
            fprintf(fp," %12.8lf", psil[p][k]);
        fprintf(fp,"\n");
    }
    for (p=0; p<=lmax; ++p) free(psil[p]);
    free(psil);


    /* append semicore corrections */
    if (r_semicore != 0.0)
    {
        fprintf(fp,"%f\n",r_semicore);
        for (k=0; k<nrl; ++k)
            fprintf(fp,"%12.8lf %12.8lf\n", rl[k],
                    fabs(sc_rhol[k]*over_fourpi));
        for (k=0; k<nrl; ++k)
            fprintf(fp,"%12.8lf %12.8lf\n", rl[k],
                    (sc_drhol[k]*over_fourpi));
    }
    free(sc_drhol);
    free(sc_rhol);
    free(rl);


    fclose(fp);


    if (debug)
    {
        printf("CPI pseudopotential Parameters\n\n");
        printf("atom : %s\n",atom);
        printf("Zion : %f\n",Zion);
        printf(" lmax: %d\n",lmax);
        printf(" locp: %d\n",locp_out);
        printf(" rlocal: %f\n\n",rlocal_out);
        printf(" r_semicore: %f\n",r_semicore);

    }

    /* free malloc memory */
    free(infile);
    free(outfile);
    free(full_filename);
    free(atom_out);

    fflush(stdout);
    return;

} /* main */



/********************************
 *				*
 *	     cpi_Spline		*
 *				*
 ********************************/

void cpi_Spline(x,y,n,yp1,ypn,y2,u)
float 	x[],
y[];
int	n;
float	yp1;
float	ypn;
float	y2[];
float	u[];
{
    int	i,k;
    float sig,qn,un,p;

    if (yp1 > 0.99e30)
    {
        y2[0] = 0.0;
        u[0]  = 0.0;
    }
    else
    {
        y2[0] = -0.5;
        u[0] = 3.0/(x[1]-x[0]) * ((y[1]-y[0])/(x[1]-x[0]) - yp1);
    }

    for (i=1; i<(n-1); ++i)
    {
        sig = (x[i]-x[i-1])/(x[i+1] - x[i-1]);
        p   = sig*y2[i-1] + 2.0;
        y2[i] = (sig-1.0)/p;
        u[i] = ( 6.0 *
                 ((y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]))
                 /(x[i+1]-x[i-1])
                 - sig*u[i-1]
               ) / p;
    }

    if (ypn > 0.99e30)
    {
        qn = 0.0;
        un = 0.0;
    }
    else
    {
        qn = 0.5;
        un = 3.0/(x[n-1]-x[n-2]) * (ypn - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }

    y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2] + 1.0);
    for (k=n-2; k>=0; --k)
        y2[k] = y2[k]*y2[k+1] + u[k];


} /* cpi_Spline */

/********************************
 *				*
 *	     cpi_Splint		*
 *				*
 ********************************/


float	cpi_Splint(xa,ya,y2a,n,nx,x)
float	xa[];
float	ya[];
float	y2a[];
int	n;
int	nx;
float	x;
{
    int khi,klo;
    float h,a,b;
    float y;

    khi = nx+1;
    klo = nx;

    while ( (xa[klo] > x) || ( xa[khi] < x))
    {
        /*
              printf("Error in Splint ");
              printf("%d ->  %le %le %le",klo,x,xa[klo],xa[khi]);
        */
        if (xa[klo] > x)
        {
            --klo;
            --khi;
            /*
                     printf("   <\n");
            */
        }
        if (xa[khi] < x)
        {
            ++klo;
            ++khi;
            /*
                     printf("   >\n");
            */
        }
    }
    h = xa[khi] - xa[klo];
    a = (xa[khi] - x)/h;
    b = (x - xa[klo])/h;
    y = a*ya[klo] + b*ya[khi]
        + ( (a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi] ) * (h*h)/6.0;

    return y;

} /* cpi_Splint */


