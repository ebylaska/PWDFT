/*
 $Id$
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include "typesf2c.h"
#include "get_word.h"

extern float tetercc(float xx);

extern float cpi_Splint();
extern void   cpi_Spline();


void teter_parse(debug_ptr,lmax_ptr,locp_ptr,rlocal_ptr,
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


    float   zatom,zion;      /* local psp parameters          */
    float over_fourpi;

    int      *nl;
    int      i,k,l,p,p1;
    int      Ngrid,nrl;
    float   *rgrid,*psi,*psp;
    float       *rl, *tmp, *tmp2, *sc_rho, *sc_rhol, *sc_drho, *sc_drhol,
    **psil,
    **pspl;
    float   drl,rmax;

    int      lmax,locp,lmaxp;
    float   r0,xx;
    int      n[10];
    int      pspdat,pspcode,pspxc;
    float   r2well,rcore[10],e99,e999;
    float   rchrg,fchrg,qchrg,pi;


    char   *w,*tc;
    FILE   *fp;

    char     comment[255];
    int      argc;



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

    pi          = 4.0*atan(1.0);
    over_fourpi = 1.0/(4.0*pi);



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







    /* Read TETER psp */
    fp = fopen(infile,"r+");
    w = get_word(fp);
    while ((w!=NIL) && (strcmp("<TETER>",w)!=0))
        w = get_word(fp);

    /* Error occured */
    if (w==NIL)
    {
        printf("Error: <TETER> section not found\n");
        fclose(fp);
        exit(99);
    }

    argc = to_eoln(fp);
    argc= get_line(fp,comment,255);

    fscanf(fp,"%f %f %d",&zatom,&zion,&pspdat);
    argc=to_eoln(fp);
    fscanf(fp,"%d %d %d %d %d %f",&pspcode,&pspxc,&lmax,&locp,&Ngrid,&r2well);
    lmaxp = lmax+1;
    argc=to_eoln(fp);


    for (p=0; p<=lmax; ++p)
    {
        fscanf(fp,"%d %f %f %d %f",&l,&e99,&e999,&(n[p]),&(rcore[p]));
        to_eoln(fp);
        to_eoln(fp);
    }
    fscanf(fp,"%f %f %f",&rchrg,&fchrg,&qchrg);




    psi     = (float *) malloc(Ngrid*sizeof(float));
    psp     = (float *) malloc(Ngrid*sizeof(float));
    rgrid   = (float *) malloc(Ngrid*sizeof(float));
    tmp     = (float *) malloc(Ngrid*sizeof(float));
    tmp2    = (float *) malloc(Ngrid*sizeof(float));
    sc_rho  = (float *) malloc(Ngrid*sizeof(float));
    sc_drho = (float *) malloc(Ngrid*sizeof(float));


    /* define Teter grid */
    for (i=0; i<Ngrid; ++i)
    {
        xx = ((float) i);
        xx=xx/((float) (Ngrid-1));
        xx = (xx+0.01);
        xx = xx*xx*xx*xx*xx;
        rgrid[i]=100.0*xx-1.0e-8;
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

    r0    = rgrid[280];
    rl[0] = rgrid[280];
    for (i=1; i<nrl; ++i)
    {
        rl[i] = drl*((float) i);
        xx = (rl[i] + 1.0e-8)/100.0;
        xx = pow(xx,0.2);
        xx = xx-0.01;
        xx = (Ngrid-1)*xx;
        nl[i] = rint(xx-0.5);
    }


    /* read in pseudopotentials */
    for (p=0; p<=lmax; ++p)
    {
        pspl[p] = (float *) malloc(nrl*sizeof(float));

        to_eoln(fp);
        to_eoln(fp);
        for (i=0; i<Ngrid; ++i)  fscanf(fp,"%f",&(psp[i]));


        cpi_Spline(rgrid,psp,Ngrid-4,0.0,0.0,tmp,tmp2);
        pspl[p][0] = psp[280];
        for (i=1; i<nrl; ++i)
        {
            pspl[p][i] = cpi_Splint(rgrid,psp,tmp,Ngrid-4,nl[i],rl[i]);
        }
    }

    /* read in wavefunctions */
    for (p=0; p<=lmax; ++p)
    {
        psil[p] = (float *) malloc(nrl*sizeof(float));

        to_eoln(fp);
        to_eoln(fp);
        for (i=0; i<Ngrid; ++i)  fscanf(fp,"%f",&(psi[i]));


        cpi_Spline(rgrid,psi,Ngrid-4,0.0,0.0,tmp,tmp2);
        psil[p][0] = psi[280];
        for (i=1; i<nrl; ++i)
        {
            psil[p][i] = cpi_Splint(rgrid,psi,tmp,Ngrid-4,nl[i],rl[i]);
        }
    }

    /* define semicore */
    if (rchrg>0.0)
    {
        for (i=0; i<Ngrid; ++i)
        {
            /*
              xx = rgrid[i]/(rchrg);
              gg=sin(2.0*pi*xx)/( (2.0*pi*xx)*(1.0-4.0*xx*xx)*(1.0-xx*xx) );
              gg=gg*gg;
            */
            xx = rgrid[i]/(rchrg);
            sc_rho[i] = 4*pi*fchrg*tetercc(xx);
        }

        cpi_Spline(rgrid,sc_rho,Ngrid-4,0.0,0.0,tmp,tmp2);
        sc_rhol[0] = sc_rho[280];
        for (i=1; i<nrl; ++i)
        {
            sc_rhol[i] = cpi_Splint(rgrid,sc_rho,tmp,Ngrid-4,nl[i],rl[i]);
        }

        /* define to be zero for now since it is not used */
        for (i=1; i<nrl; ++i)
        {
            sc_drhol[i] = 0.0;
        }

    }

    fclose(fp);

    free(tmp2);
    free(tmp);
    free(rgrid);
    free(psp);
    free(psi);
    free(sc_rho);
    free(sc_drho);
    free(nl);


    /* find the comment */
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



    /* write outfile */
    fp = fopen(outfile,"w+");
    fprintf(fp,"%s\n",atom_out);
    if (locp_out!=-1) locp=locp_out;
    fprintf(fp,"%f %f %d   %d %d %f\n",zion,0.0,lmax,lmax_out,locp,rlocal_out);
    for (p=0; p<=lmax; ++p)
        fprintf(fp,"%f ", rcore[p]);
    fprintf(fp,"\n");
    fprintf(fp,"%d %f\n",nrl,drl);
    fprintf(fp,"%s",comment);


    /* appending pseudopotentials */
    for (k=0; k<nrl; ++k)
    {
        fprintf(fp,"%14.9lf", rl[k]);
        for (p=0; p<=lmax; ++p)
            fprintf(fp," %12.8lf", pspl[p][k]);
        fprintf(fp,"\n");
    }

    for (p=0; p<=lmax; ++p) free(pspl[p]);
    free(pspl);

    /* appending pseudowavefunctions */
    for (k=0; k<nrl; ++k)
    {
        fprintf(fp,"%14.9lf", rl[k]);
        for (p=0; p<=lmax; ++p)
            fprintf(fp," %12.8lf", psil[p][k]);
        fprintf(fp,"\n");
    }
    for (p=0; p<=lmax; ++p) free(psil[p]);
    free(psil);

    /* append semicore corrections */
    if (rchrg != 0.0)
    {
        fprintf(fp,"%f\n",rchrg);
        for (k=0; k<nrl; ++k)
            fprintf(fp,"%14.9lf %12.8lf\n", rl[k],
                    fabs(sc_rhol[k]*over_fourpi));
        for (k=0; k<nrl; ++k)
            fprintf(fp,"%14.9lf %12.8lf\n", rl[k],
                    (sc_drhol[k]*over_fourpi));
    }
    free(sc_rhol);
    free(sc_drhol);
    free(rl);


    fclose(fp);


    if (debug)
    {
        printf("TETER pseudopotential Parameters\n\n");
        printf("atom : %s\n",atom_out);
        printf("Zatom= %f\n",zatom);
        printf("Zion = %f\n",zion);
        printf(" lmax= %d\n",lmax);
        printf(" locp= %d\n",locp);
        printf(" rlocal= %f\n\n",rlocal_out);
        printf(" rcrhg=%f  fchrg=%f  qchrg=%f\n",rchrg,fchrg,qchrg);
        printf("rcore: ");
        for (p=0; p<=lmax; ++p)
            printf("%f ", rcore[p]);
        printf("\n");
        printf(" nrl=%d drl=%f\n",nrl,drl);
        printf("comment:%s\n",comment);

        fflush(stdout);
    }

    /* free malloc memory */
    free(infile);
    free(outfile);
    free(full_filename);
    free(atom_out);

    return;

} /* main */


float tetercc(float xx)
{

    /*The c s are coefficients for Taylor expansion of the analytic form near xx=0, 1/2, and 1. */
    float   c21,c22,c23,c24;
    float   c31,c32,c33,c34;

    /*local variables */
    float pi,gg1cc,yy;

    pi = 4.0*atan(1.0);
    c21= 4.00/9.00;
    c22= -40.00/27.00;
    c23= 20.00/3.00-16.00*pi*pi/27.00;
    c24= -4160.00/243.00+160.00*pi*pi/81.00;
    c31= 1.00/36.00;
    c32= -25.00/108.00;
    c33= 485.00/432.00-pi*pi/27.00;
    c34=-4055.00/972.00+25.00*pi*pi/81.00;


    /* Cut off beyond 3/gcut=xcccrc */
    if (xx>3.000)
        gg1cc=0.000;

    /* Take care of difficult limits near x=0, 1/2, and 1 */
    else if (fabs(xx)<=1.e-9)
        gg1cc=1.00;

    else if (fabs(xx-0.500)<=1.0e-4)
        gg1cc=c21+(xx-0.500)*(c22+(xx-0.500)*(c23+(xx-0.500)*c24));

    else if (fabs(xx-1.00)<=1.e-04)
        gg1cc=c31+(xx-1.000)*(c32+(xx-1.000)*(c33+(xx-1.000)*c34));
    else
    {
        yy=sin(2.0*pi*xx)/( (2.0*pi*xx)*(1.0-4.0*xx*xx)*(1.0-xx*xx) );
        gg1cc=yy*yy;
    }

    return gg1cc;
}

