

/*
#include        <cstdio>
#include        <stdio.h>
#include        <cmath>
#include        <cstdlib>
using namespace std;
*/

#include        <iostream>
#include	"compressed_io.hpp"
//#include	"control.hpp"

#include	"Parallel.hpp"
#include	"Pneb.hpp"
#include	"psi.hpp"

/*****************************************************
 *                                                   *
 *                psi_get_header                     *
 *                                                   *
 *****************************************************/

void psi_get_header(Parallel *myparall,int *version, int nfft[], double unita[], int *ispin, int ne[], char *filename)
{
   if (myparall->is_master())
   {
      //char *fname = control_input_movecs_filename();
      openfile(4,filename,"r");
      iread(4,version,1);
      iread(4,nfft,3);
      dread(4,unita,9);
      iread(4,ispin,1);
      iread(4,ne,2);
      closefile(4);
   }
   myparall->Brdcst_iValue(0,0,version);
   myparall->Brdcst_iValues(0,0,3,nfft);
   myparall->Brdcst_Values(0,0,9,unita);
   myparall->Brdcst_iValue(0,0,ispin);
   myparall->Brdcst_iValues(0,0,2,ne);
}



/*****************************************************
 *                                                   *
 *                psi_read0                          *
 *                                                   *
 *****************************************************/
/* 
   Just reads psi and its header.  

   Note - the file must exist
  
*/
void psi_read0(Pneb *mypneb,int *version, int nfft[], 
               double unita[], int *ispin, int ne[],
               double *psi, char *filename)
{
   int occupation;

   Parallel *myparall = mypneb->d3db::parall;

   if (myparall->is_master())
   {
      openfile(4,filename,"r");
      iread(4,version,1);
      iread(4,nfft,3);
      dread(4,unita,9);
      iread(4,ispin,1);
      iread(4,ne,2);
      iread(4,&occupation,1);
   }
   myparall->Brdcst_iValue(0,0,version);
   myparall->Brdcst_iValues(0,0,3,nfft);
   myparall->Brdcst_Values(0,0,9,unita);
   myparall->Brdcst_iValue(0,0,ispin);
   myparall->Brdcst_iValues(0,0,2,ne);

   mypneb->g_read(4,psi);

   if (myparall->is_master()) closefile(4);
}

/*****************************************************
 *                                                   *
 *                psi_read                           *
 *                                                   *
 *****************************************************/
/* 
   Reads psi and check the header
*/
void psi_read(Pneb *mypneb, char *filename, double *psi2)
{
   int version,ispin,nfft[3],ne[2];
   double unita[9];

   /* read psi from file if psi_exist */
   if (psi_filefind(mypneb,filename))
   {
      std::cout << " input psi exists, reading from file: " << filename << std::endl;
      psi_read0(mypneb,&version,nfft,unita,&ispin,ne,psi2,filename);
   }

   /* generate new psi */
   else
   {
      std::cout << " generating random psi from scratch" << std::endl;
      mypneb->g_generate_random(psi2);
   }
}

/*****************************************************
 *                                                   *
 *                psi_write                          *
 *                                                   *
 *****************************************************/
void psi_write(Pneb *mypneb,int *version, int nfft[],
              double unita[], int *ispin, int ne[],
              double *psi, char *filename)
{  
   int occupation = -1;
   
   Parallel *myparall = mypneb->d3db::parall;
   
   if (myparall->is_master())
   {  
      openfile(6,filename,"w");
      iwrite(6,version,1);
      iwrite(6,nfft,3);
      dwrite(6,unita,9);
      iwrite(6,ispin,1);
      iwrite(6,ne,2);
      iwrite(6,&occupation,1);
   }
   
   mypneb->g_write(6,psi);
   
   if (myparall->is_master()) closefile(6);
}

/*****************************************************
 *                                                   *
 *                psi_filefind                       *
 *                                                   *
 *****************************************************/
bool psi_filefind(Pneb *mypneb, char *filename)
{
   int ifound = 0;
   Parallel *myparall = mypneb->d3db::parall;
   
   if (myparall->is_master())
   {  
      ifound = cfileexists(filename);
   }

   myparall->Brdcst_iValue(0,0,&ifound);

   return (ifound>0);
}




/*
void v_psi_read(Pneb *mypneb,int *version, int nfft[],
              double unita[], int *ispin, int ne[],
              double *psi)
{
   int occupation;

   Parallel *myparall = mypneb->d3db::parall;

   if (myparall->is_master())
   {
      openfile(4,control_input_v_movecs_filename(),"r");
      iread(4,version,1);
      iread(4,nfft,3);
      dread(4,unita,9);
      iread(4,ispin,1);
      iread(4,ne,2);
      iread(4,&occupation,1);
   }
   myparall->Brdcst_iValue(0,0,version);
   myparall->Brdcst_iValues(0,0,3,nfft);
   myparall->Brdcst_Values(0,0,9,unita);
   myparall->Brdcst_iValue(0,0,ispin);
   myparall->Brdcst_iValues(0,0,2,ne);

   mypneb->g_read(4,psi);

   if (myparall->is_master()) closefile(4);
}
*/


/*
void v_psi_write(Pneb *mypneb,int *version, int nfft[],
              double unita[], int *ispin, int ne[],
              double *psi)
{     
   int occupation = -1;
      
   Parallel *myparall = mypneb->d3db::parall;
      
   if (myparall->is_master())
   {  
      openfile(6,control_output_v_movecs_filename(),"w");
      iwrite(6,version,1);
      iwrite(6,nfft,3);
      dwrite(6,unita,9);
      iwrite(6,ispin,1);
      iwrite(6,ne,2);
      iwrite(6,&occupation,1);
   }
  
   mypneb->g_write(6,psi);
   
   if (myparall->is_master()) closefile(6);
}             

*/

