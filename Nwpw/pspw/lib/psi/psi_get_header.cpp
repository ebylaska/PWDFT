
#include	"Parallel.hpp"
#include	"compressed_io.hpp"
#include	"control.hpp"

void psi_get_header(Parallel *myparall,int *version, int nfft[], double unita[], int *ispin, int ne[])
{
   if (myparall->is_master())
   {
      //char *fname = control_input_movecs_filename();
      openfile(4,control_input_movecs_filename(),"r");
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
