/* rtdb.C
   Author - Eric Bylaska

*/


#include        <iostream>
#include        <cstdlib>



#include	"Int64.h"
extern "C" {
#include	"rtdb_seq.h"
}
#include	"rtdb.hpp"

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
RTDB::RTDB(Parallel *inparall, const char *filename, const char *mode)
{
   parall        = inparall;
   parallel_mode = 1;
   if (parall->is_master())
   {
      if (!rtdb_seq_open(filename,mode, &handle))
      {
         cout << "error opening " << filename << " mode=" << mode << "\n";
         return;
      }
   }
   parall->Brdcst_iValue(0,0,&handle);

}

static int matypesize(const int intype)
{
  switch (intype) {
  case rtdb_char:       /* char */
    return sizeof(char); break;
  case rtdb_int:        /* int */
    return sizeof(Int64); break;
  case rtdb_log:       /* log */
    return sizeof(int); break;
  case rtdb_long:       /* log */
    return sizeof(long); break;
  case rtdb_float:      /* float */
    return sizeof(float); break;
  case rtdb_double:     /* double */
    return sizeof(double); break;
  default:
    return sizeof(char); break;
  }
}


/********************************
 *                              *
 *       RTDB::get              *
 *                              *
 ********************************/
int RTDB::get(const char *tag, const int matype, const int nelem, void *array)
{
   int status;

   if (parall->is_master())
   {
      status = rtdb_seq_get(handle,tag,matype,nelem,array);
   }
   if (parallel_mode)
   {
      parall->Brdcst_iValue(0,0,&status);
      if (status)
      {
         int sz = nelem*matypesize(matype);
         parall->Brdcst_cValues(0,0,sz,array);
      }
   }
   return status;
}

/********************************
 *                              *
 *       RTDB::put              *
 *                              *
 ********************************/
int RTDB::put(const char *tag, const int matype, const int nelem, void *array)
{
   int status;

   if (parall->is_master())
   {
      status = rtdb_seq_put(handle,tag,matype,nelem,array);
      if (!status)
         cout << "rtdb error putting tag =" << tag << "\n";
   }
   if (parallel_mode)
      parall->Brdcst_iValue(0,0,&status);
   return status;
}

int RTDB::get_info(const char *tag, int *matype, int *nelem, char *date)
{
   int status;

   if (parall->is_master())
   {
      status = rtdb_seq_get_info(handle,tag,matype,nelem,date);
      if (!status)
         cout << "rtdb error get_info tag =" << tag << "\n";
   }
   if (parallel_mode)
   {
      parall->Brdcst_iValue(0,0,matype);
      parall->Brdcst_iValue(0,0,nelem);
      parall->Brdcst_iValue(0,0,&status);
   }
   return status;
}
