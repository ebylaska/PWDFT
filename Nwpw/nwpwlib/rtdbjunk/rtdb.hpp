#ifndef _RTDB_H_
#define _RTDB_H_
/* rtdb.h
   Author - Eric Bylaska
        this class is used defining nd parallel geometries
*/

#include "Parallel.hpp"

namespace pwdft {

/* datatypes */
#ifndef _RTDB_DATATYPE_H_
#define _RTDB_DATATYPE_H_
#define rtdb_char 0              /* char */
#define rtdb_long 2              /* long int */
#define rtdb_byte 1009           /* byte */
#define rtdb_int 1010            /* int */
#define rtdb_log 1011            /* log */
#define rtdb_float 1012          /* float */
#define rtdb_double 1013         /* double */
#define rtdb_complex 1014        /* complex*/
#define rtdb_double_complex 1015 /* double complex*/
#endif

class RTDB {

  int handle;

public:
  Parallel *parall;
  int parallel_mode;

  /* Constructors */
  RTDB(Parallel *, const char *, const char *);

  /* destructor */
  //~RTDB();

  int get(const char *, const int, const int, void *);
  int put(const char *, const int, const int, void *);
  int get_info(const char *, int *, int *, char *);
};
} // namespace pwdft
#endif
