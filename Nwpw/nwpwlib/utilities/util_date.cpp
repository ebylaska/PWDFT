/*$Id: util_date.c,v 1.8 2003-08-13 18:06:11 edo Exp $*/

#include <sys/types.h>
#include <ctime>
#include <cstring>


namespace pwdft {


char *util_date()
{
  time_t t = time((time_t *) 0);
  char *tmp = ctime(&t);
  int nlen = strlen(tmp);
  tmp[nlen-1] = '\0';
  return tmp;
}



#define Clock   clock()

void    seconds(double *tt)
{
   *tt = Clock/((double) CLOCKS_PER_SEC);
}

}
