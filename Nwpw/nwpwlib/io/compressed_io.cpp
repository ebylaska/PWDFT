/*
 $Id: compressed_io.c,v 1.3 2004-05-05 02:20:33 edo Exp $
*************************************************************************
*									*
* compressed_io.c							*
*									*
* These routines implement a simple input/output package for writing to	*
* and reading from gzipped files.					*
*                                                                       *
*									*
* Fortran usage:							*
*	call openfile(unit,'file', 'X', n)	X is either 'r' or 'w'	*
*					n is the length of the filename	*
*	call cwrite(unit,c, n)		c is n element character array	*
*	call iwrite(unit,i, n)		i is n element integer array	*
*	call dwrite(unit,d, n)		d is n element double array	*
*	call cread(unit,c, n)		c is n element character array	*
*	call iread(unit,i, n)		i is n element integer array	*
*	call dread(unit,d, n)		d is n element double array	*
*	call closefile(unit)		close the datafile		*
*									*
* Author:      Scott Kohn (skohn@chem.ucsd.edu)				*
* modified by: Eric Bylaska (ebylaska@chem.ucsd.edu)			*
*									*
*************************************************************************
*/

//#include        <stdlib.h>
//#include        <stdio.h>
//#include        <math.h>
//#include        <string>
//

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "Int64.h"

/*
 * Check if a file exist using stat() function
 * return 1 if the file exist otherwise return 0
 */
#include <sys/stat.h>

namespace pwdft {

int cfileexists(const char *filename) {
  struct stat buffer;
  int exist = stat(filename, &buffer);
  if (exist == 0)
    return 1;
  else // -1
    return 0;
}

#define MAX_UNIT 10

static FILE *fd[MAX_UNIT]; /* the file descriptor of the pipe */

#define BAIL(X)                                                                \
  {                                                                            \
    fprintf(stderr, X);                                                        \
    exit(-1);                                                                  \
  }


/*
const size_t BUFFER_SIZE = 4096; // 4 KB buffer
static char buffer[MAX_UNIT][BUFFER_SIZE];
static size_t bufferIndex[MAX_UNIT] = {0};
static bool bufferwrite[MAX_UNIT] = {false};;
*/

// Use IIFE to initialize readbuffer1 and readbuffer2
/*
static char* readbuffer1[MAX_UNIT] = []() -> char* [MAX_UNIT] {
    static char* temp[MAX_UNIT];
    for (size_t i = 0; i < MAX_UNIT; ++i) {
        temp[i] = buffer[i];
    }
    return temp;
}();

static char* readbuffer2[MAX_UNIT] = []() -> char* [MAX_UNIT] {
    static char* temp[MAX_UNIT];
    for (size_t i = 0; i < MAX_UNIT; ++i) {
        temp[i] = &buffer[i][BUFFER_SIZE/2];
    }
    return temp;
}();

*/


/*
static size_t refillBuffer(const int unit)
{
   size_t bytesRead = fread(buffer[unit], 1, BUFFER_SIZE, fd[unit]);
   if (bytesRead == 0 && ferror(fd[unit])) 
      BAIL("Failed to read from file\n");

   return bytesRead;
}


static void flushBufferWrite(const int unit) 
{
   size_t bytesWritten = fwrite(buffer[unit], 1, bufferIndex[unit], fd[unit]);
   // Handle error
   if (bytesWritten < bufferIndex[unit]) {
       BAIL("Failed to write BUFFER to file\n");
   }
}

*/



#define BUFFERED_WRITE(unit, data_ptr, elem_size, num_elems)                  \
    do {                                                                      \
        size_t dataSize = (num_elems) * (elem_size);                          \
                                                                              \
        if (dataSize > BUFFER_SIZE)                                           \
        {                                                                     \
            if (bufferIndex[unit] > 0)                                        \
            {                                                                 \
                flushBufferWrite(unit);                                       \
                bufferIndex[unit] = 0;                                        \
            }                                                                 \
            (void)fwrite((data_ptr), (elem_size), (num_elems), fd[unit]);     \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            if (bufferIndex[unit] + dataSize > BUFFER_SIZE)                   \
            {                                                                 \
                flushBufferWrite(unit);                                       \
                bufferIndex[unit] = 0;                                        \
            }                                                                 \
            std::memcpy(buffer[unit] + bufferIndex[unit], (data_ptr), dataSize); \
            bufferIndex[unit] += dataSize;                                    \
        }                                                                     \
    } while (0)




/*
*************************************************************************
*									*
* Define the Xwrite and Xread routines using the Fortran bindings.	*
*									*
*************************************************************************
*/

void cread(int unit, char *c, const int n) 
{
   //BUFFERED_READ(unit, c, sizeof(char), n);
   (void)fread(c, sizeof(char), n, fd[unit]);

}

void iread(const int unit, int *i, const int n) 
{
   // Int64 *itmp = (Int64 *) malloc((n+1)*sizeof(Int64));
   Int64 *itmp;
   itmp = new Int64[n];

   //BUFFERED_READ(unit, itmp, sizeof(Int64), n);
   (void)fread(itmp, sizeof(Int64), n, fd[unit]);

   for (int j = 0; j < n; ++j)
     i[j] = (int)itmp[j];
   // free(itmp);
   delete[] itmp;
}

void dread(const int unit, double *d, const int n) 
{
   //BUFFERED_READ(unit, d, sizeof(double), n);
   (void)fread(d, sizeof(double), n, fd[unit]);
}


void cwrite(int unit, char *c, const int n) 
{
   //BUFFERED_WRITE(unit, c, sizeof(char), n);
   (void)fwrite(c, sizeof(char), n, fd[unit]);
}

void iwrite(const int unit, const int *i, const int n) 
{
   // Int64 *itmp = (Int64 *) malloc((n+1)*sizeof(Int64));
   Int64 *itmp;
   itmp = new Int64[n];
   for (int j = 0; j < n; ++j)
     itmp[j] = (Int64)i[j];

   //BUFFERED_WRITE(unit, itmp, sizeof(Int64), n);
   (void)fwrite(itmp, sizeof(Int64), n, fd[unit]);

   // free(itmp);
   delete[] itmp;
}

void dwrite(const int unit, const double *d, const int n) 
{
   //BUFFERED_WRITE(unit, d, sizeof(double), n);
   (void)fwrite(d, sizeof(double), n, fd[unit]);
}




/*
*************************************************************************
*									*
* void openfile(char *filename, char *mode, Integer *n)			*
* void closefile()							*
*									*
* Function openfile opens a pipe to either gzip (to compress a stream)	*
* or zcat (to uncompress a stream).  Function closefile() closes the	*
* pipe stream created by openfile().					*
*									*
*************************************************************************
*/

#define FUDGE_FACTOR (8)

void openfile(const int unit, const char *filename, const char *mode) {

  if ((*mode == 'r') || (*mode == 'R')) {
    //bufferwrite[unit] = false;
    if (!(fd[unit] = fopen(filename, "rb")))
       BAIL("ERROR:  Could not open pipe from input file\n");
  } else {
    //bufferwrite[unit] = true;
    if (!(fd[unit] = fopen(filename, "wb")))
       BAIL("ERROR:  Could not open pipe to output file\n");
  }
  //bufferIndex[unit] = 0;
}

void closefile(const int unit) 
{ 
   //if (bufferwrite[unit] && (bufferIndex[unit] > 0))
   //   flushBufferWrite(unit);

   (void)fclose(fd[unit]); 
}

} // namespace pwdft
