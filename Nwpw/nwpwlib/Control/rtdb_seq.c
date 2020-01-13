/*$Id: rtdb_seq.c,v 1.23 2007/05/24 12:43:55 edo Exp $*/

#include <stdlib.h>
#include <sys/types.h>
#include <stdio.h>
#include "Int64.h"

#include "hdbm.h"
#define DBT datum
#define SIZE(a) ((a).dsize)
#define DATA(a) ((a).dptr)
#define WRAP(data, size, p) datum_wrap(data, size, p)

#include <unistd.h>

#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <sys/stat.h>
#include <fcntl.h>
#include "rtdb_seq.h"
#include "misc.h"

//#if !defined(LINUX)
//extern char *strdup(const char *);
//extern void *malloc(size_t);
//extern void free(void *);
//#endif

#define MAX_RTDB 32


static struct {			/* Keep track of active RTDBs */
  int active;
  char *filename;
  int scratch;
  hdbm db;

} rtdb[MAX_RTDB];

struct info_struct{		/* Matching info entry for each datum */
  int ma_type;
  int nelem;
  char date[26];
};

static char info_header[] = "!rtdb!"; /* Prefix for info entries */


static long file_size(const char *filename)
/*
  Return size of file in bytes or -1 if does not exist
*/
{
  FILE *file = fopen(filename, "r+");
  
  long length;

  if (!file)
    return (size_t) -1;

  (void) fseek(file, 0L, 2);
  length = ftell(file);
  (void) fclose(file);

  /* printf("file %s is %ld bytes long\n", filename, length); */

  return length;
}

void ma_print(FILE *file, const int ma_type, const int nelem, void *p)
{
  int i, nprint;

  switch (ma_type) {
  case rtdb_char:	/* char */

    (void) fprintf(file, "%.*s\n", nelem, (char *) p);
    break;

  case rtdb_int:	/* int */
    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%d ", ((int *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;


  case rtdb_long:	/* long int */

    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%ld ", ((long *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case rtdb_float:	/* float */

    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%.7e ", ((float *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case rtdb_double:	/* double */

    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%.14e ", ((double *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  default:

    (void) fprintf(file, " !! %d = unknown data type\n", ma_type);
  }
}

static int ma_nametype(const char *maname)
{
   int ma_type = rtdb_char;
   if (strcmp(maname,"char")==0)   ma_type = rtdb_char;
   if (strcmp(maname,"int")==0)    ma_type = rtdb_int;
   if (strcmp(maname,"long")==0)   ma_type = rtdb_long;
   if (strcmp(maname,"float")==0)  ma_type = rtdb_float;
   if (strcmp(maname,"double")==0) ma_type = rtdb_double;
   return ma_type; 
}
static int ma_typesize(const int intype)
{
  switch (intype) {
  case rtdb_char:	/* char */
    return sizeof(char); break;
  case rtdb_int:	/* int */
    return sizeof(Int64); break;
  case rtdb_log:	/* long */
    return sizeof(int); break;
  case rtdb_long:	/* long */
    return sizeof(rtdb_long); break;
  case rtdb_float:	/* float */
    return sizeof(float); break;
  case rtdb_double:	/* double */
    return sizeof(double); break;
  case rtdb_complex:	/* complex */
    return 2*sizeof(float); break;
  case rtdb_double_complex:	/* double complex */
    return 2*sizeof(double); break;
  default:
    return sizeof(char); break;
  }
}
static int ma_sizeof(const int intype, const int nelem, const int outtype)
{
   return (ma_typesize(intype)*nelem/ma_typesize(outtype));
}

const char *ma_typename(const int ma_type)
{
  switch (ma_type) {
  case rtdb_char:	/* char */
    return "char"; break;
  case rtdb_int:	/* int */
    return "int"; break;
  case rtdb_long:	/* long int */
    return "long"; break;
  case rtdb_float:	/* float */
    return "float"; break;
  case rtdb_double:	/* double */
    return "double"; break;
  default:
    return "invalid"; break;
  }
}

void rtdb_init_()
{
  int nw;
  
  for (nw=0; nw<MAX_RTDB; nw++)
    rtdb[nw].active = 0;

}

int rtdb_seq_open(const char *filename, const char *mode, int *handle)
/*
    Filename = path to file associated with the data base
    mode     = 'new'     Open only if it does not exist already
               'old',    Open only if it does exist already
               'unknown' Create new or open existing (preserving contents)
               'empty'   Create new or open existing (deleting contents)
               'scratch' Create new or open existing (deleting contents)
                         and automatically delete upon closing.  Also, items
                         cached in memory are not written to disk.
    handle   = returns handle by which all future references to the
               data base are made
*/
{
  int flags = O_RDWR | O_CREAT;
  int exists = access(filename, R_OK | W_OK) == 0;
  int nw;

  /* See if the data base is already open ... if so return and complain */

  for (nw=0; nw<MAX_RTDB; nw++)
    if (rtdb[nw].active)
      if (strcmp(filename, rtdb[nw].filename) == 0) {
	(void) fprintf(stderr, "rtdb_seq_open: %s is already open\n", 
		       filename);
	return 0;
      }

  /* Figure out if the file exists */

  if (strcmp(mode, "new") == 0) {
    if (exists) {
      (void) fprintf(stderr, 
		     "rtdb_seq_open: %s exists but is being opened new\n",
		     filename);
      return 0;
    }
  }
  else if (strcmp(mode, "old") == 0) {
    if (!exists) {
      (void) fprintf(stderr, 
		     "rtdb_seq_open: %s does not exist, cannot open old\n",
		     filename);
      return 0;
    }
  }
  else if (strcmp(mode, "empty") == 0) {
    flags |= O_TRUNC;
  }
  else if (strcmp(mode, "scratch") == 0) {
    flags |= O_TRUNC;
  }
  else if (strcmp(mode, "unknown") == 0)
    ;
  else {
    (void) fprintf(stderr, "rtdb_seq_open: unknown mode=%s\n", mode);
    return 0;
  }

  /* Find next free rtdb entry */

  for (nw=0; nw<MAX_RTDB && rtdb[nw].active; nw++)
    ;
  if (nw == MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_seq_open: too many data bases open, max=%d\n",
		   MAX_RTDB);
    return 0;
  }

  /* To avoid problems with incomplete data bases from crashed
     calculations, delete the data base if we can ... but DON'T
     do this until we have checked that it is not already open */

  if ((exists && (flags & O_TRUNC)) ||
      (strcmp(mode, "unknown") == 0 && file_size(filename) <= 0))
    (void) unlink(filename);

  /* Open the physical data base */


  /*  if (!hdbm_open(filename, 0, &rtdb[nw].db)) { */
  if (!hdbm_open(filename, 1, &rtdb[nw].db)) {
    (void) fprintf(stderr, "rtdb_seq_open: hdbm failed to open file %s\n",
		   filename);
    return 0;
  }

  /* Shove info into the RTDB structure and have at it */

  rtdb[nw].active = 1;
  rtdb[nw].filename = strdup(filename);
  rtdb[nw].scratch = (strcmp(mode, "scratch") == 0);

  *handle = nw;

  return 1;
}

static int check_handle(const int handle)
{
  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "check_handle: handle (%d) is out of range\n",
		   handle);
    return 0;
  }
  
  if (!rtdb[handle].active) {
    (void) fprintf(stderr, "check_handle: handle (%d) is inactive\n", handle);
    return 0;
  }

  return 1;
}

int rtdb_seq_first(const int handle, const int namelen, char *name)
/*
  Return the name of the first (user inserted) entry in the data base.
  The order is effectively random.
  handle  = handle to RTDB
  namelen = size of user provided buffer name
  name    = name of entry is returned in this buffer
*/
{

  hdbm db;
  int status;
  DBT key;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rdtb_first: handle (%d) is invalid\n", handle);
    return 0;
  }

  db = rtdb[handle].db;

  for (status=hdbm_first_key(db, &key);
       status;
       status=hdbm_next_key(db, &key)) {
    if (strncmp(info_header, (char *) DATA(key), strlen(info_header)) != 0)
      break;
    datum_free(key);
  }

  if (status) {
    if (namelen >= SIZE(key)) {
      strncpy(name, (char *) DATA(key), namelen);
      name[namelen-1] = 0;
    }
    else {
      (void) fprintf(stderr, 
		     "rtdb_seq_first: name is too small, need=%u got=%d\n",
		     (int) SIZE(key), (int) namelen);
      status = 0;
      }
    datum_free(key);
  }
  return status;
}

int rtdb_seq_delete(const int handle, const char *name)
/*
  Delete the entry from the database.
  Return
        1 if key was present and successfully deleted
	0 if key was not present, or if an error occured
  handle  = handle to RTDB
  name    = name of entry to delete
*/
{

  hdbm db = rtdb[handle].db;
  DBT key, info_key;
  int status;
  char info_buf[256];

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rdtb_delete: handle (%d) is invalid\n", handle);
    return 0;
  }

  if (strlen(name)+sizeof(info_header)+1 > sizeof(info_buf)) {
    (void) fprintf(stderr,
		   "rtdb_delete: info entry buffer needs to be %d\n",
		   (int) (strlen(name)+sizeof(info_header)+1));
    return 0;
  }
  
  (void) sprintf(info_buf, "%s%s", info_header, name);

  WRAP(info_buf, strlen(info_buf)+1, &info_key);
  WRAP(name, strlen(name)+1, &key);

  status = !hdbm_delete(db, key);
  if (status == 0)
    status = !hdbm_delete(db, info_key);

  if (status == 0)
    return 1;
  else if (status == 1)
    return 0;
  else {
    (void) fprintf(stderr, "rtdb_seq_delete: error deleting %s\n", name);
    return 0;
  }
}
    
int rtdb_seq_next(const int handle, const int namelen, char *name)
/*
  Return the name of the next (user inserted) entry in the data base.
  The order is effectively random.
  handle  = handle to RTDB
  namelen = size of user provided buffer name
  name    = name of entry is returned in this buffer
*/
{
  hdbm db = rtdb[handle].db;
  int status;
  DBT key;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rdtb_close: handle (%d) is invalid\n", handle);
    return 0;
  }

  while ((status = hdbm_next_key(db, &key))) {
    if (strncmp(info_header, (char *) DATA(key), strlen(info_header)) != 0)
      break;
    datum_free(key);
  }

  if (status) {
    if (namelen >= SIZE(key)) {
      strncpy(name, (char *) DATA(key), namelen);
    }
    else {
      (void) fprintf(stderr, 
		     "rtdb_seq_next: name too small, need=%u, got=%d\n",
		     (int) SIZE(key), (int) namelen);
      status = 0;
    }
    datum_free(key);
  }

  return status;
}


int rtdb_seq_close(const int handle, const char *mode)
/*
  Close the data base
  handle  = handle to RTDB
  mode    = 'keep'    Preserve the data base file to enable restart
            'delete'  Delete the data base file freeing all resources
  mode is overridden by opening the data base with mode='scratch'
  in which instance it is always deleted upon closing
*/
{
  hdbm db = rtdb[handle].db;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rdtb_close: handle (%d) is invalid\n", handle);
    return 0;
  }

  /* Mark as inactive even if we trap any other errors */

  rtdb[handle].active = 0;

  if (!hdbm_close(db)) {
    (void) fprintf(stderr, "rdtb_close: hdbm close(%s) returned error\n", 
		   rtdb[handle].filename);
    return 0;
  }

  /* Delete it if required */
    
  if (strcmp(mode,"delete") == 0 || rtdb[handle].scratch) {
    if (unlink(rtdb[handle].filename) != 0) {
      (void) fprintf(stderr, "rdtb_close: error in deleting %s\n",
		     rtdb[handle].filename);
      return 0;
    }
  }
  else {			/* Compress out dead space */
    if (!hdbm_file_compress(rtdb[handle].filename)) {
      (void) fprintf(stderr, "rdtb_close: hdbm compress(%s) returned error\n", 
		     rtdb[handle].filename);
      return 0;
    }
  }

  return 1;
}

static void get_time(char buf[26])
{
  time_t t = time((time_t *) 0);
  char *tmp = ctime(&t);

  (void) memcpy(buf, tmp, 26);
}

int rtdb_seq_print(const int handle, const int print_values)
/*
  Print the contents of the data base to stdout
  handle  = handle to RTDB
  print_values = boolean flag ... if true values as well as
                 keys are printed out.
*/
{
  char name[128];
  int status;
  int len;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_print: handle (%d) is invalid\n", handle);
    return 0;
  }

  printf("\n Contents of RTDB %s\n ", rtdb[handle].filename);
  len = strlen(rtdb[handle].filename) + 17;
  while(len--)
    printf("-");
  printf("\n\n");
  printf(" Entry                                   Type[nelem]  "
	 "         Date\n");
  printf(" ---------------------------  ----------------------  "
	 "------------------------\n");

  status = rtdb_seq_first(handle, sizeof(name), name);

  while (status) {
    char date[26];
    int nelem;
    int ma_type;

    if (!rtdb_seq_get_info(handle, name, &ma_type, &nelem, date)) {
      (void) fprintf(stderr, "rtdb_seq_print: get_info failed on %s\n", name);
    }
    else {			/* Fortran is better at this ! */
      printf(" %s", name);
      len = 29 - strlen(name);
      while (len-- > 0)
	printf(" ");
      printf("%15s[%d]", ma_typename(ma_type), nelem);
      if (nelem < 10)
	printf(" ");
      if (nelem < 100)
	printf(" ");
      if (nelem < 1000)
	printf(" ");
      if (nelem < 10000)
	printf(" ");
      if (nelem < 100000)
	printf(" ");
      printf(" %s\n", date);

      if (print_values) {
        void *ma_hndl;
        if (!rtdb_seq_ma_get(handle, name, &ma_type, &nelem, &ma_hndl)) {
          (void) fprintf(stderr, "rtdb_seq_print: ma_get failed on %s\n", name);
        }
        else {
          void *data;
          data = (void *) ma_hndl;
          ma_print(stdout, ma_type, nelem, data);
          free((void *) ma_hndl);
        }
      }
    }
    status = rtdb_seq_next(handle, sizeof(name), name);
  }
  printf("\n");

  return 1;
}

static int rtdb_seq_put_info(const int handle,
			 const char *name, const int ma_type, const int nelem)
/*
  Insert info about an entry into the data base replacing previous entry
  handle   = handle to RTDB
  name     = entry name (null terminated character string)
  ma_type  = MA type of the entry
  nelem    = no. of elements of the given type
*/  
{
  struct info_struct info;
  char info_buf[256];
  hdbm db = rtdb[handle].db;

  DBT key, value;
  int status;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_put_info: handle (%d) is invalid\n", handle);
    return 0;
  }

  info.ma_type = ma_type;
  info.nelem   = nelem;
  get_time(info.date);

  if (strlen(name)+sizeof(info_header)+1 > sizeof(info_buf)) {
    (void) fprintf(stderr,
		   "rtdb_seq_put_info: info entry buffer needs to be %d\n",
		   (int) (strlen(name)+sizeof(info_header)+1));
    return 0;
  }
  
  (void) sprintf(info_buf, "%s%s", info_header, name);

  WRAP(info_buf, strlen(info_buf)+1, &key);
  WRAP(&info, sizeof(info), &value);
  
  status = hdbm_replace(db, key, value);

  if (!status) {
    (void) fprintf(stderr, 
		   "rtdb_seq_put_info: put failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename);
    return 0;
  }

  return 1;
}

int rtdb_seq_get_info(const int handle,
		  const char *name, int *ma_type, int *nelem, char date[26])
/*
  Get info about an entry from the data base
  handle   = handle to RTDB
  name     = entry name (null terminated character string)
  ma_type  = returns MA type of the entry
  nelem    = returns no. of elements of the given type
  date     = returns date of insertion (null terminated character string)
*/  
{
  struct info_struct info;
  char info_buf[256];
  hdbm db = rtdb[handle].db;

  DBT key, value;
  int status;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_get_info: handle (%d) is invalid\n", 
		   handle);
    return 0;
  }

  if (strlen(name)+sizeof(info_header)+1 > sizeof(info_buf)) {
    (void) fprintf(stderr,
		   "rtdb_seq_get_info: info entry buffer needs to be %d\n",
		   (int) (strlen(name)+sizeof(info_header)+1));
    return 0;
  }
  
  (void) sprintf(info_buf, "%s%s", info_header, name);

  WRAP(info_buf, strlen(info_buf)+1, &key);

  status = hdbm_read(db, key, &value);

  if (!status) {
      /* Entry not found ... quietly return error so that failed
	 probes are not excessively verbose */
    return 0;
  }

  if (SIZE(value) != sizeof(info)) {
    (void) fprintf(stderr, 
		   "rtdb_seq_get_info: size mismatch : info=%d, db=%d\n",
		   (int) sizeof(info), (int) SIZE(value));
    return 0;
  }

  (void) memcpy(&info, DATA(value), SIZE(value));
  datum_free(value);

  *ma_type = info.ma_type;
  *nelem   = info.nelem;
  if (info.date[24] == '\n') info.date[24] = ' ';
  (void) memcpy(date, info.date, 26);

  return 1;
}

int rtdb_seq_put(const int handle, const char *name, const int ma_type,
	     const int nelem, const void *array)
/*
  Insert an entry into the data base replacing previous entry
  handle   = handle to RTDB
  name     = entry name (null terminated character string)
  ma_type  = MA type of the entry
  nelem    = no. of elements of the given type
  array    = data to be inserted
*/  
{
  hdbm db = rtdb[handle].db;
  DBT key, value;
  int status;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_put: handle (%d) is invalid\n", handle);
    return 0;
  }

  /* Enter the data into the data base */

  WRAP(name, strlen(name)+1,&key);
  WRAP(array, ma_typesize(ma_type)*nelem,&value);

  status = hdbm_replace(db, key, value);

  if (!status) {
    (void) fprintf(stderr, "rtdb_seq_put: put failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename);
    return 0;
  }

  /* Enter the info into the data base as "!rtdb!<name>" */

  if (!rtdb_seq_put_info(handle, name, ma_type, nelem))
    return 0;

  /* This operations flushes ALL writes immediately to disk
     ... it may be slow, in which case comment the next statement
     out and be sure to explicitly close the databse */

  status = hdbm_file_flush(db);  


  if (!status) {
    (void) fprintf(stderr, "rtdb_seq_put: flush failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename);
    return 0;
  }

  return 1;
}

int rtdb_seq_get(const int handle, const char *name, const int ma_type,
	     const int nelem, void *array)
/*
  Get an entry from the data base
  handle   = handle to RTDB
  name     = entry name (null terminated character string)
  ma_type  = MA type of the entry which must match entry type
  nelem    = size of array in units of ma_type
  array    = user provided buffer that returns data
*/
{
  hdbm db = rtdb[handle].db;

  DBT key, value;
  int status, ma_type_from_db, nelem_from_db;
  char date[26];
  int nelem_actual;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_get: handle (%d) is invalid\n", handle);
    return 0;
  }

  /* Retrieve the info from the data base to check types etc. */

  if (!rtdb_seq_get_info(handle, name, &ma_type_from_db, &nelem_from_db, date)) {

    /* In production will want to have this handled quietly  ... be
       verbose for debug only 
    (void) fprintf(stderr, "rtdb_seq_get: get info failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename); */
    return 0;
  }

  if (ma_type_from_db != ma_type) {
    (void) fprintf(stderr, 
		   "rtdb_seq_get: type mismatch \"%s\" in %s: arg=%d, db=%d\n",
		   name, rtdb[handle].filename, ma_type, ma_type_from_db);
    return 0;
  }

  /* Get the data from the data base */

  WRAP(name, strlen(name)+1, &key);

  status = hdbm_read(db, key, &value);

  if (!status) {
    (void) fprintf(stderr, "rtdb_seq_get: get failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename);
    return 0;
  }

  /* Now check that user's buffer is big enough */

  nelem_actual = ma_sizeof(rtdb_char, SIZE(value), ma_type);

  if (nelem_actual != nelem_from_db) {
    (void) fprintf(stderr, 
		   "rtdb_seq_get: size error \"%s\" in %s: info=%d, db=%d\n",
		   name, rtdb[handle].filename, nelem_from_db,
		   nelem_actual);
    return 0;
  }

  if (nelem_actual > nelem) {
    (void) fprintf(stderr, "rtdb_seq_get: \"%s\" is %d long, buffer is %d\n",
		   name, nelem_actual, nelem);
    return 0;
  }

  (void) memcpy(array, DATA(value), SIZE(value));
  datum_free(value);

  return 1;
}

int rtdb_seq_ma_get(const int handle, const char *name, int *ma_type,
		int *nelem, void **ma_hndl)
/*
  Get an entry from the data base returning an MA handle
  handle   = handle to RTDB
  name     = entry name (null terminated character string)
  ma_type  = returns MA type of the entry
  nelem    = returns no. of elements of type ma_type in data
  ma_handle= returns MA handle to data
*/
{
  char date[26];
  int ma_handle_buf;
  int ma_type_buf;
  int nelem_buf;
  void *ma_data;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_ma_get: handle (%d) invalid\n", handle);
    return 0;
  }

  /* Retrieve the type info from the data base */

  if (!rtdb_seq_get_info(handle, name, ma_type, nelem, date)) {
    /*
    (void) fprintf(stderr, "rtdb_seq_ma_get: get info failed \"%s\" in %s\n",
		   name, rtdb[handle].filename);
                   */
    return 0;
  }

  /* Allocate MA pointers */

  ma_type_buf = (int) *ma_type;
  nelem_buf   = (int) *nelem;

  ma_data = (void *) malloc(nelem_buf*ma_typesize(ma_type_buf)); 
  *ma_hndl = (void *) ma_data;
  
  if (!rtdb_seq_get(handle, name, *ma_type, *nelem, ma_data)) {
    (void) fprintf(stderr, "rtdb_seq_ma_get: rtdb_seq_get failed %s\n", name);
    (void) free(ma_data);
    return 0;
  }

  return 1;
}

int rtdb_seq_copy(const int handle, const char *suffix)
/*
  Copy the data base
  handle  = handle to RTDB
  suffix    = new file will be called rtdb_name.suffix
*/
{

    if (!hdbm_file_copy(rtdb[handle].filename, suffix)) {
	(void) fprintf(stderr,
		       "rtdb_seq_copy: copy from %s failed\n", suffix);
      return 0;
    }
  return 1;
}
int rtdb_seq_getfname(const int handle,
		  char fname[36])
/*
  Get rtdb file name
  handle   = handle to RTDB
  fname    = returns rtdb file name (null terminated character string)
*/  
{

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_get_info: handle (%d) is invalid\n", 
		   handle);
    return 0;
  }

  (void) memcpy(fname, rtdb[handle].filename, 36);

  return 1;
}
