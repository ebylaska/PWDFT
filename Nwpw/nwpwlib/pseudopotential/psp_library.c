

#include        <stdio.h>
#include        <stdlib.h>
#include        <string.h>

#include	"psp_library.h"


static const char *libraryps_dir;

void psp_library_setdir(const char *dirname)
{
    libraryps_dir = dirname;
}

