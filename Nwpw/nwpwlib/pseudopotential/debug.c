/*
 $Id$
   get_word.c -
   Author - Eric Bylaska

*/
#include "debug.h"
#include <stdio.h>

static int debug;

int debug_print() { return debug; }

void set_debug_print(int i) { debug = i; }
