#
#
#
SUBDIRS = nwpwlib pspw

include ./config/makefile.h


all:	
	$(C++) main.cpp -L./lib -lblas -lfftpack -llapack -ln2pw -o nwpw.x
	

