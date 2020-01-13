#
#
#
#########################################################################
#									#
# Sample makefile header for running with Gnue compilers  		#
#  The makefile targets are appended to  the end of this file		#
#	 Don't change anything that comes before the targets 		#
#									#
#									#
#########################################################################

MV		= mv
RM		= rm -f
LN		= ln -s
ECHO		= echo


C++ 		= /usr/local/bin/mpicxx
CC  		= /usr/local/bin/mpicc
# CC		= mpicc
FORT 		= gfortran
AR		= ar
RANLIB		= ranlib
C++LINK		= $(C++)
CLINK		= $(CC)


ifdef INCDIR
   INCLUDES = -I$(INCDIR)
endif
ifdef LIBDIR
   ifdef LIBNAME
      MYLIBRARY = $(LIBDIR)/lib$(LIBNAME).a
   endif
endif


ARCH_FLAGS      = -DMPICH_IGNORE_CXX_SEEK
# WARNINGS        = -Wno-deprecated -Wall
OPTIMIZATION    =  -O3
# OPTIMIZATION    =  -O3 -ffast-math -funroll-loops
#DEBUG          = -g

C++FLAGS        += $(INCLUDES) $(ARCH_FLAGS) $(WARNINGS) $(OPTIMIZATION) \
                  $(XTRAFLAGS) $(DEBUG)

CFLAGS		+= $(INCLUDES) $(ARCH_FLAGS) $(OPTIMIZATION) \
                  $(XTRAFLAGS) $(DEBUG)

FFLAGS		= -O3

ARFLAGS		= ru


LDFLAGS		= $(WARNINGS) $(OPTIMIZATION) $(DEBUG)


#########################################################################
# End of the System dependent prefix
#########################################################################

all: compile_subdirs $(OBJ_OPTIMIZE)
	@for i in $(OBJ_OPTIMIZE); do \
           echo compiling $$i; \
	   $(AR) r $(MYLIBRARY) $$i; \
	   $(RANLIB) $(MYLIBRARY); \
        done 

hppcopy:  hppcopy_subdirs
	@for i in $(HPPINCLUDES); do \
           echo copying include file $$i; \
           cp $$i $(INCDIR); \
        done

clean:	clean_subdirs
	@for i in $(OBJ_OPTIMIZE); do \
           echo removing $$i; \
	   $(RM) $$i; \
        done


compile_subdirs:
	@for i in $(SUBDIRS); do \
            echo compiling subdirectory $$i; \
            $(MAKE) -C  $$i; \
            echo done compiling subdirectory $$i; \
        done


hppcopy_subdirs:
	@for i in $(SUBDIRS); do \
            echo copying include files in subdirectory $$i; \
            $(MAKE) -C  $$i hppcopy; \
            echo done copying include files in subdirectory $$i; \
        done

clean_subdirs:
	@for i in $(SUBDIRS); do \
            echo cleaning subdirectory $$i; \
            $(MAKE) -C  $$i clean; \
            echo done cleaning subdirectory $$i; \
        done


#########################################################################
#									#
# Suffixes for compiling most normal C++, C, and f77 files		#
#									#
#########################################################################
.SUFFIXES:
.SUFFIXES: .cpp .c .f .o

.cpp.o:
		$(C++) $(C++FLAGS) -c $<
		@$(ECHO)

.c.o:
		$(CC) $(CFLAGS) -c $<
		@$(ECHO)


.f.o: *.f
		$(FORT)  -c $(FFLAGS) $<
