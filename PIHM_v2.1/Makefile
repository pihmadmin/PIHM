# -----------------------------------------------------------------
# Version: 2.0
# Date: "----------" 
# -----------------------------------------------------------------
# Programmer: Mukesh Kumar @ PSU
# -----------------------------------------------------------------
# Makefile for PIHM 
#
# cvode/pihm/Makefile.  
# -----------------------------------------------------------------

SHELL = /bin/sh

srcdir       = .
builddir     = .
top_builddir = ../../
top_builddir = ../../
prefix       = /projects/sundials2.2/
exec_prefix  = ${prefix}
includedir   = ${prefix}/include
libdir       = ${exec_prefix}/lib



CPP      = /usr/bin/cc -E
CPPFLAGS = 
CC       = /usr/bin/gcc
CFLAGS   = -O0 -g -pg
#CFLAGS   = 
LDFLAGS  = 
LIBS     = -lm
SRC    = pihm.c f.c initialize.c read_alloc.c is_sm_et.c print.c update.c
 

COMPILER_PREFIX = 
LINKER_PREFIX   = 

SUNDIALS_INC_DIR = $(includedir)
SUNDIALS_LIB_DIR = $(libdir)
SUNDIALS_LIBS    = -lsundials_cvode -lsundials_nvecserial
# -lsundials_shared

# EXEC_FILES = cvdx cvdxe cvbx cvkx cvkxb cvdemd cvdemk

all:
	@(echo)
	@(echo '       make pihm     - make pihm        ')
	@(echo '       make clean    - remove all executable files')
	@(echo)

pihm:
	@echo '...Compiling PIHM ...'
	@$(CC) $(CFLAGS) -I$(SUNDIALS_INC_DIR) -I$(SUNDIALS_INC_DIR)/cvode -I$(SUNDIALS_INC_DIR)/sundials -L$(SUNDIALS_LIB_DIR) -o $(builddir)/pihm $(SRC) $(SUNDIALS_LIBS) $(LIBS)

clean:
	@rm -f *.o
	@rm -f pihm

