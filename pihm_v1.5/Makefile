# -----------------------------------------------------------------
# Makefile for PIHM   
# -----------------------------------------------------------------



SHELL = /bin/sh

srcdir       = .
builddir     = .
top_builddir = ../../
top_builddir = ../../
prefix       = /home1/gxb913/pihmXM/sundials-2.2.0
exec_prefix  = ${prefix}
includedir   = ${prefix}/include
libdir       = ${exec_prefix}/lib

NETCDF_INC_DIR   = /home1/gxb913/netcdfXM/netcdf-4.0-beta1
NETCDF_LIB_DIR   = /home1/gxb913/netcdfXM/netcdf-4.0-beta1
NETCDF_LIBS	 = -lnetcdf

CPP      = /usr/bin/cc -E
CPPFLAGS = 
CC       = /usr/bin/gcc
CFLAGS   = -g -O0
#CFLAGS   = 
LDFLAGS  = 
LIBS     = -lm
SRC    = calib.c pihm.c f.c initialize.c read_alloc.c et_is.c print.c
 

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
	@$(CC) $(CFLAGS) -I$(SUNDIALS_INC_DIR) -I$(SUNDIALS_INC_DIR)/cvode -I$(SUNDIALS_INC_DIR)/sundials -L$(SUNDIALS_LIB_DIR) -I$(NETCDF_INC_DIR)/include -L$(NETCDF_LIB_DIR)/lib -o $(builddir)/pihm $(SRC) $(SUNDIALS_LIBS) $(LIBS) $(NETCDF_LIBS)

clean:
	@rm -f *.o
	@rm -f pihm

