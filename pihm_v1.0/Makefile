#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Yizhong Qu @ PENN STATE
# Version of  : 30 AUG 2004
#--------------------------------------------------------------------------------
# Makefile for the PIHM codes.
# This makefile generates the executables for the PIHM requested.
# Type 'make' to see usage instructions and list of examples.
#--------------------------------------------------------------------------------

SHELL = /bin/sh

#--------------------------
# Top of SUNDIALS directory
#--------------------------
SUNDIALS_DIR = ../..

#---------------------
# Path to header files
#---------------------
INC_DIR	= $(SUNDIALS_DIR)/include

#----------------------
# Path to library files
#----------------------
LIB_DIR = $(SUNDIALS_DIR)/lib

#-------------
# Architecture
#-------------
ARCH = `uname -s`.`uname -m`

#======================================================
# Machine-dependent variables
#======================================================

#------------------------------
# C compiler and compiler flags
#------------------------------

CC     = gcc
CFLAGS = -Wall -ffloat-store -I$(INC_DIR) -L$(LIB_DIR)
SRC    = ihm10.c f.c initialize.c read_alloc.c et_is.c print_functions.c 

#======================================================
# Make rules
#======================================================

all:
	@(echo)
	@(echo '       make pihm     - make pihm        ')
	@(echo '       make clean    - remove all executable files')
	@(echo)

pihm:
	@echo '...Compiling PIHM ...'
	@$(CC) $(CFLAGS) -o pihm $(SRC) -lcvode.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm

clean:
	@rm -f *.o
	@rm -f pihm

#---End of makefile---
