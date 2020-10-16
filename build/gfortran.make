# GCC compiler definitions

# Fortran compiler
FC=gfortran
# C compiler
CC=gcc
# Linker
LD=$(FC)
# Archiver
AR=ar
# Include directory for modules
MODINC=
# Linker flags
LIBS=$(shell nf-config --fflags --flibs)
# External names
EXTNAME=

# Compiler flags
# Optimization level
OPT=-O2
OPENMP=
DEBUG=-Wall -g
FFLAGS=$(LIBS) -fdefault-real-8 -fconvert=big-endian $(OPT) $(OPENMP) $(DEBUG)
CFLAGS=$(OPT) $(OPENMP) $(DEBUG)

# Linker flags
LDFLAGS=$(LIBS) $(OPENMP) $(DEBUG)

# Archiver flags
ARFLAGS=-r

DIRECTIVE_FLAGS="-DLEVITUS2X -DTRC -DTKE -DTKEADV -DIDLAGE"
