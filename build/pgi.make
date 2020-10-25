# PGI compiler definitions (NOTE: Experimental!)

# Fortran compiler
FC=pgfortran
# C compiler
CC=pgcc
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
DEBUG=-pg
FFLAGS=$(LIBS) $(OPT) $(OPENMP) $(DEBUG)
CFLAGS=$(OPT) $(OPENMP) $(DEBUG)

# Linker flags
LDFLAGS=$(LIBS) $(OPENMP) $(DEBUG)

# Archiver flags
ARFLAGS=-r

DIRECTIVE_FLAGS="-DLEVITUS2X -DTRC -DTKE -DTKEADV -DIDLAGE"

