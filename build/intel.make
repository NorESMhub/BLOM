# Intel compiler definitions

# Fortran compiler
FC=mpiifort
# C compiler
CC=mpiicc
# Linker
LD=$(FC)
# Archiver
AR=ar
# Include directory for modules
MODINC=
# Linker flags
LIBS=-lnetcdf -lnetcdff -lpnetcdf
# External names
EXTNAME=

# Compiler flags
# Optimization level
OPT=-O2
OPENMP=
DEBUG=-pg
FFLAGS=-real-size 64 -mkl=cluster -fp-model source -qno-opt-dynamic-align -convert big_endian -assume byterecl -ftz $(OPT) $(OPENMP) $(DEBUG)
CFLAGS=-fp-model precise $(OPENMP)

# Linker flags
LDFLAGS=$(LIBS) $(OPENMP) $(DEBUG)

# Archiver flags
ARFLAGS=-r

DIRECTIVE_FLAGS="-DMPI -DLEVITUS2X -DTRC -DTKE -DTKEADV -DIDLAGE -DPNETCDF"
