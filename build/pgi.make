# Intel compiler definitions

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
LIBS=-L$HOME/netCDF/lib/ -lnetcdf -lnetcdff 
# External names
EXTNAME=

# Compiler flags
# Optimization level
OPT=-fast
OPENMP=
DEBUG=
FFLAGS=-r8 -Kieee -byteswapio -Mrecursive -mcmodel=medium -Mflushz $(OPT) $(OPENMP) $(DEBUG) -I$HOME/netCDF/include/
CFLAGS=-Kieee $(OPENMP) -mcmodel=medium

# Linker flags
LDFLAGS=$(LIBS) $(OPENMP) $(DEBUG)

# Archiver flags
ARFLAGS=-r

#DIRECTIVE_FLAGS="-DMPI -DLEVITUS2X -DTRC -DTKE -DTKEADV -DIDLAGE -DPNETCDF"
DIRECTIVE_FLAGS="-DLEVITUS2X -DTRC -DTKE -DTKEADV -DIDLAGE"