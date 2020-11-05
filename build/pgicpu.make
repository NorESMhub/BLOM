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
LIBS=-L/scratch/project_2003520/netCDF/lib -lnetcdf -lnetcdff  
# External names
EXTNAME=

# Compiler flags
# Optimization level
OPT=-O2 -Kieee  #-O2 -acc -ta=tesla:managed -Minfo=all
OPENMP=
DEBUG=-g
#FFLAGS=-r8 -Kieee -byteswapio -Mrecursive -mcmodel=medium -Mflushz $(OPT) $(OPENMP) $(DEBUG) -I$HOME/netCDF/include/
FFLAGS=-r8 -byteswapio -Mflushz -mcmodel=medium $(OPT) $(OPENMP) $(DEBUG) -I/scratch/project_2003520/netCDF/include
CFLAGS= -Kieee -mcmodel=medium $(OPENMP) 

# Linker flags
LDFLAGS=$(LIBS) $(OPENMP) $(DEBUG) $(OPT) #-acc -ta=tesla:managed

# Archiver flags
ARFLAGS=-r

#DIRECTIVE_FLAGS="-DMPI -DLEVITUS2X -DTRC -DTKE -DTKEADV -DIDLAGE -DPNETCDF"
DIRECTIVE_FLAGS="-DLEVITUS2X -DTRC -DTKE -DTKEADV -DIDLAGE"
