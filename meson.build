project('BLOM', ['fortran', 'c'],
  version: '1.0.0',
  license: 'LGPL',
  meson_version: '>= 0.54.0',
  default_options: ['warning_level=2'])

# Setup additional flags for compilation
fcc = meson.get_compiler('fortran')
if fcc.get_id() == 'gcc'
  add_project_arguments(['-fdefault-real-8', '-fconvert=big-endian'],
    language: 'fortran')
elif fcc.get_id() == 'intel'
  add_project_arguments(['-r8', '-convert', 'big_endian', '-assume', 'byterecl'],
    language: 'fortran')
else
  warning('Unknown Fortran compiler ("' + fcc.get_id() + '"), no default flags specified')
endif

# Define global dependencies
netcdf = dependency('netcdf', language: 'fortran', version: '>=4.4.4')
omp = dependency('openmp', required: get_option('openmp'))
quadmath = fcc.find_library('quadmath')

# Directories with header files that is needed when building
phy_inc = include_directories('phy')

# Handle building of 'dimensions.F' in the 'bld' directory
subdir('bld')

# List of dependencies (created so that options can add dependencies)
deps = [netcdf, omp, quadmath]

# List of directories to include, usually with C headers
includes = [phy_inc]

# List of all sources in the project (should be added to in subfolders)
sources = [dimensions]

# Handle options and add necessary flags
if get_option('levitus2x')
  add_project_arguments('-DLEVITUS2X', language: 'fortran')
endif

if get_option('tke')
  add_project_arguments('-DTKE', language: 'fortran')
endif

if get_option('tkeadv')
  add_project_arguments('-DTKEADV', language: 'fortran')
endif

if get_option('idlage')
  add_project_arguments('-DIDLAGE', language: 'fortran')
  subdir('idlage')
endif

if get_option('mpi')
  add_project_arguments('-DMPI', language: 'fortran')
  deps += dependency('mpi', language: 'fortran')
endif

if get_option('parallel_netcdf')
  add_project_arguments('-DPNETCDF', language: 'fortran')
  deps += dependency('pnetcdf', version: '>=1.11.0')
  if not get_option('mpi')
    warning('MPI must be enabled for parallel NetCDF to work')
  endif
endif

# Process the following subdirectories which contain other 'meson.build' files
# that will add files to build into 'sources'
subdir('ben02')
subdir('cesm')
subdir('drivers')
#subdir('hamocc')
subdir('phy')
subdir('trc')

# Create BLOM executable
executable('blom', sources,
  include_directories: includes,
  dependencies: deps,
  link_language: 'fortran')
