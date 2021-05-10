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
  if get_option('buildtype').startswith('debug')
    add_project_arguments(['-fbacktrace', '-fcheck=all', '-finit-real=snan',
                           '-ffpe-trap=invalid,zero,overflow'],
                          language: 'fortran')
  endif
  # If MPI is enabled we need to allow for argument mismatch with gfortran to
  # make it able to compile
  if get_option('mpi') and fcc.has_argument('-fallow-argument-mismatch')
    add_project_arguments('-fallow-argument-mismatch',
                          language: 'fortran')
  endif
elif fcc.get_id() == 'intel'
  add_project_arguments(['-r8', '-convert', 'big_endian', '-assume',
                         'byterecl'], language: 'fortran')
  if get_option('buildtype').startswith('debug')
    add_project_arguments(['-traceback', '-check all'.split(),
                           '-init=snan,arrays', '-fpe0','-ftrapuv'],
                          language: 'fortran')
  endif
  if get_option('buildtype') == 'release'
    add_project_arguments('-fp-model precise'.split(),language: 'fortran')
  endif
  if get_option('processors') == 1 and get_option('grid') == 'channel'
    add_project_arguments('-mcmodel=medium', language: 'fortran')
  endif
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
sources = []

# Process the following subdirectories which contain other 'meson.build' files
# that will add files to build into 'sources'
subdir('phy')
subdir('drivers')
subdir('ben02')
subdir('cesm')
subdir('fuk95')
subdir('channel')
subdir('single_column')

# Handle options and add necessary flags and subfolders with source files

if get_option('iage') or get_option('turbclo').length() > 0 or get_option('ecosys')
  add_project_arguments('-DTRC', language: 'fortran')
  subdir('trc')
endif

if get_option('turbclo').length() > 0
  if not (get_option('turbclo').contains('oneeq') or get_option('turbclo').contains('twoeq'))
    error('For turbulent closure, either twoeq or oneeq must be provided as options!')
  endif
  if get_option('turbclo').contains('oneeq') and get_option('turbclo').contains('twoeq')
    error('For turbulent closure, do not use both twoeq and oneeq as options!')
  endif
  if get_option('turbclo').contains('oneeq')
    add_project_arguments('-DTKE', language: 'fortran')
  endif
  if get_option('turbclo').contains('twoeq')
    add_project_arguments('-DTKE', '-DGLS', language: 'fortran')
  endif
  if get_option('turbclo').contains('advection')
    add_project_arguments('-DTKEADV', language: 'fortran')
  endif
  if get_option('turbclo').contains('isodif')
    add_project_arguments('-DTKEIDF', language: 'fortran')
  endif
endif

if get_option('iage')
  add_project_arguments('-DIDLAGE', language: 'fortran')
  subdir('idlage')
endif

if get_option('ecosys')
  add_project_arguments('-DHAMOCC','-DRESTART_BGC', '-DWLIN', language: 'fortran')
  if get_option('hamocc_cfc')
    add_project_arguments('-DCFC', language: 'fortran')
  endif
  if get_option('hamocc_nattrc')
    add_project_arguments('-DnatDIC', language: 'fortran')
  endif
  if get_option('hamocc_sedbypass')
    add_project_arguments('-Dsedbypass', language: 'fortran')
  endif
  if get_option('hamocc_ciso')
    if not get_option('hamocc_sedbypass')
      error('hamocc_ciso=true requires hamocc_sedbypass=true!')
    endif
    add_project_arguments('-Dcisonew', language: 'fortran')
  endif
  subdir('hamocc')
endif

if get_option('levitus2x')
  add_project_arguments('-DLEVITUS2X', language: 'fortran')
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

# Handle Unit Testing
subdir('tests')

# Create BLOM executable
executable('blom', sources, dimensions,
           include_directories: includes,
           dependencies: deps,
           link_language: 'fortran')
