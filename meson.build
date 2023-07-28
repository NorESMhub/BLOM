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
    add_project_arguments('-O2 -fp-model precise'.split(),language: 'fortran')
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
subdir('pkgs/')

# Handle options and add necessary flags and subfolders with source files

if get_option('mks')
  add_project_arguments('-DMKS', language: 'fortran')
endif

turbclo = get_option('turbclo')
if turbclo.length() > 0 and get_option('vcoord') == 'cntiso_hybrid'
  message('Setting turbclo = [] for vcoord == \'cntiso_hybrid\'')
  turbclo = []
endif

if get_option('iage') or turbclo.length() > 0 or get_option('ecosys')
  add_project_arguments('-DTRC', language: 'fortran')
  subdir('trc')
endif

if turbclo.length() > 0
  if not (turbclo.contains('oneeq') or turbclo.contains('twoeq'))
    error('For turbulent closure, either twoeq or oneeq must be provided as options!')
  endif
  if turbclo.contains('oneeq') and turbclo.contains('twoeq')
    error('For turbulent closure, do not use both twoeq and oneeq as options!')
  endif
  if turbclo.contains('oneeq')
    add_project_arguments('-DTKE', language: 'fortran')
  endif
  if turbclo.contains('twoeq')
    add_project_arguments('-DTKE', '-DGLS', language: 'fortran')
  endif
  if turbclo.contains('advection')
    add_project_arguments('-DTKEADV', language: 'fortran')
  endif
  if turbclo.contains('isodif')
    add_project_arguments('-DTKEIDF', language: 'fortran')
  endif
endif

if get_option('iage')
  add_project_arguments('-DIDLAGE', language: 'fortran')
  subdir('idlage')
endif

if get_option('ecosys')
  add_project_arguments('-DHAMOCC', '-DWLIN', language: 'fortran')
  if get_option('hamocc_agg')
    add_project_arguments('-DAGG', language: 'fortran')
  endif
  if get_option('hamocc_boxatm')
    add_project_arguments('-DBOXATM', language: 'fortran')
  endif
  if get_option('hamocc_bromo')
    add_project_arguments('-DBROMO', language: 'fortran')
  endif
  if get_option('hamocc_cfc')
    add_project_arguments('-DCFC', language: 'fortran')
  endif
  if get_option('hamocc_dmsph')
    add_project_arguments('-DDMSPH', language: 'fortran')
  endif
  if get_option('hamocc_nattrc')
    add_project_arguments('-DnatDIC', language: 'fortran')
  endif
  if get_option('hamocc_fb_bgc_oce')
    add_project_arguments('-DFB_BGC_OCE', language: 'fortran')
  endif
  if get_option('hamocc_sedbypass')
    add_project_arguments('-Dsedbypass', language: 'fortran')
  endif
  if get_option('hamocc_ciso')
    if not get_option('hamocc_sedbypass')
      error('hamocc_ciso == true requires hamocc_sedbypass == true!')
    endif
    add_project_arguments('-Dcisonew', language: 'fortran')
  endif
  if get_option('hamocc_debug_timestep')
    add_project_arguments('-DPBGC_OCNP_TIMESTEP','-DPBGC_CK_TIMESTEP', language: 'fortran')
  endif
  if get_option('hamocc_pbgc_ocnp_timestep')
    add_project_arguments('-DPBGC_OCNP_TIMESTEP', language: 'fortran')
  endif
  if get_option('hamocc_pbgc_ck_timestep')
    add_project_arguments('-DPBGC_CK_TIMESTEP', language: 'fortran')
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
