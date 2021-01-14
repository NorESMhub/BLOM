project('BLOM', ['fortran', 'c'],
  version: '0.1.0',
  license: 'LGPL',
  meson_version: '>= 0.54.0',
  default_options: ['warning_level=2'])

# Setup additional flags for compilation
fcc = meson.get_compiler('fortran')
if fcc.get_id() == 'gcc'
  add_global_arguments(['-fdefault-real-8', '-fconvert=big-endian'],
    language: 'fortran')
elif fcc.get_id() == 'intel'
  add_global_arguments(['-r8', '-convert big_endian', '-assume byterecl'],
    language: 'fortran')
else
  warning('Unknown Fortran compiler ("' + fcc.get_id() + '"), no default flags specified')
endif

# Generate 'dimensions.F' based on desired patch set
blom_dims = find_program('bld/blom_dimensions')
dim_path = '..' / get_option('patch')
kdm_path = get_option('patch') / 'kdm'
if host_machine.system() in ['linux']
  dim_kdm = run_command('cat', kdm_path)
elif host_machine.system() == 'windows'
  dim_kdm = run_command('type', kdm_path)
else
  error('Could not read "kdm" content due to unknown OS (' + host_machine.system() + ')')
endif
if dim_kdm.returncode() != 0
  error('No "kdm" file found for patch set "' + get_option('patch') + '"')
endif
dim_src = configure_file(
  output: 'dimensions.F',
  command: [blom_dims,
    '-n', get_option('processors').to_string(),
    '-k', dim_kdm.stdout().strip(),
    '-d', dim_path])

# Define global dependencies that multiple subdirectories will utilize
netcdf = dependency('netcdf', language: 'fortran', version: '>=4.4.5')
omp = dependency('openmp', required: get_option('openmp'))
quadmath = fcc.find_library('quadmath')

# List of dependencies (used so that options can add dependencies)
deps = [netcdf, omp, quadmath]

# Directories with header files
cesm_inc = include_directories('cesm')
phy_inc = include_directories('phy')
trc_inc = include_directories('trc')

# List of directories to include
includes = [cesm_inc, phy_inc, trc_inc]

# Handle options and add necessary flags
if get_option('levitus2x')
  add_global_arguments('-DLEVITUS2X', language: 'fortran')
endif

if get_option('tke')
  add_global_arguments('-DTKE', language: 'fortran')
endif

if get_option('tkeadv')
  add_global_arguments('-DTKEADV', language: 'fortran')
endif

if get_option('idlage')
  add_global_arguments('-DIDLAGE', language: 'fortran')
endif

if get_option('mpi')
  add_global_arguments('-DMPI', language: 'fortran')
  deps += dependency('mpi', language: 'fortran')
endif

if get_option('parallel_netcdf')
  add_global_arguments('-DPNETCDF', language: 'fortran')
  deps += dependency('pnetcdf', version: '>=1.11.0')
  if not get_option('mpi')
    warning('MPI must be enabled for parallel NetCDF to work')
  endif
endif

# List of all sources in the project (should be added to in subfolders)
sources = [dim_src]

# Process the following subdirectories which contain other 'meson.build' files
# that will fill out 'sources'
subdir('ben02')
subdir('cesm')
subdir('drivers')
subdir('hamocc')
subdir('idlage')
subdir('phy')
subdir('trc')

# Create BLOM executable
executable('blom', sources,
  include_directories: includes,
  dependencies: deps,
  link_language: 'fortran')
