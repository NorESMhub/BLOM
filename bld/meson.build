# Generate 'dimensions.F' based on desired grid and processor count
blom_dims = find_program('blom_dimensions')
grid = get_option('grid')
processors = get_option('processors')
patch_path = meson.source_root() / 'bld' / grid
if   get_option('vcoord') == 'isopyc_bulkml'
  kdm_file = 'kdm.isopyc_bulkml'
elif get_option('vcoord') == 'cntiso_hybrid'
  kdm_file = 'kdm.cntiso_hybrid'
endif
kdm_path = grid / kdm_file
# More systems which support 'cat' can be added here, once tested
if host_machine.system() in ['linux', 'darwin']
  dim_kdm = run_command('cat', kdm_path)
elif host_machine.system() == 'windows'
  dim_kdm = run_command('type', kdm_path)
else
  error('Could not read ' + kdm_file + ' content due to unknown OS (' + host_machine.system() + ')')
endif
if dim_kdm.returncode() != 0
  error('No ' + kdm_file + ' file found for grid "' + grid + '"')
endif
dimensions = configure_file(
  output: 'dimensions.F',
  command: [blom_dims,
    '-n', processors.to_string(),
    '-k', dim_kdm.stdout().strip(),
    '-d', patch_path])

# Get grid type parameter "nreg"
blom_nreg = find_program('blom_nreg')
blom_nreg_runresult = run_command(blom_nreg, '-n', processors.to_string(), '-d', patch_path)
if blom_nreg_runresult.returncode() != 0
  error('Command "blom_nreg" failed with message\n' + blom_nreg_runresult.stderr().strip())
endif
nreg = blom_nreg_runresult.stdout().strip()
message('nreg: ', nreg)
