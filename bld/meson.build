# Generate 'dimensions.F' based on desired grid and processor count
blom_dims = find_program('blom_dimensions')
patch_path = meson.source_root() / 'bld' / get_option('grid')
kdm_path = get_option('grid') / 'kdm'
# More systems which support 'cat' can be added here, once tested
if host_machine.system() in ['linux', 'darwin']
  dim_kdm = run_command('cat', kdm_path)
elif host_machine.system() == 'windows'
  dim_kdm = run_command('type', kdm_path)
else
  error('Could not read "kdm" content due to unknown OS (' + host_machine.system() + ')')
endif
if dim_kdm.returncode() != 0
  error('No "kdm" file found for grid "' + get_option('grid') + '"')
endif
dimensions = configure_file(
  output: 'dimensions.F',
  command: [blom_dims,
    '-n', get_option('processors').to_string(),
    '-k', dim_kdm.stdout().strip(),
    '-d', patch_path])

