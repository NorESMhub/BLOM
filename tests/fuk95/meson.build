dimension_fuk95 = configure_file(
  output: 'dimensions.F',
  command: [blom_dims,
    '-n', '1',
    '-k', '12',
    '-d', meson.source_root() / 'bld'/ 'fuk95'])

fuk95 = executable('fuk95_blom',
                   sources, dimension_fuk95,
                   include_directories: includes,
                   dependencies: deps,
                   build_by_default: false,
                   link_language: 'fortran')
test('run fuk95', fuk95)

