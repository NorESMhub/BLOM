dimension_fuk95 = configure_file(
  output: 'dimensions.F',
  command: [blom_dims,
    '-n', '1',
    '-k', '12',
    '-d', meson.source_root() / 'bld'/ 'fuk95'])

configure_file(input: 'limits',
               output: 'limits',
               copy: true)

fuk95 = executable('fuk95_blom',
                   sources, dimension_fuk95,
                   include_directories: includes,
                   dependencies: deps,
                   build_by_default: false,
                   link_language: 'fortran')

test('run fuk95', fuk95,
     is_parallel: false,
     should_fail: false,
     workdir: meson.current_build_dir(),
     suite: 'fuk95',
     timeout: 900)

