dimension_single_column = configure_file(
  output: 'dimensions.F',
  command: [blom_dims,
    '-n', '1',
    '-k', '53',
    '-d', meson.source_root() / 'bld'/ 'single_column'])

single_column = executable('single_column_blom',
                           sources, dimension_single_column,
                           include_directories: includes,
                           dependencies: deps,
                           build_by_default: false,
                           link_language: 'fortran')
test('run single_column', single_column)
