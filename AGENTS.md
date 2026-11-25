# AGENTS.md

Guidance for AI coding agents working in the BLOM repository.
The contributor rules are in [CONTRIBUTING.md](CONTRIBUTING.md) and the
Fortran style is in [doc/coding-guidelines.rst](doc/coding-guidelines.rst);
this file summarizes what matters most when making changes.

This file describes conventions that apply to the whole project. Do not add
notes about individual features or changes to it; document those in the code
or in the BLOM wiki instead. If this file disagrees with the code or
configuration, follow the code and tell the user rather than editing this file.

## Project overview

BLOM (Bergen Layered Ocean Model) is the ocean component
of the Norwegian Earth System Model (NorESM). It includes the ocean
biogeochemistry model iHAMOCC. It is written in Fortran, built with Meson, and
run either stand-alone or coupled through CIME (MCT or NUOPC).

Key locations:

- `phy/` – ocean physics (MICOM/BLOM core).
- `hamocc/` – iHAMOCC ocean biogeochemistry.
- `trc/` – tracers.
- `ben02/`, `fuk95/`, `channel/`, `single_column/`, `idlage/`, `noforcing/` –
  forcing options and idealized configurations.
- `drivers/` – coupling to the host system (`mct/`, `nuopc/`, `nocoupler/`).
- `cime_config/`, `bld/` – CIME build and namelist machinery and grid
  configurations.
- `pkgs/` – external packages (git submodules).
- `tests/` – stand-alone test cases run with `meson test`.
- `meson.build`, `meson_options.txt` – build system and build options.

## Build and test

```bash
meson setup builddir -Dmpi=false -Dopenmp=disabled --buildtype=debugoptimized
meson compile -C builddir
ulimit -s 65532 && meson test -C builddir
```

- CI (`.github/workflows/ci.yml`) builds with gfortran for MPI on/off and
  OpenMP enabled/disabled, and runs the tests without MPI. Check that all
  four build combinations still compile when you change code.
- Do not commit build output (`builddir/`, `*.mod`, `*.o`).
- New source files must be added to the `meson.build` of their directory.

## Design

- Keep changes small and focused; do not mix functional changes with
  reformatting or refactoring of unrelated code.
- Prefer readable, maintainable code over efficient code.
- Code in `phy/` and `hamocc/` is shared between the MCT and NUOPC coupling
  systems and between stand-alone and coupled runs. Keep coupler-specific code
  in `drivers/`.
- Avoid new `#ifdef` directives; see the coding guidelines.
- Do not edit the external packages in `pkgs/` or the grid files and
  generated files in `bld/` unless that is the goal of the change.
- Do not copy code from other projects unless its license is compatible with
  the LGPL (see `COPYING.LESSER`).

## Code style

- Follow [doc/coding-guidelines.rst](doc/coding-guidelines.rst): free-form
  Fortran (no features beyond Fortran 2008), `implicit none`, declared `intent`
  for all dummy arguments, `use ... only`, `real(kind=wp)`, 3-space
  indentation, lines of at most 90 characters, lower-case keywords and names.
- Match the style, naming, and comment density of the surrounding code.
- Comment each new variable declaration with what it is, including units.
- Code must be safe with OpenMP: check which variables are shared or private
  in the parallel regions you touch.

## Tests

- Run `meson test` before finishing a change. Tests are stand-alone cases in
  `tests/` (e.g. `fuk95`, `single_column`).
- Add or extend a test when it is practical. Do not change existing tests or
  reference values to make new code pass: that may mean that existing
  features break. Flag such changes to the user.
- Passing tests do not show that results are scientifically correct; see
  below.

## Documentation

- Update the documentation when behavior visible to users changes. User
  documentation lives in the
  [NorESM documentation](https://noresm-docs.readthedocs.io/en/latest/);
  BLOM-specific developer information lives in the
  [BLOM wiki](https://github.com/NorESMhub/BLOM/wiki). Tell the user which
  pages need updating.
- For code that computes scientific results, add comments with references to
  papers and formula numbers.

## Backward compatibility and results

BLOM is used in NorESM production runs, so changes to namelists, input and
restart files, output fields, and CIME configuration must stay backward
compatible where possible. Flag any breaking change clearly to the user. Flag
any change that alters model results (also bit-for-bit changes, e.g. from
reordering operations or changing OpenMP or MPI decomposition) and tell the
user which output to compare with a run of the previous version. Changes in
`bld/` and `cime_config/` may require matching changes in NorESM.

## Pull requests

- Follow the AI rules in [CONTRIBUTING.md](CONTRIBUTING.md). In particular: do
  not open pull requests autonomously, do not write pull request descriptions,
  issue descriptions, or comments on behalf of the contributor, and remind the
  contributor to disclose which AI tools were used and for what.
- Development happens on feature branches that are merged into `master`
  through pull requests. Do not change the release branches (`release-*`)
  unless asked.
