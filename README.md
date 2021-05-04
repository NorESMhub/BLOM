# BLOM: Bergen Layered Ocean Model

This is the source code of BLOM and includes the ocean biogeochemistry
model iHAMOCC. BLOM is the ocean component of the Norwegian Earth System
Model (<https://github.com/NorESMhub/NorESM>).

## BLOM documentation

BLOM documetation is integrated in the general NorESM documentation on ReadTheDocs (<https://noresm-docs.readthedocs.io/en/noresm2/>).
- [Running OMIP-type experiments](https://noresm-docs.readthedocs.io/en/noresm2/configurations/omips.html#blom)
- [BLOM model description](https://noresm-docs.readthedocs.io/en/noresm2/model-description/ocn_model.html)
- [iHAMOCC model description](https://noresm-docs.readthedocs.io/en/noresm2/model-description/ocn_model.html)

### Building the code with meson
When compiling BLOM with NorESM, the NorESM build system should be used. A stand-alone
BLOM executable can be built by using the meson build system.
To build the code ensure that [`Meson`](https://mesonbuild.com/) is available.
The following will build the default version of BLOM _without_ `MPI`.

```bash
$ meson setup builddir --buildtype=debugoptimized
$ meson compile -C builddir
```

The executable `blom` file will then be stored in the `./builddir` directory.

The code has several different configuration options setup in
[`meson_options.txt`](./meson_options.txt) which can be specified either when
setting up or through [`meson
configure`](https://mesonbuild.com/Commands.html#configure).

To change the configuration, after `meson setup`, use `meson configure`. E.g. to
enable `MPI` and parallel `NetCDF` one could reconfigure the build by

```bash
$ meson configure builddir -D mpi=true -D parallel_netcdf=true
$ meson compile -C builddir
```

The same configuration when setting up would be
```bash
$ meson setup builddir --buildtype=debugoptimized -D mpi=true -D parallle_netcdf=true
$ meson compile -C builddir
```

#### Changing compiler
To change the compiler one must define the `FC` (Fortran compiler) and `CC` (C
compiler) environment variables. As an example, the following changes the
compiler suite to the Intel compilers.

```bash
$ CC=icc FC=ifort meson setup builddir --buildtype=debugoptimized
$ meson compiler -C builddir
```

After the first line all subsequent compiles (and changes with `meson
configure`) will utilize the Intel compiler to build and link BLOM.

### Running tests
After successfully building the code it can be a good idea to test that the code
behaves as expected and changes to the code does not affect the output.

Tests can be run with the following:

```bash
$ meson test -C builddir
```

The previous command will run all the test suites defined for BLOM. To run tests
quicker one can select a few tests to run or just a single test suite. To list
the available tests run `meson test -C builddir --list`. One can then run a
single test with:

```bash
$ meson test -C builddir "run single_column"
```

### Working with the BLOM git repository

The [BLOM wiki](https://github.com/NorESMhub/BLOM/wiki) includes instructions on how to contribute to the BLOM/iHAMOCC model system, and how to work with the BLOM git repository with your own fork on gitHub.

## License

BLOM is licensed under the GNU Lesser General Public License - see the
[COPYING](COPYING) and [COPYING.LESSER](COPYING.LESSER) files for
details.
