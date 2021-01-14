# BLOM: Bergen Layered Ocean Model

This is the source code of BLOM and includes the ocean biogeochemistry
model iHAMOCC. BLOM is the ocean component of the Norwegian Earth System
Model (<https://github.com/NorESMhub/NorESM>).

## BLOM documentation

BLOM documetation is integrated in the general NorESM documentation on ReadTheDocs (<https://noresm-docs.readthedocs.io/en/noresm2/>).
- [Running OMIP-type experiments](https://noresm-docs.readthedocs.io/en/noresm2/configurations/omips.html#blom)
- [BLOM model description](https://noresm-docs.readthedocs.io/en/noresm2/model-description/ocn_model.html)
- [iHAMOCC model description](https://noresm-docs.readthedocs.io/en/noresm2/model-description/ocn_model.html)

### Building the code
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

### Working with the BLOM git repository

The [BLOM wiki](https://github.com/NorESMhub/BLOM/wiki) includes instructions on how to contribute to the BLOM/iHAMOCC model system, and how to work with the BLOM git repository with your own fork on gitHub.

## License

BLOM is licensed under the GNU Lesser General Public License - see the
[COPYING](COPYING) and [COPYING.LESSER](COPYING.LESSER) files for
details.
