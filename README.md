# BLOM: Bergen Layered Ocean Model

This is the source code of BLOM and includes the ocean biogeochemistry
model iHAMOCC. BLOM is the ocean component of the Norwegian Earth System
Model (<https://github.com/NorESMhub/NorESM>).

## Structure of the BLOM repository

The BLOM repository contains source code corresponding to several versions of NorESM, which are contained in
permanent release branches. The structure of the branches and naming conventions are documented in the
[discussions item #164](https://github.com/NorESMhub/BLOM/discussions/164). In general, tags on the format
`v#.#.#` correspond to a relase version of BLOM, whereas tags on the format `dev#.#.#.#` correspond to 
development tags. Currently the following branches are actively maintained:

| branch      | note                          |
|-------------|-------------------------------|
| master      | main development branch       |
| release-1.6 | release version for NorESM2.3 |
| release-1.5 | release version for NorESM2.1 |
| release-1.4 | release version for NorESM2.0 |


## BLOM documentation

Since BLOM is mainly used in connection with the NorESM system, the BLOM user documetation has been integrated
into the general NorESM documentation on ReadTheDocs (<https://noresm-docs.readthedocs.io/en/latest/>).
- [Running OMIP-type experiments](https://noresm-docs.readthedocs.io/en/latest/configurations/omips.html#blom)
- [BLOM model description](https://noresm-docs.readthedocs.io/en/latest/model-description/ocn_model.html)
- [iHAMOCC model description](https://noresm-docs.readthedocs.io/en/latest/model-description/ocn_model.html)

The [BLOM wiki](https://github.com/NorESMhub/BLOM/wiki) contains information about
BLOM-specific topics that is not considered relevant for the general NorESM documentation:
- working with the BLOM git repository on gitHub
- running BLOM/iHAMOCC stand-alone test cases
- details about model structure

### Building a stand-alone BLOM executable with meson
When compiling BLOM with NorESM, the NorESM build system should be used. A stand-alone
BLOM executable can be built by using the meson build system.
To build the code ensure that [`Meson`](https://mesonbuild.com/) is available.
The following will build the default version of BLOM _without_ `MPI`.

```bash
$ meson setup builddir --buildtype=debugoptimized
$ meson compile -C builddir
```

The executable `blom` file will then be stored in the `./builddir` directory.

See the [BLOM/iHAMOCC stand-alone](https://github.com/NorESMhub/BLOM/wiki/Run-BLOM-iHAMOCC-stand-alone)
wiki page for further instructions on how to configure the meson build system.

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

## Contribute to BLOM/iHAMOCC development

The [CONTRIBUTING.md](CONTRIBUTING.md) file includes instructions on how to contribute
to the BLOM/iHAMOCC model system. The [BLOM wiki](https://github.com/NorESMhub/BLOM/wiki) 
includes more detailed instructions on how to work with the BLOM git repository with your
own fork on gitHub.

## License

BLOM is licensed under the GNU Lesser General Public License - see the
[COPYING](COPYING) and [COPYING.LESSER](COPYING.LESSER) files for
details.
