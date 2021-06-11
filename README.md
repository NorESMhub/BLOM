# BLOM: Bergen Layered Ocean Model

This is the source code of BLOM and includes the ocean biogeochemistry
model iHAMOCC. BLOM is the ocean component of the Norwegian Earth System
Model (<https://github.com/NorESMhub/NorESM>).

## BLOM documentation

BLOM documetation is integrated in the general NorESM documentation on ReadTheDocs (<https://noresm-docs.readthedocs.io/en/latest/>).
- [Running OMIP-type experiments](https://noresm-docs.readthedocs.io/en/latest/configurations/omips.html#blom)
- [BLOM model description](https://noresm-docs.readthedocs.io/en/latest/model-description/ocn_model.html)
- [iHAMOCC model description](https://noresm-docs.readthedocs.io/en/latest/model-description/ocn_model.html)

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
to the BLOM/iHAMOCC model system. The [BLOM wiki](https://github.com/NorESMhub/BLOM/wiki), 
includes more detailed instructions on how to work with the BLOM git repository with your
own fork on gitHub.

## License

BLOM is licensed under the GNU Lesser General Public License - see the
[COPYING](COPYING) and [COPYING.LESSER](COPYING.LESSER) files for
details.
