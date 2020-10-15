# Instructions
## Building
Start by setting the environment on the current machine, if building on your
own machine you can skip this step, however, ensure that all dependencies are
available.

```sh
# When working on Fram
$ source setintel_FRAM
# When working on Betzy
$ source setintel_BETZY
```

Then run the `depend` rule to update the automatically generated `Makefile`.
```sh
$ make depend
```

Then run `make` to build the program
```sh
$ make -j
# The '-j' enables multi-threaded building which should decrease compile times
```

To clean up any temporary files run `make clean`.

### Change compiler
To build with a different compiler than `intel` use the command line option
`COMPILER_TARGET=` to change which compiler is utilized.

```sh
$ make COMPILER_TARGET=gfortran depend
$ make -j COMPILER_TARGET=gfortran
```

### Building without `mpi`
To build without `mpi` support remove `-DMPI` and `-DPNETCDF` from
`DIRECTIVE_FLAGS` in the corresponding `<compiler>.make` file. Then ensure that
`IOTYPE` in the `limits` file (see further down for explanation of this file) is
set to `IOTYPE = 0`.

## Running
To run it on the national infrastructure
```sh
$ cd ../BLOM_channel_new01/
$ ./run.sh_normal_betzy #it runs model for 1 months on 64 processors; we set 4 nodes as short queue time limit was bit lower 
```

To run locally it's required to copy a few files first before running the
resulting executable.

```sh
$ cp ../BLOM_channel_new01/limits_BLOM_channel_new01 limits
$ cp ../regions_and_indices/*.nc .
$ cp ../regions_and_indices/*.dat .
$ ./blom
```
