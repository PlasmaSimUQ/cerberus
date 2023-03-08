# Cerberus
A solver for the ideal five-moment multi-fluid plasma equations using the [AMReX framework](https://github.com/AMReX-Codes/amrex).

**NOTE:** for detailed instructions, see [Getting started](GETTING_STARTED.md). This code is very much in development. To report bugs, [create an issue](https://github.com/PlasmaSimUQ/cerberus/issues/new/choose).
If you would like help running the code, or getting started, please [ask a question](https://github.com/PlasmaSimUQ/cerberus/discussions/5)

## Dependencies
A GNU/Linux operating system with the following packages is required to build this code:
- *C* and *C++* compiler package, *Boost*, GNU *Readline*.
- *LAPACK*
- *OpenMPI*
- *GNU make* and various other *GNU* tools such as *sed*, *git*, etc.
- __*Optional:__ Build EAR (bear)* to generate the `compile_commands.json` index file for IDE features when developing.
- __*Optional:__ D compiler and runtime*, if you wish to use the [UQ Gas Dynamics Toolkit](https://github.com/gdtk-uq/gdtk)

## Building
If you wish to change build options, copy [Make.local.template](Make.local.template) to `Make.local`. The latter is untracked, so you can make changes there without polluting the git history.
The commented values in `Make.local.template` are the defaults, unless noted otherwise.

To run the build,
``` sh
make -j$(nproc)
```

Alternatively, specify build flags on the commandline,
``` sh
make -j$(nproc) USE_EB=FALSE AMREX_PARTICLES=FALSE -j9
```

Run `source env` from the root of the repo to initialise GDTk environment variables. This is necessary for running with the GDTk gas models.

### Development builds
To compile for development, and automatically generate the `compile_commands.json` index,
ensure that the executable `bear` is available at the commandline.
Then execute either `make develop -j8` or `make develop_eilmer -j8`.

This will set `DEBUG=TRUE` and generate a `compile_commands.json` file, so you can utilise it with an IDE or editor that supports the *clangd* or *ccls* language servers.
See [Getting started](GETTING_STARTED.md#development-setup) for a little more information.


## Further documentation
Head over to [Getting started](GETTING_STARTED.md).
