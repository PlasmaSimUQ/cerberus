# Getting started
This is a general guide for those who have a basic familiarity with GNU/Linux operating systems.

**Windows users:** follow the guide to set up the [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/setup/environment).
This will allow you run a linux system with (close to) native performance provided you use WSL2
(generally the default for new installations unless you have hardware virtualization disabled in your motherboard UEFI).

## Install the dependencies
To install the dependencies on ubuntu
``` sh
sudo apt install bear build-essential liblapack-dev libopenmpi-dev libreadline-dev libboost-dev
```
**Windows users:** add *python3* to that list, as it doesn't come installed by default on *WSL*.

Similarly for *archlinux*:
``` sh
sudo pacman -Syu base-devel git lapack openmpi readline boost
```

### Eilmer gas dynamics ([UQ Gas Dynamics Toolkit](https://github.com/gdtk-uq/gdtk))
If you with to use the GDTk chemistry models, then you will need a *D* compiler package.
```sh
sudo apt install ldc
```

*Archlinux*
``` sh
sudo pacman -S ldc
```

## Compiling the code
Follow these instructions to download and build a copy of the code.
If you have any problems or weird error messages, then [ask a question](https://github.com/PlasmaSimUQ/cerberus/discussions/5).
Communication through the discussion board for most issues is preferable, so that anyone running into the same problem can be informed.

### Getting a copy of *Cerberus*
**NOTE: Do not use the "Download Zip" option in github.**

Clone *Cerberus* and `cd` into the directory:
```sh
git clone https://github.com/plasmasimuq/cerberus && cd cerberus
```

If you have an [*ssh* key set up for *github*](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)
``` sh
git clone git@github.com:plasmasimuq/cerberus && cd cerberus
```

### Building
For most simulation use cases, you shouldn't need to change any build options.
If you intend to use some of the more experimental options, or the external *UQ Gas Dynamics Toolkit*,
then copy [`Make.local.template`](Make.local.template) as `Make.local` in the root directory:
``` sh
cp Make.local.template Make.local
```
Edit `Make.local`, the build system will include it when you run `make`.

To build, run the following command in the *cerberus* root directory:
```sh
make cerberus -j$(nproc)
```

Alternatively, you can specify build flags in the command line,
``` sh
make cerberus USE_EB=FALSE AMREX_PARTICLES=FALSE -j$(nproc)
```

Look inside [`Make.local.template`](Make.local.template) for descriptions of the build flags.
Commented values are the defaults unless noted otherwise.

**NOTE:** the `-j` option specifies the number of threads for parallel compilation and should be set to an appropriate value.
In this case, we set it to the output of the `nproc` command.

Once the build process has completed, you can run the executable `MFP.2d.STUFF.ex` directly from the build directory.

#### UQ Gas Dynamics Toolkit (Eilmer)
To expand the capabilities of the code it is possible to use the *Eilmer* gas models implemented in the
[UQ Gas Dynamics Toolkit](https://github.com/gdtk-uq/gdtk).
This requires the download, compilation, and installation of the Eilmer gas library (found at `src/gas` in the Eilmer source) as well as the D-language compiler and runtime (see instructions [here](https://gdtk.uqcloud.net/)).
GDTk will be installed automatically when calling the make command with `EILMER_GAS=TRUE`. Alternatively, call `make eilmer` to just build GDTk.
This will automatically update the `gdtk` submodule to the required version and then compile libraries and executables.

**NOTE:** if building with Eilmer gas dynamics, then you will need to keep the executable in the original location, as it looks for the gas library with a relative path.

## Running simulations
Once the build process has completed, if you're using the GDTk gas models, then run:
```sh
source env
```
from the root of the repository.

See `Exec/testing/` for example simulations. We are still working on the documentation for these and the run scripts may not always work.

The executable can be run in parallel using `mpirun` and expects an inputs file that defines the simulation.
This inputs file is made up of two parts; a section using the AMReX style configuration which defines the
simulation infrastructure, and a Lua script which defines the physics side of the simulation including
initial conditions and physical reference parameters.

To run a simulation simply call the executable along with an inputs file. For example, using the Isotope-RMI test case:
```sh
cd Exec/testing/Wedge
mpirun ../../../MFP.2d.gnu.MPI.EB.PARTICLES.ex Wedge.inputs
```
All output files will be saved to the location where the run command was issued.

The accompanying `run` scripts may be used to automatically launch the test cases, but they are not guaranteed to be up-to-date.

### Visualisation
Visualisation can be carried out using  [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/)
or by using the provided *python* scripts in [`vis`](vis/).

To use visit, generate `movie.visit` after running the simulation.
Carrying on from our wedge example:
```sh
../visit.sh .
visit -o movie.visit
```

This is assuming that you have *VisIt* installed and available in your `PATH`.

### Using the gas and reaction models
An example is provided in [`tests/Eilmer-Gas`](tests/Eilmer-Gas).
This example shows how to prepare the files required by Eilmer using the `prep-gas` and `prep-chem` routines that are bundled in the Eilmer install.
The gas model specified for a state must be of type `eilmer` and specify the name of the file that defines the gas model as generated by `prep-gas`.
To define species within a state, the names and charges of each must be specified along with the mass fraction *alpha*.

The specific heats and masses of the species are calculated from the gas model (make sure to define appropriate reference quantities). 

Only `n-1` *alpha* values need to be defined, for `n` species, as the `n`-th value is calculated by assuming that all *alpha* values sum to unity.
The reactions scheme is defined as an action, where the species named within the reactions scheme are assumed to be held by the associated states. Species can be defined across multiple states.

### Advanced: Python debugging
It is possible to visualise the raw data while a simulation is running for debugging purposes.
You will require the python library headers as well as *numpy* (with headers) and *matplotlib*

## Updating Cerberus and submodules
Update cerberus with `git pull`. This should automatically update submodules to the latest compatible version as well.

To manually update submodules,
``` sh
git submodule update --init --recursive
```

## Documentation
We are still working on documenting this code. For now, see the notes in (doc/)[`doc`], and rely on the test inputs files
to learn how to run simulations.

**If you have worked with the code, please consider adding some notes to the wiki (we are happy to grant write access)**

## Developing
Pull requests are welcome. If you intend to make additions to the code and don't have a C++ development environment set up, we provide recommended extensions for the `Visual Studio Code` editor.
When you open this repository in `vscode`, it will ask if you want to install the recommended extensions.
Use the development build command to get full autocomplete, jump to reference/definition, etc.

We welcome all pull requests submissions.
