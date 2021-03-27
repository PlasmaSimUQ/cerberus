# Cerberus

A solver for the ideal five-moment multi-fluid plasma equations using the AMReX framework.

## Getting Started

The following repositories should all go into the same folder to make life easier when defining relative paths in the build procedure.

Download Cerberus:
```
git clone https://bitbucket.org/plasmasimuq/cerberus
```

Checkout  your branch of choice:
```
git checkout develop
```

Initialise the amrex submodule:
```
git submodule update --init
```
This will automatically pull the required version of amrex in `./cerberus/amrex`

### Prerequisites

```
sudo apt install gfortran liblapack-dev libopenmpi-dev g++ libreadline-dev libboost-dev
```

### Building

Choose 2D or 3D in `cerberus/Exec/Make.MFP` by setting `DIM : = 2` or `DIM:=3`

If on your PC then:

```
cd cerberus/Exec/local
make -f GNUmakefile -j4
``` 

Note that the `-j4` option is for multi-thread compilation, set whatever number here that you like, or remove this option for single thread use.

At this point you should have an executeable such as `MFP2d.gnu.MPI.ex` in ` cerberus/Exec/local`.

For a `DEBUG` version simply call 

```
make -f GNUmakefile DEBUG=TRUE -j4
```

### Running

The executeable can be run in parallel using `mpirun` and expects an `*.inputs` file that defines the simulation. This inputs file is made up of two parts, a section using the AMReX style configuration which defines  the simulation infrastructure, and a Lua script which defines the physics side of the simulation including initial conditions and physical parameters. Examples may be found in `cerberus/Exec`. 

To run a simulation simply call the executeable along with an inputs file. For example, using the Isotope-RMI test case:

```
cd Exec/testing/Isotope-RMI
mpirun -np 4 ../../local/MFP.2d.gnu.MPI.EB.PARTICLES.ex IRMI.inputs
```

Note the relative path to the executable. All output files will be saved to the location where the run command was issued. Also the `-np 4` option for `mpirun` dictates how many processes the code will be run on.

### Visualisation

Visualisation can be carried out using  [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) or by using the provided python scripts in `cerberus/vis`.
