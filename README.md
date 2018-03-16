CDMATH-CoreFlows
================

CDMATH-CoreFlows is an open source C++/Python library intended at solving PDE systems
arising from the thermalhydraulics of two phase flows in power plant boilers. It
is a simple environment meant at students and researchers to test new numerical
methods on general geometries with unstructured meshes. It is developped by
CEA Saclay since 2014 and proposes a few
basic models and finite volume numerical methods. Some of the main objectives
are the study of

- Numerical schemes for compressible flows at low Mach numbers
- Well balanced schemes for stiff source terms (heat source, phase change, pressure losses)
- Flow inversion and counter-current two phase flows
- Schemes that preserve the phasic volume fraction α ∈ [0, 1]
- Convergence of finite volume methods
- New preconditioners for implicit methods for two phase flows
- The coupling of fluid models or multiphysics coupling (eg thermal hydraulics and neutronics or thermal hydraulics and solid thermics)

CDMATH-CoreFlows relies on the numerical toolbox of the project [CDMATH]((http://cdmath.jimdo.com)) for the handling of meshes and fields, and on the library [PETSC](https://www.mcs.anl.gov/petsc/) for the handling of large sparse matrices.
You will need the packages 'doxygen' if you want to generate de documentation and 'swig' if you want to use python scripts.

The following instructions detail the compilation and installation of both CDMATH and CDMATH-CoreFlows in a linux terminal.

Download and compilation of CDMATH
----------------------------------

In order to compile 'CDMATH' you will need the packages 'cmake', 'gcc', 'gfortran', 'hdf5' plus 'numpy' and 'swig' if you intend to use CoreFlows via python scripts.
First create and access a working directory :
- `mkdir -p ~/workspace/cdmath `
- `cd ~/workspace/cdmath `

Then create build and install repositories:
- `mkdir cdmath_build `
- `mkdir cdmath_install `

Download 'CDMATH' sources and unziping them 
- `wget https://github.com/mndjinga/CDMATH/archive/master.zip`
- `unzip master.zip `
This latter command results in the creation of a directory `~/workspace/cdmath/CDMATH-CoreFlows-master` containing the source files

Go to the build directory
- `cd cdmath_build `

Then run the commands
- `cmake ../cdmath-master/ -DCMAKE_INSTALL_PREFIX=../cdmath_install -DCMAKE_BUILD_TYPE=Release -DCDMATH_WITH_PYTHON=ON -DCDMATH_WITH_PETSC=ON`
- `make install`

By default, CDMATH will compile a new sequential installation of PETSc. If an installation of PETSc (version 3.4 or later) is already available in the system, it is possible to save time by first setting the environment variables PETSC_DIR and PETSC_ARCH to the appropriate values as can be found in petscconf.h, and then run the above cmake command.

Download and compilation of CoreFlows
---------------------------------------------
First create and access a working directory :
- `mkdir -p ~/workspace/CDMATH-CoreFlows `
- `cd ~/workspace/CDMATH-CoreFlows `
Now create build and install repositories:
- `mkdir CDMATH-CoreFlows_build `
- `mkdir CDMATH-CoreFlows_install `

Download CDMATH-CoreFlows source files in zipped form
- `wget https://github.com/mndjinga/CDMATH-CoreFlows/archive/master.zip `
or clone the git repository to a folder CoreFlows-master
-`git clone https://github.com/mndjinga/CDMATH-CoreFlows.git CoreFlows-master` 

Unzip the source file
- `unzip master.zip`
This latter command results in the creation of a directory `~/workspace/CDMATH-CoreFlows/CDMATH-CoreFlows-master`  containing the source files.

In the following steps we assume that 'PETSC' (version 3.4 or more recent) and 'CDMATH' are installed on your computer.
You need to set the following variables 
- `CDMATH_DIR`, the path to your CDMATH installation, for example  `~/workspace/cdmath/cdmath_install//share/petsc-3.8.0 `
- `PETSC_DIR`, the path to your PETSC installation. If Petsc was installed by CDMATH then PETSC_DIR can be defined as `~/workspace/cdmath/cdmath_install`
- `PETSC_ARCH`, the type of installation used (usually arch-linux2-c-opt or linux-gnu-c-opt)
In order to do so, type in you linux terminal
- `export CDMATH_DIR=~/workspace/cdmath/cdmath_install`
- `export PETSC_DIR=/path/to/my/petsc/installation`
- `export PETSC_ARCH=my-petsc-arch`

then create build and install repositories next to CoreFlows-master :
- `mkdir CoreFlows_build CoreFlows_install`

Go to the build directory
- `cd CoreFlows_build `

Then run the command
- `../CDMATH-CoreFlows-master/configure  --prefix=../CDMATH-CoreFlows_install/ --with-petsc-dir=$PETSC_DIR --with-petsc-arch=$PETSC_ARCH --with-cdmath-dir=$CDMATH_DIR --with-python --with-doc`
- `make install doc`

You can add the following optional commands
- `--with-gui`, if you want to use CDMATH-CoreFlows as a Salomé module (you will need to use a Salomé shell)
- `--with-debug`, if you want to use CDMATH-CoreFlows in debug mode instead of the default optimised mode

Use of CDMATH-CoreFlows
-----------------------
First load CDMATH-CoreFlows environment from the CoreFlows-master directory
- `source ~/workspace/CDMATH-CoreFlows/CDMATH-CoreFlows_install/env_CoreFlows.sh `

If you use C language: edit the file CoreFlows-master/CoreFlows_src/main.cxx then in a terminal type
- `cd ~/workspace/CDMATH-CoreFlows/CDMATH-CoreFlows_build  `
- `make`
- `make install`
Then you can run the simulation in any directory with the command line
- `$CoreFlows `

If you use python language: edit your own python file `my_file.py` following for example the pattern of the file `CoreFlows-master/CoreFlows_src/main.py`. Then in a terminal type
- `python my_file.py `

If you use the graphic interface, you need to run a Salomé Unix shell 
- `./salome shell`
and type the command line
- `runSalome -mCOREFLOWS`
then click on new study to open CoreFlows interface

The complete documentation is available in the directory `~/workspace/CDMATH-CoreFlows/CDMATH-CoreFlows_install/share/doc/`
