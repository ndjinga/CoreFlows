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

CDMATH-CoreFlows relies on the library 'CDMATH' for the handling of meshes and fields, and on the library 'PETSC' for the handling of large sparse matrices.
You will need the packages 'doxygen' if you want to generate de documentation and 'swig' if you want to use python scripts.

Download and compilation of PETSC
---------------------------------
In order to install PETSC you need first to download and compile the sources of PETSC, version 3.4 or later. 
The sources of PETSC 3.7.2 can be downloaded with the command
-  `wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.7.2.tar.gz`
Unzip the file and move into the newly created directory using the commands
- `tar xzf petsc-3.7.2.tar.gz `
- `cd petsc-3.7.2`
Then compile petsc using on the commands
- `./configure  --prefix=~/workspace/petsc-3.7.2 --with-mpi=0 --download-f2cblaslapack=1`
- `make PETSC_DIR=~/workspace/petsc-3.7.2 PETSC_ARCH=arch-linux2-c-opt all`

For the moment, CDMATH-CoreFlows is a sequential library (no parallelism).

Download and compilation of CDMATH
----------------------------------
In order to download 'CDMATH' either unzip the following file to a directory cdmath_src
- `wget https://github.com/PROJECT-CDMATH/CDMATH/archive/master.zip`

or
- `git clone https://github.com/PROJECT-CDMATH/CDMATH.git cdmath_src`

In order to compile 'CDMATH' you will need the packages 'cmake ', and 'hdf5' plus 'numpy' and 'swig' if you intend to use cdmath functions in your python scripts. 
First create build and install repositories:

- `mkdir ~/workspace/cdmath `
- `cd ~/workspace/cdmath `
- `mkdir cdmath_build `
- `mkdir cdmath_install `

Go to the build directory
- `cd cdmath_build `

Then run the commands
- `cmake ../cdmath_src/ -DCMAKE_INSTALL_PREFIX=../cdmath_install -DCMAKE_BUILD_TYPE=Release -DCDMATH_WITH_PYTHON=ON `
- `make install -j`



Download and compilation of CoreFlows
---------------------------------------------
In order to download CDMATH-CoreFlows source files, unzip the following file to a directory CoreFlows-master
- `wget https://github.com/mndjinga/CoreFlows/archive/master.zip `
- `unzip master.zip`

The following steps assume that 'PETSC' (version 3.4 or more recent) and 'CDMATH' are installed on your computer

In CoreFlows-master, open the file CoreFlows.sh and set the variables :
- `PETSC_DIR`, the path to your PETSC installation
- `PETSC_ARCH`, the type of installation used (usually arch-linux2-c-opt or linux-gnu-c-opt)
- `CDMATH_INSTALL`, the path to your CDMATH installation
- `CoreFlows_ROOT`, the path to the CoreFlows-master directory
- `CoreFlows_INSTALL`, the path to the installation directory (default value is CoreFlows_ROOT)
- `CoreFlows_PYTHON`, set to 'ON' if you intend to use python scripts, 'OFF' otherwise
- `CoreFlows_DOC`, set to 'ON' if you want the doxygen documentation to be generated, 'OFF' otherwise
- `CoreFlows_GUI`, set to 'ON' if you want to use CDMATH-CoreFlows as a Salomé module (you will need to use a Salomé shell), 'OFF' otherwise


Once the file CoreFlows.sh has been edited, you can open a terminal, go to CoreFlows-master and type
- `./CoreFlows_install.sh` to compile the library CDMATH-CoreFlows

Use of CDMATH-CoreFlows
-----------------------
First load CDMATH-CoreFlows environment from the CoreFlows-master directory
- `source CoreFlows-master/CoreFlows.sh `

If you use C language: edit the file CoreFlows-master/CoreFlows_src/main.cxx then in a terminal type
- `cd CoreFlows-master/CoreFlows_build  `
- `make -j `
- `make install -j `

Then you can run the simulation in any directory with the command line
- `$CoreFlows `

If you use python language: edit your own python file `my_file.py` following for example the pattern of the file `CoreFlows-master/CoreFlows_src/main.py`. Then in a terminal type
- `python my_file.py `

The complete documentation is available in the directory `$CoreFlows_ROOT/CoreFlows_install/share/doc/`
