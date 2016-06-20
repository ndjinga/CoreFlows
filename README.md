CoreFlows
=========

CoreFlows is an open source C++/Python library intended at solving PDE systems
arising from the thermalhydraulics of two phase flows in power plant boilers. It
is a simple environment meant at students and researchers to test new numerical
methods on general geometries with unstructured meshes. It is developped at
CEA Saclay by Michael Ndjinga and his students since 2014 and proposes a few
basic models and finite volume numerical methods. Some of the main objectives
are the study of

- Numerical schemes for compressible flows at low Mach numbers
- Well balanced schemes for stiff source terms (heat source, phase change, pressure losses)
- Flow inversion and counter-current two phase flows
- Schemes that preserve the phasic volume fraction α ∈ [0, 1]
- Convergence of finite volume methods
- New preconditioners for implicit methods for two phase flows
- The coupling of fluid models or multiphysics coupling (eg thermal hydraulics and neutronics or thermal hydraulics and solid thermics)

CoreFlows relies on the library '\a cdmath' for the handling of meshes and fields, and on the library '\a petsc 3.4.5' for the handling of large sparse matrices.
You will need the packages '\a doxygen' if you want to generate de documentation and '\a swig' if you want to use python scripts.

Download and compilation of Petsc
---------------------------------
In order to install petsc you need first to download and compile the sources of Petsc 3.4.5. 
The sources of Petsc 3.4.5 can be downloaded from
-  `http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.5.tar.gz`

Then type on the commands
\verbatim 
 ./configure  --prefix=~/workspace/petsc-3.4.5 --with-mpi=0 --download-f-blas-lapack=1
 make all
 make install
\endverbatim

For the moment, CoreFlows is a sequential library (no parallelism).

Download and compilation of CDMATH
----------------------------------
In order to download '\a cdmath' either unzip the following file to a directory cdmath_src
- `https://github.com/PROJECT-CDMATH/CDMATH/archive/master.zip`

or
\verbatim 
git clone https://github.com/PROJECT-CDMATH/CDMATH.git cdmath_src
\endverbatim 

In order to compile '\a cdmath' you will need the packages '\a cmake ', and '\a hdf5' plus '\a numpy' and '\a swig' if you intend to use cdmath functions in your python scripts. 
First create build and install repositories:

\verbatim 
 mkdir ~/workspace/cdmath
 cd ~/workspace/cdmath
 mkdir cdmath_build
 mkdir cdmath_install
\endverbatim

Go to the build directory
\verbatim 
 cd cdmath_build
\endverbatim 

Then run the commands
\verbatim 
 cmake ../cdmath_src/ -DCMAKE_INSTALL_PREFIX=../cdmath_install -DCMAKE_BUILD_TYPE=Release -DCDMATH_WITH_PYTHON=ON
 make
 make install
\endverbatim 



Download and compilation of CoreFlows
---------------------------------------------
In order to download CoreFlows source files, either  unzip the following file to a directory CoreFlows-master (not yet available)
- `https://github.com/PROJECT-CoreFlows/CoreFlows/archive/master.zip `
or copy the directory /data/tmpletr/ndjinga/Logiciels/CoreFlows to a local directory CoreFlows-master

The following steps assume that '\a petsc' (version 3.4.5) and '\a cdmath' are installed on your computer

In CoreFlows-master, open the file CoreFlows.sh and set the variables :
- \a PETSC_DIR, the path to your PETSC installation
- \a PETSC_ARCH, the type of installation used (usually arch-linux2-c-opt or arch-linux2-c-debug)
- \a CDMATH_INSTALL, the path to your CDMATH installation
- \a CoreFlows_ROOT, the path to the CoreFlows-master directory
- \a CoreFlows_INSTALL, the path to the installation directory (default value is CoreFlows_ROOT)
- \a CoreFlows_PYTHON, set to 'ON' if you intend to use python script, 'OFF' otherwise
- \a CoreFlows_Doc, set to 'ON' if you want the doxygen documentation to be generated, 'OFF' otherwise


Once the file CoreFlows.sh has been edited, you can open a terminal, go to CoreFlows-master and type
- \a ./CoreFlows_install.sh to compile the library CoreFlows

Use of CoreFlows
----------------
First load CoreFlows environment from the CoreFlows-master directory
\verbatim
source CoreFlows-master/CoreFlows.sh
\endverbatim

If you use C language: edit the file CoreFlows-master/CoreFlows_src/main.cxx then in a terminal type
\verbatim
cd CoreFlows-master/CoreFlows_build 
make -j
make install -j
\endverbatim
Then you can run the simulation in any directory with the command line
\verbatim
$CoreFlows
\endverbatim

If you use python language: edit your own python file my_file.py following for example the pattern of the file CoreFlows-master/CoreFlows_src/main.py. Then in a terminal type
\verbatim
python my_file.py
\endverbatim

The complete documentation is available in the directory $CoreFlows_ROOT/CoreFlows_install/share/doc/
