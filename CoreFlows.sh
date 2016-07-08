
# CDMATH_DIR, PETSC_DIR, PETSC_ARCH, CoreFlows_ROOT : Paths to be set by the user 
export CDMATH_DIR=/export/home/ndjinga/Logiciels/CDMATH/INSTALL
export PETSC_DIR=/export/home/ndjinga/Logiciels/Petsc/petsc-3.7.2
export PETSC_ARCH=arch-linux2-c-opt # usually for fedora: arch-linux2-c-opt or arch-linux2-c-debug, for ubuntu: linux-gnu-c-opt or linux-gnu-c-debug
export CoreFlows_ROOT=/export/home/ndjinga/Logiciels/CoreFlows
export CoreFlows_INSTALL=$CoreFlows_ROOT/CoreFlows_install

#Compilation options (PYTHON, Doc, GUI) to be set by the user
export CoreFlows_PYTHON='ON'   # To generate the SWIG module "Python = ON or OFF "
export CoreFlows_DOC='ON'      # To generate the Documentation  "Doc = ON or OFF "
export CoreFlows_GUI='OFF'      # To generate the Graphic user interface  "GUI = ON or OFF "

#------------------------------------------------------------------------------------------------------------------- 
export CoreFlows=$CoreFlows_INSTALL/bin/Executable/CoreFlowsMainExe
export LD_LIBRARY_PATH=$CDMATH_DIR/lib/:$PETSC_DIR/$PETSC_ARCH/lib:$CoreFlows_INSTALL/lib:${LD_LIBRARY_PATH}
export PYTHONPATH=$CoreFlows_INSTALL/lib:$CoreFlows_INSTALL/lib/CoreFlows_Python:$CoreFlows_INSTALL/bin/CoreFlows_Python:$CoreFlows_INSTALL/lib/python2.7/site-packages/salome:$CDMATH_DIR/lib/:$CDMATH_DIR/lib/cdmath:$CDMATH_DIR/bin/cdmath:${PYTHONPATH}
export CoreFlowsGUI=$CoreFlows_INSTALL/bin/salome/CoreFlows_Standalone.py
export COREFLOWS_ROOT_DIR=$CoreFlows_INSTALL
