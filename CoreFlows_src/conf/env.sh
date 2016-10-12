#!/bin/bash

export CoreFlows_INSTALL=@CMAKE_INSTALL_PREFIX@
export CDMATH_DIR=@CDMATH_DIR@
export PETSC_DIR=@PETSC_DIR@
export PETSC_ARCH="@PETSC_ARCH@"



# CDMATH_DIR, PETSC_DIR, PETSC_ARCH, CoreFlows_ROOT : Paths to be set by the user 
export CDMATH_DIR=${HOME}/taff/coreflows/cdmath_install

#------------------------------------------------------------------------------------------------------------------- 
export CoreFlows=$CoreFlows_INSTALL/bin/Executable/CoreFlowsMainExe
export LD_LIBRARY_PATH=$CDMATH_DIR/lib/:$PETSC_DIR/$PETSC_ARCH/lib:$CoreFlows_INSTALL/lib:${LD_LIBRARY_PATH}
export PYTHONPATH=$CoreFlows_INSTALL/lib:$CoreFlows_INSTALL/lib/CoreFlows_Python:$CoreFlows_INSTALL/bin/CoreFlows_Python:$CoreFlows_INSTALL/lib/python2.7/site-packages/salome:$CDMATH_DIR/lib/:$CDMATH_DIR/lib/cdmath:$CDMATH_DIR/bin/cdmath:${PYTHONPATH}
export CoreFlowsGUI=$CoreFlows_INSTALL/bin/salome/CoreFlows_Standalone.py
