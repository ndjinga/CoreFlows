# Graphic interface---------------------------------------------------------------------------------------------#
														#
if (COREFLOWS_WITH_GUI)													    			 							    								#
# Find KERNEL 
# ==============
SET(KERNEL_ROOT_DIR $ENV{KERNEL_ROOT_DIR} CACHE PATH "Path to the Salome KERNEL")
IF( EXISTS ${KERNEL_ROOT_DIR} )
  LIST(APPEND CMAKE_MODULE_PATH "${KERNEL_ROOT_DIR}/salome_adm/cmake_files")
  INCLUDE(SalomeMacros)
  FIND_PACKAGE(SalomeKERNEL REQUIRED)
ELSE( EXISTS ${KERNEL_ROOT_DIR} )
  MESSAGE(FATAL_ERROR "We absolutely need a Salome KERNEL, please define KERNEL_ROOT_DIR!")
ENDIF( EXISTS ${KERNEL_ROOT_DIR} )

SET(_all_SCRIPTS
    main_tests.py
    )

# Find SALOME GUI 
# ==============
SET(GUI_ROOT_DIR $ENV{GUI_ROOT_DIR} CACHE PATH "Path to the Salome GUI")
IF(EXISTS ${GUI_ROOT_DIR})
  LIST(APPEND CMAKE_MODULE_PATH "${GUI_ROOT_DIR}/adm_local/cmake_files")
  FIND_PACKAGE(SalomeGUI)
ELSE(EXISTS ${GUI_ROOT_DIR})
  MESSAGE(FATAL_ERROR "We absolutely need a Salome GUI, please define GUI_ROOT_DIR!")
ENDIF(EXISTS ${GUI_ROOT_DIR})

# Platform setup
# ==============
INCLUDE(SalomeSetupPlatform)   # From KERNEL
    
# Prerequisites
# =============

# Mandatory products
IF(  $ENV{PARAVIEW_VERSION} STRLESS "5")
         FIND_PACKAGE(SalomePython REQUIRED)
ELSE()
         FIND_PACKAGE(SalomePythonInterp REQUIRED)
         FIND_PACKAGE(SalomePythonLibs REQUIRED)
ENDIF()

SALOME_INSTALL_SCRIPTS("${_all_SCRIPTS}" ${SALOME_INSTALL_SCRIPT_PYTHON}/tests)

SALOME_ACCUMULATE_ENVIRONMENT(PYTHONPATH ${PROJECT_BINARY_DIR}/src) 
SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)

ADD_TEST(PythonTests ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/main_tests.py)
SET_TESTS_PROPERTIES(PythonTests PROPERTIES ENVIRONMENT "${tests_env}")
						    								#
endif ()                                                                                                        #
#---------------------------------------------------------------------------------------------------------------#



############# python tests

function(CreatePythonTest pythonFile)
    add_test(${pythonFile} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${pythonFile})
    SET_PROPERTY(TEST ${pythonFile} PROPERTY ENVIRONMENT LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CDMATH_DIR}/lib:${PETSC_DIR}/${PETSC_ARCH}/lib)
    SET_PROPERTY(TEST ${pythonFile} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/CoreFlows_Python:${CMAKE_INSTALL_PREFIX}/bin/CoreFlows_Python:${CDMATH_DIR}/lib:${CDMATH_DIR}/lib/cdmath:${CDMATH_DIR}/bin/cdmath)
endfunction(CreatePythonTest)

# copy all *.py tests files before build
file(GLOB pythonTestFiles "*.py" )
file(COPY ${pythonTestFiles} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} 
    FILE_PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE
)


CreatePythonTest(DiffusionEquation_1DHeatedRod.py)
CreatePythonTest(DriftModel_1DBoilingAssembly.py)
CreatePythonTest(DriftModel_1DBoilingChannel.py)
CreatePythonTest(DriftModel_1DDepressurisation.py)
CreatePythonTest(DriftModel_1DPorosityJump.py)
CreatePythonTest(DriftModel_1DPressureLoss.py)
CreatePythonTest(DriftModel_1DRiemannProblem.py)
CreatePythonTest(DriftModel_1DVidangeReservoir.py)
CreatePythonTest(DriftModel_2BranchesBoilingChannels.py)
CreatePythonTest(DriftModel_2DInclinedChannelGravity.py)
CreatePythonTest(DriftModel_2DInclinedChannelGravityTriangles.py)
CreatePythonTest(DriftModel_2DInclinedChannelGravityBarriers.py)
CreatePythonTest(DriftModel_2DBoilingChannelBarrier.py)
CreatePythonTest(DriftModel_2DInclinedBoilingChannelBarrier.py)
CreatePythonTest(DriftModel_2DInclinedBoilingChannel.py)
CreatePythonTest(DriftModel_2DPorosityJump.py)
CreatePythonTest(DriftModel_2DPressureLoss.py)
CreatePythonTest(DriftModel_2DVidangeReservoir.py)
CreatePythonTest(DriftModel_2DVidangeReservoirUnstructured.py)
CreatePythonTest(DriftModel_3DBoilingChannelBarrier.py)
CreatePythonTest(FiveEqsTwoFluid_1DBoilingAssembly.py)
CreatePythonTest(FiveEqsTwoFluid_1DBoilingChannel.py)
CreatePythonTest(FiveEqsTwoFluid_1DVidangeReservoir.py)
CreatePythonTest(FiveEqsTwoFluid_2DInclinedBoilingChannel.py)
CreatePythonTest(FiveEqsTwoFluid_2DInclinedSedimentation.py)
CreatePythonTest(FiveEqsTwoFluid_2DVidangeReservoir.py)
CreatePythonTest(IsothermalTwoFluid_1DSedimentation.py)
CreatePythonTest(IsothermalTwoFluid_1DVidangeReservoir.py)
CreatePythonTest(IsothermalTwoFluid_2DVidangeReservoir.py)
CreatePythonTest(SinglePhase_1DDepressurisation.py)
CreatePythonTest(SinglePhase_1DHeatedAssembly.py)
CreatePythonTest(SinglePhase_1DHeatedChannel.py)
CreatePythonTest(SinglePhase_1DRiemannProblem.py)
CreatePythonTest(SinglePhase_1DWaterHammer.py)
CreatePythonTest(SinglePhase_2BranchesHeatedChannels.py)
CreatePythonTest(SinglePhase_2DHeatedChannelInclined.py)
CreatePythonTest(SinglePhase_2DLidDrivenCavity.py)
CreatePythonTest(SinglePhase_2DLidDrivenCavity_unstructured.py)
CreatePythonTest(SinglePhase_2DSphericalExplosion_unstructured.py)
CreatePythonTest(SinglePhase_2DThermalDiffusion.py)
CreatePythonTest(SinglePhase_2DVidangeReservoir.py)
CreatePythonTest(SinglePhase_2DWallHeatedChannel_ChangeSect.py)
CreatePythonTest(SinglePhase_3DHeatDrivenCavity.py)
CreatePythonTest(TransportEquation_1DHeatedChannel.py)