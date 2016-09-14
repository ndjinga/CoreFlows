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
