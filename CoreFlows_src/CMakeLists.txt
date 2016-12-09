cmake_minimum_required (VERSION 2.8.10)

# Project name
project (CoreFlows C CXX)
set (COREFLOWS_VERSION_MAJOR 0)
set (COREFLOWS_VERSION_MINOR 1)

option (CMAKE_BUILD_TYPE "Debug or Release compilation" Release)
option (COREFLOWS_WITH_DOCUMENTATION "Generate documentation for CoreFlows use." OFF)
option (COREFLOWS_WITH_PYTHON "Compile Python interface for COREFLOWS." OFF)
option (COREFLOWS_WITH_GUI "Compile Graphic user Interface for COREFLOWS." OFF)
option (COREFLOWS_WITH_PACKAGE "Generate RPM, Debian and tarball packages." OFF)

enable_testing()											   #
# Debug mode activates testing and profiling ---------------------------------------------------------------#
if (CMAKE_BUILD_TYPE STREQUAL Debug)										   #
													   #
 SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg --coverage")								   #
 SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg --coverage")						   #
 SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg --coverage")					   #
													   #
endif ()												   #
													   #
#----------------------------------------------------------------------------------------------------------#

# PETSC ----------------------------------------------------------------------------------------------------#
list (APPEND CMAKE_MODULE_PATH "${CoreFlows_SOURCE_DIR}/cmake_files")                                       #  
find_package (PETSc 3.4 REQUIRED)
petsc_get_version ()
  if (${PETSC_VERSION} VERSION_GREATER 3.5)
    add_definitions(-DPETSC_VERSION_GREATER_3_5)
  endif ()
set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, petsc-dev (>= 3.4), petsc-dev (<< 3.5)")
set (CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES}, petsc >= 3.4, petsc < 3.5") # This is not fully working yet.
                                                                                                            #
#-----------------------------------------------------------------------------------------------------------#



# CDMATH ---------------------------------------------------------------------------------------------------#
find_package (CDMATH REQUIRED)                                                                              #
#-----------------------------------------------------------------------------------------------------------#

# Base directories
set (CoreFlows_SRC ${CoreFlows_SOURCE_DIR}/Models ) 
set (CoreFlows_EXAMPLES ${CoreFlows_SOURCE_DIR}/examples)
add_subdirectory (${CoreFlows_SRC})
add_subdirectory (${CoreFlows_EXAMPLES})


# Swig interface -------------------------------------------------------------------------------------------#
                                                                                                            #
if (COREFLOWS_WITH_PYTHON)                                                                                  #
   find_package (PYTHON REQUIRED)                                                                           #
   find_package (SWIG REQUIRED)                                                                             #
   set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, python-dev, python-numpy, swig")     #
   set (CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES}, python-devel, numpy, swig")              #
   set (CoreFlows_SWIG_DIR ${CoreFlows_SOURCE_DIR}/swig)                                                    #
   add_subdirectory (${CoreFlows_SWIG_DIR})                                                                 #
endif ()                                                                                                    #
                                                                                                            #
#-----------------------------------------------------------------------------------------------------------#

# Documentation --------------------------------------------------------------------------------------------#
                                                                                                            #
if (COREFLOWS_WITH_DOCUMENTATION)                                                                           #
  find_package (Doxygen)                                                                                    #
  set (CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, doxygen, graphviz, mscgen")           #
  set (CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES}, doxygen, graphviz, mscgen")               #
  if (NOT DOXYGEN_FOUND)                                                                                    #
    message (FATAL_ERROR "Doxygen is needed to build the documentation. Please install it correctly.")      #
  endif (NOT DOXYGEN_FOUND)                                                                                 #
  configure_file (Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)                              #
  FILE(GLOB_RECURSE CoreFlows_Documentation ${CoreFlows_SOURCE_DIR}/Documentation *.dox)                    #
  SET(Documentation_SRC                                                                                     #
    ${PROJECT_BINARY_DIR}/Doxyfile                                                                          #
    ${CoreFlows_Documentation}                                                                              #
      )                                                                                                     #
  add_custom_target (doc COMMAND ${DOXYGEN_EXECUTABLE}                                                  #
  SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)			                                            #
  install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html DESTINATION share/doc/CoreFlows OPTIONAL)                  #
  INSTALL(FILES ${CoreFlows_SOURCE_DIR}/Documentation/CoreFlows.pdf DESTINATION share/doc/CoreFlows OPTIONAL)	    #
endif ()                                                                                                    #
                                                                                                            #
# ------------- --------------------------------------------------------------------------------------------#

# Graphic interface-----------------------------------------------------------------------------------------#
     													    #
if (COREFLOWS_WITH_GUI)											    							       														    #
 set (CoreFlows_GUI_DIR ${CoreFlows_SOURCE_DIR}/gui)                                                        #
 add_subdirectory(${CoreFlows_GUI_DIR})			

endif ()                                                                                                    #
													    #
# ------------- --------------------------------------------------------------------------------------------#

#--------------------- COMPILATION MAIN --------------------------------------------------------------------#
													    #
INCLUDE_DIRECTORIES(						   					    #
													    #
  ${PETSC_INCLUDES}											    #
  ${CDMATH_INCLUDES}											    #
  ${CoreFlows_SRC}/inc											    #    
  )													    #
													    #
SET(CoreFlowsMain_SOURCES										    #
    ${CoreFlows_SRC}/src/Fluide.cxx									    #
    ${CoreFlows_SRC}/src/DiffusionEquation.cxx								    #
    ${CoreFlows_SRC}/src/ProblemFluid.cxx								    #
    ${CoreFlows_SRC}/src/IsothermalTwoFluid.cxx								    #
    ${CoreFlows_SRC}/src/utilitaire_algebre.cxx								    #
    ${CoreFlows_SRC}/src/DriftModel.cxx									    #
    ${CoreFlows_SRC}/src/SinglePhase.cxx								    #
    ${CoreFlows_SRC}/src/ProblemCoreFlows.cxx								    #
    ${CoreFlows_SRC}/src/FiveEqsTwoFluid.cxx								    #
    ${CoreFlows_SRC}/src/TransportEquation.cxx								    #
  )													    #
													    #
SET(MAIN_SOURCES											    #
  main.cxx												    #
  )													    #
SET(_extra_lib_CoreFlows ${CDMATH_LIBRARIES} ${PETSC_LIBRARIES} )					    #
ADD_EXECUTABLE(CoreFlowsMainExe ${MAIN_SOURCES})							    #
SET_TARGET_PROPERTIES(CoreFlowsMainExe PROPERTIES COMPILE_FLAGS "")					    #
TARGET_LINK_LIBRARIES(CoreFlowsMainExe CoreFlows  ${_extra_lib_CoreFlows})				    #
INSTALL(TARGETS CoreFlowsMainExe DESTINATION bin/Executable)						    #
													    #

ADD_EXECUTABLE(TestEOS.exe    ${CoreFlows_SRC}/src/testEOS.cxx ${CoreFlows_SRC}/src/Fluide.cxx )
TARGET_LINK_LIBRARIES(TestEOS.exe  ${_extra_lib_CoreFlows})

configure_file(
    ${PROJECT_SOURCE_DIR}/conf/env.sh
    ${CMAKE_INSTALL_PREFIX}/env.sh
    @ONLY
)
													    
													    
# ------------- --------------------------------------------------------------------------------------------#


# Packaging
if (COREFLOWS_WITH_PACKAGE)
  include (CPackLists.txt)
endif ()
