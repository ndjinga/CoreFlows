project(test)

add_subdirectory (${CoreFlows_EXAMPLES}/Python)

INCLUDE_DIRECTORIES(
  ${PETSC_INCLUDES} 
  ${CDMATH_INCLUDES} 
  ${CoreFlows_SRC}/inc 
)

SET(CoreFlowsExample_SOURCES
    
    DiffusionEquation_1DHeatedRod.cxx
    DriftModel_1DBoilingChannel.cxx
    DriftModel_1DPressureLoss.cxx
    DriftModel_1DDepressurisation.cxx
    DriftModel_1DRiemannProblem.cxx
    DriftModel_2DInclinedBoilingChannel.cxx
    DriftModel_3DCanalCloison.cxx
    FiveEqsTwoFluid_1DBoilingChannel.cxx
    FiveEqsTwoFluid_1DDepressurisation.cxx
    FiveEqsTwoFluid_1DRiemannProblem.cxx
    FiveEqsTwoFluid_2DInclinedBoilingChannel.cxx
    FiveEqsTwoFluid_2DInclinedSedimentation.cxx
    IsothermalTwoFluid_1DDepressurisation.cxx
    IsothermalTwoFluid_1DRiemannProblem.cxx
    IsothermalTwoFluid_1DSedimentation.cxx
    IsothermalTwoFluid_2DInclinedSedimentation.cxx
    IsothermalTwoFluid_2DVidangeReservoir.cxx
    SinglePhase_1DDepressurisation.cxx
    SinglePhase_1DHeatedChannel.cxx
    SinglePhase_1DRiemannProblem.cxx
    SinglePhase_2DHeatDrivenCavity.cxx
    SinglePhase_2DHeatDrivenCavity_unstructured.cxx
    SinglePhase_2DLidDrivenCavity.cxx
    SinglePhase_2DLidDrivenCavity_unstructured.cxx
    SinglePhase_2DSphericalExplosion_unstructured.cxx
    SinglePhase_2DWallHeatedChannel_ChangeSect.cxx
    SinglePhase_2DWallHeatedChannel.cxx
    SinglePhase_2DHeatedChannelInclined.cxx
    SinglePhase_3DHeatDrivenCavity.cxx
    SinglePhase_HeatedWire_2Branches.cxx
    TransportEquation_1DHeatedChannel.cxx
    CoupledTransportDiffusionEquations_1DHeatedChannel.cxx
 
)

SET(_extra_lib_CoreFlows CoreFlows ${PETSC_LIBRARIES} ${CDMATH_LIBRARIES})

#ADD_LIBRARY(CoreFlowsExample SHARED ${CoreFlowsExample_SOURCES})
#SET(_extra_lib_CoreFlows CoreFlows ${PETSC_LIBRARIES} ${CDMATH_LIBRARIES})
#TARGET_LINK_LIBRARIES(CoreFlowsExample ${_extra_lib_CoreFlows})
#INSTALL(TARGETS CoreFlowsExample DESTINATION lib)

SET(MAIN_SOURCES
  main_tests.cxx
  )

#ADD_EXECUTABLE(CoreFlowsExampleExe ${MAIN_SOURCES})
#SET_TARGET_PROPERTIES(CoreFlowsExampleExe PROPERTIES COMPILE_FLAGS "")
#TARGET_LINK_LIBRARIES(CoreFlowsExampleExe ${_extra_lib_CoreFlows} CoreFlowsExample)
#ADD_TEST(NAME CoreFlowsCPPTests COMMAND CoreFlowsExampleExe)

if(CMAKE_COMPILER_IS_GNUCXX)
    if (CMAKE_BUILD_TYPE=Debug)
    include(CodeCoverage)
    setup_target_for_coverage(CoreFlows_coverage CoreFlowsCPPTests coverage)
    endif()
endif()

#INSTALL(TARGETS CoreFlowsExampleExe DESTINATION bin/Executable)
#-- Copier les examples dans le CoreFlows_install
#INSTALL(FILES ${CoreFlowsExample_SOURCES} DESTINATION share/examples)
#
#file(GLOB NICE_EXAMPLES_TO_INSTALL resources Python )
#
#install(DIRECTORY ${NICE_EXAMPLES_TO_INSTALL} DESTINATION share/examples)

#####################################

function(CreateTestExec SourceTestFile libList)
     get_filename_component( FILE_BASENAME ${SourceTestFile} NAME_WE) # <path>/testxxx.c --> testxxx
     set( EXECNAME "${FILE_BASENAME}.exe" )                     # testxxx          --> testxxx.exe
     add_executable(${EXECNAME} ${SourceTestFile})                    # compilation of the testxxx.exe 
     set_target_properties(${EXECNAME} PROPERTIES COMPILE_FLAGS "")
     target_link_libraries(${EXECNAME} ${libList})              # provide required lib for testxxx.exe 
     add_test(${FILE_BASENAME} ${EXECNAME} "./${EXECNAME}")     # adding a ctest Test
     message("setting test ${EXECNAME}")
endfunction(CreateTestExec)


#add_executable(test_template.exe test_template.cxx)
#set_target_properties(test_template.exe PROPERTIES COMPILE_FLAGS "")
#target_link_libraries(test_template.exe ${_extra_lib_CoreFlows})

set( libs_for_tests ${_extra_lib_CoreFlows} )

CreateTestExec(test_template.cxx  "${libs_for_tests}" )

CreateTestExec(CoupledTransportDiffusionEquations_1DHeatedChannel.cxx  "${libs_for_tests}" )
CreateTestExec(DiffusionEquation_1DHeatedRod.cxx  "${libs_for_tests}" )
CreateTestExec(DriftModel_1DBoilingAssembly.cxx  "${libs_for_tests}" )
CreateTestExec(DriftModel_1DBoilingChannel.cxx  "${libs_for_tests}" )
CreateTestExec(DriftModel_1DDepressurisation.cxx  "${libs_for_tests}" )
CreateTestExec(DriftModel_1DPorosityJump.cxx  "${libs_for_tests}" )
CreateTestExec(DriftModel_1DPressureLoss.cxx  "${libs_for_tests}" )
CreateTestExec(DriftModel_1DRiemannProblem.cxx  "${libs_for_tests}" )
CreateTestExec(DriftModel_2DInclinedBoilingChannel.cxx  "${libs_for_tests}" )
CreateTestExec(DriftModel_3DCanalCloison.cxx  "${libs_for_tests}" )
CreateTestExec(FiveEqsTwoFluid_1DBoilingChannel.cxx  "${libs_for_tests}" )
CreateTestExec(FiveEqsTwoFluid_1DDepressurisation.cxx  "${libs_for_tests}" )
CreateTestExec(FiveEqsTwoFluid_1DRiemannProblem.cxx  "${libs_for_tests}" )
CreateTestExec(FiveEqsTwoFluid_2DInclinedBoilingChannel.cxx  "${libs_for_tests}" )
CreateTestExec(FiveEqsTwoFluid_2DInclinedSedimentation.cxx  "${libs_for_tests}" )
CreateTestExec(IsothermalTwoFluid_1DDepressurisation.cxx  "${libs_for_tests}" )
CreateTestExec(IsothermalTwoFluid_1DRiemannProblem.cxx  "${libs_for_tests}" )
CreateTestExec(IsothermalTwoFluid_1DSedimentation.cxx  "${libs_for_tests}" )
CreateTestExec(IsothermalTwoFluid_2DInclinedSedimentation.cxx  "${libs_for_tests}" )
CreateTestExec(IsothermalTwoFluid_2DVidangeReservoir.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_1DDepressurisation.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_1DHeatedChannel.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_1DPorosityJump.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_1DRiemannProblem.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_2DHeatDrivenCavity.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_2DHeatDrivenCavity_unstructured.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_2DHeatedChannelInclined.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_2DLidDrivenCavity.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_2DLidDrivenCavity_unstructured.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_2DSphericalExplosion_unstructured.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_2DWallHeatedChannel_ChangeSect.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_2DWallHeatedChannel.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_3DHeatDrivenCavity.cxx  "${libs_for_tests}" )
CreateTestExec(SinglePhase_HeatedWire_2Branches.cxx  "${libs_for_tests}" )
CreateTestExec(TransportEquation_1DHeatedChannel.cxx  "${libs_for_tests}" )

# copy tests resources (med files etc.) into the build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#TODO: install for tests and test data
