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

ADD_LIBRARY(CoreFlowsExample SHARED ${CoreFlowsExample_SOURCES})
SET(_extra_lib_CoreFlows CoreFlows ${PETSC_LIBRARIES} ${CDMATH_LIBRARIES})
TARGET_LINK_LIBRARIES(CoreFlowsExample ${_extra_lib_CoreFlows})
INSTALL(TARGETS CoreFlowsExample DESTINATION lib)

SET(MAIN_SOURCES
  main_tests.cxx
  )

ADD_EXECUTABLE(CoreFlowsExampleExe ${MAIN_SOURCES})
SET_TARGET_PROPERTIES(CoreFlowsExampleExe PROPERTIES COMPILE_FLAGS "")
TARGET_LINK_LIBRARIES(CoreFlowsExampleExe ${_extra_lib_CoreFlows} CoreFlowsExample)
ADD_TEST(NAME CoreFlowsCPPTests COMMAND CoreFlowsExampleExe)

if(CMAKE_COMPILER_IS_GNUCXX)
    if (CMAKE_BUILD_TYPE=Debug)
    include(CodeCoverage)
    setup_target_for_coverage(CoreFlows_coverage CoreFlowsCPPTests coverage)
    endif()
endif()

INSTALL(TARGETS CoreFlowsExampleExe DESTINATION bin/Executable)
#-- Copier les examples dans le CoreFlows_install
INSTALL(FILES ${CoreFlowsExample_SOURCES} DESTINATION share/examples)

file(GLOB NICE_EXAMPLES_TO_INSTALL resources Python )

install(DIRECTORY ${NICE_EXAMPLES_TO_INSTALL} DESTINATION share/examples)


