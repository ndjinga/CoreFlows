
SET(SCRIPT
  ../validationStationaryDiffusionEquation.py
  convergence_StationaryDiffusion_2DFE_Dirichlet_DelaunayTriangles.py
  )

FILE(GLOB MY_MESHES_AND_PICTURES  ${CDMATH_DIR}/share/meshes/2DTriangles/squareWithTriangles_* )

file(COPY ${SCRIPT} ${MY_MESHES_AND_PICTURES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if (COREFLOWS_WITH__WITH_PYTHON )

    ADD_TEST(convergence_StationaryDiffusion_2DFE_Dirichlet_DelaunayTriangles ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/convergence_StationaryDiffusion_2DFE_Dirichlet_DelaunayTriangles.py)

    SET_PROPERTY(TEST convergence_StationaryDiffusion_2DFE_Dirichlet_DelaunayTriangles PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/medcoupling:${CMAKE_INSTALL_PREFIX}/lib/med:${PETSC_DIR}/${PETSC_ARCH}/lib:${HDF5_ROOT_DIR}/lib:${PV_LIB_DIR}")

    SET_PROPERTY(TEST convergence_StationaryDiffusion_2DFE_Dirichlet_DelaunayTriangles APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_INSTALL_PREFIX}/lib/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath:${CMAKE_INSTALL_PREFIX}/bin/cdmath/postprocessing:${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON}:${PV_PYTHON_DIR}")

endif (COREFLOWS_WITH__WITH_PYTHON )


