
############# python convergence tests

function(CreatePythonConvergenceTest pythonFile)
    file(GLOB pythonTestFiles "*.py" )
    file(COPY ${pythonTestFiles} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE )
    add_test(${pythonFile} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${pythonFile})
    SET_PROPERTY(TEST ${pythonFile} PROPERTY ENVIRONMENT LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CDMATH_DIR}/lib:${CDMATH_DIR}/lib/med:${CDMATH_DIR}/lib/medcoupling:${PETSC_DIR}/${PETSC_ARCH}/lib)
    SET_PROPERTY(TEST ${pythonFile} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib/CoreFlows_Python:${CMAKE_INSTALL_PREFIX}/bin/CoreFlows_Python:${CDMATH_DIR}/lib:${CDMATH_DIR}/lib/cdmath:${CDMATH_DIR}/bin/cdmath:${CDMATH_DIR}/lib/med:${CDMATH_DIR}/lib/medcoupling)
endfunction(CreatePythonConvergenceTest)



CreatePythonTest(DiffusionEquation_1DHeatedRod.py)

