project(test)

file(GLOB MESHES_TO_INSTALL resources )
install(DIRECTORY ${MESHES_TO_INSTALL} DESTINATION share/examples)

add_subdirectory (C)

if (COREFLOWS_WITH_PYTHON )

  add_subdirectory (Python)

endif (COREFLOWS_WITH_PYTHON )


