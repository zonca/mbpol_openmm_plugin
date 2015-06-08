set(WRAP_FILE MBPolPluginWrapper.cpp)
set(MODULE_NAME mbpolplugin)

# Execute SWIG to generate source code for the Python module.

add_custom_command(
    OUTPUT "${WRAP_FILE}"
    COMMAND "${SWIG_EXECUTABLE}"
        -python -c++
        -o "${WRAP_FILE}"
        "-I${OPENMM_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/mbpolplugin.i"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/mbpolplugin.i"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

# Compile the Python module.

add_custom_target(PythonInstall DEPENDS "${WRAP_FILE}")
set(MBPOLPLUGIN_HEADER_DIR "${CMAKE_SOURCE_DIR}/openmmapi/include")
set(MBPOLPLUGIN_LIBRARY_DIR "${CMAKE_BINARY_DIR}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mbpol.py ${CMAKE_CURRENT_BINARY_DIR}/mbpol.py COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mbpol.xml ${CMAKE_CURRENT_BINARY_DIR}/mbpol.xml COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST.in ${CMAKE_CURRENT_BINARY_DIR}/MANIFEST.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bin/mbpol_builder ${CMAKE_CURRENT_BINARY_DIR}/bin/mbpol_builder COPYONLY)
add_custom_command(TARGET PythonInstall
    COMMAND "${PYTHON_EXECUTABLE}" setup.py build
    COMMAND "${PYTHON_EXECUTABLE}" setup.py install
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
