set(WRAP_FILE MBPolPluginWrapper.cpp)
set(MODULE_NAME mbpolplugin)

# Execute SWIG to generate source code for the Python module.

add_custom_command(
    OUTPUT "${WRAP_FILE}"
    COMMAND "${SWIG_EXECUTABLE}"
        -python -c++
        -o "${WRAP_FILE}"
        "-I${OPENMM_DIR}/include"
        "-DOPENMM_MAJOR_VERSION=${OPENMM_MAJOR_VERSION}"
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

#copy python test files
file(GLOB PYTESTS ${CMAKE_CURRENT_SOURCE_DIR}/ tests/*.py)
foreach(PYTEST ${PYTESTS})
    get_filename_component(FILENAME ${PYTEST} NAME)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests/${FILENAME} ${CMAKE_CURRENT_BINARY_DIR}/tests/${FILENAME})
endforeach()
#copy pdb files for python tests
file(GLOB PDBS ${CMAKE_CURRENT_SOURCE_DIR}/ tests/pdb_files/*.pdb)
foreach(PDB ${PDBS})
    get_filename_component(FILENAME ${PDB} NAME)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests/pdb_files/${FILENAME} ${CMAKE_CURRENT_BINARY_DIR}/tests/pdb_files/${FILENAME})
endforeach()
#copy py.test infile
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests/pytest.ini ${CMAKE_CURRENT_BINARY_DIR}/tests/pytest.ini)
