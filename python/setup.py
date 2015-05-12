from setuptools import setup, Extension
import os
import sys

openmm_dir = '@OPENMM_DIR@'
mbpolplugin_header_dir = '@MBPOLPLUGIN_HEADER_DIR@'
mbpolplugin_library_dir = '@MBPOLPLUGIN_LIBRARY_DIR@'

extension = Extension(name='_mbpolplugin',
                      sources=['MBPolPluginWrapper.cpp'],
                      libraries=['OpenMM', 'OpenMMMBPol'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), mbpolplugin_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), mbpolplugin_library_dir]
                     )

setup(name='mbpol',
      version='1.0',
      py_modules=['mbpol', 'mbpolplugin'],
      ext_modules=[extension],
      scripts = ['bin/run_mbpol'],
      data_files=[('', ['mbpol.xml'])],
)
