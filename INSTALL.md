## Requirements

* `OpenMM` version 6.2 or later. 

## How to install

* Download the last release from Github <https://github.com/paesanilab/mbpol_openmm_plugin/releases>
* Create the `build_mbpol` folder outside of the source folder
* Configure the build by entering the `build_mbpol` folder and running `ccmake -i ../mbpol_openmm_plugin` (`ccmake` with 2 time `c` is a console-based GUI for `cmake`, in debian/ubuntu it is in the `cmake-curses-gui` package)
* Press `c` to start the configuration process
* Set:
  * `OPENMM_BUILD_MBPOL_CUDA_LIB`  `OFF` (not supported yet)
  * `CMAKE_INSTALL_PREFIX` and `OPENMM_DIR` should contain the path to the installed `OpenMM`, by default `/usr/local/openmm`.
  * `CMAKE_BUILD_TYPE` `Debug` (Otherwise the compiler takes a very long time to compile the large polynomials)
* Press `c` again to configure
* Press `g` to generate the configuration and exit
* Run `make` to compile the C++ library
* Run `(sudo) make install` to install the `mbpol` dynamic library to the
  `openmm` folder
* Run `make PythonInstall` to install the Python wrapper, it requires
  Python and `swig`, the best is to use the Anaconda Python distribution
* Add the OpenMM lib folder to the dynamic libraries path, generally add to `.bashrc`: `export LD_LIBRARY_PATH=/usr/local/openmm/lib:/usr/local/openmm/lib/plugins:$LD_LIBRARY_PATH` and restart `bash`

## Run unit tests

* run `make test` to run the unit test suite

## Run a test simulation

* Enter in the `python` sub-folder in the source folder
* Run `python water14.py` to run a test simulation for a cluster of 14 water molecules
