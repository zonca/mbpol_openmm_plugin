
## How to install the binary distribution for GNU/Linux-64bit with conda

Only for GNU/Linux 64bit a binary distribution is available through `conda`.

* Install the Anaconda Python distribution from <http://continuum.io/downloads>. Either Python 3 or Python 2.7 would work, Python 3 is recommended.
* Add `conda` channels:

        conda config --add  -c https://conda.anaconda.org/omnia
        conda config --add  -c https://conda.anaconda.org/paesanilab
* Install OpenMM and OpenMM-MBPol with `conda`: `conda install openmm mbpol`

## How to install from source

### Requirements

* Install OpenMM 6.2 or later, either a binary or source release, from <https://simtk.org/project/xml/downloads.xml?group_id=161>
* Install Python 3 (2.7 also works), Anaconda is recommended (<http://continuum.io/downloads>)
* Install `swig`, `conda install swig` with Anaconda (do not mix Anaconda Python with `swig` installed with other methods) or install the `swig` package with your package manager.

### Installation of OpenMM-MBPol

* Download the last release from Github <https://github.com/paesanilab/mbpol_openmm_plugin/releases> and uncompress it
* Create a `build_mbpol` folder outside of the just uncompressed source folder (`mbpol_openmm_pluginXXXX`)
* Configure the build by entering the `build_mbpol` folder and running `ccmake -i ../mbpol_openmm_pluginXXXX` (`ccmake` with 2 time `c` is a console-based GUI for `cmake`, in debian/ubuntu it is in the `cmake-curses-gui` package)
* Press `c` to start the configuration process
* Set:
  * `MBPOL_BUILD_PYTHON_WRAPPERS`  `ON` in order to build the Python wrappers (necessary to use `mbpol_builder`)
  * `CMAKE_INSTALL_PREFIX` and `OPENMM_DIR` should contain the path to the installed `OpenMM`, by default both `/usr/local/openmm`.
  * `CMAKE_BUILD_TYPE` `Debug` (Otherwise the compiler takes a very long time to compile the large polynomials)
  * `OPENMM_MAJOR_VERSION` and `OPENMM_MINOR_VERSION` based on the version of `OpenMM` you are compiling for. If you installed the OpenMM Python wrapper, you can print the version running: `python -c 'from simtk import openmm; print(openmm.version.short_version)'`.
* Press `c` again to configure
* Press `g` to generate the configuration and exit
* Run `make` to compile the C++ library
* Run `(sudo) make install` to install the `mbpol` dynamic library to the
  `openmm` folder
* Run `make PythonInstall` to install the Python wrapper, it requires
  Python and `swig`, the best is to use the Anaconda Python distribution
* Add the OpenMM lib folder to the dynamic libraries path, generally add to `.bashrc`: `export LD_LIBRARY_PATH=/usr/local/openmm/lib:/usr/local/openmm/lib/plugins:$LD_LIBRARY_PATH` and restart `bash`
* You can run `make test` to run the C++ unit test suite

## After install

Once the installation process is completed, the OpenMM and the MBPol plugins libraries are available to your Python installation, from here you can check `README.md` on how to run an example simulation or create a sample MBPol Python script with `mbpol_builder`.
