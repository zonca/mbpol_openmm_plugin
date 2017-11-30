## How to build a conda package for a OpenMM plugin

Start a Vagrant virtual machine, in this folder, run:

    vagrant up
    vagrant ssh

This will start a VM based on CentOS 6.6 as we want a old
`glibc` so that our binary can run in older distributions as well.

The VM automatically runs the `setup_centos_vm.sh` to install requirements and conda.

Finally execute the `build_conda_package_vagrant.sh` script line by line, making sure
to update version numbers.
`conda-build` works only on the root environment, so we install also `miniconda3` separately.

so we have to set the correct miniconda path in `build.sh`, check mbpol and openmm versions in `meta.yaml` and then repeat the build and upload process for python 2 and 3.

For uploading, add `-u paesanilab` to the anaconda upload command to upload to the Organization
