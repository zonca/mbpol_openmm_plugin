## How to build a conda package for a OpenMM plugin

Start a Vagrant virtual machine, in this folder, run:

    vagrant up
    vagrant ssh

This will start a VM based on CentOS 6.6 as we want a old
`glibc` so that our binary can run in older distributions as well.

The VM automatically runs the `setup_centos_vm.sh` to install requirements and conda.

Finally execute the `build_conda_package_vagrant.sh` script line by line, making sure
to update version numbers.
