## How to build a conda package for a OpenMM plugin

Start a Vagrant virtual machine, in this folder, run:

    vagrant up
    vagrant ssh

This will start a VM based on CentOS 6.6 as we want a old
`glibc` so that our binary can run in older distributions as well.

Then within the VM we can enter the folder which is shared with the host FIXME and run:

    bash setup_centos_vm.sh

Finally execute the `build_conda_package_vagrant.sh` script line by line, making sure
to update version numbers.
