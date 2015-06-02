# Prepare a vagrant CentOS 6.5 VM for building OpenMM
# Needs latest version of vagrant to auto-download the chef package
#vagrant init chef/centos-6.5
#vagrant up
#vagrant ssh

# Download and enable the EPEL RedHat EL extras repository
echo "********** Enabling EPEL RedHat EL extras repository..."
mkdir ~/Software
cd Software
sudo yum install wget -y
wget http://dl.fedoraproject.org/pub/epel/6/i386/epel-release-6-8.noarch.rpm
sudo rpm -i epel-release-6-8.noarch.rpm

# Install things needed for virtualbox guest additions.
echo "********** Removing old kernel package..."
sudo yum remove kernel-2.6.32-431.el6.x86_64
echo "********** Installing kernel headers and dkms for virtualbox guest additions..."
sudo yum install kernel-devel dkms -y

echo "********** Updating yum..."
sudo yum update -y

# Several of these come from the EPEL repo
echo "********** Installing lots of packages via yum..."
sudo yum install tar clang cmake graphviz perl flex bison rpm-build texlive texlive-latex ghostscript gcc gcc-c++ git vim emacs swig zip sphinx python-sphinx -y
# Note: changed from clang-3.4 to clang because the package has apparently been renamed.  KAB Oct 2 2014.
sudo yum install -y screen

### # Probably can't use RHEL6 version of doxygen because it's very old.
### echo "********** Compiling recent doxygen..."
### cd ~/Software
### wget http://ftp.stack.nl/pub/users/dimitri/doxygen-1.8.8.src.tar.gz
### sudo yum remove doxygen -y # Remove yum version!  Necessary as otherwise might not overwrite.
### rpmbuild -ta doxygen-1.8.8.src.tar.gz --nodeps  # Use nodeps because we will eventually use a non-standard texlive installation with no RPM
### sudo yum install -y ~/rpmbuild/RPMS/x86_64/doxygen-1.8.8-1.x86_64.rpm
### echo "exclude=doxygen" | sudo tee --append /etc/yum.conf  # The hand-built RPM package has the wrong versioning scheme and by default will be overwritten by a yum update.  This prevents overwriting.
### doxygen --version  # Should be 1.8.8
### # Although we needed to have the redhat texlive installed to build the rpm, but we can delete the redhat texlive now.
### sudo yum remove texlive texlive-latex -y  # Get rid of the system texlive in preparation for latest version.
### 
### 
### # We have to install a modern texlive 2014 distro, since the yum-installable version is missing vital components.
### echo "********** Installing texlive 2014..."
### wget http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz
### tar zxf install-tl-unx.tar.gz
### cd install-tl-*
### sudo ./install-tl -profile /vagrant/texlive.profile
### export PATH=/usr/local/texlive/2014/bin/x86_64-linux:$PATH  # texlive updates bashrc to put tex on the path, but we need to update the current shell session.
### sleep 2
### 
### # Make sure texlive install worked, as it often dies.  Only retry once, though.
### if which tex >/dev/null; then
###     echo Found texlive
### else
###     echo No texlive, resuming installation
###     sudo ./install-tl -profile /vagrant/texlive.profile
### fi

cd ..




### # Install fortran
### echo "********** Installing fortran..."
### sudo yum install gcc-gfortran -y # Used for ambermini
### 
### sudo yum install -y lapack-devel  # Used for cvxopt.  This also grabs BLAS dependency.
### 
### sudo yum clean headers
### sudo yum clean packages

### # Install CUDA6.5 for RHEL6
### echo "********** Installing CUDA..."
### cd ~/Software
### wget http://developer.download.nvidia.com/compute/cuda/repos/rhel6/x86_64/cuda-repo-rhel6-6.5-14.x86_64.rpm
### sudo rpm -i  cuda-repo-rhel6-6.5-14.x86_64.rpm
### sudo yum clean expire-cache
### sudo yum install cuda -y
### # NOTE: NVIDIA may push new MAJOR release versions of CUDA without warning.
### # This is even *before* doing the below update.  Beware.
### 
### echo "********** Forcing a second yum update in case CUDA has patches..."
### sudo yum update -y  # Force a second update, in case CUDA has necessary patches.

# Install Conda
echo "********** Installing conda..."
cd ~/Software
MINICONDA=Miniconda-latest-Linux-x86_64.sh
MINICONDA_MD5=$(curl -s http://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget http://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
echo "Miniconda MD5 mismatch"
exit 1
fi
bash $MINICONDA -b

# So there is a bug in some versions of anaconda where the path to swig files is HARDCODED.  Below is workaround.  See https://github.com/ContinuumIO/anaconda-issues/issues/48
sudo ln -s  ~/miniconda/ /opt/anaconda1anaconda2anaconda3

echo "********** Installing conda/binstar channels and packages..."
export PATH=$HOME/miniconda/bin:$PATH
conda config --add channels http://conda.binstar.org/omnia
conda install --yes fftw3f jinja2 swig sphinx conda-build cmake binstar pip

echo "********** Installing conda/binstar channels and packages..."
conda create -n py3 --yes python=3.4 anaconda fftw3f jinja2 swig sphinx conda-build cmake binstar pip

# Add conda to the path.
echo "********** Adding paths"
cd ~
echo "export PATH=$HOME/miniconda/bin:/usr/local/texlive/2014/bin/x86_64-linux:$PATH" >> $HOME/.bashrc
echo "" >> $HOME/.bashrc

### # Install additional packages via pip.
### echo "********** Installing packages via pip..."
### $HOME/miniconda/bin/pip install sphinxcontrib-bibtex

