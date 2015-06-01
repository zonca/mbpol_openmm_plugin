export CC="clang++"

conda install --yes -c https://conda.binstar.org/omnia openmm

git clone https://github.com/paesanilab/mbpol_openmm_plugin.git
cd mbpol_openmm_plugin/
git checkout tags/v1.0.0-alpha

# For CI build:
conda build devtools/conda-recipe

# For release build:
cd ../
git clone https://github.com/omnia-md/conda-recipes.git
conda build conda-recipes/openmm

# To upload the file, do something the following command but with the package version changed:

binstar upload -u omnia /home/vagrant/miniconda/conda-bld/linux-64/openmm-6.1-py27_0.tar.bz2
