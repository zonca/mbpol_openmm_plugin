export CC="clang++"

conda install --yes -c https://conda.binstar.org/omnia openmm

git clone https://github.com/paesanilab/mbpol_openmm_plugin.git
cd mbpol_openmm_plugin/
git checkout tags/v1.0.0-alpha

cd conda
conda build conda-recipe

# To upload the file, do something the following command but with the package version changed:

binstar upload -u paesanilab /home/vagrant/miniconda/conda-bld/linux-64/
