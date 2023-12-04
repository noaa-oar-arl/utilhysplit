try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

setup(name='utilhysplit',
      version='1.0.0',
      url='https://github.com/noaa-oar-arl/utilhysplit',
      license='MIT',
      include_package_data=True,
      author='Alice M. Crawford',
      author_email='alice.crawford@noaa.gov',
      maintainer='Alice Crawford',
      maintainer_email='alice.crawford@noaa.gov',
      packages=find_packages(),
      package_data={
          '': [
              'data/*.txt', 'data/*.dat', 'data/*.hdf', 'data/*.ncf',
              'data/*.jpg', 'data/*.png'
          ]
      },
      keywords=[
          'model', 'verification', 'hysplit',
          'evaluation', 'volcat'
      ],
      description='Utilities for use with HYSPLIT and VOLCAT',
      install_requires=[
          'pandas', 'netcdf4', 'xarray', 'matplotlib', 'lxml', 'scipy',
          'seaborn', 'cartopy', 'monetio', 'datetime', 'numpy', 'shapely'
      ])
# extra_requires={'xesmf;platform_system!="Windows"'})

# if notebook or ipython is not installed in the environment, then
# if you run notebook or ipython it will start up but not be using that environment.
# 10/5/2023 venv
#conda install xarray   2023.6.0
#conda install shapely
#conda install netcdf4  1.6.2
#conda install seaborn
#conda install lxml
#conda install cartopy
#conda install notebook
#conda install ipython
#conda install scipy
#conda install -c conda-forge pyresample - needed to use nearest_ij in monetio

# pyresample causes the following changes.
#The following NEW packages will be INSTALLED:

#  configobj          conda-forge/noarch::configobj-5.0.6-py_0
#  llvm-openmp        conda-forge/linux-64::llvm-openmp-12.0.1-h4bd325d_1
#  pykdtree           conda-forge/linux-64::pykdtree-1.3.6-py311h4c7f6c3_2
#  pyresample         conda-forge/linux-64::pyresample-1.25.1-py311h8b32b4d_2
#  python_abi         conda-forge/linux-64::python_abi-3.11-2_cp311

#The following packages will be REMOVED:

#  libgomp-11.2.0-h1234567_1

#The following packages will be UPDATED:

#  libgcc-ng          pkgs/main::libgcc-ng-11.2.0-h1234567_1 --> conda-forge::libgcc-ng-12.2.0-h65d4601_19
#  libstdcxx-ng       pkgs/main::libstdcxx-ng-11.2.0-h12345~ --> conda-forge::libstdcxx-ng-12.2.0-h46fd767_19

#The following packages will be SUPERSEDED by a higher-priority channel:
#
#  _libgcc_mutex           pkgs/main::_libgcc_mutex-0.1-main --> conda-forge::_libgcc_mutex-0.1-conda_forge
#  _openmp_mutex          pkgs/main::_openmp_mutex-5.1-1_gnu --> conda-forge::_openmp_mutex-4.5-2_kmp_llvm
#  ca-certificates    pkgs/main::ca-certificates-2023.08.22~ --> conda-forge::ca-certificates-2022.9.24-ha878542_0
#  certifi            pkgs/main/linux-64::certifi-2023.7.22~ --> conda-forge/noarch::certifi-2022.9.24-pyhd8ed1ab_0




# installation of scikit learn seems to have caused the import of netcdf files with
# time fields to stop working. It may have caused xarray to revert to an older version which didn't work.





