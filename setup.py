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
          'seaborn', 'cartopy', 'monetio', 'datetime', 'numpy'
      ])
# extra_requires={'xesmf;platform_system!="Windows"'})
