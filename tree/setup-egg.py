from setuptools import setup, find_packages
import os, sys
from os.path import join as pj
 
packagename = 'tree'
namespace = 'openalea'
build_prefix = "build-scons"

# Scons build directory
scons_parameters=["build_prefix="+build_prefix]


# dependencies
install_requires = []
setup_requires = ['openalea.deploy']
if("win" in sys.platform):
    setup_requires += ['MinGW','boostpython','openalea.sconsx']
    install_requires += ['boostpython']
else:
    install_requires = []
    setup_requires = []


if __name__ == '__main__':
    
    setup(name='VPlants.Tree',
          version='0.1',
          author='C. Pradal',
          description='Tree structure header',
          url='http://openalea.gforge.inria.fr',
          license='GPL',
 
          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters  
          scons_parameters=scons_parameters,
        
          # Packages
          #packages=
          #package_dir=
      
          # Add package platform libraries if any
          include_package_data=True,
          zip_safe = False,

          # Specific options of openalea.deploy
          #lib_dirs = {'lib' : pj(build_prefix, 'lib'),},
          inc_dirs = { 'include' : pj(build_prefix, 'include') },
          

          # Dependencies
          setup_requires = setup_requires,
          install_requires = install_requires,
          dependency_links = ['http://openalea.gforge.inria.fr/pi'],
          )

