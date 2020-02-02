#!/usr/bin/env python3

DESCRIPTION = " Pipeline for ZTF forced photometry "
LONG_DESCRIPTION = """ Pipeline for ZTF forced photometry """

DISTNAME = 'ztffps'
AUTHOR = 'Simeon Reusch'
MAINTAINER = 'Simeon Reusch' 
MAINTAINER_EMAIL = 'simeon.reusch@desy.de'
URL = 'https://github.com/simeonreusch/ztffps'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://github.com/MickaelRigault/ztffps/tarball/0.1'
VERSION = '0.1'

try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup
    _has_setuptools = False

def check_dependencies():
    install_requires = []
    try:
        import ztfquery
    except ImportError:
        install_requires.append('ztfquery')

    try:
        import pandas
    except ImportError:
        install_requires.append('pandas')

    try:
        import ztflc
    except ImportError:
        install_requires.append('ztflc')

    try:
        import sqlalchemy
    except ImportError:
        install_requires.append('sqlalchemy')
        
   return install_requires

if __name__ == "__main__":

    install_requires = check_dependencies()

    if _has_setuptools:
        packages = find_packages()
        print(packages)
    else:
        # This should be updated if new submodules are added
        packages = ['ztffps']
    
        
    setup(name=DISTNAME,
          author=AUTHOR,
          author_email=MAINTAINER_EMAIL,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          long_description=LONG_DESCRIPTION,
          license=LICENSE,
          url=URL,
          version=VERSION,
          download_url=DOWNLOAD_URL,
          install_requires=install_requires,
          packages=packages,
          package_data={'ztfquery': ['data/*']},
          classifiers=[
              'Intended Audience :: Science/Research',
              'Programming Language :: Python :: 3.7',
              'License :: OSI Approved :: BSD License',
              'Topic :: Scientific/Engineering :: Astronomy',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS'],
      )