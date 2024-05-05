###############################################################################
"""
    Last edited on August 15, 2019

    @author: matz

    comments: setup.py file for the nwpy module/package; copied and
    modified from:
    https://github.com/mwaskom/seaborn/blob/master/setup.py
    https://github.com/mitchnegus/seaborn/blob/master/setup.py
    
    For package organization see
    https://jeffknupp.com/blog/2013/08/16/open-sourcing-a-python-project-the-right-way/
    https://github.com/kennethreitz/samplemod
    https://docs.python.org/2.4/dist/node11.html
    https://pypi.org/project/python_boilerplate_template/
    
    
"""
###############################################################################
import os
import sys
import pkg_resources
import nwpy
from distutils.core import setup
###############################################################################


DESCRIPTION = "nwpy: Evaluate nuclear fuel cycle waste streams with python"
LONG_DESCRIPTION = """nwpy is a project to characterize nuclear fuel cycle
waste streams and evaluate waste management performance metrics
"""
DISTNAME = 'nwpy'
MAINTAINER = 'Milos Atz'
MAINTAINER_EMAIL = 'milos.atz@berkeley.edu'
URL = 'https://github.com/MilosAtz/nwpy'
LICENSE = 'MIT'
DOWNLOAD_URL = 'https://github.com/MilosAtz/nwpy'# github
VERSION = '0.0.3'


###############################################################################


def check_dependencies():
    install_requires = []
    # Just make sure dependencies exist, I haven't rigorously
    # tested what are the minimal versions that will work
    # (help on that would be awesome)
    # Built-in (standard library):
    #   os, collections, datetime, decimal, math, subprocess, imp
    try:
        import numpy
    except ImportError:
        install_requires.append('numpy')
    try:
        import scipy
    except ImportError:
        install_requires.append('scipy')
    try:
        import matplotlib
    except ImportError:
        install_requires.append('matplotlib')
    try:
        import pandas
    except ImportError:
        install_requires.append('pandas')
    try:
        import pytest
    except ImportError:
        install_requires.append('pytest')
    return install_requires


###############################################################################


if __name__ == "__main__":
    install_requires = check_dependencies()
    setup(name=DISTNAME,
          author=MAINTAINER,
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
          packages=['nwpy',
                    'nwpy.repository_area',
                    # 'nwpy.attractiveness',
                    'nwpy.criticality'],
          package_dir={'nwpy': 'nwpy',
                       'nwpy.repository_area': 'nwpy/repository_area',
                       # 'nwpy.attractiveness': 'nwpy/attractiveness',
                       'nwpy.criticality': 'nwpy/criticality'
                       },
          package_data={'nwpy': ['data/*.py',
                                 'data/*.csv',
                                 'data/fc/*.fc',
                                 'data/iso/*.csv',
                                 'data/sep/*.sep',
                                 'data/sep/*.py',
                                 'data/load/*.py',
                                 'tests/*.py',
                                 'tests/testdata/*.inp',
                                 'tests/testdata/*.plt',
                                 'tests/testdata/fc/*.fc',
                                 'tests/testdata/iso/*.csv',
                                 'tests/testdata/sep/*.sep',
                                 'tests/testdata/sep/*.py',
                                 'tests/testdata/load/*.py'
                                 ],
                        'nwpy.repository_area':['data/*.py',
                                                'data/waste/*.csv',
                                                'tests/*.py',
                                                'tests/testdata/*.csv'
                                                ],
                        # 'nwpy.attractiveness':['tests/*.py'
                        },
          include_package_data=True,
          classifiers=['Intended Audience :: Science/Research',
                       'Programming Language :: Python :: 2.7',
                       'License :: OSI Approved :: MIT License',
                       'Topic :: Scientific/Engineering :: Nuclear Energy',
                       'Topic :: Nuclear Waste :: Nuclear Fuel Cycle',
                       'Operating System :: Unix',
                       'Operating System :: MacOS'],
          )

          
###############################################################################