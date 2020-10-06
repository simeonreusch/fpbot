DESCRIPTION = "ZTF Forced Photometry Pipeline"
LONG_DESCRIPTION = """Provides a Forced Photometry Pipeline for ZTF based on ztfquery and ztflc, needs IPAC as well as Marshal or AMPEL access."""

DISTNAME = "ztffps"
AUTHOR = "Simeon Reusch"
MAINTAINER = "Simeon Reusch"
MAINTAINER_EMAIL = "simeon.reusch@desy.de"
URL = "https://github.com/simeonreusch/modelSED/"
LICENSE = "BSD (3-clause)"
DOWNLOAD_URL = "https://github.com/simeonreusch/ztffps/archive/v1.0.2.tar.gz"
VERSION = "1.0.2"

try:
    from setuptools import setup, find_packages

    _has_setuptools = True
except ImportError:
    from distutils.core import setup

    _has_setuptools = False


def check_dependencies():
    install_requires = []

    # Make sure dependencies exist. This is ongoing
    deps = [
        astropy,
        numpy,
        sncosmo,
        extinction,
        pandas,
        matplotlib,
        scipy,
        slackclient,
        sqlalchemy,
        requests,
        lxml,
        html5lib,
        bs4,
        psycopg2,
        iminuit,
        sfdmap,
        pymongo,
        fastapi,
        uvicorn,
        keyring,
        ztflc,
        ztfquery,
    ]

    for dep in deps:
        try:
            import dep
        except ImportError:
            install_requires.append(str(dep))

    # try:
    #     import astropy
    # except ImportError:
    #     install_requires.append("astropy")
    # try:
    #     import numpy
    # except ImportError:
    #     install_requires.append("numpy")
    # try:
    #     import sncosmo
    # except ImportError:
    #     install_requires.append("sncosmo")
    # try:
    #     import extinction
    # except ImportError:
    #     install_requires.append("extinction")
    # try:
    #     import pandas
    # except ImportError:
    #     install_requires.append("pandas")
    # try:
    #     import matplotlib
    # except ImportError:
    #     install_requires.append("matplotlib")
    # try:
    #     import scipy
    # except ImportError:
    #     install_requires.append("scipy")
    # try:
    #     import slackclient
    # except ImportError:
    #     install_requires.append("slackclient")
    # try:
    #     import sqlalchemy
    # except ImportError:
    #     install_requires.append("sqlalchemy")
    # try:
    #     import requests
    # except ImportError:
    #     install_requires.append("requests")
    # try:
    #     import lxml
    # except ImportError:
    #     install_requires.append("lxml")
    # try:
    #     import html5lib
    # except ImportError:
    #     install_requires.append("html5lib")
    # try:
    #     import bs4
    # except ImportError:
    #     install_requires.append("bs4")
    # try:
    #     import psycopg2
    # except ImportError:
    #     install_requires.append("psycopg2")
    # try:
    #     import iminuit
    # except ImportError:
    #     install_requires.append("iminuit")
    # try:
    #     import sfdmap
    # except ImportError:
    #     install_requires.append("sfdmap")
    # try:
    #     import pymongo
    # except ImportError:
    #     install_requires.append("pymongo")
    # try:
    #     import fastapi
    # except ImportError:
    #     install_requires.append("fastapi")
    # try:
    #     import uvicorn
    # except ImportError:
    #     install_requires.append("uvicorn")
    # try:
    #     import keyring
    # except ImportError:
    #     install_requires.append("keyring")
    # try:
    #     import ztflc
    # except ImportError:
    #     install_requires.append("ztflc")
    # try:
    #     import ztfquery
    # except ImportError:
    #     install_requires.append("ztfquery")

    print(install_requires)

    return install_requires


if __name__ == "__main__":

    install_requires = check_dependencies()

    if _has_setuptools:
        packages = find_packages()
        print(packages)
    else:
        packages = ["ztffps"]

    setup(
        name=DISTNAME,
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
        scripts=["forcedphotometry"],
        packages=packages,
        package_data={
            "ztffps": [
                "data/*.dat",
                "data/*.jpg",
            ]
        },
        classifiers=[
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 3.7",
            "License :: OSI Approved :: BSD License",
            "Topic :: Scientific/Engineering :: Astronomy",
            "Operating System :: POSIX",
            "Operating System :: Unix",
            "Operating System :: MacOS",
        ],
    )
