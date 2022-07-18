import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="fpbot",
    version="1.0.4",
    author="Simeon Reusch",
    author_email="simeon.reusch@desy.de",
    description="ZTF Forced Photometry Pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="BSD (3-clause)",
    keywords="astrophysics astronomy ZTF",
    url="https://github.com/simeonreusch/fpbot",
    packages=setuptools.find_packages(),
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8.0",
    package_data={
        "fpbot": [
            "data/*.dat",
            "data/*.jpg",
        ]
    },
    scripts=[
        "forcedphotometry",
        "fps_slackbot",
        "fps_start_slackbot",
        "fps_slackbot_spawn_screen_session",
    ],
    install_requires=[
        "astropy",
        "numpy",
        "pandas",
        "matplotlib",
        "coveralls",
        "scipy",
        "slackclient",
        "sqlalchemy",
        "requests",
        "backoff",
        "lxml",
        "html5lib",
        "bs4",
        "iminuit",
        "pymongo",
        "keyring",
        "ztflc",
        "ztfquery",
        "tqdm",
    ],
)
