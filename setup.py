import os
from setuptools import setup, find_packages


def get_version():
    with open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "zamp",
            "zamp.VERSION",
        )
    ) as f:
        return f.readline().strip()


def get_description():
    with open("README.md", "r") as fh:
        long_description = fh.read()
    return long_description


def get_data_files():
    data_files = [(".", ["README.md"])]
    return data_files


CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="zamp",
    packages=find_packages(),
    url="https://github.com/metagenlab/zAMP",
    python_requires=">=3.10",
    description="Snakemake pipeline designed for convenient, reproducible and scalable amplicon-based metagenomics",
    long_description=get_description(),
    long_description_content_type="text/markdown",
    version=get_version(),
    author="Valentin Scherz",
    author_email="valentin.scherz@chuv.ch",
    maintainer="Farid Chaabane",
    maintainer_email="farid.chaabane@chuv.ch",
    data_files=get_data_files(),
    py_modules=["zamp"],
    install_requires=[
        "snakemake>=8.10.6",
        "Click>=8.1.3",
        "attrmap>=0.0.7",
        "snaketool-utils>=0.0.5",
        "metasnek>=0.0.8",
        "biopython>=1.83",
    ],
    entry_points={"console_scripts": ["zamp=zamp.__main__:main"]},
    include_package_data=True,
)
