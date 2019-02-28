import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="zdock",
    version="0.0.1",
    author="Arjan van der Velde",
    author_email="vandervelde.ag@gmail.com",
    description="Parser for ZDOCK and M-ZDOCK output files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/weng-lab/libpdb.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

