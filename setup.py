import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sbm_utilities",
    version="0.0.1",
    author="Ana Robles, Martin Floor",
    author_email="anaroblesmar@gmail.com, martinfloor@gmail.com",
    description="A Python package to deploy SBM simulations with SBMOpenMM",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    install_requires=[],
)
