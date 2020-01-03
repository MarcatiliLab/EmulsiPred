import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    install_requires = f.read().strip().split('\n')

setuptools.setup(
    name="EmulsiPred-pkg-pamar",  # Name of package
    version="0.0.2",
    author="Paolo Marcatili, Tobias Olsen, Egon Hansen",
    author_email="pamar@dtu.dk",
    description="A package to predict emulsifying potential of peptides",
    long_description=long_description,  # README.md file as description
    long_description_content_type="text/markdown",
    url="https://github.com/MarcatiliLab/sampleproject",
    packages=setuptools.find_packages(),
    install_requires=install_requires,  # Install requirements extracted from requirements.txt
    include_package_data=True,  # Allow to include other files than .py in package
    package_data={
        '': ['NormalizationValues/*.csv']
    },  # Define which additional files should be included in package
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)