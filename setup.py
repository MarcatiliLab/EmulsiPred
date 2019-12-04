import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="EmulsiPred-pkg-pamar", # Replace with your own username
    version="0.0.1",
    author="Paolo Marcatili, Tobias Olsen, Egon Hansen",
    author_email="pamar@dtu.dk",
    description="A package to predict emulsifying potential of peptides",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MarcatiliLab/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)