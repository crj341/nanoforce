import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="nanoforce",
    version="0.0.8",
    author="Chris Jones",
    author_email="crj341@student.bham.ac.uk",
    description="Package to import and analyse AFM force curves produced using Nanoscope 6",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/crj341/nanoforce",
    packages=setuptools.find_packages(),
    install_requires=[
          'numpy',
          'plotly',
          'scipy',
          'easygui',
      ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)