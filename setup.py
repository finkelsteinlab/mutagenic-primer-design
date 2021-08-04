# Writing the setup script: https://pythonhosted.org/an_example_pypi_project/setuptools.html
import os
from setuptools import setup, find_packages

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "spikedisplay",
    version = "1.0.0",
    author = "John P. Cooper",
    author_email = "",
    description = (""),
    license = "BSD",
    keywords = "tag1 tag2 tag3",
    url = "https://github.com/johnpcooper",
    packages=find_packages(),
    package_data={'spikedisplay': ['*.csv', '*.json']},
    include_package_data=True,
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
)
