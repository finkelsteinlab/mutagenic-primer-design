# `spikedisplay`

Oligo library generation for deep mutational scanning and associated analysis

## Installation

In command line, navigate to the directory where you want to install the project and clone the repo

```sh
cd /path/to/directory
git clone https://github.com/johnpcooper/spikedisplay.git
```

Install `virtualenv` and create an environment inside the `spikedisplay` project directory

```sh
cd spikedisplay
pip install virtualenv
python -m venv .spikedisplay
```

Activate the environment

```sh
# Linux
source .spikedisplay/bin/activate
```

Install required packages

```sh
pip install -r requirements.txt
```

Install `spikedisplay` package
```sh
python setup.py install
```

Setup `jupyter notebook` compatability

```sh
# Add the .spikedisplay environment to the kernels available in jupyter
python -m ipykernel install --user --name=.spikedisplay
```

## Primer library generation notebook

Navigate to the spikedisplay/notebooks directory and open the `Primer_Library_Generation` notebook

```sh
cd notebooks
jupyter notebook Primer_Library_Generation.ipynb
```