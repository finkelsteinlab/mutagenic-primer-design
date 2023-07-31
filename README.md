# `mutagenic-primer-design`

Pooled oligo library generation for primer extension-based deep mutational scanning

## Installation

In command line, navigate to the directory where you want to install the project and clone the repo

```sh
# In this example we're cloning the repo into
# home directory
cd ~
git clone https://github.com/johnpcooper/mutagenic-primer-design.git
```

Install `virtualenv` and create an environment inside the `mutagenic-primer-design` project directory

```sh
cd mutagenic-primer-design
pip install virtualenv
python -m venv .spikedisplay
```

Activate the environment

```sh
source .spikedisplay/bin/activate
```

Install required packages to the environment

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

This [notebook](https://github.com/johnpcooper/mutagenic-primer-design/blob/main/notebooks/Primer_Library_Generation.ipynb) contains an example workflow generating a codon-optimized primer library for PCR-based saturating mutagensis the receptor-binding domain of the SARS-CoV-2 spike protein. Library generation can be applied to any open reading frame of interest by changing the input sequence file (in this case located at notebooks/Spike_RBD.txt)

Navigate to the spikedisplay/notebooks directory and open the `Primer_Library_Generation.ipynb` notebook

```sh
# Activate environment
cd ~/mutagenic-primer-design
source .spikedisplay/bin/activate
cd notebooks
jupyter notebook Primer_Library_Generation.ipynb
```

Run the notebook cell by cell and adjust sequence inputs and other parameters where desired