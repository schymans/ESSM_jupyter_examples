# Jupyter notebooks using ESSM
Author: Stan Schymanski, Luxembourg Institute of Science and Technology

ESSM stands for Environmental Science using Symbolic Math (see http://essm.rtfd.io), which is a Python package facilitating the use of sympy for mathematical derivations and computations in environmental sciences.

The repository contains another repo as a submodule (folder ESSM_plotting, from https://github.com/schymans/ESSM_plotting). After cloning, you therefore need to execute in the base folder:
```
git submodule init
git submodule update
```
At any stage, you can update all submodules in your repo with a single command:
```
git submodule update --remote --merge
```
To see all submodules, look into the file .gitmodules, e.g. by executing:
```
cat .gitmodules
```
## Contents

For a notebook showcasing general features of ESSM, look at [api_documentation_master.ipynb](https://github.com/schymans/ESSM_jupyter_examples/blob/master/api_documentation_master.ipynb).

The notebook [E_PM_eqs_essm.ipynb](https://github.com/schymans/ESSM_jupyter_examples/blob/master/E_PM_eqs_essm.ipynb) provides equations and examples related to modelling the energy balance of a plant leaf and its exchange of water and energy with the atmosphere. The definitions and equations can be imported into other notebooks by executing: `from E_PM_eqs_defs.py import *`

The notebooks [Priestley-Taylor_short.ipynb](https://github.com/schymans/ESSM_jupyter_examples/blob/master/Priestley-Taylor_short.ipynb) and [Priestley-Taylor_equation1.ipynb](https://github.com/schymans/ESSM_jupyter_examples/blob/master/Priestley-Taylor_equation1.ipynb) illustrate the benefit of ESSM's automatic dimensional consistency testing, and how to analyse and repair equations with inconsistent units.

## Running the notebooks locally
The best approach is to create a dedicated conda/mamba environment with required dependencies. The repo contains a file (environment.yml) with all required information. To create a conda/mamba environment based on this file, use a terminal, and in the cloned folder, execute:
```
mamba env create -n essm -f environment.yml
mamba activate fortress
jupyter lab
```
This will create a mamba environment called essm with the dependencies listed in environment.yml, activate the environment and then start jupyter.

