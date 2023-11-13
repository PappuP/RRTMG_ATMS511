# RRTMG In-class Lab 

Hello! This the assignment for learning how to use RRTMG, which is the global climate model version of 
the Rapid Radiative Transfer Model (RRTM). While the original model is written in Fortran, `climlab` provides a 
Python wrapper for RRTMG, which we'll be using in this in-class lab session. 

# Prior to In-class Session
## A. Installing climlab

We recommend creating a new environment or using an existing non-base environment for this. 

1. To create a new environment called `rrtmg`:
```
conda create --name rrtmg
```

(Alternative) To copy an environment you already have to make a new one: 
```
conda create --clone [existing environment] -n [new environment] 
```

2. To activate your environment:
```
conda activate rrtmg
```

3. Installing climlab in the new environment:
```
conda install pip 
pip install ipykernel # create a new kernel to run the notebook in 
conda install -c conda-forge climlab matplotlib numba # Download climlab
conda install climlab-rrtmg --channel conda-forge # Download the RRTMG specific package of climlab
```

To do this all in one line:
```
conda create --name rrtmg python=3.10 climlab climlab-rrtmg numpy=1.23 numba matplotlib ipykernel --channel conda-forge
```

This should be fairly straightforward. All dependencies should be installed automatically if they're not already present in your directory.

## B. Running the code
You can download all the files or ``git clone`` this repo into your home directory by the following line of code, and it should be ready to run!
```
git clone https://github.com/PappuP/RRTMG_ATMS511.git
```

. If you'd like to set up one in Keeling but haven't before, we suggest using [Max Grover's tutorial](https://github.com/mgrover1/keeling-crash-course).

## C. Debugging: Change the kernel in notebook
If notebook crashes while running,  the steps below are necessary.

1. Install ``nb_conda_kernels`` in your ``base`` environment (i.e. after you typed `bash` + `enter` on Keeling). 

```
conda install -c conda-forge nb_conda_kernels
```

Once this is installed, kernels in all environment that has ipykernel installed will be automatically shown in any notebook running from the base environment (i.e. jupyter notebooks launched in the Keeling login node). 

2. Open Jupyter Notebooks from your base environment and follow the steps below
In the top pannel of the notebook find ``Kernel>>change kernel>>select your new env kernel``

<img src="Supplementary/kernal.png" alt="Change the kernel in notebook" title="Change the kernel in notebook">

# Resources
- [Climlab](https://climlab.readthedocs.io/en/latest/)
- [Conda cheat sheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf)
