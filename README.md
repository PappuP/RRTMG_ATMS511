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
```


This should be fairly straightforward. All dependencies should be installed automatically if they're not already present in your directory.

# Running the code

You can download all the files or ``git clone`` this repo into your home directory by the following line of code, and it should be ready to run!
```
git clone https://github.com/PappuP/RRTMG_ATMS511.git
```

. If you'd like to set up one in Keeling but haven't before, we suggest using [Max Grover's tutorial](https://github.com/mgrover1/keeling-crash-course).
# Resources
- [Climlab](https://climlab.readthedocs.io/en/latest/)
- [Conda cheat sheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf)
