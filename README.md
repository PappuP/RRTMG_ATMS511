# RRTMG Assignment

Hello! This the assignment for learning how to use RRTMG, which is the global climate model version of 
the Rapid Radiative Transfer Model (RRTM). While the original model is written in Fortran, climlab provides a 
Python wrapper for RRTMG, which we'll be using.

# Installing climlab

We recommend creating a new environment or using an existing non-base environment for this. To create a new environment:
`` conda create --name [new env name] ``

To copy an environment you already have to make a new one: 
``  conda create --clone [existing environment] -n [new environment] ``

To activate your environment:
`` conda activate [env name]``

Installing climlab:
``conda install -c conda-forge climlab``

This should be fairly straightforward. All dependencies should be installed automatically if they're not already present in your directory.

# Running the code

You can download all the files or ``git clone`` this repo into your home directory, and it should be ready to run!

We'll be using Jupyter notebooks. If you'd like to set up one in Keeling but haven't before, we suggest using (Max Grover's tutorial)[https://github.com/mgrover1/keeling-crash-course]

# Resources
- (Climlab)[https://climlab.readthedocs.io/en/latest/]
- (Conda cheat sheet)[https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf]
- 
