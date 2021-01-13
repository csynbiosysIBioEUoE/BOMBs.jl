# BOMBS.jl

[![Build Status](https://travis-ci.com/DavidGomezC/BOMBS.jl.svg?branch=master)](https://travis-ci.com/DavidGomezC/BOMBS.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/DavidGomezC/BOMBS.jl?svg=true)](https://ci.appveyor.com/project/DavidGomezC/BOMBS-jl)
[![Coverage](https://codecov.io/gh/DavidGomezC/BOMBS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/DavidGomezC/BOMBS.jl)
[![Coverage](https://coveralls.io/repos/github/DavidGomezC/BOMBS.jl/badge.svg?branch=master)](https://coveralls.io/github/DavidGomezC/BOMBS.jl?branch=master)

                            ___________   __________   ___  ___   ___________   __________
                           /  _____   /  / ______  /  /   \/   \  \   _____  \  \  _______\
                          /  /____/  /  / /     / /  /          \  \  \____\  \  \ \________
                         /  _____   /  / /     / /  /            \  \   _____  \  \________ \
                        /  /    /  /  / /     / /  /   /\____/\   \  \  \    \  \          \ \
                       /  /____/  /  / /_____/ /  /   /        \   \  \  \____\  \   _______\ \
                      /__________/  /_________/  /__ /          \___\  \__________\  \_________\


BOMBS is a simulation and optimisation package for Julia (http://julialang.org/) to make ODE model simulations and design of experiments simpler in the field of Biology (but potentially applied to others). The package is designed to allow the user simulate an ODE model with external inputs, generate pseudo-data, estimate model parameters (Maximum Likelihood Estimation and Bayesian Inference) and design optimal experiments for model selection and calibration in a simple and general way. Advanced programming skills are not required, just make sure to know how to use dictionaries in Julia (https://docs.julialang.org/en/v1/base/collections/#Dictionaries). 

## Installation 
```julia
using Pkg
Pkg.add("BOMBS")
```
or latest master directly from github: 
```julia
Pkg.clone("https://github.com/DavidGomezC/BOMBS.jl")
```
## Package Sections
  ### 1.- Model Generation
  Read more [here](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/1_GenerateModel.ipynb)
    [TEXT here SHOW](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/1_GenerateModel.ipynb)
  ### 2.- Model Simulation
  
  ### 3.- Pseudo-Data Generation
  
  ### 4.- Maximum Likelihood Estimation (MLE)
  
  ### 5.- Bayesian Inference of Parameters (Stan)
  
  ### 6.- Optimal Experimental Design for Model Selection
  
  ### 7.- Optimal Experimental Design for Model Calibration

## Documentation 
All source code for the package is included in the directory src and tests sets in the test directory. \
For deppendencies and version specification, please look at the file Project.toml. \
A set of example jupyter ntoebooks and files can be found in the directory Examples. \
Instructions on how to install cmdstan (at least one way to do it, the one that worked for me) can be found in the file InstallStanInJulia.pdf under the directory InstallCmdstanInfo. \
For more information about what functions are included in each section of the package and some basic information about them have a look at the file BOMBS_Functions_Documentation.pdf under the directory FunctionDocs. 

## References
  1. 
  2. 
  3. 
  4. 
  5. 
  6. 
  7. 
  8. 
  9. 
  
## Usefull links to some of the packages used in BOMBS
  **1. BayesianOptimization.jl:** https://github.com/jbrea/BayesianOptimization.jl \
  **2. BlackBoxOptim.jl:** https://github.com/robertfeldt/BlackBoxOptim.jl \
  **3. Calculus.jl:** https://github.com/JuliaMath/Calculus.jl \
  **4. CmdStan.jl:** https://github.com/StanJulia/CmdStan.jl \
  **5. DifferentialEquations.jl:** https://github.com/SciML/DifferentialEquations.jl \
  **6. GaussianProcesses.jl:** https://github.com/STOR-i/GaussianProcesses.jl \
  **7. GaussianMixtures.jl:** https://github.com/davidavdav/GaussianMixtures.jl \
  **8. ScikitLearnBase.jl:** https://github.com/DavidGomezC/BOMBS.jl/blob/main/src/BOMBS.jl














