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
  For more information about this section and how to use it have a look at the [Notebook1](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/1_GenerateModel.ipynb)
    
  ### 2.- Model Simulation
  For more information about this section and how to use it have a look at the [Notebook2](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/2_SimulateModel.ipynb)
  
  ### 3.- Pseudo-Data Generation
  For more information about this section and how to use it have a look at the [Notebook3](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/3_GeneratePseudoData.ipynb)
  
  ### 4.- Maximum Likelihood Estimation (MLE)
  For more information about this section and how to use it have a look at the [Notebook4](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/4_MaximumLikelihoodEstimation.ipynb)
  
  ### 5.- Bayesian Inference of Parameters (Stan)
  For more information about this section and how to use it have a look at the [Notebook5](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/5_BayesianInferenceStan.ipynb)
  
  ### 6.- Optimal Experimental Design for Model Selection
  For more information about this section and how to use it have a look at the [Notebook6](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/6_OEDModelSelection.ipynb)
  
  ### 7.- Optimal Experimental Design for Model Calibration
  For more information about this section and how to use it have a look at the [Notebook7](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/7_OEDModelCalibration.ipynb)
  
## Documentation 
All source code for the package is included in the directory src and tests sets in the test directory. \
For deppendencies and version specification, please look at the file Project.toml. \
A set of example jupyter ntoebooks and files can be found in the directory Examples. \
Instructions on how to install cmdstan (at least one way to do it, the one that worked for me) can be found in the file InstallStanInJulia.pdf under the directory InstallCmdstanInfo. \
For more information about what functions are included in each section of the package and some basic information about them have a look at the file BOMBS_Functions_Documentation.pdf under the directory FunctionDocs. 

## References
  **1. D. G. Cabeza, L. Bandiera, E. Balsa-Canto and F. Menolascina, "Information content analysis reveals desirable aspects of in vivo experiments of a synthetic circuit,"** 2019 IEEE Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB), Siena, Italy, 2019, pp. 1-8, doi: 10.1109/CIBCB.2019.8791449. \
  
  **2. Bandiera, Lucia & Cabeza, D. & Balsa-Canto, Eva & Menolascina, Filippo. (2019). Bayesian model selection in synthetic biology: factor levels and observation functions.** IFAC-PapersOnLine. 52. 24-31. 10.1016/j.ifacol.2019.12.231.  \
  
  **3. Bandiera L, Gomez-Cabeza D, Gilman J, Balsa-Canto E, Menolascina F. Optimally Designed Model Selection for Synthetic Biology.** ACS Synth Biol. 2020 Nov 20;9(11):3134-3144. doi: 10.1021/acssynbio.0c00393. Epub 2020 Nov 5. PMID: 33152239. \
  
  **4. Stan Development Team. 2020. Stan Modeling Language Users Guide and Reference Manual**, 2.25. https://mc-stan.org \
  
  **5. Balsa-Canto E, Henriques D, Gábor A, Banga JR. AMIGO2, a toolbox for dynamic modeling, optimization and control in systems biology.** Bioinformatics. 2016 Nov 1;32(21):3357-3359. doi: 10.1093/bioinformatics/btw411. Epub 2016 Jul 4. PMID: 27378288; PMCID: PMC5079478. \
  
  **6. M. F. Huber, T. Bailey, H. Durrant-Whyte and U. D. Hanebeck, "On entropy approximation for Gaussian mixture random vectors,"** 2008 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems, Seoul, 2008, pp. 181-188, doi: 10.1109/MFI.2008.4648062. \
  
  **7. OED model calibration paper**

## Usefull links to some of the packages used in BOMBS
  **1. BayesianOptimization.jl:** https://github.com/jbrea/BayesianOptimization.jl \
  **2. BlackBoxOptim.jl:** https://github.com/robertfeldt/BlackBoxOptim.jl \
  **3. Calculus.jl:** https://github.com/JuliaMath/Calculus.jl \
  **4. CmdStan.jl:** https://github.com/StanJulia/CmdStan.jl \
  **5. DifferentialEquations.jl:** https://github.com/SciML/DifferentialEquations.jl \
  **6. GaussianProcesses.jl:** https://github.com/STOR-i/GaussianProcesses.jl \
  **7. GaussianMixtures.jl:** https://github.com/davidavdav/GaussianMixtures.jl \
  **8. ScikitLearnBase.jl:** https://github.com/DavidGomezC/BOMBS.jl/blob/main/src/BOMBS.jl














