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


BOMBS (**B**ayesian **O**ptimisation of **M**odels for **B**iological **S**ystems) is a simulation and optimisation package for Julia (http://julialang.org/) to make ODE model simulations and design of experiments simpler in the field of Biology (but potentially applied to others). The package is designed to allow the user simulate an ODE model with external inputs, generate pseudo-data, estimate model parameters (Maximum Likelihood Estimation and Bayesian Inference) and design optimal experiments for model selection and calibration in a simple and general way. Advanced programming skills are not required, just make sure to know how to use dictionaries in Julia (https://docs.julialang.org/en/v1/base/collections/#Dictionaries). 

## Installation 
```julia
using Pkg
Pkg.add("BOMBS") # Note that this installation mode is still not available!
```
or latest master directly from github: 
```julia
Pkg.clone("https://github.com/DavidGomezC/BOMBS.jl")
```
## Package Sections
  ### 1.- Model Generation
  The model generation section allows you to generate the necessary Julia functions for simulation of your ODEs with time-varying inputs in an easy manner. You only need to provide some information about the model inside a dictionary structure and that is it!
  
  The way on how to introduce the model information is inspired in the Matlab toolbox AMIGO2 [**[1]**](https://sites.google.com/site/amigo2toolbox/), a really nice (frequentist) toolbox worth to check if you have a Matlab licence. 
  
  For more information about this section and how to use it have a look at [Notebook1](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/1_GenerateModel.ipynb) or the brief function documentation of the section from [this document](https://github.com/DavidGomezC/BOMBS.jl/blob/main/FunctionDocs/BOMBS_Functions_Documentation.pdf).
    
  ### 2.- Model Simulation
  The model simulation section allows you to simulate your ODEs with time-varying inputs by providing the model generated in section one and introducing the specifications of the experiment(s) in a dictionary. Alternatively, you can provide CSV files from which the experimental details will be extracted.  
  
  For more information about this section and how to use it have a look at [Notebook2](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/2_SimulateModel.ipynb) or the brief function documentation of the section from [this document](https://github.com/DavidGomezC/BOMBS.jl/blob/main/FunctionDocs/BOMBS_Functions_Documentation.pdf).
  
  ### 3.- Pseudo-Data Generation
  The pseudo-data generation section works in a similar way than section two, but not also provides you with the simulations for your model but also generates pseudo-data for the observables you select. The only noise option (for now) is additive heteroscedastic noise for the observable where you can choose the percentage of noise used from the simulation. \
  Thus, the general formulation would be: \
  <img src="https://render.githubusercontent.com/render/math?math=\hat{y}_t = y_t + \epsilon_t"> \
  where \
  <img src="https://render.githubusercontent.com/render/math?math=\epsilon_t = \mathcal{N}(0,y_t*pr)"> \
  t indicates a time-point, y the observable and pr the percentage of y to be used for the noise. 
  
  For more information about this section and how to use it have a look at [Notebook3](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/3_GeneratePseudoData.ipynb) or the brief function documentation of the section from [this document](https://github.com/DavidGomezC/BOMBS.jl/blob/main/FunctionDocs/BOMBS_Functions_Documentation.pdf).
  
  ### 4.- Maximum Likelihood Estimation (MLE)
  The MLE section allows you to quickly estimate the parameters of your model using a (univariate or multivariate) Gaussian distribution form for your likelihood (as a distance measure for the optimisation) between simulations and data. If more than one observable or experiment is present, the average of these will be considered. \
  For now (expansion might happen in the future), only global optimisation using genetic algorithms (thanks to the package BlackBoxOptim.jl) is present, but multiple instances of the optimisation can be run in parallel if the user desires it. However, if you desire to use a differnt optimiser, a new cost function script gets generated each time and you can use it, but note that you will need to apply modifications due to the use of global variables.\
  Additionally, you can use the (univariate or multivariate) Gaussian log-likelihood distributions to perform a cross-validation step with additional data. \
  As in for the other sections, you only need to provide a dictionary with the specifications of the model, experiment and data making the process easy and quick. 
  
  For more information about this section and how to use it have a look at [Notebook4](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/4_MaximumLikelihoodEstimation.ipynb) or the brief function documentation of the section from [this document](https://github.com/DavidGomezC/BOMBS.jl/blob/main/FunctionDocs/BOMBS_Functions_Documentation.pdf).
  
  ### 5.- Bayesian Inference of Parameters (Stan)
  For more information about this section and how to use it have a look at [Notebook5](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/5_BayesianInferenceStan.ipynb) or the brief function documentation of the section from [this document](https://github.com/DavidGomezC/BOMBS.jl/blob/main/FunctionDocs/BOMBS_Functions_Documentation.pdf).
  
  ### 6.- Optimal Experimental Design for Model Selection
  For more information about this section and how to use it have a look at [Notebook6](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/6_OEDModelSelection.ipynb) or the brief function documentation of the section from [this document](https://github.com/DavidGomezC/BOMBS.jl/blob/main/FunctionDocs/BOMBS_Functions_Documentation.pdf).
  
  ### 7.- Optimal Experimental Design for Model Calibration
  For more information about this section and how to use it have a look at [Notebook7](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/7_OEDModelCalibration.ipynb) or the brief function documentation of the section from [this document](https://github.com/DavidGomezC/BOMBS.jl/blob/main/FunctionDocs/BOMBS_Functions_Documentation.pdf).
  
## Documentation 
All source code for the package is included in the directory src and tests sets in the test directory. \
For deppendencies and version specification, please look at the file Project.toml. \
A set of example jupyter ntoebooks and files can be found in the directory Examples. \
Instructions on how to install cmdstan (at least one way to do it, the one that worked for me) can be found in the file InstallStanInJulia.pdf under the directory InstallCmdstanInfo. \
For more information about what functions are included in each section of the package and some basic information about them have a look at the file BOMBS_Functions_Documentation.pdf under the directory FunctionDocs. 

## References
  **1. Balsa-Canto E, Henriques D, Gábor A, Banga JR. AMIGO2, a toolbox for dynamic modeling, optimization and control in systems biology.** Bioinformatics. 2016 Nov 1;32(21):3357-3359. doi: 10.1093/bioinformatics/btw411. Epub 2016 Jul 4. PMID: 27378288; PMCID: PMC5079478. \
  
  **2. D. G. Cabeza, L. Bandiera, E. Balsa-Canto and F. Menolascina, "Information content analysis reveals desirable aspects of in vivo experiments of a synthetic circuit,"** 2019 IEEE Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB), Siena, Italy, 2019, pp. 1-8, doi: 10.1109/CIBCB.2019.8791449. \
  
  **3. Bandiera, Lucia & Cabeza, D. & Balsa-Canto, Eva & Menolascina, Filippo. (2019). Bayesian model selection in synthetic biology: factor levels and observation functions.** IFAC-PapersOnLine. 52. 24-31. 10.1016/j.ifacol.2019.12.231.  \
  
  **4. Bandiera L, Gomez-Cabeza D, Gilman J, Balsa-Canto E, Menolascina F. Optimally Designed Model Selection for Synthetic Biology.** ACS Synth Biol. 2020 Nov 20;9(11):3134-3144. doi: 10.1021/acssynbio.0c00393. Epub 2020 Nov 5. PMID: 33152239. \
  
  **5. Stan Development Team. 2020. Stan Modeling Language Users Guide and Reference Manual**, 2.25. https://mc-stan.org \
  
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














