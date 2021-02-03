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
  
  The inspiration on how to introduce the model information to the package comes from the Matlab toolbox AMIGO2 [**[1]**](https://sites.google.com/site/amigo2toolbox/), a really nice (frequentist) toolbox worth to check if you have a Matlab licence! 
  
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
  For now (expansion might happen in the future), only global optimisation using genetic algorithms (thanks to the package [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl)) is present, but multiple instances of the optimisation can be run in parallel if the user desires it. However, if you desire to use a different optimiser, a new cost function script gets generated each time and you can use it, but note that you will need to apply modifications due to the use of global variables.
  
  Additionally, you can use the (univariate or multivariate) Gaussian log-likelihood distributions to perform a cross-validation step with additional data. \
  As in for the other sections, you only need to provide a dictionary with the specifications of the model, experiment and data making the process easy and quick. 
  
  For more information about this section and how to use it have a look at [Notebook4](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/4_MaximumLikelihoodEstimation.ipynb) or the brief function documentation of the section from [this document](https://github.com/DavidGomezC/BOMBS.jl/blob/main/FunctionDocs/BOMBS_Functions_Documentation.pdf).
  
  ### 5.- Bayesian Inference of Parameters (Stan)
  The Bayesian inference section allows you to quickly and automatically perform Bayesian inference (multi-experimental inference is considered) of your model parameters using Stan [**[2]**](https://mc-stan.org/), a wonderful tool worth to check. To be able to perform the inference within the package, you will need to have CmdStan installed in your computer. Please, refer to [this document](https://github.com/DavidGomezC/BOMBS.jl/blob/main/InstallCmdstanInfo/InstallStanInJulia.pdf) to see how I installed it (there might be other ways, just make sure that it works and that you apply the mentioned modifications).
  
  This section of the package will automatically generate the data structure and Stan model script to perform the inference in a similar way than shown in references **[[3](https://ieeexplore.ieee.org/document/8791449),[4](https://www.sciencedirect.com/science/article/pii/S2405896319321123),[5](https://pubs.acs.org/doi/abs/10.1021/acssynbio.0c00393)]** (note that that code is in R, here I translated everything to Julia), allowing as much freedom as possible in the definition of your priors. The main intention of the section is to automatically generate the data structure and the main code to perform inference of your model considering event input experiments (external inputs that change across an experiment) without you having to spend a lot of time coding and debugging it (believe me, it takes a long time). Alternatively, you can use this section to generate the aforementioned script and then use it in a different environment for Stan (such as RStan or PyStan). 
  
  Additionally, this section allows you to compute an entropy approximation of your priors samples or definitions and posterior samples as we did in **[[3](https://ieeexplore.ieee.org/document/8791449)]** following the work in **[[6](https://ieeexplore.ieee.org/document/4648062)]** but in a general way so it can be quickly and easily applied to your model. 
  
  For more information about this section and how to use it have a look at [Notebook5](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/5_BayesianInferenceStan.ipynb) or the brief function documentation of the section from [this document](https://github.com/DavidGomezC/BOMBS.jl/blob/main/FunctionDocs/BOMBS_Functions_Documentation.pdf).
  
  ### 6.- Optimal Experimental Design for Model Selection
  The OED for model selection section allows you to design an experiment to optimally discriminate between 2 competing models. For now, the only design available is the optimisation of the concentration for your external inputs where the rest of the experimental details (sampling times, input switching times, observables, etc.) have to be provided. However, a lot of work was been done to allow maximum flexibility on how to optimise your inputs. This is; allow some inputs to be fixed across all the experiment, set a fixed value for some of your inputs to a fixed value to not be optimised, set for your inputs to have the same optimised value in some of the selected steps and any possible combination between these. Average between different observables will be considered. 
  
  The strategy selected for the optimisation is the one described in [**[5]**](https://pubs.acs.org/doi/abs/10.1021/acssynbio.0c00393), where if you only provide one sample for the parameters of your models the frequentist option will take place (Euclidean distance between simulations) and if multiple samples are given for the two models (sorry, but you cannot give only one sample for one model and multiple for the other) the Bayesian option will take place (Bhattacharyya distance between simulations). For the optimisation, everything is set to use the nice package [BayesianOptimization.jl](https://github.com/jbrea/BayesianOptimization.jl) (more options will probably be added in the future, but for now with a little modification of the code you can implement different optimisers since a script for the cost/utility function will be generated in each run). For this, default settings will be used, but do not worry, a nice print when using the package will show you where to go (find a function called settingsBayesOpt in the script [ModelOEDSelection.jl](https://github.com/DavidGomezC/BOMBS.jl/blob/main/src/ModelOEDSelection.jl)) in the source code to modify it, and just in case, a backup file with the default settings is also present. 
  
  For more information about this section and how to use it have a look at [Notebook6](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/6_OEDModelSelection.ipynb) or the brief function documentation of the section from [this document](https://github.com/DavidGomezC/BOMBS.jl/blob/main/FunctionDocs/BOMBS_Functions_Documentation.pdf).
  
  ### 7.- Optimal Experimental Design for Model Calibration
  The OED for model calibration section allows you to design an experiment to optimally infer your model parameters by finding the experiment where your model is most uncertain about (this is, where the posterior predictive distribution for the simulation has maximum uncertainty). For now, the only design available is the optimisation of the concentration for your external inputs where the rest of the experimental details (sampling times, input switching times, observables, etc.) have to be provided. However, a lot of work was been done to allow maximum flexibility on how to optimise your inputs. This is; allow some inputs to be fixed across all the experiment, set a fixed value for some of your inputs to a fixed value to not be optimised, set for your inputs to have the same optimised value in some of the selected steps and any possible combination between these. Average between different observables will be considered. 
  
  However, note that in this case the utility function can be defined in two ways: 1) using distance between percentiles, 2) using Entropy approximation for each time point. The former assesses model uncertainty by computing the Euclidean distance between the 0.5 and 99.5 percentiles for a given observable. The latter loops across each time point, fitting a series of distributions (Beta, Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform) for each and selecting the one with the best likelihood. Then, the closed Entropy formulation for the resultant distribution will be computed and the sum across all time points entropy will be taken. 
  
  The strategy selected for the optimisation is the one described in [**[7]**](), where only a Bayesian approach is available since we are finding the experiment that yields the most uncertain simulation using your posterior parameters, a frequentist and quick alternative cannot be implemented here using the same strategy (FIM or confidence intervals could be used here, but for that you can have a look at the Matlab toolbox [AMIGO2](https://sites.google.com/site/amigo2toolbox/)). For the optimisation, everything is set to use the nice package [BayesianOptimization.jl](https://github.com/jbrea/BayesianOptimization.jl) (more options will probably be added in the future, but for now with a little modification of the code you can implement different optimisers since a script for the cost/utility function will be generated in each run). For this, default settings will be used, but do not worry, a nice print when using the package will show you where to go (find a function called settingsBayesOpt in the script [ModelOEDSelection.jl](https://github.com/DavidGomezC/BOMBS.jl/blob/main/src/ModelOEDSelection.jl)) in the source code to modify it, and just in case, a backup file with the default settings is also present. 
  
  For more information about this section and how to use it have a look at [Notebook7](https://github.com/DavidGomezC/BOMBS.jl/blob/main/Examples/7_OEDModelCalibration.ipynb) or the brief function documentation of the section from [this document](https://github.com/DavidGomezC/BOMBS.jl/blob/main/FunctionDocs/BOMBS_Functions_Documentation.pdf).
  
## Documentation 
All source code for the package is included in the directory src and tests sets in the test directory. \
For dependencies and version specification, please look at the file Project.toml. \
A set of example jupyter notebooks and files can be found in the directory Examples. \
Instructions on how to install cmdstan (at least one way to do it, the one that worked for me) can be found in the file InstallStanInJulia.pdf under the directory InstallCmdstanInfo. \
For more information about what functions are included in each section of the package and some basic information about them have a look at the file BOMBS_Functions_Documentation.pdf under the directory FunctionDocs. 

## References
  **1. Balsa-Canto E, Henriques D, Gábor A, Banga JR. AMIGO2, a toolbox for dynamic modeling, optimization and control in systems biology.** Bioinformatics. 2016 Nov 1;32(21):3357-3359. doi: 10.1093/bioinformatics/btw411. Epub 2016 Jul 4. PMID: 27378288; PMCID: PMC5079478. 
  
  **2. Stan Development Team. 2020. Stan Modeling Language Users Guide and Reference Manual**, 2.25. https://mc-stan.org 
  
  **3. D. G. Cabeza, L. Bandiera, E. Balsa-Canto and F. Menolascina, "Information content analysis reveals desirable aspects of in vivo experiments of a synthetic circuit,"** 2019 IEEE Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB), Siena, Italy, 2019, pp. 1-8, doi: 10.1109/CIBCB.2019.8791449. 
  
  **4. Bandiera, Lucia & Cabeza, D. & Balsa-Canto, Eva & Menolascina, Filippo. (2019). Bayesian model selection in synthetic biology: factor levels and observation functions.** IFAC-PapersOnLine. 52. 24-31. 10.1016/j.ifacol.2019.12.231.  
  
  **5. Bandiera L, Gomez-Cabeza D, Gilman J, Balsa-Canto E, Menolascina F. Optimally Designed Model Selection for Synthetic Biology.** ACS Synth Biol. 2020 Nov 20;9(11):3134-3144. doi: 10.1021/acssynbio.0c00393. Epub 2020 Nov 5. PMID: 33152239. 
  
  **6. M. F. Huber, T. Bailey, H. Durrant-Whyte and U. D. Hanebeck, "On entropy approximation for Gaussian mixture random vectors,"** 2008 IEEE International Conference on Multisensor Fusion and Integration for Intelligent Systems, Seoul, 2008, pp. 181-188, doi: 10.1109/MFI.2008.4648062. 
  
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














