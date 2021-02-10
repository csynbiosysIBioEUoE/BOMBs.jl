## This function prints some information on what should be the contents the user has to give in the different sections of the package.
# To select from which section you want the information from, just type the keywords:
#        Model, Simulation/s, Pseudo-Data/Pseudodata, MLE, Inference/Stan


function infoAll(woo=[])

    if typeof(woo) != String
        println("Sorry, but the only allowed input to this function is a string :(");

        print(string("

Please, if you want information about one of the package structures, type:
    1) ","\"","model","\""," for the model generation section
    2) ","\"","simulation","\"",", ","\"","simulation","\""," or ","\"","simul","\""," for the model simulation section
    3) ","\"","pseudo-data","\"",", ","\"","pseudodata","\"",", ","\"","pseudo data","\""," for the model pseudo-data generation section
    4) ","\"","mle","\"",", ","\"","likelihood","\""," for the maximum likelihood estimation section
    5) ","\"","inference","\"",", ","\"","stan","\"",", ","\"","stan inference","\"",", ","\"","staninference","\""," for the Bayesian parameter inference section
    6) ","\"","oedms","\"",", ","\"","model selection","\"",", ","\"","modelselection","\"",", ","\"","oed model selection","\""," for the optimal experimental design for model selection section
    7) ","\"","oedmc","\"",", ","\"","model calibration","\"",", ","\"","modelcalibration","\"",", ","\"","oed model calibration","\""," for the optimal experimental design for model calibration section

        "))

        return
    end

    woohoo = lowercase(woo);

    if woohoo == "model"
        print(string(raw"

                CALL defModStruct()

model_def[","\"","NameF","\"","] = [];
        # String containing the name of the model. Scripts and results will be stored using this name

model_def[","\"","nStat","\"","] = [];
        # Integer indicating the total number of steps of the model.

model_def[","\"","nPar","\"","] = [];
        # Integer indicating the total number of parameters of the model.

model_def[","\"","nInp","\"","] = [];
        # Integer indicating the total number of stimuli (inducers) of the model. If the model has no inputs, set it to 0.

model_def[","\"","stName","\"","] = [];
        # Vector of strings indicating the name of all the states of the model (without a d letter in front).

model_def[","\"","parName","\"","] = [];
        # Vector of strings indicating the name of all the parameters of the model.

model_def[","\"","inpName","\"","] = [];
        # Vector of strings indicating the name of all the stimuli (inducers) of the model. If the model has no inputs
        # just give an empty vector.

model_def[","\"","eqns","\"","] = [];
        # Vector of strings containing all the equations for the model (left and right-hand sides).
        # If an equation represents a state, the left-hand side has to be one of the strings contained in
        # stName but with a d in front (example: Prot -> dProt = ...).
        # Equations that are not states of ODEs are also allowed. Same as Julia expressions (println, for, if, etc.)
        # If you want to include a condition (if) for a state variable, each one of the if statement has to be
        # written in a separate string. Be careful with this, since then your Stan model will not work. To make it
        # work you have to write it in Stan language, but then the Julia code will not work. We recommend that if so,
        # first generate the stan code and then apply all the necessary modifications there. However, if you do not care
        # for the Julia code and just want the Stan code go on. Just know that in this case, if there is a condition in
        # one of the states, in stan you need to type the full if statement in one same string.
        # Note that models without external inputs are also supported except for the optimal experimental design sections.

model_def[","\"","Y0eqs","\"","] = [];
        # Vector of strings containing the steady-state equations of the model if desired (if not, just leave
        # it as an empty vector). These equations will be used to compute y0 assuming steady-state reached
        # before the experiment.
        # All the states must appear (identified in the left-hand side, but this time without the d in front),
        # however other equations are allowed.
        # If some  element of the equation requires an experimental value for the calculation,
        # please add exp at the beginning of the state (example: Cmrna -> expCmrna).
        # Please, do not use the name alp for anything, since this is reserved.

model_def[","\"","Y0Sim","\"","] = [];
        # If analytical solution of the model at steady state is not accurate enough and you want to add
        # an Over-Night simulation before the actual experiment simulation.
        # Allowed vales are true, false, ","\"","Yes","\"",", ","\"","yes","\"",", ","\"","No","\"",", ","\"","no","\"",".
        # Default value is false.
        # Time-scale for the simulation is assumed in minutes (1440 min). If this wants to be changed, what should be
        # introduced here is a number for the time conversion (e.g. 1/60 if to convert to days or 60 if to convert to seconds)

model_def[","\"","tols","\"","] = [];
        # Vector of 2 floats containing the relative and absolute tolerances for the solver (in this order).
        # If left empty, 1e-6 will be assumed for both.

model_def[","\"","solver","\"","] = [];
        # IVP solver to solve the ODEs. If nothing specified, the default will be Tsit5().
        # For more info check https://diffeq.sciml.ai/v2.0/tutorials/ode_example.html

"));


    elseif woohoo == "simulation" || woohoo == "simulations" || woohoo == "simul"
        print(string(raw"

                CALL defSimulStruct()
        MAIN STRUCTURE


simul_def[","\"","Nexp","\"","] = [];
        # Integer indicating the number of experiments to be simulated

simul_def[","\"","finalTime","\"","] = [];
        # Vector of final times for each simulation (initial time will always be assumed as 0,
        # so please consider that).

simul_def[","\"","switchT","\"","] = [];
        # Array with the switching times of the inducer in the simulation (time 0 and final
        # time need to be considered)

simul_def[","\"","y0","\"","] = [];
        # Array (single simulation) or matrix (multiple simulations) of Y0s for the simulations for each
        # experiment. If you are computing the steady-state this vector might not be used, however, you
        # still need to introduce it with some random numbers.

simul_def[","\"","preInd","\"","] = [];
        # Vector of numbers with the values for the stimuli (inducer) in the over-night. It might be the case
        # that this entry is not required since only the y0 vector is considered for the initial point of
        # the simulation. However, you still need to introduce a random value for it to avoid future issues.

simul_def[","\"","uInd","\"","] = [];
        # Array containing the values for the stimuli at each step for each experiment. In each sub-array, columns
        # indicate a step and rows indicate an inducer (if multiple ones are considered).

simul_def[","\"","theta","\"","] = [];
        # Vector/Matrix with the parameter samples or directory and file location of CSV file with them.

simul_def[","\"","tsamps","\"","] = [];
        # Array of sampling time vectors for the experiments.

simul_def[","\"","plot","\"","] = [];
        # Boolean or yes/no string to save the resulting simulations in the results directory (false will be
        # considered as default).

simul_def[","\"","flag","\"","] = [];
        # String to attach a unique flag to the generated scripts and result files so they are not overwritten.
        # If empty, nothing will be added.

"))


print(string(raw"

            CALL defSimulStructFiles()

        STRUCTURE IF CSV FILES ARE USED

simul_def[","\"","ObservablesFile","\"","] = [];
        # Vector of strings containing the name of the files that have the information about the sampling
        # times and y0 values.

simul_def[","\"","EventInputsFile","\"","] = [];
        # Vector of strings containing the name of all the files that have the information about the stimuli
        # and events of the experiment.

simul_def[","\"","theta","\"","] = [];
        # Vector/Matrix with the parameter samples or directory and file location of CSV file with them.

simul_def[","\"","MainDir","\"","] = [];
        # Main directory path for the files in ObervablesFile and EventInputsFile.

simul_def[","\"","plot","\"","] = [];
        # Boolean or yes/no string to save the resulting simulations in the results directory (false will be
        # considered as default).

simul_def[","\"","flag","\"","] = [];
        # String to attach a unique flag to the generated scripts and result files so they are not overwritten.
        # If empty, nothing will be added.

"))


    elseif woohoo == "pseudo-data" || woohoo == "pseudodata" || woohoo == "pseudo data"
        print(string(raw"

            CALL defPseudoDatStruct()

        MAIN STRUCTURE

pseudo_def[","\"","Nexp","\"","] = [];
        # Integer indicating the number of experiments to be simulated

pseudo_def[","\"","finalTime","\"","] = [];
        # Vector of final times for each simulation (initial time will always be assumed as 0,
        # so please consider that).

pseudo_def[","\"","switchT","\"","] = [];
        # Array with the switching times of the inducer in the simulation (time 0 and final
        # time need to be considered)

pseudo_def[","\"","y0","\"","] = [];
        # Array (single simulation) or matrix (multiple simulations) of Y0s for the simulations for each
        # experiment. If you are computing the steady-state this vector might not be used, however, you
        # still need to introduce it with some random numbers.

pseudo_def[","\"","preInd","\"","] = [];
        # Vector of numbers with the values for the stimuli (inducer) in the over-night. It might be the case
        # that this entry is not required since only the y0 vector is considered for the initial point of
        # the simulation. However, you still need to introduce a random value for it to avoid future issues.

pseudo_def[","\"","uInd","\"","] = [];
        # Array containing the values for the stimuli at each step for each experiment. In each sub-array, columns
        # indicate a step and rows indicate an inducer (if multiple ones are considered).

pseudo_def[","\"","theta","\"","] = [];
        # Vector/Matrix with the parameter samples or directory and file location of CSV file with them.

pseudo_def[","\"","tsamps","\"","] = [];
        # Array of sampling time vectors for the experiments.

pseudo_def[","\"","plot","\"","] = [];
        # Boolean or yes/no string to save the resulting simulations in the results directory (false will be
        # considered as default).

pseudo_def[","\"","flag","\"","] = [];
        # String to attach a unique flag to the generated scripts and result files so they are not overwritten.
        # If empty, nothing will be added.

pseudo_def[","\"","Obs","\"","] = [];
        # States of the model that are observables. This is either a vector of strings, a vector of integers
        # indicating which entries from model_def[","\"","stName","\"","] are observables. If a vector of
        # strings is given, these could also be an expression combining states (Only +,-,*,/ and ^ will be
        # considered).

pseudo_def[","\"","Noise","\"","] = [];
        # Percentage of heteroscedastic noise (introduced as value from 0 to 1). If empty 10% will be
        # assumed. This has to be a vector of noise values for each observable.

"))

        print(string(raw"

            CALL defPseudoDatStructFiles()

        STRUCTURE IF CSV FILES ARE USED

pseudo_def[","\"","ObservablesFile","\"","] = [];
        # Vector of strings containing the name of the files that have the information about the sampling
        # times and y0 values.

pseudo_def[","\"","EventInputsFile","\"","] = [];
        # Vector of strings containing the name of all the files that have the information about the stimuli
        # and events of the experiment.

pseudo_def[","\"","theta","\"","] = [];
        # Vector/Matrix with the parameter samples or directory and file location of CSV file with them.

pseudo_def[","\"","MainDir","\"","] = [];
        # Main directory path for the files in ObervablesFile and EventInputsFile.

pseudo_def[","\"","plot","\"","] = [];
        # Boolean or yes/no string to save the resulting simulations in the results directory (false will be
        # considered as default).

pseudo_def[","\"","flag","\"","] = [];
        # String to attach a unique flag to the generated scripts and result files so they are not overwritten.
        # If empty, nothing will be added.

pseudo_def[","\"","Obs","\"","] = [];
        # States of the model that are observables. This is either a vector of strings, a vector of integers
        # indicating which entries from model_def[","\"","stName","\"","] are observables. If a vector of
        # strings is given, these could also be an expression combining states (Only +,-,*,/ and ^ will be
        # considered).

pseudo_def[","\"","Noise","\"","] = [];
        # Percentage of heteroscedastic noise (introduced as value from 0 to 1). If empty 10% will be
        # assumed. This has to be a vector of noise values for each observable.

"))

    elseif woohoo == "mle" || woohoo == "likelihood"
        print(string(raw"

                CALL defMLEStruct()

mle_def[","\"","Nexp","\"","] = [];
        # Integer indicating the number of experiments to be simulated

mle_def[","\"","finalTime","\"","] = [];
        # Vector of final times for each simulation (initial time will always be assumed as 0,
        # so please consider that).

mle_def[","\"","switchT","\"","] = [];
        # Array with the switching times of the inducer in the simulation (time 0 and final
        # time need to be considered)

mle_def[","\"","y0","\"","] = [];
        # Array (single simulation) of Y0s for the simulations for the
        # experiment. If you are computing the steady-state this vector might not be used, however, you
        # still need to introduce it with some random numbers.

mle_def[","\"","preInd","\"","] = [];
        # Vector of numbers with the values for the stimuli (inducer) in the over-night. It might be the case
        # that this entry is not required since only the y0 vector is considered for the initial point of
        # the simulation. However, you still need to introduce a random value for it to avoid future issues.

mle_def[","\"","uInd","\"","] = [];
        # Array containing the values for the stimuli at each step for each experiment. In each sub-array, columns
        # indicate a step and rows indicate an inducer (if multiple ones are considered).

mle_def[","\"","tsamps","\"","] = [];
        # Array of sampling time vectors for the experiments.

mle_def[","\"","plot","\"","] = [];
        # Boolean or yes/no string to save the resulting simulations in the results directory (false will be
        # considered as default).

mle_def[","\"","flag","\"","] = [];
        # String to attach a unique flag to the generated scripts and result files so they are not overwritten.
        # If empty, nothing will be added.

mle_def[","\"","thetaMAX","\"","] = [];
        # Vector containing the maximum bounds for theta (no files can be introduced)

mle_def[","\"","thetaMIN","\"","] = [];
        # Vector containing the minimum bounds for theta (no files can be introduced)

mle_def[","\"","runs","\"","] = [];
        # Integer indicating how many runs of Optimisation will be done. You will get as many theta vectors
        # as runs selected.

mle_def[","\"","parallel","\"","] = [];
        # Boolean or yes/no string indicating if the different runs want to be done in parallel (true) or
        # series (false). Default is false.


    # For the two following fields, you can introduce a string pointing to the observable files (same strings in)
    #     both fields having the same structure as the ones generated in the PseudoData section. If multiple theta are
    #     considered in the file, then the covariance matrix will be taken.
mle_def[","\"","DataMean","\"","] = [];
        # Array containing the vector of means for each experiment.

mle_def[","\"","DataError","\"","] = [];
        # Array containing the vector or matrices (covariance included) of errors for the data for each experiment.

                                        IMPORTANT!!!!

        # Whilst each entry (experiment) of DataMean can be a matrix where each column is an observable of the
        # system, for DataError this is not the case. Each entry of the array (experiment) will have as many
        # entries as observables, where it would be a vector of errors or a matrix. This is done this way to
        # generalise the presence of both options.

mle_def[","\"","Obs","\"","] = [];
        # States of the model that are observables. This is either a vector of strings, a vector of integers
        # indicating which entries from model_def[","\"","stName","\"","] are observables. If a vector of
        # strings is given, these could also be an expression combining states (Only +,-,*,/ and ^ will be
        # considered).

mle_def[","\"","OPTsolver","\"","] = [];
        # For now we only use the package BlackBoxOptim, so any of their options can be used. The default
        # is adaptive_de_rand_1_bin_radiuslimited. Please, introduce it as a string.

mle_def[","\"","MaxTime","\"","] = [];
        # Integer indicating the maximum number of time allowed for the optimisation as a stop criterion.
        # If this is selected MaxFuncEvals has to be empty (only 1 stop criterion). If both are empty,
        # the default will be to do 1000 function evaluations

mle_def[","\"","MaxFuncEvals","\"","] = [];
        # Integer indicating the maximum number of function evaluations allowed for the optimisation as a
        # stop criterion. If this is selected MaxTime has to be empty (only 1 stop criterion). If both are
        # empty, the default will be to do 1000 function evaluations

"))


    elseif woohoo == "inference" || woohoo == "stan" || woohoo == "stan inference" || woohoo == "staninference"
        print(string(raw"

            CALL defBayInfStruct()

        MAIN STRUCTURE

bayinf_def[","\"","Priors","\"","] = []; # Five options for this:
        # 1) A 2*N array containing the bounds for the parameters. Order of parameters will be assumed as
            # the one introduced in model_def[","\"","parName","\"","]. As a prior a truncated Normal covering 2
            # standard deviations in the bounds given will be generated as prior.
        # 2) Path to a CSV file containing samples for the parameters. Fitting of the samples to
            # different type of distributions will be done to generate the priors. Order of parameters
            # will be assumed as the one introduced in model_def[","\"","parName","\"","]. You can also
            # introduce a 2D array of floats with the samples (Array{Float64,2}).
        # 3) Dictionary with fields pars, transpars and pridis defining the
            # parameters, subsequent desired transformations and prior distributions.
            # In the field pars, parameters
            # have to be defined with the same names and order as in model_def[","\"","parName","\"","],
            # otherwise the script will not proceed.
        # 4) An empty array. If this is the case, the stan model will be generated with nothing in
            # the parameters and transformed parameters section. The path to the stan file will
            # be given so the user can fill these sections.
        # 5) Path to a Stan model file if you already have one. For example, if first you introduce an
            # empty array, then you can fill the parameters section and run from that script instead
            # of having to copy the things in here.

bayinf_def[","\"","Data","\"","] = []; # Two options:
        # 1) A dictionary containing 3 fields. Observables, for the path to the files containing the
            # experimental data. Inputs, for the path to the files containing the input profiles for the
            # experiments. Obs, a string vector containing the observables of the experiments.
            # The format of the files has to be the same one as the ones generated in the
            # pseudo-data section. For each of the 2 entries, more than one file can be given if a
            # multi-experimental inference wants to be done. Obs file also needs to be given.
            # A y0 in a field with the same name has to be given for each experiment. The field is
            # compulsory, even if not used.
            # WARNING: For now, this option does not include extraction of covariance matrix for data.
        # 2) A dictionary containing the same structure as the simul_def plus the fields DataMean and
            # DataError containing the experimental data. You can call the function defBayInfDataStruct()
            # to obtain the empty structure of the dictionary.

bayinf_def[","\"","StanSettings","\"","] = []; # Two options:
        # 1) A dictionary with the basic fields from a stan run. The structure of the dictionary can be
            # extracted calling the function defBasicStanSettingsStruct().
        # 2) An empty array. If this is the case, no inference will be done after calling the main
            # function. Instead, you will be given a StanModel file and data structure and you will
            # have to use this to run a call of Stan by yourself. An example will be provided.

bayinf_def[","\"","flag","\"","] = [];
            # String to attach a unique flag to the generated files so it is not overwritten.
            # If empty, nothing will be added.

bayinf_def[","\"","plot","\"","] = [];
            # true, flase, ","\"","Yes","\"",", ","\"","yes","\"",", ","\"","No","\"",", ","\"","no","\""," or []
            # indicating if plots with the results will be generated. Default is false.

bayinf_def[","\"","runInf","\"","] = [];
            # true, flase, ","\"","Yes","\"",", ","\"","yes","\"",", ","\"","No","\"",", ","\"","no","\""," or []
            # indicating if inference wants to be performed or only stan structure is given.
            # Default is true (if all the necessary elements are present)

bayinf_def[","\"","MultiNormFit","\"","] = [];
            # true, false, ","\"","yes","\"",", ","\"","no","\"",", ","\"","Yes","\"",", ","\"","No","\""," or []
            # indicating that if samples are given to fit, the prior will be fitted as a MultiNormal distribution
            # or not. the default is false. If some parameter is better fit with a Log-Normal distribution, this
            # will be used with a distribution reparameterisation to begin to be included in the multinormal.
            # If for some parameter a Uniform distribution is better, this parameter will be excluded from
            # the multinormal and defined as a separate parameter. This can be used as an example of how to set
            # your Stan Model to use multi_normal priors.

"))

print(string(raw"

            CALL defBayInfDataStruct()

        DATA DICTIONARY STRUCTURE

data_def[","\"","Nexp","\"","] = [];
        # Integer indicating the number of experiments to be simulated

data_def[","\"","finalTime","\"","] = [];
        # Vector of final times for each simulation (initial time will always be assumed as 0,
        # so please consider that).

data_def[","\"","switchT","\"","] = [];
        # Array with the switching times of the inducer in the simulation (time 0 and final
        # time need to be considered)

data_def[","\"","y0","\"","] = [];
        # Array (single simulation) or matrix (multiple simulations) of Y0s for the simulations for each
        # experiment. If you are computing the steady-state this vector might not be used, however, you
        # still need to introduce it with some random numbers.

data_def[","\"","preInd","\"","] = [];
        # Vector of numbers with the values for the stimuli (inducer) in the over-night. It might be the case
        # that this entry is not required since only the y0 vector is considered for the initial point of
        # the simulation. However, you still need to introduce a random value for it to avoid future issues.

data_def[","\"","uInd","\"","] = [];
        # Array containing the values for the stimuli at each step for each experiment. In each sub-array, columns
        # indicate a step and rows indicate an inducer (if multiple ones are considered).

data_def[","\"","tsamps","\"","] = [];
        # Array of sampling time vectors for the experiments.

data_def[","\"","DataMean","\"","] = [];
        # Array containing the vector of means for each experiment.

data_def[","\"","DataError","\"","] = [];
        # Array containing the vector or matrices (covariance included) of errors for the data for each experiment.

                                        IMPORTANT!!!!

        # Whilst each entry (experiment) of DataMean can be a matrix where each column is an observable of the
        # system, for DataError this is not the case. Each entry of the array (experiment) will have as many
        # entries as observables, where it would be a vector of errors or a matrix. This is done this way to
        # generalise the presence of both options.

data_def[","\"","Obs","\"","] = [];
        # States of the model that are observables. This is either a vector of strings, a vector of integers
        # indicating which entries from model_def[","\"","stName","\"","] are observables. If a vector of
        # strings is given, these could also be an expression combining states (Only +,-,*,/ and ^ will be
        # considered).

"))


print(string(raw"

            CALL defBayInfDataFromFilesStruct()

        DATA DICTIONARY STRUCTURE IF FILES ARE GIVEN

data_def[","\"","Observables","\"","] = [];
        # Vector of strings containing the name of the files that have the information about the sampling
        # times and y0 values.

data_def[","\"","Inputs","\"","] = [];
        # Vector of strings containing the name of all the files that have the information about the stimuli
        # and events of the experiment.

data_def[","\"","Obs","\"","] = [];
        # States of the model that are observables. This is either a vector of strings, a vector of integers
        # indicating which entries from model_def[","\"","stName","\"","] are observables. If a vector of
        # strings is given, these could also be an expression combining states (Only +,-,*,/ and ^ will be
        # considered).

data_def[","\"","y0","\"","] = [];
        # Array (single simulation) or matrix (multiple simulations) of Y0s for the simulations for each
        # experiment. If you are computing the steady-state this vector might not be used, however, you
        # still need to introduce it with some random numbers.

"))

print(string(raw"

            CALL defBasicStanSettingsStruct()

        STAN SETTINGS DICTIONARY STRUCTURE

stan_def[","\"","cmdstan_home","\"","] = [];
        # String containing the full path to the installation directory for cmdstan

stan_def[","\"","nchains","\"","] = [];
        # Number of chains for the inference (integer)

stan_def[","\"","nsamples","\"","] = [];
        # Number of post-warmup samples for each chain (integer)

stan_def[","\"","nwarmup","\"","] = [];
        # Number of warm-up samples for each chain (integer)

stan_def[","\"","printsummary","\"","] = [];
        # Print sumary of inference at the end. This can be either true or false. Default will be true.

stan_def[","\"","init","\"","] = [];
        # Initial point for the parameters in each chain. See output of MLE results. This field can be empty.
        # Please introduce the parameter values in their true parameter range.

stan_def[","\"","maxdepth","\"","] = [];
        #  Maximum tree-depth. This field can be empty. Check Stan documentation for more information.

stan_def[","\"","adaptdelta","\"","] = [];
        # Delta value between 0 and 1. This field can be empty.  Check Stan documentation for more information.

stan_def[","\"","jitter","\"","] = [];
        # Jitter value between 0 and 1. This field can be empty.  Check Stan documentation for more information.

"));



    elseif woohoo == "oedms" || woohoo == "model selection" || woohoo == "modelselection" || replace("woohoo", " "=>"") == "oedmodelselection"
print(string(raw"

                CALL defODEModelSelectStruct()

oedms_def[","\"","Model_1","\"","] = [];
        # Model structure for Model 1. Dict. See Model Generation Section.
        # Note that the order of the inputs is the one defined in this model. If there is a chance that
        # Model_2 has an input that does not exist in this model, these will be appended at the end
        # of the ones in Model_1.

oedms_def[","\"","Model_2","\"","] = [];
        # Model structure for Model 2. Dict. See Model Generation Section.

oedms_def[","\"","Obs","\"","] = [];
        # States of the model that are observables. This is a vector of strings. These could also be an
        # expression combining states (Only +,-,*,/ and ^ will be considered).
        # Note that due to the use of the covariance matrix in some computations an addition of a small
        # number (0.1) is done in the diagonal elements (for now) to ensure this to be positive definite.
        # If the observable(s) of your model are normalised or in a range between 0 and 1 we recommend to
        # re-scale this to a larger range (0 to 100 for example) so the variances considered are as close
        # as possible to the real ones. This can be done by just multiplying the observable by a constant
        # (Obs*100 for example) in each entry of the vector.

oedms_def[","\"","Theta_M1","\"","] = [];
        # Theta vector (frequentist OED or model with 1 parameter) or matrix (Bayesian OED) for model 1.
        # Path to file with samples is also allowed.

oedms_def[","\"","Theta_M2","\"","] = [];
        # Theta vector (frequentist OED or model with 2 parameter) or matrix (Bayesian OED) for model 2.
        # Path to file with samples is also allowed.

oedms_def[","\"","y0_M1","\"","] = [];
        # For Model 1:
        # Array (single simulation) or matrix (multiple simulations) of Y0s for the simulations for each
        # experiment. If you are computing the steady-state this vector might not be used, however, you
        # still need to introduce it with some random numbers.

oedms_def[","\"","y0_M2","\"","] = [];
        # For Model 2:
        # Array (single simulation) or matrix (multiple simulations) of Y0s for the simulations for each
        # experiment. If you are computing the steady-state this vector might not be used, however, you
        # still need to introduce it with some random numbers.

oedms_def[","\"","preInd_M1","\"","] = [];
        # For Model 1:
        # Vector of numbers with the values for the stimuli (inducer) in the over-night. It might be the case
        # that this entry is not required since only the y0 vector is considered for the initial point of
        # the simulation. However, you still need to introduce a random value for it to avoid future issues.

oedms_def[","\"","preInd_M2","\"","] = [];
        # For Model 2:
        # Vector of numbers with the values for the stimuli (inducer) in the over-night. It might be the case
        # that this entry is not required since only the y0 vector is considered for the initial point of
        # the simulation. However, you still need to introduce a random value for it to avoid future issues.

oedms_def[","\"","finalTime","\"","] = [];
        # Vector of final times for each simulation (initial time will always be assumed as 0,
        # so please consider that).

oedms_def[","\"","switchT","\"","] = [];
        # Array with the switching times of the inducer in the simulation (time 0 and final
        # time need to be considered)

oedms_def[","\"","tsamps","\"","] = [];
        # Array of sampling time vectors for the experiments.

oedms_def[","\"","fixedInp","\"","] = [];
        # If more than 1 inducer for the system exists, but only 1 input can be dynamic, give a vector
        # of strings indicating which inputs are going to be optimised but as a constant input instead
        # of dynamic. If none, just give an empty vector.
        # If there is only 1 stimuli in the model, this field will be ignored.

oedms_def[","\"","fixedStep","\"","] = [];
        # If you want any of the steps to be fixed to a value. This has to be an empty array if none is fixed
        # or an array of tuples, where each tuple is a step to be fixed. The first entry of the tuple is the
        # index of the step (as an Integer), and the second and array of values for each inducer. Note that
        # the fixed inputs will be ignored, so do not take them into account here.
        # The example type should be Array{Array{Int,1},1} or an empty array ([]) if not used.

oedms_def[","\"","equalStep","\"","] = [];
        # If you want a series of steps to have the same optimised value (for example if you want to design a
        # pulse experiment) you can introduce inside this array different arrays with the indexes of the steps
        # that will have the same value. The values introduced in each array need to be integers.

oedms_def[","\"","plot","\"","] = [];
        # true, flase, ","\"","Yes","\"",", ","\"","yes","\"",", ","\"","No","\"",", ","\"","no","\""," or []
        # indicating if plots with the results will be generated. Default is false.

oedms_def[","\"","flag","\"","] = [];
        # String to attach a unique flag to the generated files so it is not overwritten.
        # If empty, nothing will be added.

            # The order of uUpper and uLower will be taken as the order of inducers defined in Model 1.
            # If Model 2 has an aditional input, this will be added after, so consider this when introducing
            # the bounds for them.
oedms_def[","\"","uUpper","\"","] = [];
        # Vector indicating the upper bounds for the inducers

oedms_def[","\"","uLower","\"","] = [];
        # Vector indicating the lower bounds for the inducers

oedms_def[","\"","maxiter","\"","] = [];
        # Maximum number of iterations for the Bayesian Optimisation. If nothing is introduced a default
        # of 100 iterations will be taken

"))




    elseif woohoo == "oedmc" || woohoo == "model calibration" || woohoo == "modelcalibration" || replace("woohoo", " "=>"") == "oedmodelcalibration"
print(string(raw"

                CALL defODEModelCalibrStruct()

oedmc_def[","\"","Model","\"","] = [];
        # Dict with Model. See Model Generation Section.

oedmc_def[","\"","Obs","\"","] = [];
        # States of the model that are observables. This is a vector of strings.These could also be an
        # expression combining states (Only +,-,*,/ and ^ will be considered).

oedmc_def[","\"","Theta","\"","] = [];
        # Theta matrix (Bayesian OED) for the model. No single vectors will be allowed.
        # Path to file is also allowed.

oedmc_def[","\"","y0","\"","] = [];
        # Array (single simulation) or matrix (multiple simulations) of Y0s for the simulations for each
        # experiment. If you are computing the steady-state this vector might not be used, however, you
        # still need to introduce it with some random numbers.

oedmc_def[","\"","preInd","\"","] = [];
        # Vector of numbers with the values for the stimuli (inducer) in the over-night. It might be the case
        # that this entry is not required since only the y0 vector is considered for the initial point of
        # the simulation. However, you still need to introduce a random value for it to avoid future issues.

oedmc_def[","\"","finalTime","\"","] = [];
        # Vector of final times for each simulation (initial time will always be assumed as 0,
        # so please consider that).

oedmc_def[","\"","switchT","\"","] = [];
        # Array with the switching times of the inducer in the simulation (time 0 and final
        # time need to be considered)

oedmc_def[","\"","tsamps","\"","] = [];
        # Array of sampling time vectors for the experiments.

oedmc_def[","\"","fixedInp","\"","] = [];
        # If more than 1 inducer for the system exists, but only 1 input can be dynamic, give a vector
        # of strings indicating which inputs are going to be optimised but as a constant input instead
        # of dynamic. If none, just give an empty vector.
        # If there is only 1 stimuli in the model, this field will be ignored.

oedmc_def[","\"","fixedStep","\"","] = [];
        # If you want any of the steps to be fixed to a value. This has to be an empty array if none is fixed
        # or an array of tuples, where each tuple is a step to be fixed. The first entry of the tuple is the
        # index of the step (as an Integer), and the second and array of values for each inducer. Note that
        # the fixed inputs will be ignored, so do not take them into account here.

oedmc_def[","\"","equalStep","\"","] = [];
        # If you want a series of steps to have the same optimised value (for example if you want to design a
        # pulse experiment) you can introduce inside this array different arrays with the indexes of the steps
        # that will have the same value. The values introduced in each array need to be integers.

oedmc_def[","\"","plot","\"","] = [];
        # true, flase, ","\"","Yes","\"",", ","\"","yes","\"",", ","\"","No","\"",", ","\"","no","\""," or []
        # indicating if plots with the results will be generated. Default is false.

oedmc_def[","\"","flag","\"","] = [];
        # String to attach a unique flag to the generated files so it is not overwritten.
        # If empty, nothing will be added.

oedmc_def[","\"","uUpper","\"","] = [];
        # Vector indicating the upper bounds for the inducers

oedmc_def[","\"","uLower","\"","] = [];
        # Vector indicating the lower bounds for the inducers

oedmc_def[","\"","maxiter","\"","] = [];
        # Maximum number of iterations for the Bayesian Optimisation. If nothing is introduced a default
        # of 100 iterations will be taken

oedmc_def[","\"","util","\"","] = [];
        # String indicating entropy or perc (or percentile) as the core of the utility function to
        # compute the uncertainty of the model simulations. The default will be to use percentiles.

"))

    else
        print(string("

Please, if you want information about one of the package structures, type:
    1) ","\"","model","\""," for the model generation section
    2) ","\"","simulation","\"",", ","\"","simulation","\""," or ","\"","simul","\""," for the model simulation section
    3) ","\"","pseudo-data","\"",", ","\"","pseudodata","\"",", ","\"","pseudo data","\""," for the model pseudo-data generation section
    4) ","\"","mle","\"",", ","\"","likelihood","\""," for the maximum likelihood estimation section
    5) ","\"","inference","\"",", ","\"","stan","\"",", ","\"","stan inference","\"",", ","\"","staninference","\""," for the Bayesian parameter inference section
    6) ","\"","oedms","\"",", ","\"","model selection","\"",", ","\"","modelselection","\"",", ","\"","oed model selection","\""," for the optimal experimental design for model selection section
    7) ","\"","oedmc","\"",", ","\"","model calibration","\"",", ","\"","modelcalibration","\"",", ","\"","oed model calibration","\""," for the optimal experimental design for model calibration section

        "))
    end

end
