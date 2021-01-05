## Inputs structure function
function defMLEStruct()

    mle_def = Dict();

    mle_def["Nexp"] = []; # Integer indicating the number of experiments to be simulated
    mle_def["finalTime"] = []; # -> Vector of final times for each simulation (initial time will allways be asumed as 0, so please consider that)
    mle_def["switchT"] = []; # Array with the switching times of the inducer in the simulation (time 0 and final time need to be considered)
    mle_def["y0"] = []; # Array of Y0s for the simulations for each experiment. If you are computing the steady state this vector might not be used, however you still need to introduce it
    mle_def["preInd"] = []; # Level of the inducer in the ON. This entry can be empty if no Steady State equations are given and only the Y0 given by the user is considered.
    mle_def["uInd"] = []; # Array with the values for the inducer for each experiment at each step
    # mle_def["theta"] = []; # Array with the parameter samples or directory and file location of CSV file with them
    mle_def["tsamps"] = []; # Array of Sampling times vectors
    mle_def["plot"] = []; # Bollean or yes/no string to save the resulting optimisation plots in the results directory
    mle_def["flag"] = []; # String to attach a unique flag to the generated files so it is not overwritten. If empty, nothing will be added.


    mle_def["thetaMAX"] = []; # Vector containing the maximum bounds for theta (no files can be introduced)
    mle_def["thetaMIN"] = []; # Vector containing the minimum bounds for theta (no files can be introduced)
    mle_def["runs"] = []; # Integer indicating how many runs of Optimisation will be done

    # Optimiser does not allow an initial guess for now, but apparently they are working on it. Will just comment this for now.
#     mle_def["thetaGUESS"] = []; # Initial guess for theta for each run. This can be empty, so a random initial theta will be generated for each run


    mle_def["parallel"] = []; # Bollean or yes/no string indicating if the different runs whant to be done in parallel (true) or series (false). Default is false.

    # For the two following fields, you can introduce a string pointing to the observable files (same strings in)
    #     both fields having the same structure as the ones generated in the PseudoData section. If multiple theta are
    #     considered in the file, then the covariance matrix will be taken.
    mle_def["DataMean"] = []; # Array containin ght evector of means for each experiment.
    mle_def["DataError"] = []; # Array containing the vector or matrices (covariance included) of errors for the data for each experiment

    mle_def["Obs"] = []; # States that are observable. This is either a vecotr of strings (that has the same order asmodel_def["stName"] and model_def["eqns"]), a vector of integers insidating which states are observable or the string "All" if all states are observables. This can also be an experssion combining states (Only +,-,*,\ and ^ will be considered)

    mle_def["OPTsolver"] = []; # For now we only use the package BlackBoxOptim, so any of their options can be used. The default is adaptive_de_rand_1_bin_radiuslimited.
    mle_def["MaxTime"] = []; # Maximum number of time allowed for the optimisation as a stop criterion. If this is selected MaxFuncEvals has to be empty (only 1 stop criterion). If both are empty, the default will be to do 1000 function evaluations
    mle_def["MaxFuncEvals"] = []; # Maximum number of function evaluations allowed for the optimisation as a stop criterion. If this is selected MaxTime has to be empty (only 1 stop criterion). If both are empty, the default will be to do 1000 function evaluations


    return(mle_def)
end


## Function to extract the fields from a simulation structure if that is given (this can also be done with the pseudodata file)
function SimToMle(mle_def, simul_def)

    for (i,j) in simul_def
        if haskey(mle_def, i)
            mle_def[i] = simul_def[i];
        end
    end

    return(mle_def)
end

## Function to check the structure introduced by the user
function checkStructMLE(model_def, mle_def)

    # Check taht all the dictionary entries are correct ---- , "thetaGUESS"
    entries = ["Nexp", "finalTime", "switchT", "y0", "preInd", "uInd", "tsamps", "plot", "flag",
                "thetaMAX", "thetaMIN", "runs", "parallel", "DataMean", "DataError", "Obs",
                "OPTsolver", "MaxTime", "MaxFuncEvals"]
    if symdiff(entries,keys(mle_def))!=[] && symdiff(entries,keys(mle_def)) != ["savepath", "savename"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is soemthign wrong...")
        println(symdiff(entries,keys(mle_def)))
        return
    end

    # Check one to be sure there is no empty entry
    nee = ["Nexp", "finalTime", "switchT", "y0", "uInd", "tsamps", "thetaMAX", "thetaMIN", "runs","DataMean","DataError","Obs"]; # No empty entries
    for i in 1:length(nee)
        if mle_def[nee[i]] == []
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Plese, check the field ", nee[i], ", this cannot be empty"))
            return
        end
    end

    # Checks to see if the content of the fields is the expected (also consier if single value entries have been given as a number or a vector)
    if ((typeof(mle_def["Nexp"][1])!=Int) || length(mle_def["Nexp"])!=1) && (typeof(mle_def["Nexp"])!=Int)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Nexp! This should be an integer. ")
        return
    elseif (typeof(mle_def["finalTime"]) != Array{Int,1}) && (typeof(mle_def["finalTime"]) != Array{Float64,1}) && (typeof(mle_def["finalTime"]) != Array{Int32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field finalTime! This should be a vector of final times. ")
        return
    elseif (typeof(mle_def["switchT"]) != Array{Array{Int,1},1}) && (typeof(mle_def["switchT"]) != Array{Array{Float64,1},1}) && (typeof(mle_def["switchT"]) != Array{Array{Float32,1},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field switchT! This should be an array of arrays of numbers. ")
        return
    elseif (typeof(mle_def["y0"]) != Array{Array{Int,1},1}) && (typeof(mle_def["y0"]) != Array{Array{Float64,1},1}) && (typeof(mle_def["y0"]) != Array{Array{Float32,1},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field y0! This should be an array of arrays of numbers. ")
        return
    elseif (typeof(mle_def["preInd"]) != Array{Array{Int,1},1}) && (typeof(mle_def["preInd"]) != Array{Array{Float64,1},1}) && (typeof(mle_def["preInd"]) != Array{Array{Float32,1},1}) && mle_def["preInd"] != []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field preInd! This should be an array of arrays of numbers. ")
        return
    elseif (typeof(mle_def["tsamps"]) != Array{Array{Int,1},1}) && (typeof(mle_def["tsamps"]) != Array{Array{Float64,1},1}) && (typeof(mle_def["tsamps"]) != Array{Array{Float32,1},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field tsamps! This should be an array of arrays of numbers. ")
        return
    elseif (mle_def["plot"] != [])
        try
            if (typeof(mle_def["plot"]) != Bool) && (typeof(mle_def["plot"]) != Array{Bool,1}) &&
            (typeof(mle_def["plot"]) != String) && (typeof(mle_def["plot"]) != Array{String,1})
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check the field plot! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
                return
            end
        catch
             if (typeof(mle_def["plot"]) != Array{Bool,1}) &&
            (typeof(mle_def["plot"]) != Array{String,1})
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check the field plot! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
                return
            end
        end
    elseif (typeof(mle_def["flag"]) != Array{String,1}) && (typeof(mle_def["flag"]) != String) && ((mle_def["flag"]) != [])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field flag! This should be a vector, matrix or string. ")
        return
    elseif (typeof(mle_def["thetaMAX"]) != Array{Float64,1}) && (typeof(mle_def["thetaMAX"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field thetaMAX! This should be a vector. ")
        return
    elseif (typeof(mle_def["thetaMIN"]) != Array{Float64,1}) && (typeof(mle_def["thetaMIN"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field thetaMIN! This should be a vector. ")
        return
    elseif ((typeof(mle_def["runs"])!=Array{Int,1}) || length(mle_def["runs"])!=1) && (typeof(mle_def["runs"])!=Int)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Nexp! This should be an integer. ")
        return
#     elseif (typeof(mle_def["thetaGUESS"]) != Array{Float64,1}) && (typeof(mle_def["thetaGUESS"]) != Array{Float32,1}) &&
#             (typeof(mle_def["thetaGUESS"]) != Array{Float64,2}) && (typeof(mle_def["thetaGUESS"]) != Array{Float32,2}) &&
#             (typeof(mle_def["thetaGUESS"]) != Array{String,1}) && (typeof(mle_def["thetaGUESS"]) != String)
#         println("-------------------------- Process STOPPED!!! --------------------------")
#         println("Please, check the field thetaGUESS! This should be a vector, matrix or string. ")
#         return
    elseif (mle_def["parallel"] != [])
        try
            if (typeof(mle_def["parallel"]) != Bool) && (typeof(mle_def["parallel"]) != Array{Bool,1}) &&
            (typeof(mle_def["parallel"]) != String) && (typeof(mle_def["parallel"]) != Array{String,1})
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check the field parallel! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
                return
            end
        catch
             if (typeof(mle_def["parallel"]) != Array{Bool,1}) &&
            (typeof(mle_def["parallel"]) != Array{String,1})
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check the field parallel! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
                return
            end
        end
    elseif (typeof(mle_def["Obs"]) != Array{String,1}) && (typeof(mle_def["Obs"]) != Array{Int,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Obs! This should be a vector of strings or integers. ")
        return
    # elseif (typeof(mle_def["DataMean"]) != Array{String,1}) && (typeof(mle_def["DataMean"]) != Array{Array{Float64,1},1})
    #     && (typeof(mle_def["DataMean"]) != Array{Array{Float32,1},1}) && (typeof(mle_def["DataMean"]) != Array{Array{Float32,2},1}) &&
    #     (typeof(mle_def["DataMean"]) != Array{Array{Float64,2},1})
    #     println("-------------------------- Process STOPPED!!! --------------------------")
    #     println("Please, check the field DataMean! This should be an array of vectors, matrices or an array of strings. ")
    #     return
    # elseif (typeof(mle_def["DataError"]) != Array{String,1}) && (typeof(mle_def["DataError"]) != Array{Array{Float64,1},1}) &&
    #     (typeof(mle_def["DataError"]) != Array{Array{Float32,1},1})&& (typeof(mle_def["DataError"]) != Array{Array{Float64,2},1}) &&
    #     (typeof(mle_def["DataError"]) != Array{Array{Float64,2},1})
    #     println("-------------------------- Process STOPPED!!! --------------------------")
    #     println("Please, check the field DataError! This should be an array of vectors, matrices or an array of strings. ")
    #     return
    elseif (mle_def["OPTsolver"]!=[]) && (typeof(mle_def["OPTsolver"])!=String) && (typeof(mle_def["OPTsolver"])!=Array{String,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field OPTsolver! This should be a string containing the solver or an empty vector. ")
        return
    elseif (mle_def["MaxTime"]!=[]) && (typeof(mle_def["MaxTime"])!=Int)&& (typeof(mle_def["MaxTime"])!=Array{Int,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field MaxTime! This should be an integer or an empty vector. ")
        return
    elseif (mle_def["MaxFuncEvals"]!=[]) && (typeof(mle_def["MaxFuncEvals"])!=Int)&& (typeof(mle_def["MaxFuncEvals"])!=Array{Int,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field MaxFuncEvals! This should be an integer or an empty vector. ")
        return
    end

    # Extract necessary elements to ease generalisation
    if typeof(mle_def["Nexp"]) == Array{Int,1}
        mle_def["Nexp"] = mle_def["Nexp"][1];
    end

    if mle_def["plot"]==[true] || mle_def["plot"]==["Yes"] || mle_def["plot"]==["yes"] || mle_def["plot"]=="Yes" || mle_def["plot"]=="yes"
        mle_def["plot"]=true
    elseif mle_def["plot"]==[false] || mle_def["plot"]==["No"] || mle_def["plot"]==["no"] || mle_def["plot"]=="No" || mle_def["plot"]=="no" || mle_def["plot"]==[]
        mle_def["plot"]=false
    end

#     if (typeof(mle_def["thetaGUESS"]) == Array{String,1})
#         mle_def["thetaGUESS"] = mle_def["thetaGUESS"][1];
#     end

    if typeof(mle_def["flag"]) == Array{String,1}
        mle_def["flag"] = mle_def["flag"][1];
    elseif typeof(mle_def["flag"]) == String
        mle_def["flag"] = mle_def["flag"];
    elseif ((mle_def["flag"]) == [])
        mle_def["flag"] = "";
    end

    if typeof(mle_def["runs"]) == Array{Int,1}
        mle_def["runs"] = mle_def["runs"][1];
    end

    if mle_def["parallel"]==[true] || mle_def["parallel"]==["Yes"] || mle_def["parallel"]==["yes"] || mle_def["parallel"]=="Yes" || mle_def["parallel"]=="yes"
        mle_def["parallel"]=true;
    elseif mle_def["parallel"]==[false] || mle_def["parallel"]==["No"] || mle_def["parallel"]==["no"] || mle_def["parallel"]=="No" || mle_def["parallel"]=="no" || mle_def["parallel"]==[]
        mle_def["parallel"]=false;
    end

    if (mle_def["MaxTime"]!=[]) && (mle_def["MaxFuncEvals"]!=[])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, only introduce the value for one stopping criteria for the optimiser or leave both empty!")
        return
    elseif (mle_def["MaxTime"]==[]) && (mle_def["MaxFuncEvals"]==[])
        mle_def["MaxFuncEvals"] = 1000;
    elseif (typeof(mle_def["MaxTime"])==Array{Int,1})
        mle_def["MaxTime"] = mle_def["MaxTime"][1];
    elseif (typeof(mle_def["MaxFuncEvals"])==Array{Int,1})
        mle_def["MaxFuncEvals"] = mle_def["MaxFuncEvals"][1];
    end

    if (typeof(mle_def["OPTsolver"])==Array{String,1})
        mle_def["OPTsolver"] = mle_def["OPTsolver"][1];
    elseif mle_def["OPTsolver"] == []
        mle_def["OPTsolver"] = "adaptive_de_rand_1_bin_radiuslimited";
    end

    # Check that all the contents make sense
    if (mle_def["Nexp"] != length(mle_def["finalTime"])) || (mle_def["Nexp"] != length(mle_def["switchT"])) ||
        (mle_def["Nexp"] != length(mle_def["y0"])) || (mle_def["Nexp"] != length(mle_def["uInd"])) ||
        (mle_def["Nexp"] != length(mle_def["tsamps"]))
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, some field does not match the number of experiments selected. Check that the number of")
        println("entries for finalTime, switchT, y0, uInd or tsamps matches the number of experiments in Nexp.")
        return
    end

    for i in 1:mle_def["Nexp"]
        if mle_def["tsamps"][i][end] > mle_def["finalTime"][i]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check finalTime. You have selected a sampling point past it.")
            return
        end

        if length(mle_def["uInd"][i]) != (length(mle_def["switchT"][i])-1)
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check uInd and switchT. Number of steps does not match the number of values for the inputs.")
            return
        end

        if length(mle_def["y0"][i]) != model_def["nStat"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but Y0 does not match the number or states in experiment ", i))
            return
        end
    end

    if (model_def["Y0eqs"] != []) && mle_def["preInd"]==[]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("You specified computation of Y0 as steady state but you have not introduced an inducer value for it!")
        return
    end


#     if (typeof(mle_def["thetaGUESS"]) == Array{Float64,1}) || (typeof(mle_def["thetaGUESS"]) == Array{Float32,1})
#         if length(mle_def["thetaGUESS"]) != model_def["nPar"]
#             println("-------------------------- Process STOPPED!!! --------------------------")
#             println("Number of parameters introduced for thetaGUESS does not match the specidied")
#             return
#         end
#     elseif (typeof(mle_def["thetaGUESS"]) == Array{Float64,2}) || (typeof(mle_def["thetaGUESS"]) == Array{Float32,2})
#         if size(mle_def["thetaGUESS"])[1] != model_def["nPar"] && size(mle_def["thetaGUESS"])[2] != model_def["nPar"]
#             println("-------------------------- Process STOPPED!!! --------------------------")
#             println("Number of parameters introduced for thetaGUESS does not match the specidied")
#             return
#         end
#     elseif (typeof(mle_def["thetaGUESS"]) == String)
#         if mle_def["thetaGUESS"][end-3:end] != ".csv"
#             println("-------------------------- Process STOPPED!!! --------------------------")
#             println("Please, add the .csv termination to the theta file! And no spaces after!")
#             return
#         end
#         if isfile(mle_def["thetaGUESS"])
#             mle_def["thetaGUESS"] = Matrix(CSV.read(mle_def["thetaGUESS"]));
#         else
#             println("-------------------------- Process STOPPED!!! --------------------------")
#             println("Sorry, but the file path you introduced does not exist!")
#             return
#         end
#     end

    if (model_def["nPar"] != length(mle_def["thetaMAX"])) || (model_def["nPar"] != length(mle_def["thetaMIN"]))
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, the dimensions of the bound vectors for theta does not match the number of parameters in the model!")
        return
    end

#     if mle_def["thetaGUESS"] != []
#         if length(mle_def["thetaGUESS"])/mle_def["runs"] != model_def["nPar"]
#             println("-------------------------- Process STOPPED!!! --------------------------")
#             println("Sorry, the dimensions of the thetaGUESS introduced does not match the number of runs of optimisation!")
#             return
#         end
#     end

    if length(mle_def["DataMean"]) != length(mle_def["DataError"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, the dimensions the data means and errors do not match!")
        return
    end

    if length(mle_def["DataMean"]) != mle_def["Nexp"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, the dimensions the data and number of experiments indicated do not match!")
        return
    end

    if (typeof(mle_def["DataMean"]) == Array{String,1} && typeof(mle_def["DataError"]) != Array{String,1} ) ||
        (typeof(mle_def["DataMean"]) != Array{String,1} && typeof(mle_def["DataError"]) == Array{String,1} )
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, if files are introduced for the data please introduce the same in DataMean and DataError!")
        return
    end

    optsolvs = ["separable_nes","dxnes","xnes","adaptive_de_rand_1_bin","adaptive_de_rand_1_bin_radiuslimited",
        "de_rand_1_bin","de_rand_1_bin_radiuslimited","de_rand_2_bin","de_rand_2_bin_radiuslimited",
        "generating_set_search","probabilistic_descent","resampling_memetic_search",
        "resampling_inheritance_memetic_search","simultaneous_perturbation_stochastic_approximation","random_search"];
    if mle_def["OPTsolver"] in optsolvs
        nothing
    else
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, select a correct optimisation solver method from the package BlackBoxOptim.jl!")
        println("The list is: ")
        println(optsolvs)
        return
    end

    # Load data if a file is given
    if typeof(mle_def["DataMean"]) == Array{String,1}

        dat = Array{Any,1}(undef,mle_def["Nexp"]);
        err = Array{Any,1}(undef,mle_def["Nexp"]);

        for i in 1:mle_def["Nexp"]
            if mle_def["DataMean"][i][end-3:end] != ".csv"
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, add the .csv termination to the data files! And no spaces after!")
                return
            end
            if isfile(mle_def["DataMean"][i])
                try
                    tmpD = Matrix(CSV.read(mle_def["DataMean"][i]));
                    if size(tmpD)[2] != 3*length(mle_def["Obs"])
                        tmpdat = zeros(size(tmpD)[1],convert(Int,size(tmpD)[2]/(3*length(mle_def["Obs"]))), length(mle_def["Obs"]));
                        co = 2;
                        for k in 1:convert(Int,size(tmpD)[2]/(3*length(mle_def["Obs"])))
                            for j in 1:length(mle_def["Obs"])
                                tmpdat[:,k] = tmpD[:,co,j];
                                co +=3;
                            end
                        end
                        dat[i] = mean(tmpdat, dims=2);
                        tmco = Array{Any}(undef,length(mle_def["Obs"]))
                        for j in 1:length(mle_def["Obs"])
                            tmco[j] = cov(tmpdat[:,:,j]');
                        end
                        err[i] = tmco;
                    else
                        tmpdat = zeros(size(tmpD)[1], length(mle_def["Obs"]));
                        tmperr = Array{Any}(undef,length(mle_def["Obs"]));;

                        co = 2;
                        for j in 1:length(mle_def["Obs"])
                            tmpdat[:,j] = tmpD[:,co];
                            tmperr[j] = tmpD[:,co+1];
                            co +=3;
                        end
                        dat[i] = tmpdat;
                        err[i] = tmperr;

                    end
                catch
                    println("-------------------------- Process STOPPED!!! --------------------------")
                    println("Sorry, but there is an issue with the data files introduced!")
                    return
                end

            else
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Sorry, but the file path you introduced does not exist!")
                return
            end
        end

        mle_def["DataMean"] = dat;
        mle_def["DataError"] = err;

    end

    # Check that the dimensions of the data make sense
    for i in 1:mle_def["Nexp"]

        if (size(mle_def["DataMean"][i])[2] != length(mle_def["Obs"])) || (length(mle_def["DataError"][i]) != length(mle_def["Obs"]))
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but the data dimensions do not match the number of observables!")
            return
        end

        if (size(mle_def["DataMean"][i])[1] != length(mle_def["tsamps"][i]))
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but the number of data points do not match the number of sampling points!")
            return
        end

        for j in 1:length(mle_def["Obs"])
            if size(mle_def["DataError"][i][j])[1] != length(mle_def["tsamps"][i])
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Sorry, but the number of error points do not match the number of sampling points!")
                return
            end

            if length(size(mle_def["DataError"][i][j])) == 2
                if (size(mle_def["DataError"][i][j]))[1] != (size(mle_def["DataError"][i][j]))[2]
                    println("-------------------------- Process STOPPED!!! --------------------------")
                    println("Sorry, but it seems that you have introduced a covariance matrix for the data but this is not simetric!")
                    return
                end
            end
        end

    end

    # Check warning in case the matrix introduced is simetric
#     if length(size(mle_def["thetaGUESS"])) == 2
#         if size(mle_def["thetaGUESS"])[1] == size(mle_def["thetaGUESS"])[2]
#             println("-------------------------- WARNING --------------------------")
#             println("Sorry, but the number of rows and columns of the theta matrix is the same, so the checks on correct ")
#             println("correct orientation will not work. Please make sure that the dimensions follow: ")
#             println("theta[samples, parameters]")
#         end
#     end

    # Check that no step is no smaller than 2 unit of time
    for i in 1:mle_def["Nexp"]
        for j in 1:length(mle_def["uInd"][i])
            if (mle_def["switchT"][i][j+1]-mle_def["switchT"][i][j])<=4
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but 2 of the steps in experiment ", i, " are too close. This package cannot "))
                println("handle this for now. We are working on it!")
                return
            end
        end
    end

    # Check that the inputs for the experiment are given as collumns
    try
        if model_def["nInp"]>1
            for i in 1:mle_def["Nexp"]
                if size(mle_def["uInd"][i])[2] != model_def["nInp"];
                    mle_def["uInd"][i] = mle_def["uInd"][i]';
                end
            end
        end
    catch
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but there is some issue with the contents of the field uInd."))
        return
    end

    if mle_def["Obs"] == Array{Int,1}
        if (maximum(mle_def["Obs"]) > model_def["nStat"]) || (length(mle_def["Obs"]) > model_def["nStat"])
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but there is some issue with the contents of the field Obs. It seems that "))
            println("     you have selected more observables than states or an index higher than the number of observables")
            return
        end
    end

    if typeof(mle_def["Obs"]) == Array{String,1}
        if convert(Bool,sum(occursin.("+", mle_def["Obs"]))) || convert(Bool,sum(occursin.("-", mle_def["Obs"]))) ||
           convert(Bool,sum(occursin.("/", mle_def["Obs"]))) ||
           convert(Bool,sum(occursin.("*", mle_def["Obs"]))) || convert(Bool,sum(occursin.("^", mle_def["Obs"])))
            nothing
        else
            if sum(occursin.(model_def["stName"], mle_def["Obs"])) == 0
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but there is some issue with the contents of the field Obs."))
                println("It seems that the observable(s) selected do not match any state")
                return
            else
                mle_def["Obs"] = sort([findall(x->x==mle_def["Obs"][i], model_def["stName"])[1] for i in 1:length(mle_def["Obs"])])
            end
        end

    elseif typeof(mle_def["Obs"]) == Array{Int,1}
        mle_def["Obs"] = sort(mle_def["Obs"]);
    end

    return(mle_def)
end

function selectObsSim_te(simul, Obs, stName)

    simObs = zeros(size(simul)[1], length(Obs), size(simul)[3]);
    if typeof(Obs)==Array{Int,1}
        simObs = simul[:,(Obs),:];
    else
        coun = 1;
        for j in 1:length(Obs)
           if !(occursin.("+", Obs[j])) && !(occursin.("-", Obs[j])) &&
           !(occursin.("/", Obs[j])) && !(occursin.("*", Obs[j])) &&
            !(occursin.("^", Obs[j]))
                tind = findall(x->x==Obs[j], stName)[1];
                simObs[:,coun,:] = simul[:,tind,:];
                coun +=1;
            else #This is the complicated stuff where we need to evaluate the function written as a string
                if (occursin(".+", Obs[j]))
                    Obs[j] = replace(Obs[j], ".+"=>" .+ ")
                elseif (occursin("+", Obs[j]))
                    Obs[j] = replace(Obs[j], "+"=>" .+ ")
                end

                if (occursin(".-", Obs[j]))
                    Obs[j] = replace(Obs[j], ".-"=>" .- ")
                elseif (occursin("-", Obs[j]))
                    Obs[j] = replace(Obs[j], "-"=>" .- ")
                end

                if (occursin("./", Obs[j]))
                    Obs[j] = replace(Obs[j], "./"=>" ./ ")
                elseif (occursin("/", Obs[j]))
                    Obs[j] = replace(Obs[j], "/"=>" ./ ")
                end

                if (occursin(".*", Obs[j]))
                    Obs[j] = replace(Obs[j], ".*"=>" .* ")
                elseif (occursin("*", Obs[j]))
                    Obs[j] = replace(Obs[j], "*"=>" .* ")
                end

                if (occursin(".^", Obs[j]))
                    Obs[j] = replace(Obs[j], ".^"=>" .^ ")
                elseif (occursin("^", Obs[j]))
                    Obs[j] = replace(Obs[j], "^"=>" .^ ")
                end

                tmp1 = Obs[j];
                for k in 1:length(stName)
                    if occursin(stName[k], Obs[j])
                        tmp1 = replace(tmp1, stName[k]=>string(simul[:,k,:]))
                    end
                end
                simObs[:,coun,:] = eval(Meta.parse(tmp1));
                coun +=1;
            end
        end
    end
    return(simObs)
end

function restructInputs_te(model_def, mle_def, expp)
    r = convert.(Int, 1:model_def["nInp"]:(length(mle_def["uInd"][expp])));
    inputs = zeros(convert.(Int,length(mle_def["uInd"][expp])));
    for j in 1:convert.(Int,length(mle_def["uInd"][expp])/model_def["nInp"])
        for k in 0:(model_def["nInp"]-1)
            inputs[r[j]+k] = mle_def["uInd"][expp][j,(k+1)];
        end
    end
    return(inputs)
end

# Univariate
function UVloglike(dats, mes, errs)

    # Check that the dimensions of all the vectors are the same
    if (size(dats) != size(mes)) || (size(dats) != size(errs))
        println("Please, check the vectors dimensions!")
        return
    end

    # Computation of log likelihood
    nd = length(dats); # Number of time points
    llk_nd = zeros(nd);

    for i in 1:nd
        tmp1 = log(2*pi);
        tmp2 = log(errs[i]^2);
        tmp3 = ((mes[i]-dats[i])^2)/(errs[i]^2);

        llk_nd[i] = (-0.5)*(tmp1+tmp2+tmp3);

        # Check to avoid Nans (should not be required)
        if isnan(llk_nd[i])
            llk_nd[i] = 0;
        end
    end
    llk = sum(llk_nd);

    return(llk);
end

# Multivariate
function MVloglike(dats, mes, errs)

    # Check that the dimensions of all the vectors are the same
    if (length(dats) != length(mes))
        println("Please, check the vectors dimensions!")
        return
    end

    # Constant term
    tmp1 = length(dats)*log(2*pi);

    # Determinant term
    correcmat = Diagonal(ones(length(dats))).*0.1;
    Em1 = det(errs.+correcmat)
    if Em1 == convert(Float64, Inf)
        Em1 = 1e300
    elseif Em1 == convert(Float64, -Inf)
        Em1 = 1e-300
    end
    tmp2 = log(Em1)

    # Distance term
    dis = dats.-mes;
    inn = inv(errs.+correcmat);

    if size(dis)[1] == 1 # Check for correct orientation of vectors
        t1 = dis;
        t2 = dis';
    else
        t1 = dis';
        t2 = dis;
    end
    tmp3 = t1*inn*t2;
    tmp3 = tmp3[1];


    # LogLikelihood
    llk = -0.5*(tmp1+tmp2+tmp3)

    return(llk);

end

## Plot results function
function plotMLEResults(mle_res,model_def,mle_def)

    cudi = pwd(); #@__DIR__;

    # Convergence curve plots
    x1 = [mle_res["convCurv"][1][i][1] for i in 1:length(mle_res["convCurv"][1])];
    y1 = [mle_res["convCurv"][1][i][2] for i in 1:length(mle_res["convCurv"][1])];
    pc = plot(x1,y1, ylabel = "CFU", xlabel = "Function Evaluations",  yscale = :log10, linetype=:step, label = "")
    for j in 2:length(mle_res["convCurv"])
        x = [mle_res["convCurv"][j][i][1] for i in 1:length(mle_res["convCurv"][j])];
        y = [mle_res["convCurv"][j][i][2] for i in 1:length(mle_res["convCurv"][j])];
        pc = plot!(x,y, ylabel = "CFU", xlabel = "Function Evaluations",  yscale = :log10, linetype=:step, label = "")
    end
    savefig(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\Plot_MLEConvergence_", mle_def["flag"], ".png"));

    # Simulation Results
    simul_def = defSimulStruct()
    simul_def = SimToMle(simul_def, mle_def)
    simul_def["plot"] = false
    simul_def["theta"] = convert(Array, mle_res["Theta"]);
    simul_def["flag"] = "MLEsimulations1";

    simuls, model_def, simul_def = simulateODEs(model_def, simul_def);
    simul_def["theta"] = convert(Array, mle_res["BestTheta"]);
    simul_def["flag"] = "MLEsimulations2";
    simulsBest, ~, ~ = simulateODEs(model_def, simul_def);

    for i in 1:mle_def["Nexp"]

        tit = "";
        yl1 = "";
        for k in 1:length(mle_def["Obs"])
            if typeof(mle_def["Obs"]) == Array{Int,1};
                tit = hcat(tit, string(model_def["stName"][mle_def["Obs"][k]]));
            else
                tit = hcat(tit, string(mle_def["Obs"][k]));
            end

            yl1 = hcat(yl1, "y");
        end
        tit = tit[:,2:end];
        yl1 = yl1[:,2:end];

        titu = "";
        yl2 = "";
        for k in 1:model_def["nInp"]
            titu = hcat(titu, string(model_def["inpName"][k]));
            yl2 = hcat(yl2, "u");
        end
        titu = titu[:,2:end];
        yl2 = yl2[:,2:end];

        tuu = hcat(tit, titu);
        yuu = hcat(yl1, yl2);

        SimObs = selectObsSim_te(simuls[string("Exp_",i)], mle_def["Obs"],model_def["stName"]);
        SimObsBest = selectObsSim_te(simulsBest[string("Exp_",i)], mle_def["Obs"],model_def["stName"]);

        errorss = zeros(size(SimObs)[1],length(mle_def["Obs"]));
        if length(size(mle_def["DataError"][i][1])) == 1
            for k in 1:length(mle_def["Obs"])
                errorss[:,k] = mle_def["DataError"][i][k];
            end
        else
            for k in 1:length(mle_def["Obs"])
                errorss[:,k] = [mle_def["DataError"][i][k][f,f] for f in 1:size(mle_def["DataError"][i][k])[1]];
            end
        end

        # Elements for the inducers
        tu = Array{Array{Any,1}}(undef, length(mle_def["Obs"])+model_def["nInp"])
        [tu[k] = [] for k in 1:length(mle_def["Obs"])]
        su = Array{Array{Any,1}}(undef, length(mle_def["Obs"])+model_def["nInp"])
        [su[k] = [] for k in 1:length(mle_def["Obs"])]
        for k in 1:model_def["nInp"]
            tu[k+length(mle_def["Obs"])] = round.(mle_def["switchT"][i]);
            su[k+length(mle_def["Obs"])] = vcat(mle_def["uInd"][i][:,(k)], mle_def["uInd"][i][end,(k)])
        end

        pl = plot(round.(simul_def["tsamps"][i]), mle_def["DataMean"][i],yerror = errorss,
        layout=length(mle_def["Obs"])+model_def["nInp"], label = "", title = tit, xlabel = "time", ylabel = yuu,
                color = "gray", size = [2000, 1200]);

        for p in 1:size(SimObs)[3]
            pl = plot!(round.(simul_def["tsamps"][i]),SimObs[:,:,p], layout=length(mle_def["Obs"]), label = "", color = "blue");
        end

        pl = plot!(round.(simul_def["tsamps"][i]),SimObsBest[:,:,1], layout=length(mle_def["Obs"]), label = "Best", color = "red")

        pl = plot!(tu, su, layout=length(mle_def["Obs"])+model_def["nInp"], label = "", linetype=:step, title = tuu);

        savefig(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\PlotMLEResults_Exp", i,"_", mle_def["flag"], ".png"));

    end
end

## Cross validation def structure function
# Inputs structure function
function defCrossValMLEStruct()

    cvmle_def = Dict()

    cvmle_def["Nexp"] = []; # Integer indicating the number of experiments to be simulated
    cvmle_def["finalTime"] = []; # -> Vector of final times for each simulation (initial time will allways be asumed as 0, so please consider that)
    cvmle_def["switchT"] = []; # Array with the switching times of the inducer in the simulation (time 0 and final time need to be considered)
    cvmle_def["y0"] = []; # Array of Y0s for the simulations for each experiment. If you are computing the steady state this vector might not be used, however you still need to introduce it
    cvmle_def["preInd"] = []; # Level of the inducer in the ON. This entry can be empty if no Steady State equations are given and only the Y0 given by the user is considered.
    cvmle_def["uInd"] = []; # Array with the values for the inducer for each experiment at each step
    cvmle_def["theta"] = []; # Array with the parameter samples or directory and file location of CSV file with them
    cvmle_def["tsamps"] = []; # Array of Sampling times vectors
    cvmle_def["plot"] = []; # Bollean or yes/no string to save the resulting optimisation plots in the results directory
    cvmle_def["flag"] = []; # String to attach a unique flag to the generated files so it is not overwritten. If empty, nothing will be added.

    # For the two following fields, you can introduce a string pointing to the observable files (same strings in)
    #     both fields having the same structure as the ones generated in the PseudoData section. If multiple theta are
    #     considered in the file, then the covariance matrix will be taken.
    cvmle_def["DataMean"] = []; # Array containin ght evector of means for each experiment.
    cvmle_def["DataError"] = []; # Array containing the vector or matrices (covariance included) of errors for the data for each experiment

    cvmle_def["Obs"] = []; # States that are observable. This is either a vecotr of strings (that has the same order asmodel_def["stName"] and model_def["eqns"]), a vector of integers insidating which states are observable or the string "All" if all states are observables. This can also be an experssion combining states (Only +,-,*,\ and ^ will be considered)

    return(cvmle_def)
end

## Cross-Validation check structure function
# Function to check the structure introduced by the user
function checkStructCrossValMLE(model_def, cvmle_def)

    # Check taht all the dictionary entries are correct ---- , "thetaGUESS"
    entries = ["Nexp", "finalTime", "switchT", "y0", "preInd", "uInd", "tsamps", "plot", "flag",
                "theta", "DataMean", "DataError", "Obs"]

    if symdiff(entries,keys(cvmle_def))!=[] && symdiff(entries,keys(cvmle_def)) != ["savepath", "savename"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is soemthign wrong...")
        println(symdiff(entries,keys(cvmle_def)))
        return
    end

    # Check one to be sure there is no empty entry
    nee = ["Nexp", "finalTime", "switchT", "y0", "uInd", "tsamps", "theta","DataMean","DataError","Obs"]; # No empty entries
    for i in 1:length(nee)
        if cvmle_def[nee[i]] == []
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Plese, check the field ", nee[i], ", this cannot be empty"))
            return
        end
    end

    # Checks to see if the content of the fields is the expected (also consier if single value entries have been given as a number or a vector)
    if ((typeof(cvmle_def["Nexp"][1])!=Int) || length(cvmle_def["Nexp"])!=1) && (typeof(cvmle_def["Nexp"])!=Int)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Nexp! This should be an integer. ")
        return
    elseif (typeof(cvmle_def["finalTime"]) != Array{Int,1}) && (typeof(cvmle_def["finalTime"]) != Array{Float64,1}) && (typeof(cvmle_def["finalTime"]) != Array{Int32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field finalTime! This should be a vector of final times. ")
        return
    elseif (typeof(cvmle_def["switchT"]) != Array{Array{Int,1},1}) && (typeof(cvmle_def["switchT"]) != Array{Array{Float64,1},1}) && (typeof(cvmle_def["switchT"]) != Array{Array{Float32,1},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field switchT! This should be an array of arrays of numbers. ")
        return
    elseif (typeof(cvmle_def["y0"]) != Array{Array{Int,1},1}) && (typeof(cvmle_def["y0"]) != Array{Array{Float64,1},1}) && (typeof(cvmle_def["y0"]) != Array{Array{Float32,1},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field y0! This should be an array of arrays of numbers. ")
        return
    elseif (typeof(cvmle_def["preInd"]) != Array{Array{Int,1},1}) && (typeof(cvmle_def["preInd"]) != Array{Array{Float64,1},1}) && (typeof(cvmle_def["preInd"]) != Array{Array{Float32,1},1}) && cvmle_def["preInd"] != []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field preInd! This should be an array of arrays of numbers. ")
        return
    elseif (typeof(cvmle_def["tsamps"]) != Array{Array{Int,1},1}) && (typeof(cvmle_def["tsamps"]) != Array{Array{Float64,1},1}) && (typeof(cvmle_def["tsamps"]) != Array{Array{Float32,1},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field tsamps! This should be an array of arrays of numbers. ")
        return
    elseif (cvmle_def["plot"] != [])
        try
            if (typeof(cvmle_def["plot"]) != Bool) && (typeof(cvmle_def["plot"]) != Array{Bool,1}) &&
            (typeof(cvmle_def["plot"]) != String) && (typeof(cvmle_def["plot"]) != Array{String,1})
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check the field plot! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
                return
            end
        catch
             if (typeof(cvmle_def["plot"]) != Array{Bool,1}) &&
            (typeof(cvmle_def["plot"]) != Array{String,1})
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check the field plot! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
                return
            end
        end
    elseif (typeof(cvmle_def["flag"]) != Array{String,1}) && (typeof(cvmle_def["flag"]) != String) && ((cvmle_def["flag"]) != [])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field flag! This should be a vector, matrix or string. ")
        return
    elseif (typeof(cvmle_def["theta"]) != Array{Float64,1}) && (typeof(cvmle_def["theta"]) != Array{Float32,1}) &&
            (typeof(cvmle_def["theta"]) != Array{Float64,2}) && (typeof(cvmle_def["theta"]) != Array{Float32,2}) &&
            (typeof(cvmle_def["theta"]) != Array{String,1}) && (typeof(cvmle_def["theta"]) != String)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field theta! This should be a matrix or string. If a single vector is given")
        println("Cross Validation between different samples cannot be performed!")
        return
    elseif (typeof(cvmle_def["Obs"]) != Array{String,1}) && (typeof(cvmle_def["Obs"]) != Array{Int,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Obs! This should be a vector of strings or integers. ")
        return
    # elseif (typeof(cvmle_def["DataMean"]) != Array{String,1}) && (typeof(cvmle_def["DataMean"]) != Array{Array{Float64,1},1})
    #     && (typeof(cvmle_def["DataMean"]) != Array{Array{Float32,1},1}) && (typeof(cvmle_def["DataMean"]) != Array{Array{Float32,2},1}) &&
    #     (typeof(cvmle_def["DataMean"]) != Array{Array{Float64,2},1})
    #     println("-------------------------- Process STOPPED!!! --------------------------")
    #     println("Please, check the field DataMean! This should be an array of vectors, matrices or an array of strings. ")
    #     return
    # elseif (typeof(cvmle_def["DataError"]) != Array{String,1}) && (typeof(cvmle_def["DataError"]) != Array{Array{Float64,1},1}) &&
    #     (typeof(cvmle_def["DataError"]) != Array{Array{Float32,1},1})&& (typeof(cvmle_def["DataError"]) != Array{Array{Float64,2},1}) &&
    #     (typeof(cvmle_def["DataError"]) != Array{Array{Float64,2},1})
    #     println("-------------------------- Process STOPPED!!! --------------------------")
    #     println("Please, check the field DataError! This should be an array of vectors, matrices or an array of strings. ")
    #     return
    end

    # Extract necessary elements to ease generalisation
    if typeof(cvmle_def["Nexp"]) == Array{Int,1}
        cvmle_def["Nexp"] = cvmle_def["Nexp"][1];
    end

    if cvmle_def["plot"]==[true] || cvmle_def["plot"]==["Yes"] || cvmle_def["plot"]==["yes"] || cvmle_def["plot"]=="Yes" || cvmle_def["plot"]=="yes"
        cvmle_def["plot"]=true
    elseif cvmle_def["plot"]==[false] || cvmle_def["plot"]==["No"] || cvmle_def["plot"]==["no"] || cvmle_def["plot"]=="No" || cvmle_def["plot"]=="no" || cvmle_def["plot"]==[]
        cvmle_def["plot"]=false
    end

    if (typeof(cvmle_def["theta"]) == Array{String,1})
        cvmle_def["theta"] = cvmle_def["theta"][1];
    end

    if typeof(cvmle_def["flag"]) == Array{String,1}
        cvmle_def["flag"] = cvmle_def["flag"][1];
    elseif typeof(cvmle_def["flag"]) == String
        cvmle_def["flag"] = cvmle_def["flag"];
    elseif ((cvmle_def["flag"]) == [])
        cvmle_def["flag"] = "";
    end

    # Check that all the contents make sense
    if (cvmle_def["Nexp"] != length(cvmle_def["finalTime"])) || (cvmle_def["Nexp"] != length(cvmle_def["switchT"])) ||
        (cvmle_def["Nexp"] != length(cvmle_def["y0"])) || (cvmle_def["Nexp"] != length(cvmle_def["uInd"])) ||
        (cvmle_def["Nexp"] != length(cvmle_def["tsamps"]))
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, some field does not match the number of experiments selected. Check that the number of")
        println("entries for finalTime, switchT, y0, uInd or tsamps matches the number of experiments in Nexp.")
        return
    end

    for i in 1:cvmle_def["Nexp"]
        if cvmle_def["tsamps"][i][end] > cvmle_def["finalTime"][i]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check finalTime. You have selected a sampling point past it.")
            return
        end

        if length(cvmle_def["uInd"][i]) != (length(cvmle_def["switchT"][i])-1)
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check uInd and switchT. Number of steps does not match the number of values for the inputs.")
            return
        end

        if length(cvmle_def["y0"][i]) != model_def["nStat"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but Y0 does not match the number or states in experiment ", i))
            return
        end
    end

    if (model_def["Y0eqs"] != []) && cvmle_def["preInd"]==[]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("You specified computation of Y0 as steady state but you have not introduced an inducer value for it!")
        return
    end


    if (typeof(cvmle_def["theta"]) == Array{Float64,1}) || (typeof(cvmle_def["theta"]) == Array{Float32,1})
        if length(cvmle_def["theta"]) != model_def["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Number of parameters introduced for theta does not match the specidied")
            return
        end
    elseif (typeof(cvmle_def["theta"]) == Array{Float64,2}) || (typeof(cvmle_def["theta"]) == Array{Float32,2})
        if size(cvmle_def["theta"])[1] != model_def["nPar"] && size(cvmle_def["theta"])[2] != model_def["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Number of parameters introduced for theta does not match the specidied")
            return
        end
    elseif (typeof(cvmle_def["theta"]) == String)
        if cvmle_def["theta"][end-3:end] != ".csv"
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, add the .csv termination to the theta file! And no spaces after!")
            return
        end
        if isfile(cvmle_def["theta"])
            cvmle_def["theta"] = Matrix(CSV.read(cvmle_def["theta"]));
        else
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but the file path you introduced does not exist!")
            return
        end
    end

    if length(cvmle_def["DataMean"]) != length(cvmle_def["DataError"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, the dimensions the data means and errors do not match!")
        return
    end

    if length(cvmle_def["DataMean"]) != cvmle_def["Nexp"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, the dimensions the data and number of experiments indicated do not match!")
        return
    end

    if (typeof(cvmle_def["DataMean"]) == Array{String,1} && typeof(cvmle_def["DataError"]) != Array{String,1} ) ||
        (typeof(cvmle_def["DataMean"]) != Array{String,1} && typeof(cvmle_def["DataError"]) == Array{String,1} )
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, if files are introduced for the data please introduce the same in DataMean and DataError!")
        return
    end

    # Load data if a file is given
    if typeof(cvmle_def["DataMean"]) == Array{String,1}

        dat = Array{Any,1}(undef,cvmle_def["Nexp"]);
        err = Array{Any,1}(undef,cvmle_def["Nexp"]);

        for i in 1:cvmle_def["Nexp"]
            if cvmle_def["DataMean"][i][end-3:end] != ".csv"
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, add the .csv termination to the data files! And no spaces after!")
                return
            end
            if isfile(cvmle_def["DataMean"][i])
                try
                    tmpD = Matrix(CSV.read(cvmle_def["DataMean"][i]));
                    if size(tmpD)[2] != 3*length(cvmle_def["Obs"])
                        tmpdat = zeros(size(tmpD)[1],convert(Int,size(tmpD)[2]/(3*length(cvmle_def["Obs"]))), length(cvmle_def["Obs"]));
                        co = 2;
                        for k in 1:convert(Int,size(tmpD)[2]/(3*length(cvmle_def["Obs"])))
                            for j in 1:length(cvmle_def["Obs"])
                                tmpdat[:,k] = tmpD[:,co,j];
                                co +=3;
                            end
                        end
                        dat[i] = mean(tmpdat, dims=2);
                        tmco = Array{Any}(undef,length(cvmle_def["Obs"]))
                        for j in 1:length(cvmle_def["Obs"])
                            tmco[j] = cov(tmpdat[:,:,j]');
                        end
                        err[i] = tmco;
                    else
                        tmpdat = zeros(size(tmpD)[1], length(cvmle_def["Obs"]));
                        tmperr = Array{Any}(undef,length(cvmle_def["Obs"]));;

                        co = 2;
                        for j in 1:length(cvmle_def["Obs"])
                            tmpdat[:,j] = tmpD[:,co];
                            tmperr[j] = tmpD[:,co+1];
                            co +=3;
                        end
                        dat[i] = tmpdat;
                        err[i] = tmperr;

                    end
                catch
                    println("-------------------------- Process STOPPED!!! --------------------------")
                    println("Sorry, but there is an issue with the data files introduced!")
                    return
                end

            else
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Sorry, but the file path you introduced does not exist!")
                return
            end
        end

        cvmle_def["DataMean"] = dat;
        cvmle_def["DataError"] = err;

    end

    # Check that the dimensions of the data make sense
    for i in 1:cvmle_def["Nexp"]

        if (size(cvmle_def["DataMean"][i])[2] != length(cvmle_def["Obs"])) || (length(cvmle_def["DataError"][i]) != length(cvmle_def["Obs"]))
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but the data dimensions do not match the number of observables!")
            return
        end

        if (size(cvmle_def["DataMean"][i])[1] != length(cvmle_def["tsamps"][i]))
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but the number of data points do not match the number of sampling points!")
            return
        end

        for j in 1:length(cvmle_def["Obs"])
            if size(cvmle_def["DataError"][i][j])[1] != length(cvmle_def["tsamps"][i])
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Sorry, but the number of error points do not match the number of sampling points!")
                return
            end

            if length(size(cvmle_def["DataError"][i][j])) == 2
                if (size(cvmle_def["DataError"][i][j]))[1] != (size(cvmle_def["DataError"][i][j]))[2]
                    println("-------------------------- Process STOPPED!!! --------------------------")
                    println("Sorry, but it seems that you have introduced a covariance matrix for the data but this is not simetric!")
                    return
                end
            end
        end

    end

#     Check warning in case the matrix introduced is simetric
    if length(size(cvmle_def["theta"])) == 2
        if size(cvmle_def["theta"])[1] == size(cvmle_def["theta"])[2]
            println("-------------------------- WARNING --------------------------")
            println("Sorry, but the number of rows and columns of the theta matrix is the same, so the checks on correct ")
            println("correct orientation will not work. Please make sure that the dimensions follow: ")
            println("theta[samples, parameters]")
        end
    end

    # Check that no step is no smaller than 2 unit of time
    for i in 1:cvmle_def["Nexp"]
        for j in 1:length(cvmle_def["uInd"][i])
            if (cvmle_def["switchT"][i][j+1]-cvmle_def["switchT"][i][j])<=4
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but 2 of the steps in experiment ", i, " are too close. This package cannot "))
                println("handle this for now. We are working on it!")
                return
            end
        end
    end

    # Check that the inputs for the experiment are given as collumns
    try
        if model_def["nInp"]>1
            for i in 1:cvmle_def["Nexp"]
                if size(cvmle_def["uInd"][i])[2] != model_def["nInp"];
                    cvmle_def["uInd"][i] = cvmle_def["uInd"][i]';
                end
            end
        end
    catch
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but there is some issue with the contents of the field uInd."))
        return
    end

    if cvmle_def["Obs"] == Array{Int,1}
        if (maximum(cvmle_def["Obs"]) > model_def["nStat"]) || (length(cvmle_def["Obs"]) > model_def["nStat"])
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but there is some issue with the contents of the field Obs. It seems that "))
            println("     you have selected more observables than states or an index higher than the number of observables")
            return
        end
    end

    if typeof(cvmle_def["Obs"]) == Array{String,1}
        if convert(Bool,sum(occursin.("+", cvmle_def["Obs"]))) || convert(Bool,sum(occursin.("-", cvmle_def["Obs"]))) ||
           convert(Bool,sum(occursin.("/", cvmle_def["Obs"]))) ||
           convert(Bool,sum(occursin.("*", cvmle_def["Obs"]))) || convert(Bool,sum(occursin.("^", cvmle_def["Obs"])))
            nothing
        else
            if sum(occursin.(model_def["stName"], cvmle_def["Obs"])) == 0
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but there is some issue with the contents of the field Obs."))
                println("It seems that the observable(s) selected do not match any state")
                return
            else
                cvmle_def["Obs"] = sort([findall(x->x==cvmle_def["Obs"][i], model_def["stName"])[1] for i in 1:length(cvmle_def["Obs"])])
            end
        end

    elseif typeof(cvmle_def["Obs"]) == Array{Int,1}
        cvmle_def["Obs"] = sort(cvmle_def["Obs"]);
    end

    return(cvmle_def)
end

## Plot of Cross-Validation RESULTS
function plotCrossValMLEResults(cvmle_res,model_def,cvmle_def)

    cudi = pwd(); #@__DIR__;

    for i in 1:cvmle_def["Nexp"]

        tit = "";
        yl1 = "";
        for k in 1:length(cvmle_def["Obs"])
            if typeof(cvmle_def["Obs"]) == Array{Int,1};
                tit = hcat(tit, string(model_def["stName"][cvmle_def["Obs"][k]]));
            else
                tit = hcat(tit, string(cvmle_def["Obs"][k]));
            end

            yl1 = hcat(yl1, "y");
        end
        tit = tit[:,2:end];
        yl1 = yl1[:,2:end];

        titu = "";
        yl2 = "";
        for k in 1:model_def["nInp"]
            titu = hcat(titu, string(model_def["inpName"][k]));
            yl2 = hcat(yl2, "u");
        end
        titu = titu[:,2:end];
        yl2 = yl2[:,2:end];

        tuu = hcat(tit, titu);
        yuu = hcat(yl1, yl2);


        errorss = zeros(size(cvmle_res["BestSimObservables"][i])[1],length(cvmle_def["Obs"]));
        if length(size(cvmle_def["DataError"][i][1])) == 1
            for k in 1:length(cvmle_def["Obs"])
                errorss[:,k] = cvmle_def["DataError"][i][k];
            end
        else
            for k in 1:length(cvmle_def["Obs"])
                errorss[:,k] = [cvmle_def["DataError"][i][k][f,f] for f in 1:size(cvmle_def["DataError"][i][k])[1]];
            end
        end


        # Elements for the inducers
        tu = Array{Array{Any,1}}(undef, length(cvmle_def["Obs"])+model_def["nInp"])
        [tu[k] = [] for k in 1:length(cvmle_def["Obs"])]
        su = Array{Array{Any,1}}(undef, length(cvmle_def["Obs"])+model_def["nInp"])
        [su[k] = [] for k in 1:length(cvmle_def["Obs"])]
        for k in 1:model_def["nInp"]
            tu[k+length(cvmle_def["Obs"])] = round.(cvmle_def["switchT"][i]);
            su[k+length(cvmle_def["Obs"])] = vcat(cvmle_def["uInd"][i][:,(k)], cvmle_def["uInd"][i][end,(k)])
        end

        pl = plot(round.(simul_def["tsamps"][i]), cvmle_def["DataMean"][i],yerror = errorss,
            layout=length(cvmle_def["Obs"])+model_def["nInp"], label = "", title = tit, xlabel = "time", ylabel = yuu,
                    color = "gray", size = [2000, 1200]);
        for p in 1:size(cvmle_res["SimObservables"][string("ExpObs_",i)])[3]
            pl = plot!(round.(simul_def["tsamps"][i]),cvmle_res["SimObservables"][string("ExpObs_",i)][:,:,p], layout=length(cvmle_def["Obs"]), label = "", color = "blue");
        end

        pl = plot!(round.(simul_def["tsamps"][i]),cvmle_res["BestSimObservables"][i][:,:,1], layout=length(cvmle_def["Obs"]), label = "Best", color = "red")

        pl = plot!(tu, su, layout=length(cvmle_def["Obs"])+model_def["nInp"], label = "", linetype=:step, title = tuu);

        savefig(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\PlotCrossValMLEResults_Exp", i,"_", cvmle_def["flag"], ".png"));


    end
end

## CrossValidation function for MLE results main function
# Main cross validation function
function CrossValMLE(model_def, cvmle_def)

    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
    end

    cvmle_def = checkStructCrossValMLE(model_def, cvmle_def);

    # Simulation Results
    simul_def = defSimulStruct()
    simul_def = SimToMle(simul_def, cvmle_def)
    simul_def["plot"] = false
    simul_def["flag"] = "MLEsimulations1";

    global model_def2 = model_def;
    global simul_def2 = simul_def;

    simuls, model_def, simul_def = simulateODEs(model_def2, simul_def2);

    SimObs = Dict();
    for k in 1:cvmle_def["Nexp"]
        SimObs[string("ExpObs_", k)] = selectObsSim_te(simuls[string("Exp_",k)], mle_def["Obs"],model_def["stName"]);
    end

    dataM = cvmle_def["DataMean"];
    dataE = cvmle_def["DataError"];

    LLK_Exps = zeros(size(SimObs[string("ExpObs_", i)])[3], cvmle_def["Nexp"]);
    for i in 1:cvmle_def["Nexp"]
        llkobs = zeros(size(SimObs[string("ExpObs_", i)])[3],length(cvmle_def["Obs"]));
        for j in 1:length(cvmle_def["Obs"])
            for k in 1:size(SimObs[string("ExpObs_", i)])[3]
                if length(size(dataE[i][j])) == 1
                    llkobs[k,j] = UVloglike(dataM[i][:,j], SimObs[string("ExpObs_", i)][:,j,k], dataE[i][j]);
                else
                    llkobs[k,j] = MVloglike(dataM[i][:,j], SimObs[string("ExpObs_", i)][:,j,k], dataE[i][j]);
                end
            end
        end
        LLK_Exps[:,i] = sum(llkobs, dims=2);
    end

    LLK = -sum(LLK_Exps, dims=2);

    tpin = findfirst(x->x==minimum(LLK), LLK);

    cvmle_res = Dict();
    cvmle_res["BestTheta"] = cvmle_def["theta"][:,tpin];
    cvmle_res["Costs"] = LLK;
    cvmle_res["Simulations"] = simuls;
    cvmle_res["SimObservables"] = SimObs;
    cvmle_res["BestSimulations"] = [simuls[string("Exp_",i)][:,:,tpin] for i in 1:cvmle_def["Nexp"]];
    cvmle_res["BestSimObservables"] = [SimObs[string("ExpObs_",i)][:,:,tpin] for i in 1:cvmle_def["Nexp"]];

    cvmle_def["savepath"] = string(cudi, "\\Results\\", model_def["NameF"],"_",today());
    cvmle_def["savename"] = string(model_def["NameF"],"_",today(), "_MLEresults_",cvmle_def["flag"],".jld");

    save(string(cvmle_def["savepath"], "\\", cvmle_def["savename"]), "CrossValMLEresults", cvmle_res, "model_def", model_def, "cvmle_def", cvmle_def);

    println("")
    println("----------------------------------------- RESULTS -----------------------------------------")
    println("Cross Validation of MLE results are saved in the directory: ")
    println(string("                 ", cvmle_def["savepath"]))
    println(string("Under the name ",cvmle_def["savename"]))
    println("--------------------------------------------------------------------------------------")
    println("")

    if cvmle_def["plot"] == true
        plotCrossValMLEResults(cvmle_res,model_def,cvmle_def)
        println("")
        println("----------------------------------------- PLOTS -----------------------------------------")
        println("Simulation PLOTS are saved in the directory: ")
        println(string("                 ", cvmle_def["savepath"]))
        println(string("Under the names PlotCrossValMLEResults_Exp(i)_",cvmle_def["flag"],".png"))
        println("--------------------------------------------------------------------------------------")
        println("")
    end

    return cvmle_res, model_def, cvmle_def

end


## Main Function
function MLEtheta(model_def, mle_def)

    # Generate results directory
    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\MLEscripts"))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\MLEscripts"))
    end

    # Check structure introduced by the user
    mle_def = checkStructMLE(model_def, mle_def);

    # Check if multiple runs are going to be run in parallel
    paral = "";
    cores = "";
    if mle_def["parallel"] == true
        paral = "@everywhere ";
        if length(Sys.cpu_info()) > nworkers()
            cores = string("addprocs(length(Sys.cpu_info())-",nworkers(),")");
        else
            cores = "";
        end
    end

    # Packages needed for the simulations
    Head = string("
    using Distributed
    ",join(cores),"

    ",join(paral),"using DifferentialEquations
    ",join(paral),"using OrdinaryDiffEq
    using DiffEqBase
    ",join(paral),"using Sundials
    using Plots
    # using DiffEqGPU
    ",join(paral),"using ODEInterfaceDiffEq
    using DataFrames
    using CSV
    ",join(paral),"using Statistics
    ",join(paral),"using LinearAlgebra
    using JLD
    ",join(paral),"using StatsBase
    ",join(paral),"using Random
    using SharedArrays
    ",join(paral),"using BlackBoxOptim
    ",join(paral),"using BlackBoxOptim: num_func_evals
        ");

    # Include file with the model
    Head2 = string("
    ",join(paral),"include(","\"",replace(string(cudi,"\\ModelsFunctions\\", model_def["NameF"], "_Model.jl"), "\\"=>"\\\\"),"\"",")
        ");

    # Generate wither univariate or multivariate log-likelihood function deppending on the data introdued

    if length(size(mle_def["DataError"][1][1])) != 1 # Multivariate case
    llkfy = string("
    ",join(paral),"function loglike(dats, mes, errs)

        # Check that the dimensions of all the vectors are the same
        if (length(dats) != length(mes))
            println(","\"","Please, check the vectors dimensions!","\"",")
            return
        end

        # Constant term
        tmp1 = length(dats)*log(2*pi);

        # Determinant term
        correcmat = Diagonal(ones(length(dats))).*0.1;
        Em1 = det(errs.+correcmat)
        if Em1 == convert(Float64, Inf)
            Em1 = 1e300
        elseif Em1 == convert(Float64, -Inf)
            Em1 = 1e-300
        end
        tmp2 = log(Em1)

        # Distance term
        dis = dats.-mes;
        inn = inv(errs.+correcmat);

        if size(dis)[1] == 1 # Check for correct orientation of vectors
            t1 = dis;
            t2 = dis';
        else
            t1 = dis';
            t2 = dis;
        end
        tmp3 = t1*inn*t2;
        tmp3 = tmp3[1];


        # LogLikelihood
        llk = -0.5*(tmp1+tmp2+tmp3)

        return(llk);

    end

        ");
    else # Univariate case
            # Define log-likelihood deppending on the case
    llkfy = string("
    ",join(paral),"function loglike(dats, mes, errs)

        # Check that the dimensions of all the vectors are the same
        if (size(dats) != size(mes)) || (size(dats) != size(errs))
            println(","\"","Please, check the vectors dimensions!","\"",")
            return
        end

        # Computation of log likelihood
        nd = length(dats); # Number of time points
        llk_nd = zeros(nd);

        for i in 1:nd
            tmp1 = log(2*pi);
            tmp2 = log(errs[i]^2);
            tmp3 = ((mes[i]-dats[i])^2)/(errs[i]^2);

            llk_nd[i] = (-0.5)*(tmp1+tmp2+tmp3);

            # Check to avoid Nans (should not be required)
            if isnan(llk_nd[i])
                llk_nd[i] = 0;
            end
        end
        llk = sum(llk_nd);

        return(llk);
    end

    ");
    end;

    # Function to extract the observables from the simulations
    selectObsS = string("

    ",join(paral),"function selectObsSim(simul)

        simObs = zeros(size(simul)[1], length(Obs));
        if typeof(Obs)==Array{Int,1}
            simObs = simul[:,(Obs),1];
        else
            coun = 1;
            for j in 1:length(Obs)
               if !(occursin.(","\"","+","\"",", Obs[j])) && !(occursin.(","\"","-","\"",", Obs[j])) &&
               !(occursin.(","\"","/","\"",", Obs[j])) && !(occursin.(","\"","*","\"",", Obs[j])) &&
                !(occursin.(","\"","^","\"",", Obs[j]))
                    tind = findall(x->x==Obs[j], stName)[1];
                    simObs[:,coun] = simul[:,tind,1];
                    coun +=1;
                else #This is the complicated stuff where we need to evaluate the function written as a string
                    tmp1 = Obs[j];
                    for k in 1:length(stName)
                        if occursin(stName[k], Obs[j])
                            tmp1 = replace(tmp1, stName[k]=>string(simul[:,k,1]))
                        end
                    end
                    simObs[:,coun] = eval(Meta.parse(tmp1));
                    coun +=1;
                end
            end
        end
        return(simObs)
    end


    ");


    if mle_def["parallel"] == true
        parfun = string("
    ",join(paral),"function parOptim(opts)
        ress = convert(SharedArray{Float64,2}, zeros(",length(mle_def["thetaMIN"]),",length(opts)));
        arr = @sync @distributed (vcat) for i in 1:length(opts)
            resT = bboptimize(opts[string(","\"","Opt_","\"",", i)]);
            ress[:,i] = resT.archive_output.best_candidate;

            chr = (0,i);
            tmpcv = vcat(chr, convcur[i]);
            tmpcv
        end;
        return ress, arr;
    end

    ");

    else
        parfun = string("
    function parOptim(opts)
        ress = convert(SharedArray{Float64,2}, zeros(",length(mle_def["thetaMIN"]),",length(opts)));
        for i in 1:length(opts)
            resT = bboptimize(opts[string(","\"","Opt_","\"",", i)]);
            ress[:,i] = resT.archive_output.best_candidate;
        end;
        return(ress);
    end
    ");
    end;

    # Generation of objective function and function to run optimisation in parallel
    objec = string("

    ",join(paral),"function objectiveMLE",model_def["NameF"],"(theta)

        LLK_Exps = zeros(1, Nexp);
        simobs = Array{Any,1}(undef,Nexp);
        for i in 1:Nexp
            simul = ",join(model_def["NameF"]),"_SolveAll(tim[i], theta, esp[i], inp[i], ini[i], smp[i], prr[i])
            simobs[i] = selectObsSim(simul);
            llkobs = zeros(1,length(Obs));
            for j in 1:length(Obs)
                llkobs[1,j] = loglike(dataM[i][:,j], simobs[i][:,j], dataE[i][j]);
            end
            LLK_Exps[1,i] = sum(llkobs);
        end
        return(-sum(LLK_Exps))
    end

    ");

    # Function to re-structure inputs to desired form
    inpre = string("

    ",join(paral),"function restructInputs(model_def, mle_def, expp)
        r = convert.(Int, 1:model_def[","\"","nInp","\"","]:(length(mle_def[","\"","uInd","\"","][expp])));
        inputs = zeros(convert.(Int,length(mle_def[","\"","uInd","\"","][expp])));
        for j in 1:convert.(Int,length(mle_def[","\"","uInd","\"","][expp])/model_def[","\"","nInp","\"","])
            for k in 0:(model_def[","\"","nInp","\"","]-1)
                inputs[r[j]+k] = mle_def[","\"","uInd","\"","][expp][j,(k+1)];
            end
        end
        return(inputs)
    end

    ");

    # Check stop criterion for the optimiser
    stopopt = "";
    if mle_def["MaxTime"] != []
        stopopt = string("MaxTime = ", mle_def["MaxTime"]);
    elseif mle_def["MaxFuncEvals"] != []
        stopopt = string("MaxFuncEvals = ", mle_def["MaxFuncEvals"]);
    end



    # We need to spawn the structure mle_def in all cores
    spaw1 = "";
    spaw2 = "";
    spaw3 = "";
    if mle_def["parallel"] == true
        for k in 1:length(Sys.cpu_info())
            spaw1 = vcat(spaw1, string("        @spawnat ",k," mle_def \n"));
            spaw2 = vcat(spaw2, string("        @spawnat ",k," model_def \n"));
            spaw3 = vcat(spaw3, string("        @spawnat ",k," ts \n", "        @spawnat ",k," sp \n", "        @spawnat ",k," ivss \n",
            "        @spawnat ",k," pre \n", "        @spawnat ",k," samps \n", "        @spawnat ",k," inputs \n"));
        end
    end

    # Needed to extract the convergence curve from the parallel version in correct order.
    if mle_def["parallel"] == true
        funrun = string("tet, convcur2 = parOptim(opts);

        convcur3 = Array{Any}(undef,chains)
        i1 = 1;
        i2 = 1;
        k = 1;
        for h in 2:length(convcur2)
            if convcur2[h][1] < convcur2[h-1][1]
                i2 = h-1;
                convcur3[k] = convcur2[i1:i2];
                i1 = h;
                k += 1;
            elseif h == length(convcur2)
                i2 = h;
                convcur3[k] = convcur2[i1:i2];
            end
        end
        ",join(paral),"convcur = Array{Any}(undef,chains);
        for k in 1:chains
            tmpch = convcur3[k][1][2];
            convcur[tmpch] = convcur3[k][2:end];
        end


        ");
    else
        funrun = "tet= parOptim(opts);";
    end;


    # Main function for optimisation
    mainfun = string("

    function RunMLE", model_def["NameF"],"(model_def, mle_def)

",join(spaw1),"
",join(spaw2),"

        ",join(paral),"global Nexp = mle_def[","\"","Nexp","\"","];
        ",join(paral),"global Obs = mle_def[","\"","Obs","\"","];
        ",join(paral),"global stName = model_def[","\"","stName","\"","];

        ",join(paral),"global dataM = mle_def[","\"","DataMean","\"","];
        ",join(paral),"global dataE = mle_def[","\"","DataError","\"","];

        ",join(paral),"ts = Array{Any,1}(undef,mle_def[","\"","Nexp","\"","]);
        ",join(paral),"sp = Array{Any,1}(undef,mle_def[","\"","Nexp","\"","]);
        ",join(paral),"ivss = Array{Any,1}(undef,mle_def[","\"","Nexp","\"","]);
        ",join(paral),"pre = Array{Any,1}(undef,mle_def[","\"","Nexp","\"","]);
        ",join(paral),"samps = Array{Any,1}(undef,mle_def[","\"","Nexp","\"","]);
        ",join(paral),"inputs = Array{Any,1}(undef,mle_def[","\"","Nexp","\"","]);

        for i in 1:mle_def[","\"","Nexp","\"","]
            global ts[i] = collect(0.0:round(mle_def[","\"","finalTime","\"","][i]))';
            global sp[i] = [convert(Int,v) for v in (round.(mle_def[","\"","switchT","\"","][i])')];
            global ivss[i] = mle_def[","\"","y0","\"","][i];
            if model_def[","\"","Y0eqs","\"","] != [] # ON inputs
                global pre[i] = mle_def[","\"","preInd","\"","][i];
            else
                global pre[i] = [];
            end
            global samps[i] = convert.(Int, round.(mle_def[","\"","tsamps","\"","][i]));
            global inputs[i] = restructInputs(model_def, mle_def, i);
        end

",join(spaw3),"

        ",join(paral),"global tim = ts;
        ",join(paral),"global esp = sp;
        ",join(paral),"global inp = inputs;
        ",join(paral),"global ini = ivss;
        ",join(paral),"global prr = pre;
        ",join(paral),"global smp = samps;


        ",join(paral),"lb = mle_def[","\"","thetaMIN","\"","];
        ",join(paral),"up = mle_def[","\"","thetaMAX","\"","];

        ",join(paral),"rang = [Array{Tuple{Int, Int}}(undef, length(lb))];
        ",join(paral),"rang = [(lb[i], up[i]) for i in 1:length(lb)];

        ",join(paral),"chains = mle_def[","\"","runs","\"","];
        ",join(paral),"opts = Dict()


        ",join(paral),"convcur = Array{Any,1}(undef, chains);
        for i in 1:chains
            convcur[i] = Array{Tuple{Int, Float64},1}();
            callback = oc -> push!(convcur[i], (num_func_evals(oc), best_fitness(oc)))
            opts[string(","\"","Opt_","\"",", i)] = bbsetup(objectiveMLE",model_def["NameF"],"; SearchRange = rang, NumDimensions = length(lb),
                Method = :",mle_def["OPTsolver"],", ",join(stopopt),", TraceMode = :silent, CallbackFunction = callback, CallbackInterval = 0.0);
        end


        println(","\"","----------------------------------------- OPTIMISATION STARTS! -----------------------------------------","\"",")

        @time begin
            ",funrun,"
        end

        println(","\"","----------------------------------------- OPTIMISATION ENDED -----------------------------------------","\"",")


        names = model_def[","\"","parName","\"","];

        alltog = Array{Dict{String,Any},1}(undef,chains)
        for i in 1:(chains)
            global tmp = Dict{String, Any}()
            for j in 1:length(names)
                tmp[names[j]] = tet[j,i];
            end
            alltog[i] = tmp
        end

        mle_res = Dict();
        mle_res[","\"","Theta","\"","] = convert(Array,tet); # Best thetas results
        mle_res[","\"","convCurv","\"","] = convcur; # Convergence curves
        mle_res[","\"","StanDict","\"","] = alltog; # Best thetas in the structure Stan needs as initial guess.

        bcfv = zeros(length(convcur))
        for k in 1:length(convcur)
            bcfv[k] = convcur[k][end][2]
        end
        mle_res[","\"","BestCFV","\"","] = bcfv; # Best cost function values for each run

        tpin = findfirst(x->x==minimum(bcfv), bcfv);
        mle_res[","\"","BestTheta","\"","] = tet[:,tpin]; # Best theta value from all runs

        return(mle_res)

    end
    ");

    # Put all together and generate script
    fungen = join([Head, Head2,selectObsS,parfun,llkfy,objec,inpre, mainfun]);

    open(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\MLEscripts\\", model_def["NameF"], "_MLE.jl"), "w") do io
       write(io, fungen);
    end;


    println("")
    println("----------------------------------------- SCRIPTS -----------------------------------------")
    println("MLE scripts have been generated in the directory: ")
    println(string("         ", string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\MLEscripts\\", model_def["NameF"], "_MLE.jl")))
    println("--------------------------------------------------------------------------------------")
    println("")

    # Run Main script and save results
    include(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\MLEscripts\\", model_def["NameF"], "_MLE.jl"))

    if mle_def["parallel"] == true
        for k in 1:length(Sys.cpu_info())
            spaw1a = string("        @spawnat ",k," mle_def \n");
            spaw2a = string("        @spawnat ",k," model_def \n");
            p1 = Meta.parse(spaw1a)
            p2 = Meta.parse(spaw2a)

            @eval $p1
            @eval $p2
        end
    end
    global model_def2 = model_def;
    global mle_def2 = mle_def;

    MLEfun = Symbol(string("RunMLE", model_def["NameF"]));

    mle_res = @eval $MLEfun(model_def2, mle_def2);

    mle_def["savepath"] = string(cudi, "\\Results\\", model_def["NameF"],"_",today());
    mle_def["savename"] = string(model_def["NameF"],"_",today(), "_MLEresults_",mle_def["flag"],".jld");

    save(string(mle_def["savepath"], "\\", mle_def["savename"]), "MLEresults", mle_res, "model_def", model_def, "mle_def", mle_def);

    println("")
    println("----------------------------------------- RESULTS -----------------------------------------")
    println("MLE results are saved in the directory: ")
    println(string("                 ", mle_def["savepath"]))
    println(string("Under the name ",mle_def["savename"]))
    println("--------------------------------------------------------------------------------------")
    println("")

    if mle_def["plot"] == true
        plotMLEResults(mle_res,model_def,mle_def)
        println("")
        println("----------------------------------------- PLOTS -----------------------------------------")
        println("Simulation PLOTS are saved in the directory: ")
        println(string("                 ", mle_def["savepath"]))
        println(string("Under the names PlotMLEResults_Exp(i)_",mle_def["flag"],".png", " and Plot_MLEConvergence_", mle_def["flag"], ".png"))
        println("--------------------------------------------------------------------------------------")
        println("")
    end


    return(mle_res, model_def, mle_def)
end
