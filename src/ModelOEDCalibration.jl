## Define main Structure
function defODEModelCalibrStruct()
    oedmc_def = Dict();

    oedmc_def["Model"] = []; # Model structure for the model Model. Dict
    oedmc_def["Obs"] = [];# Observables of the models. Array of strings.
    oedmc_def["Theta"] = []; # Theta matrix (Bayesian OED) for the model. No single vectors will be allowed

    oedmc_def["y0"] = []; # Array of Y0s for the simulations for each experiment. If you are computing the steady state this vector might not be used, however you still need to introduce it
    oedmc_def["preInd"] = []; # Level of the inducer in the ON. This entry can be empty if no Steady State equations are given and only the Y0 given by the user is considered.
    oedmc_def["finalTime"] = []; # Vector of final time for the experiment (initial time will allways be asumed as 0, so please consider that)
    oedmc_def["switchT"] = []; # Array with the switching times of the inducer in the simulation (time 0 and final time need to be considered)
    oedmc_def["tsamps"] = []; # Vector of Sampling times
    oedmc_def["fixedInp"] = []; # If more than 1 inducer for the system exists, but only 1 input can be dynamic, give a vector of strings indicating which inputs are going to be optimised but as a constant input instead of dynamic. If none, just give an empty vector.
    oedmc_def["fixedStep"] = []; # If you want any of the steps to be fixed to a value. This has to be an empty array if none is fixed or an array of tuples, where each tuple is a step to be fixed. The first entry of the tuple is the index of the step (as an Integer), and the second and array of values for each inducer. Note that the fixed inputs will be ignored, so do not take them into account here.
    oedmc_def["equalStep"] = []; # If you want a series of steps to have the same optimised value (for example if you want to design a pulse experiment) you can introduce inside this array different arrays with the indexes of the steps that will have the same value. The values introduced in each array need to be integers.

    oedmc_def["plot"] = []; # Bollean or yes/no string to save the resulting optimisation plots in the results directory
    oedmc_def["flag"] = []; # String containing an identifier for the saved files that will be generated.

    oedmc_def["uUpper"] = []; # Vector indicating the upper bounds for the inducers
    oedmc_def["uLower"] = []; # Vector indicating the lower bounds for the inducers
    oedmc_def["maxiter"] = []; # Maximum number of iterations for the Bayesian Optimisation. If nothing is introduced a default of 100 iterations will be taken
    oedmc_def["maxtime"] = []; # Maximum time allowed for the Bayesian Optimisation. If nothing is introduced, only the maximum number of iterations will be taken into account.

    oedmc_def["util"] = []; # String indicating entropy or perc (or percentile) as the core of the utility function to compute the uncertainty of the model simulations. The default will be to use percentiles.

    return oedmc_def
end

## Function to check the section Structure
function checkStructOEDMC(oedmc_def)

    # Check that all the fields (and nothing more) are present
    entries = ["Model", "Obs", "Theta", "y0", "preInd", "finalTime", "switchT", "tsamps", "equalStep",
                "fixedInp", "fixedStep", "plot", "flag", "uUpper", "uLower", "maxiter", "maxtime", "util"];
    if symdiff(entries,keys(oedmc_def))!=[] && symdiff(entries,keys(oedmc_def)) != ["savepath", "savename"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is something wrong...")
        println(symdiff(entries,keys(oedmc_def)))
        return
    end

    # Check one to be sure there is no empty entry
    nee = ["Model", "Obs", "Theta", "y0", "preInd", "uUpper", "uLower"]; # No empty entries
    for i in 1:length(nee)
        if oedmc_def[nee[i]] == []
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Plese, check the field ", nee[i], ", this cannot be empty"))
            return
        end
    end

    # Check that the type of all the contents is correct
    if typeof(oedmc_def["Model"]) <: Dict
        oedmc_def["Model"] = checkStruct(oedmc_def["Model"]);
    elseif typeof(oedmc_def["Model"]) <: Array
        oedmc_def["Model"] = checkStruct(oedmc_def["Model"][1]);
    end;


    if (typeof(oedmc_def["Obs"]) != Array{String,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Obs! This should be a vector of strings or integers. ")
        return
    elseif (typeof(oedmc_def["Theta"]) != Array{Float64,1}) && (typeof(oedmc_def["Theta"]) != Array{Float32,1}) &&
            (typeof(oedmc_def["Theta"]) != Array{Float64,2}) && (typeof(oedmc_def["Theta"]) != Array{Float32,2}) &&
            (typeof(oedmc_def["Theta"]) != Array{String,1}) && (typeof(oedmc_def["Theta"]) != String)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Theta! This should be a vector, matrix or string. ")
        return
    elseif (typeof(oedmc_def["y0"]) != Array{Int,1}) && (typeof(oedmc_def["y0"]) != Array{Float64,1}) &&
        (typeof(oedmc_def["y0"]) != Array{Float32,1}) &&
        (typeof(oedmc_def["y0"]) != Array{Int,2}) && (typeof(oedmc_def["y0"]) != Array{Float64,2}) &&
            (typeof(oedmc_def["y0"]) != Array{Float32,2})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field y0! This should be an array of numbers. ")
        return
    elseif (typeof(oedmc_def["preInd"]) != Array{Int,1}) && (typeof(oedmc_def["preInd"]) != Array{Float64,1}) &&
        (typeof(oedmc_def["preInd"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field preInd! This should be an array of numbers. ")
        return
    elseif (typeof(oedmc_def["finalTime"]) != Array{Int,1}) && (typeof(oedmc_def["finalTime"]) != Array{Float64,1}) &&
        (typeof(oedmc_def["finalTime"]) != Array{Float32,1}) &&
        typeof(oedmc_def["finalTime"]) !=Int && typeof(oedmc_def["finalTime"]) !=Float32 &&
        typeof(oedmc_def["finalTime"]) !=Float64
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field finalTime! This should be a vector or a single number")
        return
    elseif (typeof(oedmc_def["switchT"]) != Array{Array{Int,1},1}) &&
        (typeof(oedmc_def["switchT"]) != Array{Array{Float64,1},1}) &&
        (typeof(oedmc_def["switchT"]) != Array{Array{Float32,1},1}) &&
        (typeof(oedmc_def["switchT"]) != Array{Int,1}) &&
        (typeof(oedmc_def["switchT"]) != Array{Float64,1}) &&
        (typeof(oedmc_def["switchT"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field switchT! This should be an array or an array of arrays of numbers")
        return
    elseif (typeof(oedmc_def["tsamps"]) != Array{Array{Int,1},1}) &&
        (typeof(oedmc_def["tsamps"]) != Array{Array{Float64,1},1}) &&
        (typeof(oedmc_def["tsamps"]) != Array{Array{Float32,1},1}) &&
        (typeof(oedmc_def["tsamps"]) != Array{Int,1}) &&
        (typeof(oedmc_def["tsamps"]) != Array{Float64,1}) &&
        (typeof(oedmc_def["tsamps"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field tsamps! This should be an array or an array of arrays of numbers")
        return
    elseif (typeof(oedmc_def["fixedInp"]) != Array{String,1})&& oedmc_def["fixedInp"] != []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field fixedInp! This should be a vector of strings or an empty vector. ")
        return
    elseif oedmc_def["fixedStep"] != [] && typeof(oedmc_def["fixedStep"]) != Array{Tuple{Int,Array{Int64,1}},1} && typeof(oedmc_def["fixedStep"]) != Array{Tuple{Int,Array{Float64,1}},1} &&
        typeof(oedmc_def["fixedStep"]) != Array{Tuple{Int,Array{Float32,1}},1} && typeof(oedmc_def["fixedStep"]) != Array{Tuple{Int,Array{Any,1}},1} &&
        typeof(oedmc_def["fixedStep"]) != Array{Tuple{Int,Array{T,1} where T},1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field fixedStep! This should have the structure of Array{Tuple{Int,Array{Any,1}},1}. ")
        return
    elseif oedmc_def["equalStep"] != [] && typeof(oedmc_def["equalStep"]) != Array{Array{Int,1},1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field equalStep! This should have the structure of Array{Array{Int,1},1}. ")
        println("Remember, it has to be an array where as each entry you have another array(s) with the indexes")
        println("of the steps that will be considered as the same step in the optimisation value wise.")
        return
    elseif (typeof(oedmc_def["flag"]) != Array{String,1}) && (typeof(oedmc_def["flag"]) != String) && ((oedmc_def["flag"]) != [])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field flag! This should be a vector, matrix or string. ")
        return
    elseif (typeof(oedmc_def["plot"]) != Bool) && (typeof(oedmc_def["plot"]) != Array{Bool,1}) &&
            (typeof(oedmc_def["plot"]) != String) && (typeof(oedmc_def["plot"]) != Array{String,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field plot! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
        return
    elseif (typeof(oedmc_def["uUpper"]) != Array{Int,1}) && (typeof(oedmc_def["uUpper"]) != Array{Float64,1}) &&
            (typeof(oedmc_def["uUpper"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field uUpper! This should be a vector of numbers")
        return
    elseif (typeof(oedmc_def["uLower"]) != Array{Int,1}) && (typeof(oedmc_def["uLower"]) != Array{Float64,1}) &&
            (typeof(oedmc_def["uLower"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field uUpper! This should be a vector of numbers")
        return
    elseif oedmc_def["maxiter"] != [] && typeof(oedmc_def["maxiter"]) != Array{Int,1} != typeof(oedmc_def["maxiter"]) != Int
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field maxiter! This should be an integer or an empty vector")
        return
    elseif oedmc_def["maxtime"] != [] && typeof(oedmc_def["maxtime"]) != Array{Int,1} != typeof(oedmc_def["maxtime"]) != Int &&
        oedmc_def["maxtime"] != [Inf] && oedmc_def["maxtime"] != Inf
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field maxtime! This should be an integer or an empty vector")
        return
    elseif oedmc_def["util"] != [] && typeof(oedmc_def["util"]) != String && typeof(oedmc_def["util"]) != Array{String,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field util! This should be an empty vector or a string containing the word entropy or perc.")
        return
    end

    # Extract necessary elements to ease generalisation
    if oedmc_def["plot"]==[true] || oedmc_def["plot"]==["Yes"] || oedmc_def["plot"]==["yes"] || oedmc_def["plot"]=="Yes" || oedmc_def["plot"]=="yes"
        oedmc_def["plot"]=true
    elseif oedmc_def["plot"]==[false] || oedmc_def["plot"]==["No"] || oedmc_def["plot"]==["no"] || oedmc_def["plot"]=="No" || oedmc_def["plot"]=="no" || oedmc_def["plot"]==[]
        oedmc_def["plot"]=false
    end

    if typeof(oedmc_def["flag"]) == Array{String,1}
        oedmc_def["flag"] = oedmc_def["flag"][1];
    elseif typeof(oedmc_def["flag"]) == String
        oedmc_def["flag"] = oedmc_def["flag"];
    elseif ((oedmc_def["flag"]) == [])
        oedmc_def["flag"] = "";
    end

    if oedmc_def["maxiter"] == []
        oedmc_def["maxiter"] = 100;
    elseif typeof(oedmc_def["maxiter"]) == Array{Int,1}
        oedmc_def["maxiter"] = oedmc_def["maxiter"][1];
    end

    if oedmc_def["maxtime"] == []
        oedmc_def["maxtime"] = Inf;
    elseif typeof(oedmc_def["maxtime"]) == Array{Int,1}
        oedmc_def["maxtime"] = oedmc_def["maxtime"][1];
    end

    if typeof(oedmc_def["Theta"]) == Array{String,1}
        oedmc_def["Theta"] = oedmc_def["Theta"][1];
    end

    if typeof(oedmc_def["finalTime"]) <: Array
        oedmc_def["finalTime"] = oedmc_def["finalTime"][1];
    end

    if typeof(oedmc_def["switchT"]) <: Array
        if typeof(oedmc_def["switchT"][1]) <: Array
            oedmc_def["switchT"] = oedmc_def["switchT"][1];
        end
    end

    if typeof(oedmc_def["tsamps"]) <: Array
        if typeof(oedmc_def["tsamps"][1]) <: Array
            oedmc_def["tsamps"] = oedmc_def["tsamps"][1];
        end
    end

    if typeof(oedmc_def["util"]) == Array{String, 1}
        oedmc_def["util"] = oedmc_def["util"][1];
    end

    if oedmc_def["util"] == []
        oedmc_def["util"] = "perc"
    elseif lowercase(oedmc_def["util"]) == "percentile"
        oedmc_def["util"] = "perc"
    elseif lowercase(oedmc_def["util"]) == "perc"
        oedmc_def["util"] = lowercase(oedmc_def["util"]);
    elseif lowercase(oedmc_def["util"]) == "entropy"
        oedmc_def["util"] = lowercase(oedmc_def["util"]);
    else
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field util! This should be an empty vector or a string containing the word entropy or perc according to the utility you wanna use.")
        return
    end

    # Checks on thetas and load if file path is given

    if (typeof(oedmc_def["Theta"]) == Array{Float64,1}) || (typeof(oedmc_def["Theta"]) == Array{Float32,1})
        if length(oedmc_def["Theta"]) != oedmc_def["Model"]["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Number of parameters introduced does not match the specified")
            return
        end
    elseif (typeof(oedmc_def["Theta"]) == Array{Float64,2}) || (typeof(oedmc_def["Theta"]) == Array{Float32,2})
        if size(oedmc_def["Theta"])[1] != oedmc_def["Model"]["nPar"] && size(oedmc_def["Theta"])[2] != oedmc_def["Model"]["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Number of parameters introduced does not match the specified")
            return
        end
    elseif (typeof(oedmc_def["Theta"]) == String)
        if oedmc_def["Theta"][end-3:end] != ".csv"
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, add the .csv termination to the theta file! And no spaces after!")
            return
        end
        if isfile(oedmc_def["Theta"])
            oedmc_def["Theta"] = Matrix(CSV.read(oedmc_def["Theta"]));
        else
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but the file path you introduced does not exist!")
            return
        end
    end

    # Check that there is more than 1 sample for the parameters
    if length(oedmc_def["Theta"])/oedmc_def["Model"]["nPar"] == 1
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, but you only introduced one sample for the parameters! This method only works with multiple samples. ")
        return
    end

    # Check warning in case the matrix introduced is symmetric
    if length(size(oedmc_def["Theta"])) == 2
        if size(oedmc_def["Theta"])[1] == size(oedmc_def["Theta"])[2]
            println("-------------------------- WARNING --------------------------")
            println("Sorry, but the number of rows and columns of the Theta matrix is the same, so the checks on ")
            println("correct orientation will not work. Please make sure that the dimensions follow: ")
            println("theta[samples, parameters]")
            println("-------------------------------------------------------------")
        end
    end

    # Check that there are multiple samples or single samples for both models (either frequentist or bayesian)
    if length(oedmc_def["Theta"])/oedmc_def["Model"]["nPar"] != 1
        if length(oedmc_def["Theta"])/oedmc_def["Model"]["nPar"] <50
            println("-------------------------- WARNING --------------------------")
            println("Less than 50 samples for theta are given for the Model. ")
            println("Consider using more samples for better results. ")
            println("-------------------------------------------------------------")
        end
    end

    if oedmc_def["tsamps"][end] > oedmc_def["finalTime"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check finalTime. You have selected a sampling point past it.")
        return
    end
    if (oedmc_def["Model"]["Y0eqs"] != []) && oedmc_def["preInd"]==[]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("You specified computation of Y0 as steady-state but you have not introduced an inducer value for it!")
        return
    end

    # Check that no step is no smaller than 2 unit of time
    for j in 1:length(oedmc_def["switchT"])-1
        if (oedmc_def["switchT"][j+1]-oedmc_def["switchT"][j])<=4
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but 2 of the steps in the experiment are too close. This package cannot "))
            println("handle this for now. We are working on it!")
            return
        end
    end

    if typeof(oedmc_def["Obs"]) == Array{String,1}
        if convert(Bool,sum(occursin.("+", oedmc_def["Obs"]))) || convert(Bool,sum(occursin.("-", oedmc_def["Obs"]))) ||
           convert(Bool,sum(occursin.("/", oedmc_def["Obs"]))) ||
           convert(Bool,sum(occursin.("*", oedmc_def["Obs"]))) || convert(Bool,sum(occursin.("^", oedmc_def["Obs"])))
            nothing
        else
            if sum(sum.([occursin.(oedmc_def["Model"]["stName"][k], oedmc_def["Obs"]) for k in 1:length(oedmc_def["Model"]["stName"])])) == 0
            # if sum(occursin.(oedmc_def["Model"]["stName"], oedmc_def["Obs"])) == 0
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but there is some issue with the contents of the field Obs."))
                println("It seems that the observable(s) selected do not match any state for one of the models")
                return
            end
        end
    end

    for j in 1:length(oedmc_def["Obs"])
        s1=[occursin(oedmc_def["Model"]["stName"][i], oedmc_def["Obs"][j]) for i in 1:(oedmc_def["Model"]["nStat"])]
        if sum(s1) == 0
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but you have indicated as observable something that is not a state in the Model."))
            println("The same observable has to appear in both models for the optimisation to work.")
            return
        end
    end

    if  length(oedmc_def["y0"])/oedmc_def["Model"]["nStat"] == 1
        if oedmc_def["Model"]["nStat"] != length(oedmc_def["y0"])
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but the y0 vector for the model does not have the same amount of values as states in the model."))
            println("Remember that even if it will not be used you need to provide a y0 vector.")
            return
        end
    elseif length(oedmc_def["y0"])/oedmc_def["Model"]["nStat"] > 1
        if length(oedmc_def["y0"])/oedmc_def["Model"]["nStat"] != length(oedmc_def["Theta"])/oedmc_def["Model"]["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but the y0 matrix for the model does not have the same amount of entries as parameter samples."))
            println("Remember that even if it will not be used you need to provide a y0 vector.")
            return
        end
        if size(oedmc_def["y0"])[1] == size(oedmc_def["y0"])[2]
            println("-------------------------- WARNING --------------------------")
            println(string("Sorry, but the number of rows and columns of the y0 matrix is the same, so the checks on "))
            println("correct orientation will not work. Please make sure that the dimensions follow: ")
            println("y0[samples, states]")
        end

    end

    if oedmc_def["Model"]["nInp"] != length(oedmc_def["preInd"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but the preInd vector for model 1 does not have the same amount of values as inputs in the model."))
        println("Remember that even if it will not be used you need to provide a preInd_M1 vector.")
        return
    end

    if round(oedmc_def["finalTime"]) == 0
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but it seems that the final time is 0?"))
        println("Remember that for now, the minimum time discretisation allowed is 1 so adjust the parameter scales according to it.")
        return
    end

    if length(oedmc_def["switchT"]) <= 1
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but it seems that you have not introduced a correct switchT"))
        println("Remember, these are the switching times for the inducer, including time 0 and the final time. ")
        return
    end

    # Check that no step is no smaller than 2 unit of time
    for j in 1:length(oedmc_def["tsamps"])-1
        if (round(oedmc_def["tsamps"][j+1])-round(oedmc_def["tsamps"][j]))<1
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but 2 of the sampling times are too close. "))
            println("Remember that for now, the minimum time discretisation allowed is 1 so adjust the parameter scales according to it.")
            println("Sorry for the inconvenience, I am working on it...")
            return
        end
    end

    if oedmc_def["Model"]["nInp"] == 1
        if oedmc_def["fixedInp"] != []
            println("-------------------------- WARNING --------------------------")
            println("You have selected one of the input of the models to be fixed, but there is only 1 inducer in the model.")
            println("This entry will be ignored. If you want to design a single step experiment adjust switchT accordingly.")
            println("-------------------------------------------------------------")
            oedmc_def["fixedInp"] = [];
        end
    else
        if oedmc_def["fixedInp"] != []
            s1=[(oedmc_def["Model"]["inpName"][i] in oedmc_def["fixedInp"]) for i in 1:(oedmc_def["Model"]["nInp"])]

            if sum(s1) == 0
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but you have indicated as fixedInp something that is not an input in the Model."))
                println("The same inputs have to appear in both models for the optimisation to work.")
                return
            end
        end

    end

    tmp1 = [];
    for i in 1:(oedmc_def["Model"]["nInp"])
        if oedmc_def["Model"]["inpName"][i] in oedmc_def["fixedInp"]
            tmp1 = vcat(tmp1, oedmc_def["Model"]["inpName"][i])
        end
    end

    if length(tmp1) >= oedmc_def["Model"]["nInp"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but it seems that you want fixed all the inputs."))
        println("To do this, please adjust switchT instead of trying it through fixedInp")
        return
    end

    if length(oedmc_def["uUpper"]) != length(oedmc_def["uLower"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("You have given a different number of upper and lower bounds for the inducers.")
        println("Please, check this before proceeding.")
        return
    end

    if maximum([oedmc_def["Model"]["nInp"], oedmc_def["Model"]["nInp"]]) != length(oedmc_def["uUpper"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("You have introduced less bounds than inputs the models have or vice-versa")
        println("Please, check this before proceeding.")
        return
    end

    tm1 = oedmc_def["uUpper"].>= oedmc_def["uLower"]
    if sum(tm1) < length(tm1)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Check the input bounds please, it seems that upper and lower bounds are not correct.")
        return
    end

    if oedmc_def["fixedStep"] != []
        fs = [oedmc_def["fixedStep"][j][1] for j in 1:length(oedmc_def["fixedStep"])];
        if minimum(fs) < 1
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have introduced and index for the step smaller than 1. This does not make much sense.")
            return
        elseif maximum(fs) > (length(oedmc_def["switchT"])-1)
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have introduced and index for the step larger than the number of steps of the experiment. This does not make much sense.")
            return
        elseif length(fs) > (length(oedmc_def["switchT"])-1)
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have introduced more fixed steps than steps in the experiment. This does not make much sense.")
            return
        elseif unique(fs) != fs
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have introduced a repeated fixed step for the inputs. This does not make much sense.")
            return
        end

        fi = [oedmc_def["fixedStep"][j][2] for j in 1:length(oedmc_def["fixedStep"])];
        for k in 1:length(fi)
            if length(fi[k]) != oedmc_def["Model"]["nInp"] - length(oedmc_def["fixedInp"])
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Sorry, but it seems that there is a mismatch between the number of inputs and the number of values in the second entry of the tuple for the fixed step entry.")
                println("Remember that any fixed inputs (fixedInp) are not considered in this field, so you should not include them here.")
                println("Tip: If you want a fixed input with a fixed value in the experiment, hard-code its value in the definition of the model.")
                return
            end
        end
    end

    for k in 1:length(oedmc_def["equalStep"])
        if length(oedmc_def["equalStep"][k]) == 1
            println("-------------------------- WARNING !!! --------------------------")
            println("One of the entries for equalStep only has 1 index, please avoid doing this.")
        end
        if minimum(oedmc_def["equalStep"][k]) < 1
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have not introduced and index for a step smaller than 1. ")
            return
        elseif maximum(oedmc_def["equalStep"][k]) > (length(oedmc_def["switchT"])-1)
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have not introduced and index for a step higher than the actual number of steps in the experiment.")
            return
        end
        oedmc_def["equalStep"][k] = sort(oedmc_def["equalStep"][k]);
    end

    if oedmc_def["equalStep"] != []
        unni = [oedmc_def["equalStep"][i][j] for i in 1:length(oedmc_def["equalStep"]) for j in 1:length(oedmc_def["equalStep"])];
        if length(unni) != length(unique(unni))
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have not introduced a repeated index somewhere in the definition of equalStep.")
            return
        end
    end

    if oedmc_def["equalStep"] != [] && oedmc_def["fixedStep"] != []
        fs = [oedmc_def["fixedStep"][j][1] for j in 1:length(oedmc_def["fixedStep"])];
        if typeof(oedmc_def["fixedStep"]) != Array{Tuple{Int,Array{Any,1}},1}
            for k in 1:length(oedmc_def["equalStep"])
                for h in 1:length(fs)
                    if fs[h] in oedmc_def["equalStep"][k] &&  typeof(oedmc_def["fixedStep"][h]) != Tuple{Int,Array{Any,1}}
                        println("-------------------------- Process STOPPED!!! --------------------------")
                        println("Sorry, but it seems that you have defined an equalStep that is also a fixedStep.")
                        return
                    end
                end
            end
        end
    end

    for i in 1:length(oedmc_def["fixedStep"])
        for j in 1:length(oedmc_def["fixedStep"][i][2])
            if typeof(oedmc_def["fixedStep"][i][2][j]) != Int && typeof(oedmc_def["fixedStep"][i][2][j]) != Float64 && typeof(oedmc_def["fixedStep"][i][2][j]) != Float32 &&
                (oedmc_def["fixedStep"][i][2][j]) != Any
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Sorry, but it seems you have introduced something wrong in the second entry of the tuple.")
                println("This has to be a number indicating a concentration or the variable Any if that step is not fixed for that specific inducer. ")
                return
            end
        end
    end

    return(oedmc_def)
end

## Function that will generate the necessary scripts for the optimisation
function genOptimMCFuncts(oedmc_def)

    oedmc_def = checkStructOEDMC(oedmc_def);
    if typeof(oedmc_def["Model"]) <: Dict
        oedmc_def["Model"] = checkStruct(oedmc_def["Model"]);
        oedmc_def["Model"] = GenerateModel(oedmc_def["Model"]);
    elseif typeof(oedmc_def["Model"]) <: Array
        oedmc_def["Model"] = checkStruct(oedmc_def["Model"][1]);
        oedmc_def["Model"] = GenerateModel(oedmc_def["Model"]);
    end;

    # Generate results directory
    cudi = pwd(); #@__DIR__;
    cudi2 = @__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", oedmc_def["Model"]["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", oedmc_def["Model"]["NameF"],"_",today()))
    end
    if !isdir(string(cudi, "\\Results\\", oedmc_def["Model"]["NameF"],"_",today(), "\\OEDModelCalibrationScripts"))
        mkdir(string(cudi, "\\Results\\", oedmc_def["Model"]["NameF"],"_",today(), "\\OEDModelCalibrationScripts"))
    end

    # Packages needed for the simulations
    Head = string("
        ");

    # Manage inputs order
    inpnam = oedmc_def["Model"]["inpName"];

    steps = length(oedmc_def["switchT"])-1;
    induc = oedmc_def["Model"]["nInp"];
    fiii = length(oedmc_def["fixedInp"]);

    doof = (steps*(induc-fiii))+fiii;

    ustr = Array{String,1}(undef, doof); # String vector that will contain the names for the inputs to be optimised at each step. Inducer names are in order (first all variables of on einducer, then all for the next and so on.)
                                        # If an input is fixed to a value, then it will not be considered here
                                        # If an input is fixed it will only be taken into account Once
                                        # If an input concentration is repeated in multiple steps, then it will only be considered once
    us = [oedmc_def["fixedStep"][j][1] for j in 1:length(oedmc_def["fixedStep"])]; # Fixed Steps (the actual step number)
    counn = 1; # Index for the entry of ustr that will advance only if certain conditions are meet (that that position has been already checked and a variable or number has been introduced)
                # Tecnichally this indicated the first entry of ustr for each inducer variables.
    fixin = 0; # Index to account for fixedInp. If the fixedInp is the first one in the inpam list this is the best way (that I have found) to account for it. This is because of the use of i (input index). This is because first I consider the case where there is a fisedInp and then check the rest.
                # Necessary since fixedStep does not consider fixedInps
    # THE FIRST MOUNTAINS!
    for i in 1:length(inpnam) # Loop over each inducer. i will be the index from the list inpnam
        if inpnam[i] in oedmc_def["fixedInp"] # Check if input i is in the list of fixed inputs
            ustr[counn] = inpnam[i]; # If yes, add the name of the input as the only variable to optimise for that input.
            counn += 1; # Need to move index in ustr
            fixin += 1; # Need to take into account that one input is being fixed
        else # Now check the rest of the cases
            for j in 0:steps-1 # j will be the index of the steps -1. If 6 steps, then a value between 0 and 5. Only used to check fixedStep.
                if j+1 in us # Checl if the step is present in the list of fixedStep
                    for u in 1:length(us) # To account for the case that the user might have introduced the index order not sorted, so we fing the j+1 that matches the us entry
                        if j+1 == us[u]
                            if oedmc_def["fixedStep"][u][2][i-fixin] != Any # If user set that as a fixed value do nothing, since this one will be optimised, not fixed. Coded like this to account for future possible modifications (and to not get lost).
                                nothing
                            else
                                ustr[counn+(j)] = string(inpnam[i], j+1); # Add the variable name, which will be the input name from inpnam plus the number of this step (to avoid problems with repetitions in the folowwing section)
                            end

                        end
                    end
                else
                    ustr[counn+(j)] = string(inpnam[i], j+1); # The current step is not in the list of fixedStep, so this is optimised so introduce a name for the variable.
                end
            end
            if oedmc_def["equalStep"] == [] # This section accounts for the contents of equalStep, so if there is nothing, skip.
                nothing
            else
                if typeof(oedmc_def["fixedStep"]) != Array{Tuple{Int,Array{Any,1}},1} && # This section is only used if all entries in fixed step are numbers (so no Any from the user)
                typeof(oedmc_def["fixedStep"]) != Array{Tuple{Int,Array{T,1} where T},1}
                    tmpin = ustr[counn:counn+steps-1]; # Temporary variable that takes all the entries of ustr that correspond to the current inducer i.
                    for k in 1:length(oedmc_def["equalStep"]) # Loop over each equalStep entry and add the same variable name in all the steps in that entry (steps where the inducer is going to have the same value)
                        tmpin[oedmc_def["equalStep"][k]] .= tmpin[oedmc_def["equalStep"][k][1]]
                    end
                    ustr[counn:counn+steps-1] = tmpin; # Define the corresponding ustr filling it with the temporary variable
                else
                    tmpin = ustr[counn:counn+steps-1]; # Temporary variable that takes all the entries of ustr that correspond to the current inducer i.
                    for k in 1:length(oedmc_def["equalStep"]) # Loop over each entry of equalStep to see if a variable name will be assigned or not (account for fixedStep with Any entries)
                        tmpfe = zeros(length(oedmc_def["equalStep"][k])); # Temporary variable that will contain the same entries of equalStep but only if that input for that step is not fixed (can happen due to the option of introducing Any)
                        for j in 1:length(tmpfe) # Loop over each entry of the equalStep entry that we are in now.
                            if oedmc_def["equalStep"][k][j] in us # Check if the entry is in the list of equal step and if for that inducer in this step the user specified Any (input to be optimised there) then add the step index in the temporary variable so a variable name can be assigned later. If not nothing since that value would have been fixed before.
                                for w in 1:length(us)
                                    if oedmc_def["equalStep"][k][j] == us[w] && oedmc_def["fixedStep"][w][2][i-fixin] == Any
                                        tmpfe[j] = oedmc_def["equalStep"][k][j];
                                    end
                                end
                            else
                                tmpfe[j] = oedmc_def["equalStep"][k][j]; # If the entrie is not present in fixedStep, then add the entrie in the temporary variable.
                            end
                        end
                        tmpfe = convert.(Int, filter!(e->e∉[0],tmpfe)); # Eliminate all the entries that are 0 (the ones for steps where the input considered is fixed. )
                        if tmpfe != [] # For all those entries, add the same variable name in the temporary variable tmpin.
                            tmpin[tmpfe] .= tmpin[tmpfe[1]];
                        end
                    end
                    ustr[counn:counn+steps-1] = tmpin; # Define the entries of ustr ussing the filled tmpin
                end
            end
            counn += (steps); # Needed to move to the indexes of the second inducer.
        end
    end

    ustr = [ustr[i] for i in 1:length(ustr) if isassigned(ustr, i)]; # If an entriy has still nothign assigned yet is that it will not be used (for exampel the case of having a fixedInp is not taken into acciunt in the calculation of degrees of freedom doof, so more than one entry will be generated, but only 1 will be used.)
    ustr = unique(ustr); # To avoid any repetitions (just in case)

    inpstr1 = Array{String, 1}(undef, steps*induc) # Inputs for model 1 that will be given for the simulation. This is needed to be ordered grouping the inputs for each step together (for example, if we have 3 steps and the inputs IPTG and aTc, the order has to be IPTG1,aTc1,IPTG2,aTc2,IPTG3,aTc3) due to the model generation scripts structure so this is taken into account.
    inpstr2 = Array{String, 1}(undef, steps*induc) # Same as previous but for model 2 (to account for the case where the models do not contain the same number of inducers)
    r = 1:induc:steps*induc # To ease indexing later. This goes from 1 to the end of the vectors inpstr1 and inpstr2 every total number of inducers (hence, grouping the window by all inducers in each step of the experiment)
    fs = [oedmc_def["fixedStep"][j][2] for j in 1:length(oedmc_def["fixedStep"])]; # Input values (or Any entry) for each step specified in fixedStep
    cnt = 1; # Variables not used in current version. Left for checks to see where the if statemetns go. This is for model 1 the one below for model 2.
    cnt2 = 1;

    # THE SECOND MOUNTAINS!
    for i in 1:steps # length(r) to go through each step
        fixin1 = 0; # To discard the fixed inputs (as done before). Separated for both models, but technically only using one should do.
        fixin2 = 0;
        for j in 0:induc-1 # j will be the index of the inducer considered -1. If 2 inducers, then 0 and 1 (done to ease later indexing in previous versions).
            if inpnam[j+1] in oedmc_def["fixedInp"] # Consider case that the current inducer is fixed. If not check the rest (fixedStep and equalStep)
                fixin1 += 1; # Same reason ad fixin in the FIRSTS MOUNTAINS
                fixin2 += 1;
                if inpnam[j+1] in oedmc_def["Model"]["inpName"] # The input is fixed, so add the same variable name for each step. For model 1
                    inpstr1[r[i]+j] = inpnam[j+1];
                else
                    inpstr1[r[i]+j] = "";
                end
            else
                if inpnam[j+1] in oedmc_def["Model"]["inpName"] # Check if the current inducer is present in Model 1. If so proceed.
                    if typeof(oedmc_def["fixedStep"]) != Array{Tuple{Int,Array{Any,1}},1} && # Consider the case where in fixedStep all the inducers are fixed to a value and there is no Any. Was easier to consider the case of presence or absence of Anys separately.
                    typeof(oedmc_def["fixedStep"]) != Array{Tuple{Int,Array{T,1} where T},1}
                        if i in us # If the current step is present in fixedSteps.
                            for w in 1:length(us) # This loop and inequality is to account for the case where the user does not introduce the entries in fixedStep in order (sorted by step)
                                if i == us[w]
                                    if j+1-fixin1 != 0 # To avoid the case of indexing by 0. This should not happen, but had some issues during tests so left it just in case.
                                        inpstr1[r[i]+j] = string(fs[w][j+1-fixin1]); # Hard code the input value introduced by the user for that inducer and that step.
                                        cnt += 1;
                                    end
                                end
                            end
                        else # This input for this step is not hard-coded, so add the variable name
                            inpstr1[r[i]+j] = string(inpnam[j+1],i);
                            cnt += 1;
                        end
                        if oedmc_def["equalStep"] == []  # This section accounts for the contents of equalStep, so if there is nothing, skip.
                            if !isassigned(inpstr1, r[i]+j) # If current entry is not assigned to anything (from the fixedStep contents), the add the corersponding variable name to the current step for the current input
                                inpstr1[r[i]+j] = string(inpnam[j+1],i);
                                cnt += 1
                            end
                        else
                            for k in 1:length(oedmc_def["equalStep"]) # Check for all the equal steps and add the same variable name for the current input.
                                if i in oedmc_def["equalStep"][k]
                                    inpstr1[r[i]+j] = string(inpnam[j+1],oedmc_def["equalStep"][k][1]);
                                end
                            end
                        end
                    else # This accounts for Any entries in fixedStep
                        if i in us # If the current step is present in fixedSteps.
                            for w in 1:length(us) # This loop and inequality is to account for the case where the user does not introduce the entries in fixedStep in order (sorted by step)
                                if j+1-fixin1 != 0
                                    if i == us[w] && oedmc_def["fixedStep"][w][2][j+1-fixin1] != Any # If this input for the current step is fixed to a value for the user
                                        inpstr1[r[i]+j] = string(fs[w][j+1-fixin1]); # Hard code the value introduced by the user
                                        cnt += 1;
                                    elseif i == us[w] && oedmc_def["fixedStep"][w][2][j+1-fixin1] == Any # If this input for the current step is set to Any, hence to be optimised
                                        inpstr1[r[i]+j] = string(inpnam[j+1],i); # Add the corresponting variable name.
                                        cnt += 1;
                                    end
                                end
                            end
                        else # Entry not in fixedStep, so add the corresponding variable name for the step and input.
                            inpstr1[r[i]+j] = string(inpnam[j+1],i);
                            cnt += 1;
                        end
                        if oedmc_def["equalStep"] == [] # This section accounts for the contents of equalStep, so if there is nothing, skip.
                            if !isassigned(inpstr1, r[i]+j) # If current entry is not assigned to anything (from the fixedStep contents), the add the corersponding variable name to the current step for the current input
                                inpstr1[r[i]+j] = string(inpnam[j+1],i);
                                cnt += 1;
                            end
                        else # Now account for steps that have to have the same variable name (even though it has been coded differently in the previous loop)
                            for k in 1:length(oedmc_def["equalStep"]) # Loop over every entry in equalStep (arrays with the step indexes to have the same value)
                                tmpfe = zeros(length(oedmc_def["equalStep"][k])); # Temporary variable that will have the same values than the current equalStep but without the ones where the input value is fixed to a number by the user
                                for z in 1:length(tmpfe) # Loop over each entry in the current equalStep
                                    if oedmc_def["equalStep"][k][z] in us # If the current value of the equalStep entry (step index) is present in fixedStep
                                        for w in 1:length(us) # To avoid indexing issues if user does not introduce things sorted by step index in fixedStep
                                            if 1+j-fixin1 != 0  # To avoid the case of indexing by 0. This should not happen, but had some issues during tests so left it just in case.
                                                if oedmc_def["equalStep"][k][z] == us[w] && oedmc_def["fixedStep"][w][2][1+j-fixin1] == Any # If the current step from the equalStep entry is also in fixedStep and this is set to Any (so to be optimised)
                                                    tmpfe[z] = oedmc_def["equalStep"][k][z]; # Add the step index in the temporary variable
                                                end
                                            end
                                        end
                                    else # If the step in the equalStep entry is not in fixedStep then it has to be optimsied, so add the entry.
                                        tmpfe[z] = oedmc_def["equalStep"][k][z];
                                    end
                                end
                                tmpfe = convert.(Int, filter!(e->e∉[0],tmpfe)); # Delete zeros (corresponding to fixed values in the inducer)
                                if i in tmpfe
                                    inpstr1[r[i]+j] = string(inpnam[j+1],tmpfe[1]); # For each equalStep entry add the same inducer variable name
                                end
                            end
                        end
                    end
                else
                    inpstr1[r[i]+j] = ""; # If no condition is meet, then add an empty string. This will happen if the current inducer is not present un this model.
                end
            end
        end
    end
    inpstr1 = [inpstr1[i] for i in 1:length(inpstr1) if isassigned(inpstr1, i) && inpstr1[i] != ""]; # Delete unasigned entries (should not happen but just in case) and emtpy strings (some input not present in the model)


    obsM1 = Array{String,1}(undef, length(oedmc_def["Obs"]));
    for i in 1:length(oedmc_def["Obs"])
        oedmc_def["Obs"][i] = replace(oedmc_def["Obs"][i], ".+"=>"+")
        oedmc_def["Obs"][i] = replace(oedmc_def["Obs"][i], ".-"=>"-")
        oedmc_def["Obs"][i] = replace(oedmc_def["Obs"][i], ".*"=>"*")
        oedmc_def["Obs"][i] = replace(oedmc_def["Obs"][i], "./"=>"/")
        oedmc_def["Obs"][i] = replace(oedmc_def["Obs"][i], ".^"=>"^")

        oedmc_def["Obs"][i] = replace(oedmc_def["Obs"][i], "+"=>" .+")
        oedmc_def["Obs"][i] = replace(oedmc_def["Obs"][i], "-"=>" .-")
        oedmc_def["Obs"][i] = replace(oedmc_def["Obs"][i], "*"=>" .*")
        oedmc_def["Obs"][i] = replace(oedmc_def["Obs"][i], "/"=>" ./")
        oedmc_def["Obs"][i] = replace(oedmc_def["Obs"][i], "^"=>" .^")
        for j in 1:oedmc_def["Model"]["nStat"]
            if occursin(oedmc_def["Model"]["stName"][j], oedmc_def["Obs"][i])
                tmpo1 = replace(oedmc_def["Obs"][i], oedmc_def["Model"]["stName"][j] => string("solMC[:,",j,",:]"));
                obsM1[i] = string("    Obs",i,"_MC = ",tmpo1,"; \n")
            end
        end
    end

    # Asses which utility the user has chossen and write the function
    if oedmc_def["util"] == "perc"
        ccomp = string("
    ",join([string("    LowqObs",i,"_MC = zeros(size(solMC)[1]); \n") for i in 1:length(oedmc_def["Obs"])]),"
    ",join([string("    HighqObs",i,"_MC = zeros(size(solMC)[1]); \n") for i in 1:length(oedmc_def["Obs"])]),"
        for i in 1:size(solMC)[1]
        ",join([string("    LowqObs",i,"_MC[i] = percentile(Obs",i,"_MC[i,:],0.5); \n") for i in 1:length(oedmc_def["Obs"])]),"
        ",join([string("    HighqObs",i,"_MC[i] = percentile(Obs",i,"_MC[i,:],99.5); \n") for i in 1:length(oedmc_def["Obs"])]),"
        end

        # Compute Euclidean distances
    ",join([string("    EuObs",i,"_MC = sqrt(sum((LowqObs",i,"_MC-HighqObs",i,"_MC).^2)); \n") for i in 1:length(oedmc_def["Obs"])]),"
        util = (",join([string("EuObs",i,"_MC") for i in 1:length(oedmc_def["Obs"])], "+"),")*(1/",length(oedmc_def["Obs"]),");
    ")

    elseif oedmc_def["util"] == "entropy"

        # If y0 is the same for all runs we need to ignore this time point since that would mean having a variance of 0 which will cause errors in the distribution fits.
        if oedmc_def["Model"]["Y0eqs"] == [] && length(size(oedmc_def["y0"])) == 1
            inn = 2;
        else
            inn = 1;
        end
        ccomp = string("
    ",join([string("    EntObs",i,"_MC = zeros(size(solMC)[1]); \n") for i in 1:length(oedmc_def["Obs"])]),"

        dists = [Beta, Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform];
        dists2 = [Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform];

        names = [","\"","timePoint","\"","];

        for j in ",inn,":size(solMC)[1]
    ",join([string("        fitts",i," = Dict(); \n        bestfit",i," = Dict(); \n        bestfitInd",i," = zeros(1,1); \n") for i in 1:length(oedmc_def["Obs"])]),"
    ",join([string("
            try
                fitts",i,"[names[1]] = fit.(dists, Ref(Obs",i,"_MC[j,:]));
                bestfitInd",i,"[1] = findmax(loglikelihood.(fitts",i,"[names[1]], Ref(Obs",i,"_MC[j,:])))[2];
                bestfit",i,"[names[1]] = fitts",i,"[names[1]][findmax(loglikelihood.(fitts",i,"[names[1]], Ref(Obs",i,"_MC[j,:])))[2]];
            catch
                fitts",i,"[names[1]] = fit.(dists2, Ref(Obs",i,"_MC[j,:]));
                bestfitInd",i,"[1] = findmax(loglikelihood.(fitts",i,"[names[1]], Ref(Obs",i,"_MC[j,:])))[2]+1;
                bestfit",i,"[names[1]] = fitts",i,"[names[1]][findmax(loglikelihood.(fitts",i,"[names[1]], Ref(Obs",i,"_MC[j,:])))[2]];
            end") for i in 1:length(oedmc_def["Obs"])]),"

        ",join([string("    EntObs",i,"_MC[j] = entropy(bestfit",i,"[names[1]]); \n") for i in 1:length(oedmc_def["Obs"])]),"
        end

        HES = zeros(1,",length(oedmc_def["Obs"]),");
    ",join([string("    HES[1,",i,"] = sum(EntObs",i,"_MC); \n") for i in 1:length(oedmc_def["Obs"])]),"
        util = sum(HES)/",length(oedmc_def["Obs"]),"
    ")
    else
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, but there seems to be an issue with the definition of the Utility function. ")
        return
    end;

    Util = string("

    function ", oedmc_def["Model"]["NameF"],"Utility(ins)

        # Definition of the inputs for the ODEs
        if length(ins) == 1
            ",join(ustr, ",")," = ins[1];
        else
            ",join(ustr, ",")," = ins;
        end
        inputsMC = [",replace(join(inpstr1, ","), ",,"=>","),"];

        # Solve ODEs
        solMC = ",string(join(oedmc_def["Model"]["NameF"]),"_SolveAll"),"(tsMC, pD1MC, spMC, inputsMC, ivss1MC, sampsMC, pre1MC);

        # Extracte wanted vectors (observables) with time reduction
    ",join(obsM1),"

    ",join(ccomp),"

        return(util)
    end

        ");

    fungen = join([Head, Util]);

    open(string(cudi, "\\Results\\", oedmc_def["Model"]["NameF"],"_",
        today(), "\\OEDModelCalibrationScripts\\",oedmc_def["Model"]["NameF"],"_OEDMC.jl"), "w") do io
       write(io, fungen);
    end;

    println("")
    println("----------------------------------------- SCRIPTS OED Model Calibration -----------------------------------------")
    println("Utility function script to perform OED for Model Calibration has been generated in the directory: ")
    println(string(cudi, "\\Results\\", oedmc_def["Model"]["NameF"],"_",
            today(), "\\OEDModelCalibrationScripts\\",oedmc_def["Model"]["NameF"],"_OEDMC.jl"))
    println("--------------------------------------------------------------------------------------")
    println("")

    include(string(cudi, "\\Results\\", oedmc_def["Model"]["NameF"],"_",
        today(), "\\OEDModelCalibrationScripts\\",oedmc_def["Model"]["NameF"],"_OEDMC.jl"))

    println("")
    println("---------------------------------------------- OPTIMISATION INFO ----------------------------------------------")
    println("")
    println("If you wish to modify any of the settings for the Bayesian Optimisation used (from the package ")
    println("BayesianOptimization.jl) please search the function settingsBayesOptMC located in the file ModelOEDCalibration.jl ")
    println("that is located in the directory: ")
    println(string("        ", cudi2))
    println("and change any desired setting (under your own risk). ")
    println("If any irreversible change is made by mistake, just look for the file settingsBayesOptMCBackUp.jl and copy the")
    println("contents of the function inside (which are the same as the original except for the function name ;) )")
    println("")
    println("---------------------------------------------------------------------------------------------------------------")
    println("")

    return(oedmc_def)
end

## Function to plot results
function plotOEDMCResults(oedmc_res, oedmc_def)

    # Utility convergence curve
    pu = plot(oedmc_res["ResOptim"]["modelY"], label = "UtilityVal", legend=:topleft, xlabel = "Function Evaluations", ylabel = "Utility")
    pu = plot!(oedmc_res["ConvCurv"][:,2], linetype=:step, label = "BestUtil", title = "OED Model Calibration")
    savefig(string(oedmc_def["savepath"], "\\Plot_OEDMCConvergence_", oedmc_def["flag"], ".png"));

    # Simulations
    # Manage inputs order
    inpnam = oedmc_def["Model"]["inpName"];

    tit = "";
    yl1 = "";
    for k in 1:length(oedmc_def["Obs"])
        tit = hcat(tit, string(oedmc_def["Obs"][k]));
        yl1 = hcat(yl1, "y");
    end
    tit = tit[:,2:end];
    yl1 = yl1[:,2:end];

    titu = "";
    yl2 = "";
    for k in 1:length(inpnam)
        titu = hcat(titu, string(inpnam[k]));
        yl2 = hcat(yl2, "u");
    end
    titu = titu[:,2:end];
    yl2 = yl2[:,2:end];

    tuu = hcat(tit, titu);
    yuu = hcat(yl1, yl2);

    # Elements for the inducers
    tu = Array{Array{Any,1}}(undef, length(oedmc_def["Obs"])+length(inpnam))
    [tu[k] = [] for k in 1:length(oedmc_def["Obs"])]
    su = Array{Array{Any,1}}(undef, length(oedmc_def["Obs"])+length(inpnam))
    [su[k] = [] for k in 1:length(oedmc_def["Obs"])]
    for k in 1:length(inpnam)
        if oedmc_def["fixedInp"] == []
            tu[k+length(oedmc_def["Obs"])] = round.(oedmc_def["switchT"]);
        else
            if inpnam[k] in oedmc_def["fixedInp"]
                tu[k+length(oedmc_def["Obs"])] = [round.(oedmc_def["switchT"])[1], round.(oedmc_def["switchT"])[end]];
            else
                tu[k+length(oedmc_def["Obs"])] = round.(oedmc_def["switchT"]);
            end
        end
        su[k+length(oedmc_def["Obs"])] = vcat(oedmc_res["uInpOpt"][inpnam[k]], oedmc_res["uInpOpt"][inpnam[k]][end])
    end

    pl = plot(round.(oedmc_def["tsamps"]), oedmc_res["SimulObs_MC"][:,:,1],
            layout=length(oedmc_def["Obs"])+length(inpnam), label = "", title = tit, xlabel = "time", ylabel = yuu,
                    color = "gray", size = [2000, 1200]);

    for samp in 2:size(oedmc_res["SimulObs_MC"])[3]
        pl = plot!(round.(oedmc_def["tsamps"]), oedmc_res["SimulObs_MC"][:,:,samp],
            layout=length(oedmc_def["Obs"])+length(inpnam), label = "",color = "gray", size = [2000, 1200]);
    end

    pl = plot!(tu, su, layout=length(oedmc_def["Obs"])+length(inpnam), label = "", linetype=:step, title = tuu);

    savefig(string(oedmc_def["savepath"], "\\PlotOEDMCResults_Exp1_", oedmc_def["flag"], ".png"));
end

## Main function to run Bayesian OPTIMISATION
function settingsBayesOptMC(oedmc_def)

    # Manage inputs order
    inpnam = oedmc_def["Model"]["inpName"];

    # Optimisation section
    steps = length(oedmc_def["switchT"])-1;
    induc = oedmc_def["Model"]["nInp"];
    fiii = length(oedmc_def["fixedInp"]);

    doof = (steps*(induc-fiii))+fiii;

    # THE THIRD MOUNTAINS!
    uppe = Array{Any,1}(undef, doof); # Vector of upper bounds. Length set with doof but deppending on the user settings from fixedStep and equalStep this can change.
    lowe = Array{Any,1}(undef, doof); # Vector of lower bounds. Same ^
    coun = 1; # Index for the entry of ustr that will advance only if certain conditions are meet (that that position has been already checked and a variable or number has been introduced)
                # Tecnichally this indicated the first entry of uppe and lowwe for each inducer variables.
    fs = [oedmc_def["fixedStep"][j][1] for j in 1:length(oedmc_def["fixedStep"])]; # Fixed Steps (the actual step number)
    us = [oedmc_def["fixedStep"][j][2] for j in 1:length(oedmc_def["fixedStep"])]; # Input values (or Any entry) for each step specified in fixedStep
    fixin = 0; # Index to account for fixedInp. If the fixedInp is the first one in the inpam list this is the best way (that I have found) to account for it. This is because of the use of i (input index). This is because first I consider the case where there is a fisedInp and then check the rest.

    for i in 1:length(inpnam) # Loop over the number of inducers
        if inpnam[i] in oedmc_def["fixedInp"] # If the current inducer is fixed then add only one bound and move index (coun). Fixin works as described in the other MOUNTAINS
            uppe[coun] = oedmc_def["uUpper"][i];
            lowe[coun] = oedmc_def["uLower"][i];
            coun += 1;
            fixin += 1;
        else
            try # Using coun you take all the indexes corresponding to the current inducer. The you fill it with the specified uper and lower bounds. Added the try just in case someoen runns it in a machien that does not support Float 64? (I am not sure)
                uppe[coun:coun+(steps-1)] = repeat([convert(Float64, oedmc_def["uUpper"][i])], outer = steps);
                lowe[coun:coun+(steps-1)] = repeat([convert(Float64, oedmc_def["uLower"][i])], outer = steps);
            catch
                uppe[coun:coun+(steps-1)] = repeat([convert(Float32, oedmc_def["uUpper"][i])], outer = steps);
                lowe[coun:coun+(steps-1)] = repeat([convert(Float32, oedmc_def["uLower"][i])], outer = steps);
            end

            for j in 1:length(uppe[coun:coun+(steps-1)]) # Now we loop over each entry that we filled to check if we need to delete it according to the users fixedStep and equalStep
                if j in fs # Since this is not a fixed input j also indicates the step index, so check if this is in the list of fixedStep
                    for w in 1:length(fs) # As commentetd in other loops, this for and if is to consider the case the user does not introduce fixedStep entries sorted by step.
                        if j == fs[w]
                            if us[w][i-fixin] != Any # Check if current input contained in a fixedStep is defined as Any (to be optimised) or not (has been fixed. )
                                uppe[coun-1+j] = "";# Since it was fixed (different from Any) ad an empty string there to remove later in both bounds.
                                lowe[coun-1+j] = "";
                            end
                        end
                    end
                end
                if oedmc_def["equalStep"] != [] # If equalStep is empty then we do not need to do anything else, if it has something then there are entries to remove from the bounds
                    for k in 1:length(oedmc_def["equalStep"]) # Loop over each entry of equalStep
                        tmpk = zeros(length(oedmc_def["equalStep"][k])); # Temporal variable that will contain the same values as the equalStep entry except the ones that have been fixed for that step (this will happen if the user set an Any for the current inducer)
                        for h in 1:length(oedmc_def["equalStep"][k]) # Loop over each value in the current entry
                            if oedmc_def["equalStep"][k][h] in fs # Check if the current step considered is also in fixedStep
                                for w in 1:length(fs) # Account for ht euser not introducng steps sorted
                                    if oedmc_def["equalStep"][k][h] == fs[w]
                                        if us[w][i-fixin] == Any # If the current step considered was defined as Any for the current inducer add the index, if it was fixed to a value ignore (this has been taken into account before)
                                            tmpk[h] = oedmc_def["equalStep"][k][h];
                                        end
                                    end
                                end
                            else
                                tmpk[h] = oedmc_def["equalStep"][k][h] # If the current step was not in fixedStep just add it
                            end
                        end
                        tmpk = convert.(Int, filter!(e->e∉[0],tmpk)); # Remove 0 entries (fixed values)
                        if j in tmpk && j != tmpk[1] # Now, check if the current step is in our temporary index. If it is the first entry in the vector then do nothing, if not remove (it will be the same as the first, so no bound for it)
                            uppe[coun-1+j] = "";
                            lowe[coun-1+j] = "";
                        end
                    end
                end
            end
            coun += (steps); # Add the number of steps to the varable coun so we can move to the next inducer
        end
    end

    uppe = [uppe[i] for i in 1:length(uppe) if uppe[i] != ""]; # Remove empty strings so we are left with a vector of numbers with the bounds for the inducers.
    lowe = [lowe[i] for i in 1:length(lowe) if lowe[i] != ""];


    #################################### OPTIMISER SETINGS THAT CAN BE MODIFIED START HERE ####################################

    model = ElasticGPE(convert(Int, length(uppe)),                            #
           mean = MeanConst(0.),
           kernel = SEArd(zeros(length(uppe)), 5.),
           logNoise = 0.,
           capacity = 3000);

    modeloptimizer = MAPGPOptimizer(every = 50, noisebounds = [-4, 3],       # bounds of the logNoise
                   kernbounds = [[-1*ones(length(uppe)); 0], [4*ones(length(uppe)); 10]],  # bounds of the 3 parameters GaussianProcesses.get_param_names(model.kernel)
                   maxeval = 40);

    sut = Symbol(string(oedmc_def["Model"]["NameF"],"Utility"));
    opt = BOpt((@eval $sut),
          model,
          ExpectedImprovement(),                   # type of acquisition
          modeloptimizer,
          lowe,
          uppe,                                     # lowerbounds, upperbounds
          repetitions = 1,                          # evaluate the function for each input 1 times
          maxiterations = oedmc_def["maxiter"],     # evaluate at 50 input positions
          maxduration = oedmc_def["maxtime"],
          sense = Max,                              # maximise the function
          verbosity = Progress);

    #################################### OPTIMISER SETINGS THAT CAN BE MODIFIED STOP HERE ####################################

    return opt

end

## Main function that runs the optimisation
function mainOEDMC(oedmc_def)

    oedmc_def = checkStructOEDMC(oedmc_def);

    # Generate results directory
    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", oedmc_def["Model"]["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", oedmc_def["Model"]["NameF"],"_",today()))
    end

    try
        include(string(cudi, "\\Results\\", oedmc_def["Model"]["NameF"],"_",
                today(), "\\OEDModelCalibrationScripts\\",oedmc_def["Model"]["NameF"],"_OEDMC.jl"))
    catch
        nothing
    end


    include(string(oedmc_def["Model"]["modelpath"]));

    global tsMC = collect(0.0:round(oedmc_def["finalTime"]))';
    global pD1MC = oedmc_def["Theta"];
    global spMC = [convert(Int,v) for v in (round.(oedmc_def["switchT"])')];
    global ivss1MC = oedmc_def["y0"];
    global sampsMC = oedmc_def["tsamps"];
    global pre1MC = oedmc_def["preInd"];

    # Manage inputs order
    inpnam = oedmc_def["Model"]["inpName"];

    # Optimisation section
    steps = length(spMC)-1;
    induc = oedmc_def["Model"]["nInp"];
    fiii = length(oedmc_def["fixedInp"]);

    doof = (steps*(induc-fiii))+fiii;

    opt = settingsBayesOptMC(oedmc_def);

    result = boptimize!(opt);

    # results

    # Extract results
    oedmc_res = Dict()
    oedmc_res["BestResOptim"] = result;
    oedmc_res["ResOptim"] = Dict("acquisitionoptions"=>opt.acquisitionoptions,
                    "duration"=>opt.duration,
                    "model_optimizer"=>opt.model_optimizer,
                    "model_optimum"=>opt.model_optimum,
                    "modeloptimizer"=>opt.modeloptimizer,
                    "observed_optimizer"=>opt.observed_optimizer,
                    "observed_optimum"=>opt.observed_optimum,
                    "repetitions"=>opt.repetitions,
                    "sense"=>opt.sense,
                    "upperbounds"=>opt.upperbounds,
                    "modelX"=>Array(opt.model.x),
                    "modelY"=>Array(opt.model.y));
    oedmc_res["BestUtil"] = maximum(opt.model.y);

    convc = zeros(length(opt.model.y));
    convc[1] = opt.model.y[1];
    for i in 2:length(opt.model.y)
        if opt.model.y[i] > convc[i-1];
            convc[i] = opt.model.y[i];
        else
            convc[i] = convc[i-1];
        end
    end
    oedmc_res["ConvCurv"] = hcat(collect(1:length(convc)), convc);

    # THE FOURTH MOUNTAINS!
    inducers = Dict(); # Library that will contain an entry for each of the inducers with its values ordered (including fixed values and repetitions)
    cnt2 = 1; # Index for the optimisation result (we move it by one once we have already extracted the current value for all the necessary inducers entries)
    fs = [oedmc_def["fixedStep"][j][1] for j in 1:length(oedmc_def["fixedStep"])]; # Fixed Steps (the actual step number)
    us = [oedmc_def["fixedStep"][j][2] for j in 1:length(oedmc_def["fixedStep"])]; # Input values (or Any entry) for each step specified in fixedStep
    fixin = 0; # Same use as in the other MOUNTAINS
    for i in 1:induc # Loop over each inducer
        if inpnam[i] in oedmc_def["fixedInp"] # Check if the current inducer was fixed to a single value to be optimised
            inducers[inpnam[i]] = [result.observed_optimizer[cnt2]]; # Take the value
            fixin += 1; # Add this so there will not be bad indexing in fixedStep
            cnt2 += 1; # Move the index of the result vector to one since we already extracted this
        else
            allstp = Array{Any,1}(undef, steps); # Temporary variable that will have the inducer values for each step of the experiment.
            for j in 1:steps # Loop over the number of steps of the experiment
                if j in fs # Check if current step is in the list of fixedStep
                    for w in 1:length(fs) # Consider entries not sorted by step
                        if j==fs[w]
                            if oedmc_def["fixedStep"][w][2][i-fixin] != Any # If inducer is in the fixedStep but set as Any (optimise), nothing. If it is fixed get the fixed value into the temporary variable allstp.
                                allstp[j] = oedmc_def["fixedStep"][w][2][i-fixin];
                            end
                        end
                    end
                end
                if oedmc_def["equalStep"] == [] # Case where equalStep is empty (no equal steps)
                    if typeof(oedmc_def["fixedStep"]) != Array{Tuple{Int,Array{Any,1}},1} && # Case where there is no Any in the fixedStep case. Also check if the current step is nto assigned to a value. This can happen if there is also nothing in fixedStep. In this case just fill the allstp with the current value and move index.
                        typeof(oedmc_def["fixedStep"]) != Array{Tuple{Int,Array{T,1} where T},1} &&
                        !isassigned(allstp, j)
                        allstp[j] = result.observed_optimizer[cnt2];
                        cnt2 += 1;
                    end
                else
                    for k in 1:length(oedmc_def["equalStep"]) # Loop over each entry in equalStep
                        tmpk = zeros(length(oedmc_def["equalStep"][k])); # Temporary variable to have same values as equalStep entry without the fixed steps
                        for h in 1:length(oedmc_def["equalStep"][k]) # Loop over each step in the entry
                            if oedmc_def["equalStep"][k][h] in fs # Check if current step is in fixedStep
                                for w in 1:length(fs) # Consider no sorted index in fixedStep
                                    if oedmc_def["equalStep"][k][h] == fs[w]
                                        if us[w][i-fixin] == Any # If current inducer was not fixed (Any) then consider the step, if not no.
                                            tmpk[h] = oedmc_def["equalStep"][k][h];
                                        end
                                    end
                                end
                            else
                                tmpk[h] = oedmc_def["equalStep"][k][h] # Step not in fixedStep, so just consider it.
                            end
                        end
                        tmpk = convert.(Int, filter!(e->e∉[0],tmpk)); # Remove zeros (fixed inputs)
                        if j in tmpk # Check if the current step is in the temporary list
                            if j == tmpk[1] # If it is the first entry then modify all the entries in allstep to the curren result value. (only done if fits the first step so it is only done once)
                                allstp[tmpk] .= result.observed_optimizer[cnt2];
                                cnt2 += 1;
                            end
                        end
                    end
                end
                if !isassigned(allstp, j) # If there is still a no assigned value (chouldn't happen) then fill it with the current result value
                    allstp[j] = result.observed_optimizer[cnt2];
                    cnt2 += 1;
                end
            end
            inducers[inpnam[i]] = allstp; # Add the temporary variable into the inducers list.
        end
    end

    oedmc_res["uInpOpt"] = inducers;

    OEDsFunM1 = Symbol(string(join(oedmc_def["Model"]["NameF"]),"_SolveAll"));

    global inpuMC = zeros(oedmc_def["Model"]["nInp"]*(length(oedmc_def["switchT"])-1));

    r1 = 1:(oedmc_def["Model"]["nInp"]):length(inpuMC);

    for i in 1:length(inpnam)
        if inpnam[i] in oedmc_def["Model"]["inpName"]
            if inpnam[i] in oedmc_def["fixedInp"]
                tmpin = oedmc_res["uInpOpt"][inpnam[i]];
                inpuMC[r1.+(i-1)] .= tmpin;
            else
                tmpin = oedmc_res["uInpOpt"][inpnam[i]];
                inpuMC[r1.+(i-1)] = tmpin;
            end
        end
    end

    simulMC = @eval $OEDsFunM1(tsMC, pD1MC, spMC, inpuMC, ivss1MC, sampsMC, pre1MC);

    SimObsMC = selectObsSim_te(simulMC, oedmc_def["Obs"],oedmc_def["Model"]["stName"]);

    oedmc_res["Simul_MC"] = simulMC;
    oedmc_res["SimulObs_MC"] = SimObsMC;

    oedmc_def["savepath"] = string(cudi, "\\Results\\", oedmc_def["Model"]["NameF"],"_",today());
    oedmc_def["savename"] = string("OEDModelCalibrationResults_",oedmc_def["flag"],".jld");

    JLD.save(string(oedmc_def["savepath"], "\\", oedmc_def["savename"]), "oedmc_res", oedmc_res, "oedmc_def", oedmc_def);

    println("")
    println("----------------------------------------- RESULTS -----------------------------------------")
    println("OED for Model Calibration results are saved in the directory: ")
    println(string("                 ", oedmc_def["savepath"]))
    println(string("Under the name ",oedmc_def["savename"]))
    println("--------------------------------------------------------------------------------------")
    println("")

    if oedmc_def["plot"] == true
        plotOEDMCResults(oedmc_res, oedmc_def)
        println("")
        println("----------------------------------------- PLOTS -----------------------------------------")
        println("PLOTS are saved in the directory: ")
        println(string("                 ", oedmc_def["savepath"]))
        println(string("Under the names PlotOEDMCResults_Exp1_",oedmc_def["flag"],".png", " and Plot_OEDMCConvergence_", oedmc_def["flag"], ".png"))
        println("--------------------------------------------------------------------------------------")
        println("")
    end

    return(oedmc_res, oedmc_def)
end
