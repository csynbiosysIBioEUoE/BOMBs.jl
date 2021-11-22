## Define main Structure
function defODEModelSelectStruct()

    # The structure of things user should give
    oedms_def = Dict()

    oedms_def["Model_1"] = []; # Model structure for Model 1. Dict
    oedms_def["Model_2"] = []; # Model structure for Model 2. Dict
    oedms_def["Obs"] = []; # Observables of the models. Array of strings.
    oedms_def["Theta_M1"] = []; # Theta vector (frequentist OED or model with 1 parameter) or matrix (Bayesian OED) for model 1
    oedms_def["Theta_M2"] = []; # Theta vector (frequentist OED or model with 1 parameter) or matrix (Bayesian OED) for model 2

    oedms_def["y0_M1"] = []; # Array of Y0s for the simulations for each experiment. If you are computing the steady state this vector might not be used, however you still need to introduce it
    oedms_def["y0_M2"] = [];
    oedms_def["preInd_M1"] = []; # Level of the inducer in the ON. This entry can be empty if no Steady State equations are given and only the Y0 given by the user is considered.
    oedms_def["preInd_M2"] = [];
    oedms_def["finalTime"] = []; # Vector of final time for the experiment (initial time will allways be asumed as 0, so please consider that)
    oedms_def["switchT"] = []; # Array with the switching times of the inducer in the simulation (time 0 and final time need to be considered)
    oedms_def["tsamps"] = []; # Vector of Sampling times
    oedms_def["fixedInp"] = []; # If more than 1 inducer for the system exists, but only 1 input can be dynamic, give a vector of strings indicating which inputs are going to be optimised but as a constant input instead of dynamic. If none, just give an empty vector.
    oedms_def["fixedStep"] = []; # If you want any of the steps to be fixed to a value. This has to be an empty array if none is fixed or an array of tuples, where each tuple is a step to be fixed. The first entry of the tuple is the index of the step (as an Integer), and the second and array of values for each inducer. Note that the fixed inputs will be ignored, so do not take them into account here.
    oedms_def["equalStep"] = []; # If you want a series of steps to have the same optimised value (for example if you want to design a pulse experiment) you can introduce inside this array different arrays with the indexes of the steps that will have the same value. The values introduced in each array need to be integers.

    oedms_def["plot"] = []; # Bollean or yes/no string to save the resulting optimisation plots in the results directory
    oedms_def["flag"] = [];

    # The order of uUpper and uLower will be taken as the order of inducers defined in Model 1. If Model 2 has an aditional input, this will be added after, so consider this when introducing the bounds for them.
    oedms_def["uUpper"] = []; # Vector indicating the upper bounds for the inducers
    oedms_def["uLower"] = []; # Vector indicating the lower bounds for the inducers
    oedms_def["maxiter"] = []; # Maximum number of iterations for the Bayesian Optimisation. If nothing is introduced a default of 100 iterations will be taken

    return(oedms_def)

end

## Function to check contents of Structure

function checkStructOEDMS(oedms_def)

    # Check that all the fields (and nothing more) are present
    entries = ["Model_1", "Model_2", "Obs", "Theta_M1", "Theta_M2", "y0_M1", "y0_M2", "preInd_M1", "preInd_M2",
                "finalTime", "switchT", "tsamps", "equalStep",
                "fixedInp", "fixedStep", "plot", "flag", "uUpper", "uLower", "maxiter"];
    if symdiff(entries,keys(oedms_def))!=[] && symdiff(entries,keys(oedms_def)) != ["savepath", "savename"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is something wrong...")
        println(symdiff(entries,keys(oedms_def)))
        return
    end

    # Check one to be sure there is no empty entry
    nee = ["Model_1", "Model_2", "Obs", "Theta_M1", "Theta_M2", "y0_M1", "y0_M2", "preInd_M1", "preInd_M2",
            "uUpper", "uLower"]; # No empty entries
    for i in 1:length(nee)
        if oedms_def[nee[i]] == []
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Plese, check the field ", nee[i], ", this cannot be empty"))
            return
        end
    end

    # Check Models and generate the script. Also check that they do not have the same name
    if typeof(oedms_def["Model_1"]) <: Dict
        oedms_def["Model_1"] = checkStruct(oedms_def["Model_1"]);
    elseif typeof(oedms_def["Model_1"]) <: Array
        oedms_def["Model_1"] = checkStruct(oedms_def["Model_1"][1]);
    end;

    if typeof(oedms_def["Model_2"]) <: Dict
        oedms_def["Model_2"] = checkStruct(oedms_def["Model_2"]);
    elseif typeof(oedms_def["Model_2"]) <: Array
        oedms_def["Model_2"] = checkStruct(oedms_def["Model_2"][1]);
    end;

    if oedms_def["Model_1"]["NameF"] == oedms_def["Model_2"]["NameF"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Plese, give different names to the 2 models in the field NameF, otherwise only 1 model would be generated"))
        return
    end;

    if (typeof(oedms_def["Obs"]) != Array{String,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Obs! This should be a vector of strings or integers. ")
        return
    elseif (typeof(oedms_def["Theta_M1"]) != Array{Float64,1}) && (typeof(oedms_def["Theta_M1"]) != Array{Float32,1}) &&
            (typeof(oedms_def["Theta_M1"]) != Array{Float64,2}) && (typeof(oedms_def["Theta_M1"]) != Array{Float32,2}) &&
            (typeof(oedms_def["Theta_M1"]) != Array{String,1}) && (typeof(oedms_def["Theta_M1"]) != String)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Theta_M1! This should be a vector, matrix or string. ")
        return
    elseif (typeof(oedms_def["Theta_M2"]) != Array{Float64,1}) && (typeof(oedms_def["Theta_M2"]) != Array{Float32,1}) &&
            (typeof(oedms_def["Theta_M2"]) != Array{Float64,2}) && (typeof(oedms_def["Theta_M2"]) != Array{Float32,2}) &&
            (typeof(oedms_def["Theta_M2"]) != Array{String,1}) && (typeof(oedms_def["Theta_M2"]) != String)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Theta_M2! This should be a vector, matrix or string. ")
        return
    elseif (typeof(oedms_def["y0_M1"]) != Array{Int,1}) && (typeof(oedms_def["y0_M1"]) != Array{Float64,1}) &&
        (typeof(oedms_def["y0_M1"]) != Array{Float32,1}) &&
        (typeof(oedms_def["y0_M1"]) != Array{Int,2}) && (typeof(oedms_def["y0_M1"]) != Array{Float64,2}) &&
            (typeof(oedms_def["y0_M1"]) != Array{Float32,2})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field y0_M1! This should be an array of numbers. ")
        return
    elseif (typeof(oedms_def["y0_M2"]) != Array{Int,1}) && (typeof(oedms_def["y0_M2"]) != Array{Float64,1}) &&
        (typeof(oedms_def["y0_M2"]) != Array{Float32,1}) &&
        (typeof(oedms_def["y0_M2"]) != Array{Int,2}) && (typeof(oedms_def["y0_M2"]) != Array{Float64,2}) &&
            (typeof(oedms_def["y0_M2"]) != Array{Float32,2})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field y0_M2! This should be an array of numbers. ")
        return
    elseif (typeof(oedms_def["preInd_M1"]) != Array{Int,1}) && (typeof(oedms_def["preInd_M1"]) != Array{Float64,1}) &&
        (typeof(oedms_def["preInd_M1"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field preInd_M1! This should be an array of numbers. ")
        return
    elseif (typeof(oedms_def["preInd_M2"]) != Array{Int,1}) && (typeof(oedms_def["preInd_M2"]) != Array{Float64,1}) &&
        (typeof(oedms_def["preInd_M2"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field preInd_M2! This should be an array of numbers. ")
        return
    elseif (typeof(oedms_def["finalTime"]) != Array{Int,1}) && (typeof(oedms_def["finalTime"]) != Array{Float64,1}) &&
        (typeof(oedms_def["finalTime"]) != Array{Float32,1}) &&
        typeof(oedms_def["finalTime"]) !=Int && typeof(oedms_def["finalTime"]) !=Float32 &&
        typeof(oedms_def["finalTime"]) !=Float64
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field finalTime! This should be a vector or a single number")
        return

    elseif (typeof(oedms_def["switchT"]) != Array{Array{Int,1},1}) &&
        (typeof(oedms_def["switchT"]) != Array{Array{Float64,1},1}) &&
        (typeof(oedms_def["switchT"]) != Array{Array{Float32,1},1}) &&
        (typeof(oedms_def["switchT"]) != Array{Int,1}) &&
        (typeof(oedms_def["switchT"]) != Array{Float64,1}) &&
        (typeof(oedms_def["switchT"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field switchT! This should be an array or an array of arrays of numbers")
        return
    elseif (typeof(oedms_def["tsamps"]) != Array{Array{Int,1},1}) &&
        (typeof(oedms_def["tsamps"]) != Array{Array{Float64,1},1}) &&
        (typeof(oedms_def["tsamps"]) != Array{Array{Float32,1},1}) &&
        (typeof(oedms_def["tsamps"]) != Array{Int,1}) &&
        (typeof(oedms_def["tsamps"]) != Array{Float64,1}) &&
        (typeof(oedms_def["tsamps"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field tsamps! This should be an array or an array of arrays of numbers")
        return
    elseif (typeof(oedms_def["fixedInp"]) != Array{String,1})&& oedms_def["fixedInp"] != []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field fixedInp! This should be a vector of strings or an empty vector. ")
        return
    elseif oedms_def["fixedStep"] != [] && typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{Int,1}},1} && typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{Float64,1}},1} &&
        typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{Float32,1}},1} && typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{Any,1}},1} &&
        typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{T,1} where T},1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field fixedStep! This should have the structure of Array{Tuple{Int,Array{Number,1}},1}. ")
        return
    elseif oedms_def["equalStep"] != [] && typeof(oedms_def["equalStep"]) != Array{Array{Int,1},1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field equalStep! This should have the structure of Array{Array{Int,1},1}. ")
        println("Remember, it has to be an array where as each entry you have another array(s) with the indexes")
        println("of the steps that will be considered as the same step in the optimisation value wise.")
        return
    elseif (typeof(oedms_def["flag"]) != Array{String,1}) && (typeof(oedms_def["flag"]) != String) && ((oedms_def["flag"]) != [])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field flag! This should be a vector, matrix or string. ")
        return
    elseif (typeof(oedms_def["plot"]) != Bool) && (typeof(oedms_def["plot"]) != Array{Bool,1}) &&
            (typeof(oedms_def["plot"]) != String) && (typeof(oedms_def["plot"]) != Array{String,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field plot! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
        return
    elseif (typeof(oedms_def["uUpper"]) != Array{Int,1}) && (typeof(oedms_def["uUpper"]) != Array{Float64,1}) &&
            (typeof(oedms_def["uUpper"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field uUpper! This should be a vector of numbers")
        return
    elseif (typeof(oedms_def["uLower"]) != Array{Int,1}) && (typeof(oedms_def["uLower"]) != Array{Float64,1}) &&
            (typeof(oedms_def["uLower"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field uUpper! This should be a vector of numbers")
        return
    elseif oedms_def["maxiter"] != [] && typeof(oedms_def["maxiter"]) != Array{Int,1} != typeof(oedms_def["maxiter"]) != Int
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field maxiter! This should be an integer or an empty vector")
        return
    end

    # Extract necessary elements to ease generalisation
    if oedms_def["plot"]==[true] || oedms_def["plot"]==["Yes"] || oedms_def["plot"]==["yes"] || oedms_def["plot"]=="Yes" || oedms_def["plot"]=="yes"
        oedms_def["plot"]=true
    elseif oedms_def["plot"]==[false] || oedms_def["plot"]==["No"] || oedms_def["plot"]==["no"] || oedms_def["plot"]=="No" || oedms_def["plot"]=="no" || oedms_def["plot"]==[]
        oedms_def["plot"]=false
    end

    if typeof(oedms_def["flag"]) == Array{String,1}
        oedms_def["flag"] = oedms_def["flag"][1];
    elseif typeof(oedms_def["flag"]) == String
        oedms_def["flag"] = oedms_def["flag"];
    elseif ((oedms_def["flag"]) == [])
        oedms_def["flag"] = "";
    end

    if oedms_def["maxiter"] == []
        oedms_def["maxiter"] = 100;
    elseif typeof(oedms_def["maxiter"]) == Array{Int,1}
        oedms_def["maxiter"] = oedms_def["maxiter"][1];
    end

    if typeof(oedms_def["Theta_M1"]) == Array{String,1}
        oedms_def["Theta_M1"] = oedms_def["Theta_M1"][1];
    end

    if typeof(oedms_def["Theta_M2"]) == Array{String,1}
        oedms_def["Theta_M2"] = oedms_def["Theta_M2"][1];
    end

    if typeof(oedms_def["finalTime"]) <: Array
        oedms_def["finalTime"] = oedms_def["finalTime"][1];
    end

    if typeof(oedms_def["switchT"]) <: Array
        if typeof(oedms_def["switchT"][1]) <: Array
            oedms_def["switchT"] = oedms_def["switchT"][1];
        end
    end

    if typeof(oedms_def["tsamps"]) <: Array
        if typeof(oedms_def["tsamps"][1]) <: Array
            oedms_def["tsamps"] = oedms_def["tsamps"][1];
        end
    end

        # Check that all the contents make sense

    # Checks on thetas and load if file path is given

    if (typeof(oedms_def["Theta_M1"]) == Array{Float64,1}) || (typeof(oedms_def["Theta_M1"]) == Array{Float32,1})
        if length(oedms_def["Theta_M1"]) != oedms_def["Model_1"]["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Number of parameters introduced does not match the specified")
            return
        end
    elseif (typeof(oedms_def["Theta_M1"]) == Array{Float64,2}) || (typeof(oedms_def["Theta_M1"]) == Array{Float32,2})
        if size(oedms_def["Theta_M1"])[1] != oedms_def["Model_1"]["nPar"] && size(oedms_def["Theta_M1"])[2] != oedms_def["Model_1"]["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Number of parameters introduced does not match the specified")
            return
        end
    elseif (typeof(oedms_def["Theta_M1"]) == String)
        if oedms_def["Theta_M1"][end-3:end] != ".csv"
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, add the .csv termination to the theta file! And no spaces after!")
            return
        end
        if isfile(oedms_def["Theta_M1"])
            oedms_def["Theta_M1"] = Matrix(CSV.read(oedms_def["Theta_M1"]));
        else
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but the file path you introduced does not exist!")
            return
        end
    end

    if (typeof(oedms_def["Theta_M2"]) == Array{Float64,1}) || (typeof(oedms_def["Theta_M2"]) == Array{Float32,1})
        if length(oedms_def["Theta_M2"]) != oedms_def["Model_2"]["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Number of parameters introduced does not match the specified")
            return
        end
    elseif (typeof(oedms_def["Theta_M2"]) == Array{Float64,2}) || (typeof(oedms_def["Theta_M2"]) == Array{Float32,2})
        if size(oedms_def["Theta_M2"])[1] != oedms_def["Model_2"]["nPar"] && size(oedms_def["Theta_M2"])[2] != oedms_def["Model_2"]["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Number of parameters introduced does not match the specified")
            return
        end
    elseif (typeof(oedms_def["Theta_M2"]) == String)
        if oedms_def["Theta_M2"][end-3:end] != ".csv"
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, add the .csv termination to the theta file! And no spaces after!")
            return
        end
        if isfile(oedms_def["Theta_M2"])
            oedms_def["Theta_M2"] = Matrix(CSV.read(oedms_def["Theta_M2"]));
        else
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but the file path you introduced does not exist!")
            return
        end
    end

    # Check warning in case the matrix introduced is symmetric
    if length(size(oedms_def["Theta_M1"])) == 2
        if size(oedms_def["Theta_M1"])[1] == size(oedms_def["Theta_M1"])[2]
            println("-------------------------- WARNING --------------------------")
            println("Sorry, but the number of rows and columns of the Theta_M1 matrix is the same, so the checks on ")
            println("correct orientation will not work. Please make sure that the dimensions follow: ")
            println("theta[samples, parameters]")
            println("-------------------------------------------------------------")
        end
    end

    # Check warning in case the matrix introduced is symmetric
    if length(size(oedms_def["Theta_M2"])) == 2
        if size(oedms_def["Theta_M2"])[1] == size(oedms_def["Theta_M2"])[2]
            println("-------------------------- WARNING --------------------------")
            println("Sorry, but the number of rows and columns of the Theta_M2 matrix is the same, so the checks on ")
            println("correct orientation will not work. Please make sure that the dimensions follow: ")
            println("theta[samples, parameters]")
            println("-------------------------------------------------------------")
        end
    end

    # Check that there are multiple samples or single samples for both models (either frequentist or bayesian)
    if length(oedms_def["Theta_M1"])/oedms_def["Model_1"]["nPar"] == 1 &&
        length(oedms_def["Theta_M2"])/oedms_def["Model_2"]["nPar"] == 1
        println("-------------------------- INFO --------------------------")
        println("Single samples for both models parameters are given. Frequentist scenario will be used (Euclidean Distance)")
        println("----------------------------------------------------------")
    elseif length(oedms_def["Theta_M1"])/oedms_def["Model_1"]["nPar"] != 1 &&
        length(oedms_def["Theta_M2"])/oedms_def["Model_2"]["nPar"] != 1
        if length(oedms_def["Theta_M1"])/oedms_def["Model_1"]["nPar"] <50
            println("-------------------------- WARNING --------------------------")
            println("Less than 50 samples for theta are given for Model 1. ")
            println("Consider using more samples for better results. ")
            println("-------------------------------------------------------------")
        end
        if length(oedms_def["Theta_M2"])/oedms_def["Model_2"]["nPar"] < 50
            println("-------------------------- WARNING --------------------------")
            println("Less than 50 samples for theta are given for Model 2. ")
            println("Consider using more samples for better results. ")
            println("-------------------------------------------------------------")
        end
    else
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, it seems that you gave a single vector of parameters for one model ")
        println("and various samples for the other. Please be consistent and either give a single")
        println("or else multiple samples for both models")
        return
    end


    if oedms_def["tsamps"][end] > oedms_def["finalTime"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check finalTime. You have selected a sampling point past it.")
        return
    end
    if (oedms_def["Model_1"]["Y0eqs"] != []) && oedms_def["preInd_M1"]==[]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("You specified computation of Y0 as steady state but you have not introduced an inducer value for it!")
        return
    end
    if (oedms_def["Model_2"]["Y0eqs"] != []) && oedms_def["preInd_M2"]==[]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("You specified computation of Y0 as steady state but you have not introduced an inducer value for it!")
        return
    end

    # Check that no step is no smaller than 2 unit of time
    for j in 1:length(oedms_def["switchT"])-1
        if (oedms_def["switchT"][j+1]-oedms_def["switchT"][j])<=4
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but 2 of the steps in the experiment are too close. This package cannot "))
            println("handle this for now. We are working on it!")
            return
        end
    end


    if typeof(oedms_def["Obs"]) == Array{String,1}
        if convert(Bool,sum(occursin.("+", oedms_def["Obs"]))) || convert(Bool,sum(occursin.("-", oedms_def["Obs"]))) ||
           convert(Bool,sum(occursin.("/", oedms_def["Obs"]))) ||
           convert(Bool,sum(occursin.("*", oedms_def["Obs"]))) || convert(Bool,sum(occursin.("^", oedms_def["Obs"])))
            nothing
        else
            if sum(sum.([occursin.(oedms_def["Model_1"]["stName"][k], oedms_def["Model_1"]["Obs"]) for k in 1:length(oedms_def["Model_1"]["stName"])])) == 0 ||
                sum(sum.([occursin.(oedms_def["Model_2"]["stName"][k], oedms_def["Model_2"]["Obs"]) for k in 1:length(oedms_def["Model_2"]["stName"])])) == 0
            # if sum(occursin.(oedms_def["Model_1"]["stName"], oedms_def["Obs"])) == 0 || sum(occursin.(oedms_def["Model_2"]["stName"], oedms_def["Obs"])) == 0
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but there is some issue with the contents of the field Obs."))
                println("It seems that the observable(s) selected do not match any state for one of the models")
                return
            end
        end
    end

    for j in 1:length(oedms_def["Obs"])
        s1=[occursin(oedms_def["Model_1"]["stName"][i], oedms_def["Obs"][j]) for i in 1:(oedms_def["Model_1"]["nStat"])]
        s2=[occursin(oedms_def["Model_2"]["stName"][i], oedms_def["Obs"][j]) for i in 1:(oedms_def["Model_2"]["nStat"])]

        if sum(s1) == 0
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but you have indicated as observable something that is not a state in Model 1."))
            println("The same observable has to appear in both models for the optimisation to work.")
            return
        end

        if sum(s2) == 0
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but you have indicated as observable something that is not a state in Model 2."))
            println("The same observable has to appear in both models for the optimisation to work.")
            return
        end

    end

    if  length(oedms_def["y0_M1"])/oedms_def["Model_1"]["nStat"] == 1
        if oedms_def["Model_1"]["nStat"] != length(oedms_def["y0_M1"])
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but the y0 vector for model 1 does not have the same amount of values as states in the model."))
            println("Remember that even if it will not be used you need to provide a y0 vector.")
            return
        end
    elseif length(oedms_def["y0_M1"])/oedms_def["Model_1"]["nStat"] > 1
        if length(oedms_def["y0_M1"])/oedms_def["Model_1"]["nStat"] != length(oedms_def["Theta_M1"])/oedms_def["Model_1"]["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but the y0 matrix for model 1 does not have the same amount of entries as parameter samples."))
            println("Remember that even if it will not be used you need to provide a y0 vector.")
            return
        end
        if size(oedms_def["y0_M1"])[1] == size(oedms_def["y0_M1"])[2]
            println("-------------------------- WARNING --------------------------")
            println(string("Sorry, but the number of rows and columns of the y0_M1 matrix is the same, so the checks on "))
            println("correct orientation will not work. Please make sure that the dimensions follow: ")
            println("y0_M1[samples, states]")
        end

    end

    if  length(oedms_def["y0_M2"])/oedms_def["Model_2"]["nStat"] == 1
        if oedms_def["Model_2"]["nStat"] != length(oedms_def["y0_M2"])
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but the y0 vector for model 2 does not have the same amount of values as states in the model."))
            println("Remember that even if it will not be used you need to provide a y0 vector.")
            return
        end
    elseif length(oedms_def["y0_M2"])/oedms_def["Model_2"]["nStat"] > 1
        if length(oedms_def["y0_M2"])/oedms_def["Model_2"]["nStat"] != length(oedms_def["Theta_M2"])/oedms_def["Model_2"]["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but the y0 matrix for model 2 does not have the same amount of entries as parameter samples."))
            println("Remember that even if it will not be used you need to provide a y0 vector.")
            return
        end
        if size(oedms_def["y0_M2"])[1] == size(oedms_def["y0_M2"])[2]
            println("-------------------------- WARNING --------------------------")
            println(string("Sorry, but the number of rows and columns of the y0_M2 matrix is the same, so the checks on "))
            println("correct orientation will not work. Please make sure that the dimensions follow: ")
            println("y0_M2[samples, states]")
        end

    end

    if oedms_def["Model_1"]["nInp"] != length(oedms_def["preInd_M1"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but the preInd_M! vector for model 1 does not have the same amount of values as inputs in the model."))
        println("Remember that even if it will not be used you need to provide a preInd_M1 vector.")
        return
    elseif oedms_def["Model_2"]["nInp"] != length(oedms_def["preInd_M2"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but the preInd_M2 vector for model 2 does not have the same amount of values as inputs in the model."))
        println("Remember that even if it will not be used you need to provide a preInd_M2 vector.")
        return
    end

    if round(oedms_def["finalTime"]) == 0
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but it seems that the final time is 0?"))
        println("Remember that for now, the minimum time discretisation allowed is 1 so adjust the parameter scales according to it.")
        return
    end

    if length(oedms_def["switchT"]) <= 1
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but it seems that you have not introduced a correct switchT"))
        println("Remember, this are the switching times for the inducer, including time 0 and the final time. ")
        return
    end

    # Check that no step is no smaller than 2 unit of time
    for j in 1:length(oedms_def["tsamps"])-1
        if (round(oedms_def["tsamps"][j+1])-round(oedms_def["tsamps"][j]))<=1
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but 2 of the sampling times are too close. "))
            println("Remember that for now, the minimum time discretisation allowed is 1 so adjust the parameter scales according to it.")
            println("Sorry for the inconvenience, I am working on it...")
            return
        end
    end

    if oedms_def["Model_1"]["nInp"] != oedms_def["Model_2"]["nInp"]
        println("-------------------------- WARNING --------------------------")
        println("The amount of inducers is different in the two models.")
        println("Be careful on how you configure that. ")
        println("-------------------------------------------------------------")
    end

    if oedms_def["Model_1"]["nInp"] == 1 && oedms_def["Model_2"]["nInp"] == 1
        if oedms_def["fixedInp"] != []
            println("-------------------------- WARNING --------------------------")
            println("You have selected one of the input of the models to be fixed, but there is only 1 inducer in the models.")
            println("This entry will be ignored. If you want to design a single step experiment adjust switchT accordingly.")
            println("-------------------------------------------------------------")
            oedms_def["fixedInp"] = [];
        end
    else
        if oedms_def["fixedInp"] != []
            s1=[(oedms_def["Model_1"]["inpName"][i] in oedms_def["fixedInp"]) for i in 1:(oedms_def["Model_1"]["nInp"])]
            s2=[(oedms_def["Model_2"]["inpName"][i] in oedms_def["fixedInp"]) for i in 1:(oedms_def["Model_2"]["nInp"])]

            if sum(s1)+sum(s2) == 0
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but you have indicated as fixedInp something that is not an input in Model 1 or Model 2."))
                println("The inputs have to appear in one of the models for the optimisation to work.")
                return
            end
        end
        # if sum(s1) == 0
        #     println("-------------------------- Process STOPPED!!! --------------------------")
        #     println(string("Sorry, but you have indicated as fixedInp something that is not an input in Model 1."))
        #     println("The same inputs have to appear in both models for the optimisation to work.")
        #     return
        # end
        #
        # if sum(s2) == 0
        #     println("-------------------------- Process STOPPED!!! --------------------------")
        #     println(string("Sorry, but you have indicated as fixedInp something that is not an input in Model 2."))
        #     println("The same inputs have to appear in both models for the optimisation to work.")
        #     return
        # end
    end

    tmp1 = [];
    tmp2 = [];
    for i in 1:(oedms_def["Model_1"]["nInp"])
        if oedms_def["Model_1"]["inpName"][i] in oedms_def["fixedInp"]
            tmp1 = vcat(tmp1, oedms_def["Model_1"]["inpName"][i])
        end
    end

    for i in 1:(oedms_def["Model_2"]["nInp"])
        if oedms_def["Model_2"]["inpName"][i] in oedms_def["fixedInp"]
            tmp2 = vcat(tmp2, oedms_def["Model_2"]["inpName"][i])
        end
    end

    if length(tmp1) >= oedms_def["Model_1"]["nInp"] && length(tmp2) >= oedms_def["Model_2"]["nInp"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but it seems that you want fixed all the inputs."))
        println("To do this, please adjust switchT instead of trying it through fixedInp")
        return
    elseif length(tmp1) >= oedms_def["Model_1"]["nInp"] && length(tmp2) <= oedms_def["Model_2"]["nInp"]
        println("-------------------------- WARNING --------------------------")
        println("Be careful, it seems that all the inputs for Model 1 will be fixed.")
        println("-------------------------------------------------------------")
    elseif length(tmp1) <= oedms_def["Model_1"]["nInp"] && length(tmp2) >= oedms_def["Model_2"]["nInp"]
        println("-------------------------- WARNING --------------------------")
        println("Be careful, it seems that all the inputs for Model 2 will be fixed.")
        println("-------------------------------------------------------------")
    end

    if length(oedms_def["uUpper"]) != length(oedms_def["uLower"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("You have given a different number of upper and lower bounds for the inducers.")
        println("Please, check this before proceeding.")
        return
    end

    if maximum([oedms_def["Model_1"]["nInp"], oedms_def["Model_2"]["nInp"]]) != length(oedms_def["uUpper"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("You have introduced less bounds than inputs the models have or vice-versa")
        println("Please, check this before proceeding.")
        return
    end

    tm1 = oedms_def["uUpper"].>= oedms_def["uLower"]
    if sum(tm1) < length(tm1)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Check the input bounds please, it seems that upper and lower bounds are not correct.")
        return
    end

    if length(oedms_def["switchT"]) == 2 && oedms_def["fixedStep"] != []
        if typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{Any,1}},1} # Ther might be 1 step but multiple inputs and one or more are to be fixed and the others optimised.
            println("-------------------------- WARNING !!! --------------------------")
            println("You are designing a single step experiment but you have introduced something in the field fixedStep.")
            println("Since this would result in no optimisation at all, the contents of this field will be ignored.")
            oedms_def["fixedStep"] = [];
        end
    end

    if oedms_def["fixedStep"] != []
        fs = [oedms_def["fixedStep"][j][1] for j in 1:length(oedms_def["fixedStep"])];
        if minimum(fs) < 1
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have introduced and index for the step smaller than 1. This does not make much sense.")
            return
        elseif maximum(fs) > (length(oedms_def["switchT"])-1)
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have introduced and index for the step larger than the number of steps of the experiment. This does not make much sense.")
            return
        elseif length(fs) > (length(oedms_def["switchT"])-1)
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have introduced more fixed steps than steps in the experiment. This does not make much sense.")
            return
        elseif unique(fs) != fs
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have introduced a repeated fixed step for the inputs. This does not make much sense.")
            return
        end

        fi = [oedms_def["fixedStep"][j][2] for j in 1:length(oedms_def["fixedStep"])];
        for j in 1:length(fi)
            if length(fi[j]) != maximum([oedms_def["Model_1"]["nInp"], oedms_def["Model_2"]["nInp"]]) - length(oedms_def["fixedInp"])
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Sorry, but it seems that there is a mismatch between the number of inputs and the number of values in the second entry of the tuple for the fixed step entry.")
                println("Remember that any fixed inputs (fixedInp) are not considered in this field, so you should not include them here.")
                println("Tip: If you want a fixed input with a fixed value in the experiment, hard-code its value in the definition of the model.")
                return
            end
        end
    end

    if length(oedms_def["uUpper"]) != maximum([oedms_def["Model_1"]["nInp"], oedms_def["Model_2"]["nInp"]]) ||
        length(oedms_def["uLower"]) != maximum([oedms_def["Model_1"]["nInp"], oedms_def["Model_2"]["nInp"]])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, but it seems that you have not introduced as many bounds as total amount of inducers between the two models.")
        return
    end

    for k in 1:length(oedms_def["equalStep"])
        if length(oedms_def["equalStep"][k]) == 1
            println("-------------------------- WARNING !!! --------------------------")
            println("One of the entries for equalStep only has 1 index, please avoid doing this.")
        end
        if minimum(oedms_def["equalStep"][k]) < 1
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have not introduced and index for a step smaller than 1. ")
            return
        elseif maximum(oedms_def["equalStep"][k]) > (length(oedms_def["switchT"])-1)
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have not introduced and index for a step higher than the actual number of steps in the experiment.")
            return
        end
        oedms_def["equalStep"][k] = sort(oedms_def["equalStep"][k]);
    end

    if oedms_def["equalStep"] != []
        unni = [oedms_def["equalStep"][i][j] for i in 1:length(oedms_def["equalStep"]) for j in 1:length(oedms_def["equalStep"])];
        if length(unni) != length(unique(unni))
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but it seems that you have not introduced a repeated index somewhere in the definition of equalStep.")
            return
        end
    end

    if oedms_def["equalStep"] != [] && oedms_def["fixedStep"] != []
        fs = [oedms_def["fixedStep"][j][1] for j in 1:length(oedms_def["fixedStep"])];
        if typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{Any,1}},1}
            for k in 1:length(oedms_def["equalStep"])
                for h in 1:length(fs)
                    if fs[h] in oedms_def["equalStep"][k]  &&  typeof(oedms_def["fixedStep"][h]) != Tuple{Int,Array{Any,1}}
                        println("-------------------------- Process STOPPED!!! --------------------------")
                        println("Sorry, but it seems that you have defined an equalStep that is also a fixedStep.")
                        return
                    end
                end
            end
        end
    end

    for i in 1:length(oedms_def["fixedStep"])
        for j in 1:length(oedms_def["fixedStep"][i][2])
            if typeof(oedms_def["fixedStep"][i][2][j]) != Int && typeof(oedms_def["fixedStep"][i][2][j]) != Float64 && typeof(oedms_def["fixedStep"][i][2][j]) != Float32 &&
                (oedms_def["fixedStep"][i][2][j]) != Any
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Sorry, but it seems you have introduced something wrong in the second entry of the tuple.")
                println("This has to be a number indicating a concentration or the variable Any if that step is not fixed for that specific inducer. ")
                return
            end
        end
    end

    return(oedms_def)
end

## Bhattacharyya Distance
function BhattacharyyaDist(mu1, mu2, sd1, sd2)

    E = (sd1+sd2)/2;
    Em1 = inv(E);
    dE = abs(det(E));

    t1 = mu1'-mu2';
    t2 = mu1-mu2;

    ft = (1/8)*t1*Em1*t2;
    trk = sqrt(abs(det(sd1))*abs(det(sd2)));

    if trk == 0
        trk = 1e-300;
    end

    st = dE/trk;

    bhd = ft+0.5*log(st);

    return(bhd);

end

## Euclidean Distance
function EuclideanDist(sm1, sm2)

    eud = sum((sm1.-sm2).^2)^0.5;

    return(eud)

end

## Function that will generate the necessary scripts for the optimisation

# Note for the user: covariance matrices need to be regularised to ensure that they are positive definite (some really low variances will be considered as 0 as well) and full ranked. For this we add 0.1 in all the elements of the diagonal. So, the maximum resolution of this method is a variance of 0.1 for a given time point.
# If you are working with a normalised observable note that variances might be really low and hence all of them will be considered 0.1. The method would not work anyway doe to the really low variances (issues computing determinatn and amtrix inverse) so for now it is recomended to not use normalised observables ranged between 0 and 1 (make it be between 0 and 100 minimum so the regularisation does not have a strong effect.)
# Other ways to avoid issues with computation of determinatns and inverse of covariance matrices will be investigated and inplemented in the future (reduction of the number of time points considered can also be a way, since the more spaced the sampling points are the more likely numerical issues will be avoided and the number of linearly independent colums could increase, hence increase the rank). More work will be done in this section in the future.

function genOptimMSFuncts(oedms_def)

    oedms_def = checkStructOEDMS(oedms_def);

    # Check Models and generate the script. Also check that they do not have the same name
    if typeof(oedms_def["Model_1"]) <: Dict
        oedms_def["Model_1"] = checkStruct(oedms_def["Model_1"]);
    elseif typeof(oedms_def["Model_1"]) <: Array
        oedms_def["Model_1"] = checkStruct(oedms_def["Model_1"][1]);
    end;
    oedms_def["Model_1"] = GenerateModel(oedms_def["Model_1"]);

    if typeof(oedms_def["Model_2"]) <: Dict
        oedms_def["Model_2"] = checkStruct(oedms_def["Model_2"]);
    elseif typeof(oedms_def["Model_2"]) <: Array
        oedms_def["Model_2"] = checkStruct(oedms_def["Model_2"][1]);
    end;

    # All this is to extract the correct order for the inputs (order of Model_1 and tagged after any extra inputs for model 2)
    # so it can be changed in model 2 before generating the scripts.
    if sort(oedms_def["Model_2"]["inpName"]) == sort(oedms_def["Model_1"]["inpName"]) # Check if both models have the same inputs
        oedms_def["Model_2"]["inpName"] = oedms_def["Model_1"]["inpName"]; # Put order of model 1
    else
        inpnam = oedms_def["Model_1"]["inpName"]; # Extrsct inputs of model 1 and tag the rest after. Then change the inpNam field in model 2
        if oedms_def["Model_1"]["nInp"] < oedms_def["Model_2"]["nInp"]
            tmp1 = [];
            try
                for i in 1:oedms_def["Model_2"]["nInp"]
                    if oedms_def["Model_2"]["inpName"][i] in oedms_def["Model_1"]["inpName"]
                        nothing
                    else
                        tmp1 = vcat(tmp1, oedms_def["Model_2"]["inpName"][i])
                    end
                end
                inpnam = vcat(inpnam, tmp1)
            catch
                for i in 1:oedms_def["Model_2"]["nInp"]
                    if oedms_def["Model_2"]["inpName"][i] in oedms_def["Model_1"]["inpName"]
                        nothing
                    else
                        tmp1 = hcat(tmp1, oedms_def["Model_2"]["inpName"][i])
                    end
                end
                inpnam = hcat(inpnam, tmp1)
            end
        end;
        oedms_def["Model_2"]["inpName"] = convert(Array{String,1},inpnam);
    end
    oedms_def["Model_2"] = checkStruct(oedms_def["Model_2"]); # Check structure again just in case
    oedms_def["Model_2"] = GenerateModel(oedms_def["Model_2"]);


    # Generate results directory
    cudi = pwd(); #@__DIR__;
    cudi2 = @__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"_",today()))
    end
    if !isdir(string(cudi, "\\Results\\", oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"_",today(), "\\OEDModelSelectionScripts"))
        mkdir(string(cudi, "\\Results\\", oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"_",today(), "\\OEDModelSelectionScripts"))
    end

    # Packages needed for the simulations
    Head = string("
        ");

    # Manage inputs order (see before)
    inpnam = oedms_def["Model_1"]["inpName"];
    if oedms_def["Model_1"]["nInp"] < oedms_def["Model_2"]["nInp"]
        tmp1 = [];
        try
            for i in 1:oedms_def["Model_2"]["nInp"]
                if oedms_def["Model_2"]["inpName"][i] in oedms_def["Model_1"]["inpName"]
                    nothing
                else
                    tmp1 = vcat(tmp1, oedms_def["Model_2"]["inpName"][i])
                end
            end
            inpnam = vcat(inpnam, tmp1)
        catch
            for i in 1:oedms_def["Model_2"]["nInp"]
                if oedms_def["Model_2"]["inpName"][i] in oedms_def["Model_1"]["inpName"]
                    nothing
                else
                    tmp1 = hcat(tmp1, oedms_def["Model_2"]["inpName"][i])
                end
            end
            inpnam = hcat(inpnam, tmp1)
        end
    end;

    steps = length(oedms_def["switchT"])-1; # Number of steps of the experiment
    induc = maximum([oedms_def["Model_1"]["nInp"], oedms_def["Model_2"]["nInp"]]); # Total number of inducers
    fiii = length(oedms_def["fixedInp"]); # Number of fixed inputs

    doof = (steps*(induc-fiii))+fiii; # Degrees of freedom of the problem.

    # Some of the following sections are heavily commented since otherwise it is really hard to follow the code.
    # This sections are called MOUNTAINS since there is so many loops and exceptions that if you tilt the screen it draws a series of mountains.

    ustr = Array{String,1}(undef, doof); # String vector that will contain the names for the inputs to be optimised at each step. Inducer names are in order (first all variables of on einducer, then all for the next and so on.)
                                        # If an input is fixed to a value, then it will not be considered here
                                        # If an input is fixed it will only be taken into account Once
                                        # If an input concentration is repeated in multiple steps, then it will only be considered once
    us = [oedms_def["fixedStep"][j][1] for j in 1:length(oedms_def["fixedStep"])]; # Fixed Steps (the actual step number)
    counn = 1; # Index for the entry of ustr that will advance only if certain conditions are meet (that that position has been already checked and a variable or number has been introduced)
                # Tecnichally this indicated the first entry of ustr for each inducer variables.
    fixin = 0; # Index to account for fixedInp. If the fixedInp is the first one in the inpam list this is the best way (that I have found) to account for it. This is because of the use of i (input index). This is because first I consider the case where there is a fisedInp and then check the rest.
                # Necessary since fixedStep does not consider fixedInps
    # THE FIRST MOUNTAINS!
    for i in 1:length(inpnam) # Loop over each inducer. i will be the index from the list inpnam
        if inpnam[i] in oedms_def["fixedInp"] # Check if input i is in the list of fixed inputs
            ustr[counn] = inpnam[i]; # If yes, add the name of the input as the only variable to optimise for that input.
            counn += 1; # Need to move index in ustr
            fixin += 1; # Need to take into account that one input is being fixed
        else # Now check the rest of the cases
            for j in 0:steps-1 # j will be the index of the steps -1. If 6 steps, then a value between 0 and 5. Only used to check fixedStep.
                if j+1 in us # Checl if the step is present in the list of fixedStep
                    for u in 1:length(us) # To account for the case that the user might have introduced the index order not sorted, so we fing the j+1 that matches the us entry
                        if j+1 == us[u]
                            if oedms_def["fixedStep"][u][2][i-fixin] != Any # If user set that as a fixed value do nothing, since this one will be optimised, not fixed. Coded like this to account for future possible modifications (and to not get lost).
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
            if oedms_def["equalStep"] == [] # This section accounts for the contents of equalStep, so if there is nothing, skip.
                nothing
            else
                if typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{Any,1}},1} && # This section is only used if all entries in fixed step are numbers (so no Any from the user)
                typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{T,1} where T},1}
                    tmpin = ustr[counn:counn+steps-1]; # Temporary variable that takes all the entries of ustr that correspond to the current inducer i.
                    for k in 1:length(oedms_def["equalStep"]) # Loop over each equalStep entry and add the same variable name in all the steps in that entry (steps where the inducer is going to have the same value)
                        tmpin[oedms_def["equalStep"][k]] .= tmpin[oedms_def["equalStep"][k][1]]
                    end
                    ustr[counn:counn+steps-1] = tmpin; # Define the corresponding ustr filling it with the temporary variable
                else
                    tmpin = ustr[counn:counn+steps-1]; # Temporary variable that takes all the entries of ustr that correspond to the current inducer i.
                    for k in 1:length(oedms_def["equalStep"]) # Loop over each entry of equalStep to see if a variable name will be assigned or not (account for fixedStep with Any entries)
                        tmpfe = zeros(length(oedms_def["equalStep"][k])); # Temporary variable that will contain the same entries of equalStep but only if that input for that step is not fixed (can happen due to the option of introducing Any)
                        for j in 1:length(tmpfe) # Loop over each entry of the equalStep entry that we are in now.
                            if oedms_def["equalStep"][k][j] in us # Check if the entry is in the list of equal step and if for that inducer in this step the user specified Any (input to be optimised there) then add the step index in the temporary variable so a variable name can be assigned later. If not nothing since that value would have been fixed before.
                                for w in 1:length(us)
                                    if oedms_def["equalStep"][k][j] == us[w] && oedms_def["fixedStep"][w][2][i-fixin] == Any
                                        tmpfe[j] = oedms_def["equalStep"][k][j];
                                    end
                                end
                            else
                                tmpfe[j] = oedms_def["equalStep"][k][j]; # If the entrie is not present in fixedStep, then add the entrie in the temporary variable.
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
    fs = [oedms_def["fixedStep"][j][2] for j in 1:length(oedms_def["fixedStep"])]; # Input values (or Any entry) for each step specified in fixedStep
    cnt = 1; # Variables not used in current version. Left for checks to see where the if statemetns go. This is for model 1 the one below for model 2.
    cnt2 = 1;

    # THE SECOND MOUNTAINS!
    for i in 1:steps # length(r) to go through each step
        fixin1 = 0; # To discard the fixed inputs (as done before). Separated for both models, but technically only using one should do.
        fixin2 = 0;
        for j in 0:induc-1 # j will be the index of the inducer considered -1. If 2 inducers, then 0 and 1 (done to ease later indexing in previous versions).
            if inpnam[j+1] in oedms_def["fixedInp"] # Consider case that the current inducer is fixed. If not check the rest (fixedStep and equalStep)
                fixin1 += 1; # Same reason ad fixin in the FIRSTS MOUNTAINS
                fixin2 += 1;
                if inpnam[j+1] in oedms_def["Model_1"]["inpName"] # The input is fixed, so add the same variable name for each step. For model 1
                    inpstr1[r[i]+j] = inpnam[j+1];
                else
                    inpstr1[r[i]+j] = "";
                end
                if inpnam[j+1] in oedms_def["Model_2"]["inpName"] # The input is fixed, so add the same variable name for each step. For model 2
                    inpstr2[r[i]+j] = inpnam[j+1];
                else
                    inpstr2[r[i]+j] = "";
                end
            else
                if inpnam[j+1] in oedms_def["Model_1"]["inpName"] # Check if the current inducer is present in Model 1. If so proceed.
                    if typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{Any,1}},1} && # Consider the case where in fixedStep all the inducers are fixed to a value and there is no Any. Was easier to consider the case of presence or absence of Anys separately.
                    typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{T,1} where T},1}
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
                        if oedms_def["equalStep"] == []  # This section accounts for the contents of equalStep, so if there is nothing, skip.
                            if !isassigned(inpstr1, r[i]+j) # If current entry is not assigned to anything (from the fixedStep contents), the add the corersponding variable name to the current step for the current input
                                inpstr1[r[i]+j] = string(inpnam[j+1],i);
                                cnt += 1
                            end
                        else
                            for k in 1:length(oedms_def["equalStep"]) # Check for all the equal steps and add the same variable name for the current input.
                                if i in oedms_def["equalStep"][k]
                                    inpstr1[r[i]+j] = string(inpnam[j+1],oedms_def["equalStep"][k][1]);
                                end
                            end
                        end
                    else # This accounts for Any entries in fixedStep
                        if i in us # If the current step is present in fixedSteps.
                            for w in 1:length(us) # This loop and inequality is to account for the case where the user does not introduce the entries in fixedStep in order (sorted by step)
                                if j+1-fixin1 != 0
                                    if i == us[w] && oedms_def["fixedStep"][w][2][j+1-fixin1] != Any # If this input for the current step is fixed to a value for the user
                                        inpstr1[r[i]+j] = string(fs[w][j+1-fixin1]); # Hard code the value introduced by the user
                                        cnt += 1;
                                    elseif i == us[w] && oedms_def["fixedStep"][w][2][j+1-fixin1] == Any # If this input for the current step is set to Any, hence to be optimised
                                        inpstr1[r[i]+j] = string(inpnam[j+1],i); # Add the corresponting variable name.
                                        cnt += 1;
                                    end
                                end
                            end
                        else # Entry not in fixedStep, so add the corresponding variable name for the step and input.
                            inpstr1[r[i]+j] = string(inpnam[j+1],i);
                            cnt += 1;
                        end
                        if oedms_def["equalStep"] == [] # This section accounts for the contents of equalStep, so if there is nothing, skip.
                            if !isassigned(inpstr1, r[i]+j) # If current entry is not assigned to anything (from the fixedStep contents), the add the corersponding variable name to the current step for the current input
                                inpstr1[r[i]+j] = string(inpnam[j+1],i);
                                cnt += 1;
                            end
                        else # Now account for steps that have to have the same variable name (even though it has been coded differently in the previous loop)
                            for k in 1:length(oedms_def["equalStep"]) # Loop over every entry in equalStep (arrays with the step indexes to have the same value)
                                tmpfe = zeros(length(oedms_def["equalStep"][k])); # Temporary variable that will have the same values than the current equalStep but without the ones where the input value is fixed to a number by the user
                                for z in 1:length(tmpfe) # Loop over each entry in the current equalStep
                                    if oedms_def["equalStep"][k][z] in us # If the current value of the equalStep entry (step index) is present in fixedStep
                                        for w in 1:length(us) # To avoid indexing issues if user does not introduce things sorted by step index in fixedStep
                                            if 1+j-fixin1 != 0  # To avoid the case of indexing by 0. This should not happen, but had some issues during tests so left it just in case.
                                                if oedms_def["equalStep"][k][z] == us[w] && oedms_def["fixedStep"][w][2][1+j-fixin1] == Any # If the current step from the equalStep entry is also in fixedStep and this is set to Any (so to be optimised)
                                                    tmpfe[z] = oedms_def["equalStep"][k][z]; # Add the step index in the temporary variable
                                                end
                                            end
                                        end
                                    else # If the step in the equalStep entry is not in fixedStep then it has to be optimsied, so add the entry.
                                        tmpfe[z] = oedms_def["equalStep"][k][z];
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


                # Exactly the same as before but for model 2
                if inpnam[j+1] in oedms_def["Model_2"]["inpName"]
                    if typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{Any,1}},1} &&
                    typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{T,1} where T},1}
                        if i in us
                            for w in 1:length(us)
                                if i == us[w]
                                    if j+1-fixin2 != 0
                                        inpstr2[r[i]+j] = string(fs[w][j+1-fixin2]);
                                        cnt2 += 1;
                                    end
                                end
                            end
                        else
                            inpstr2[r[i]+j] = string(inpnam[j+1],i);
                            cnt2 += 1;
                        end
                        if oedms_def["equalStep"] == []
                            if !isassigned(inpstr2, r[i]+j)
                                inpstr2[r[i]+j] = string(inpnam[j+1],i);
                                cnt2 += 1;
                            end
                        else
                            for k in 1:length(oedms_def["equalStep"])
                                if i in oedms_def["equalStep"][k]
                                    inpstr2[r[i]+j] = string(inpnam[j+1],oedms_def["equalStep"][k][1]);
                                end
                            end
                        end
                    else
                        if i in us
                            for w in 1:length(us)
                                if j+1-fixin2 != 0
                                    if i == us[w] && oedms_def["fixedStep"][w][2][j+1-fixin2] != Any
                                        inpstr2[r[i]+j] = string(fs[w][j+1-fixin2]);
                                        cnt2 += 1;
                                    elseif i == us[w] && oedms_def["fixedStep"][w][2][j+1-fixin2] == Any
                                        inpstr2[r[i]+j] = string(inpnam[j+1],i);
                                        cnt2 += 1;
                                    end
                                end
                            end
                        else
                            inpstr2[r[i]+j] = string(inpnam[j+1],i);
                        end

                        if oedms_def["equalStep"] == []
                            if !isassigned(inpstr2, r[i]+j)
                                inpstr2[r[i]+j] = string(inpnam[j+1],i);
                                cnt2 += 1;
                            end
                        else
                            for k in 1:length(oedms_def["equalStep"])
                                tmpfe = zeros(length(oedms_def["equalStep"][k]));
                                for z in 1:length(tmpfe)
                                    if oedms_def["equalStep"][k][z] in us
                                        for w in 1:length(us)
                                            if j+1-fixin2 != 0
                                                if oedms_def["equalStep"][k][z] == us[w] && oedms_def["fixedStep"][w][2][1+j-fixin2] == Any
                                                    tmpfe[z] = oedms_def["equalStep"][k][z];
                                                end
                                            end
                                        end
                                    else
                                        tmpfe[z] = oedms_def["equalStep"][k][z];
                                    end
                                end
                                tmpfe = convert.(Int, filter!(e->e∉[0],tmpfe));
                                if i in tmpfe
                                    inpstr2[r[i]+j] = string(inpnam[j+1],tmpfe[1]);
                                end
                            end
                        end

                    end
                else
                    inpstr2[r[i]+j] = "";
                end
            end
        end
    end
    inpstr1 = [inpstr1[i] for i in 1:length(inpstr1) if isassigned(inpstr1, i) && inpstr1[i] != ""]; # Delete unasigned entries (should not happen but just in case) and emtpy strings (some input not present in the model)
    inpstr2 = [inpstr2[i] for i in 1:length(inpstr2) if isassigned(inpstr2, i) && inpstr2[i] != ""];

    # If it turns out that for one model the simulation ends up having just 1 step and 1 inducer (important) then reduce the input to 1, so the simulation will be faster.
    if length(unique(inpstr1)) == 1
        inpstr1 = unique(inpstr1);
        spm1 = "[spMS[1] spMS[end]]";
    else
        spm1 = "spMS";
    end

    if length(unique(inpstr2)) == 1
        inpstr2 = unique(inpstr2);
        spm2 = "[spMS[1] spMS[end]]";
    else
        spm2 = "spMS";
    end

    obsM1 = Array{String,1}(undef, length(oedms_def["Obs"]));
    obsM2 = Array{String,1}(undef, length(oedms_def["Obs"]));
    for i in 1:length(oedms_def["Obs"])
        oedms_def["Obs"][i] = replace(oedms_def["Obs"][i], ".+"=>"+")
        oedms_def["Obs"][i] = replace(oedms_def["Obs"][i], ".-"=>"-")
        oedms_def["Obs"][i] = replace(oedms_def["Obs"][i], ".*"=>"*")
        oedms_def["Obs"][i] = replace(oedms_def["Obs"][i], "./"=>"/")
        oedms_def["Obs"][i] = replace(oedms_def["Obs"][i], ".^"=>"^")

        oedms_def["Obs"][i] = replace(oedms_def["Obs"][i], "+"=>" .+")
        oedms_def["Obs"][i] = replace(oedms_def["Obs"][i], "-"=>" .-")
        oedms_def["Obs"][i] = replace(oedms_def["Obs"][i], "*"=>" .*")
        oedms_def["Obs"][i] = replace(oedms_def["Obs"][i], "/"=>" ./")
        oedms_def["Obs"][i] = replace(oedms_def["Obs"][i], "^"=>" .^")
        for j in 1:oedms_def["Model_1"]["nStat"]
            if occursin(oedms_def["Model_1"]["stName"][j], oedms_def["Obs"][i])
                tmpo1 = replace(oedms_def["Obs"][i], oedms_def["Model_1"]["stName"][j] => string("solM1[:,",j,",:]"));
                obsM1[i] = string("    Obs",i,"_M1 = ",tmpo1,"; \n")
            end
        end

        for j in 1:oedms_def["Model_2"]["nStat"]
            if occursin(oedms_def["Model_2"]["stName"][j], oedms_def["Obs"][i])
                tmpo2 = replace(oedms_def["Obs"][i], oedms_def["Model_2"]["stName"][j] => string("solM2[:,",j,",:]"));
                obsM2[i] = string("    Obs",i,"_M2 = ",tmpo2,"; \n")
            end
        end
    end

    # Assess if frequentist or bayesian case and write the utility/cost computation step

    if length(oedms_def["Theta_M1"])/oedms_def["Model_1"]["nPar"]==1
        ccomp = string("
",join([string("    Obs",i,"_M1 = Obs",i,"_M1[:,1]; \n") for i in 1:length(oedms_def["Obs"])]),"
",join([string("    Obs",i,"_M2 = Obs",i,"_M2[:,1]; \n") for i in 1:length(oedms_def["Obs"])]),"

",join([string("    EUD",i, " = EuclideanDist(Obs",i,"_M1, Obs",i,"_M2); \n") for i in 1:length(oedms_def["Obs"])]),"

    util = mean(",join([string("EUD",i) for i in 1:length(oedms_def["Obs"])], ","),");
");
else

    ccomp = string("
",join([string("    mu",i,"_M1 = mean(Obs",i,"_M1, dims=2); \n") for i in 1:length(oedms_def["Obs"])]),"
",join([string("    sd",i,"_M1 = cov(Obs",i,"_M1, dims=2).+(0.1*Matrix((Diagonal(ones(size(Obs",i,"_M1)[1]))))); \n") for i in 1:length(oedms_def["Obs"])]),"

",join([string("    mu",i,"_M2 = mean(Obs",i,"_M2, dims=2); \n") for i in 1:length(oedms_def["Obs"])]),"
",join([string("    sd",i,"_M2 = cov(Obs",i,"_M2, dims=2).+(0.1*Matrix((Diagonal(ones(size(Obs",i,"_M2)[1]))))); \n") for i in 1:length(oedms_def["Obs"])]),"

",join([string("    BHD",i, " = BhattacharyyaDist(mu",i,"_M1[:,1], mu",i,"_M2[:,1], sd",i,"_M1, sd",i,"_M2); \n") for i in 1:length(oedms_def["Obs"])]),"

    util = mean(",join([string("BHD",i) for i in 1:length(oedms_def["Obs"])], ","),");

        ");
    end;

    Util = string("

function ", oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"Utility(ins)

    # Definition of the inputs for the ODEs
    if length(ins) == 1
        ",join(ustr, ",")," = ins[1];
    else
        ",join(ustr, ",")," = ins;
    end
    inputsM1 = [",replace(join(inpstr1, ","), ",,"=>","),"];
    inputsM2 = [",replace(join(inpstr2, ","), ",,"=>","),"];

    # Solve ODEs

    solM1 = ",string(join(oedms_def["Model_1"]["NameF"]),"_SolveAll"),"(tsMS, pD1MS, ",spm1,", inputsM1, ivss1MS, sampsMS, pre1MS);
    solM2 = ",string(join(oedms_def["Model_2"]["NameF"]),"_SolveAll"),"(tsMS, pD2MS, ",spm2,", inputsM2, ivss2MS, sampsMS, pre2MS);

    # Extracte wanted vectors (observables) with time reduction
",join(obsM1),"
",join(obsM2),"
    ",join(ccomp),"

    return(util)
end

    ");

    fungen = join([Head, Util]);

    open(string(cudi, "\\Results\\", oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"_",
        today(), "\\OEDModelSelectionScripts\\",oedms_def["Model_1"]["NameF"], "_VS_",
        oedms_def["Model_2"]["NameF"],"_OEDMS.jl"), "w") do io
       write(io, fungen);
    end;

    println("")
    println("----------------------------------------- SCRIPTS OED Model Selection -----------------------------------------")
    println("Utility function script to perform OED for Model Selection has been generated in the directory: ")
    println(string(cudi, "\\Results\\", oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"_",
            today(), "\\OEDModelSelectionScripts\\",oedms_def["Model_1"]["NameF"], "_VS_",
            oedms_def["Model_2"]["NameF"],"_OEDMS.jl"))
    println("--------------------------------------------------------------------------------------")
    println("")

    include(string(cudi, "\\Results\\", oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"_",
            today(), "\\OEDModelSelectionScripts\\",oedms_def["Model_1"]["NameF"], "_VS_",
            oedms_def["Model_2"]["NameF"],"_OEDMS.jl"))

    println("")
    println("---------------------------------------------- OPTIMISATION INFO ----------------------------------------------")
    println("")
    println("If you wish to modify any of the settings for the Bayesian Optimisation used (from the package ")
    println("BayesianOptimization.jl) please search the function settingsBayesOpt located in the file ModelOEDSelection.jl ")
    println("that is located in the directory: ")
    println(string("        ", cudi2))
    println("and change any desired setting (under your own risk). ")
    println("If any irreversible change is made by mistake, just look for the file settingsBayesOptBackUp.jl and copy the")
    println("contents of the function inside (which are the same as the original except for the function name ;) )")
    println("")
    println("---------------------------------------------------------------------------------------------------------------")
    println("")

    return(oedms_def)
end


function plotOEDMSResults(oedms_res, oedms_def)

    # Utility convergence curve
    pu = plot(oedms_res["ResOptim"]["modelY"], label = "UtilityVal", legend=:topleft, xlabel = "Function Evaluations", ylabel = "Utility")
    pu = plot!(oedms_res["ConvCurv"][:,2], linetype=:step, label = "BestUtil", title = "OED Model Selection")
    savefig(string(oedms_def["savepath"], "\\Plot_OEDMSConvergence_", oedms_def["flag"], ".png"));

    # Manage inputs order
    inpnam = oedms_def["Model_1"]["inpName"];
    if oedms_def["Model_1"]["nInp"] < oedms_def["Model_2"]["nInp"]
        tmp1 = [];
        try
            for i in 1:oedms_def["Model_2"]["nInp"]
                if oedms_def["Model_2"]["inpName"][i] in oedms_def["Model_1"]["inpName"]
                    nothing
                else
                    tmp1 = vcat(tmp1, oedms_def["Model_2"]["inpName"][i])
                end
            end
            inpnam = vcat(inpnam, tmp1)
        catch
            for i in 1:oedms_def["Model_2"]["nInp"]
                if oedms_def["Model_2"]["inpName"][i] in oedms_def["Model_1"]["inpName"]
                    nothing
                else
                    tmp1 = hcat(tmp1, oedms_def["Model_2"]["inpName"][i])
                end
            end
            inpnam = hcat(inpnam, tmp1)
        end
    end;

    tit = "";
    yl1 = "";
    for k in 1:length(oedms_def["Obs"])
        tit = hcat(tit, string(oedms_def["Obs"][k]));
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
    tu = Array{Array{Any,1}}(undef, length(oedms_def["Obs"])+length(inpnam))
    [tu[k] = [] for k in 1:length(oedms_def["Obs"])]
    su = Array{Array{Any,1}}(undef, length(oedms_def["Obs"])+length(inpnam))
    [su[k] = [] for k in 1:length(oedms_def["Obs"])]
    for k in 1:length(inpnam)
        if oedms_def["fixedInp"] == []
            tu[k+length(oedms_def["Obs"])] = round.(oedms_def["switchT"]);
        else
            if inpnam[k] in oedms_def["fixedInp"]
                tu[k+length(oedms_def["Obs"])] = [round.(oedms_def["switchT"])[1], round.(oedms_def["switchT"])[end]];
            else
                tu[k+length(oedms_def["Obs"])] = round.(oedms_def["switchT"]);
            end
        end
        su[k+length(oedms_def["Obs"])] = vcat(oedms_res["uInpOpt"][inpnam[k]], oedms_res["uInpOpt"][inpnam[k]][end])
    end


    pl = plot(round.(oedms_def["tsamps"]), oedms_res["SimulObs_M1"][:,:,1],
            layout=length(oedms_def["Obs"])+length(inpnam), label = "Model 1", title = tit, xlabel = "time", ylabel = yuu,
                    color = "gray", size = [2000, 1200]);

    pl = plot!(round.(oedms_def["tsamps"]), oedms_res["SimulObs_M2"][:,:,1],
            layout=length(oedms_def["Obs"])+length(inpnam), label = "Model 2", title = tit, xlabel = "time", ylabel = yuu,
                    color = "red", size = [2000, 1200]);

    for samp in 2:size(oedms_res["SimulObs_M1"])[3]
        pl = plot!(round.(oedms_def["tsamps"]), oedms_res["SimulObs_M1"][:,:,samp],
            layout=length(oedms_def["Obs"])+length(inpnam), label = "",color = "gray", size = [2000, 1200]);

        pl = plot!(round.(oedms_def["tsamps"]), oedms_res["SimulObs_M2"][:,:,samp],
            layout=length(oedms_def["Obs"])+length(inpnam), label = "",color = "red", size = [2000, 1200]);
    end

    pl = plot!(tu, su, layout=length(oedms_def["Obs"])+length(inpnam), label = "", linetype=:step, title = tuu);

    savefig(string(oedms_def["savepath"], "\\PlotOEDMSResults_Exp1_", oedms_def["flag"], ".png"));

end

## Main function to run Bayesian OPTIMISATION
function settingsBayesOpt(oedms_def)

    # Manage inputs order
    inpnam = oedms_def["Model_1"]["inpName"];
    if oedms_def["Model_1"]["nInp"] < oedms_def["Model_2"]["nInp"]
        tmp1 = [];
        try
            for i in 1:oedms_def["Model_2"]["nInp"]
                if oedms_def["Model_2"]["inpName"][i] in oedms_def["Model_1"]["inpName"]
                    nothing
                else
                    tmp1 = vcat(tmp1, oedms_def["Model_2"]["inpName"][i])
                end
            end
            inpnam = vcat(inpnam, tmp1)
        catch
            for i in 1:oedms_def["Model_2"]["nInp"]
                if oedms_def["Model_2"]["inpName"][i] in oedms_def["Model_1"]["inpName"]
                    nothing
                else
                    tmp1 = hcat(tmp1, oedms_def["Model_2"]["inpName"][i])
                end
            end
            inpnam = hcat(inpnam, tmp1)
        end
    end;

    # Optimisation section
    steps = length(oedms_def["switchT"])-1;
    induc = maximum([oedms_def["Model_1"]["nInp"], oedms_def["Model_2"]["nInp"]]);
    fiii = length(oedms_def["fixedInp"]);

    doof = (steps*(induc-fiii))+fiii;

    # THE THIRD MOUNTAINS!
    uppe = Array{Any,1}(undef, doof); # Vector of upper bounds. Length set with doof but deppending on the user settings from fixedStep and equalStep this can change.
    lowe = Array{Any,1}(undef, doof); # Vector of lower bounds. Same ^
    coun = 1; # Index for the entry of ustr that will advance only if certain conditions are meet (that that position has been already checked and a variable or number has been introduced)
                # Tecnichally this indicated the first entry of uppe and lowwe for each inducer variables.
    fs = [oedms_def["fixedStep"][j][1] for j in 1:length(oedms_def["fixedStep"])]; # Fixed Steps (the actual step number)
    us = [oedms_def["fixedStep"][j][2] for j in 1:length(oedms_def["fixedStep"])]; # Input values (or Any entry) for each step specified in fixedStep
    fixin = 0; # Index to account for fixedInp. If the fixedInp is the first one in the inpam list this is the best way (that I have found) to account for it. This is because of the use of i (input index). This is because first I consider the case where there is a fisedInp and then check the rest.

    for i in 1:length(inpnam) # Loop over the number of inducers
        if inpnam[i] in oedms_def["fixedInp"] # If the current inducer is fixed then add only one bound and move index (coun). Fixin works as described in the other MOUNTAINS
            uppe[coun] = oedms_def["uUpper"][i];
            lowe[coun] = oedms_def["uLower"][i];
            coun += 1;
            fixin += 1;
        else
            try # Using coun you take all the indexes corresponding to the current inducer. The you fill it with the specified uper and lower bounds. Added the try just in case someoen runns it in a machien that does not support Float 64? (I am not sure)
                uppe[coun:coun+(steps-1)] = repeat([convert(Float64, oedms_def["uUpper"][i])], outer = steps);
                lowe[coun:coun+(steps-1)] = repeat([convert(Float64, oedms_def["uLower"][i])], outer = steps);
            catch
                uppe[coun:coun+(steps-1)] = repeat([convert(Float32, oedms_def["uUpper"][i])], outer = steps);
                lowe[coun:coun+(steps-1)] = repeat([convert(Float32, oedms_def["uLower"][i])], outer = steps);
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
                if oedms_def["equalStep"] != [] # If equalStep is empty then we do not need to do anything else, if it has something then there are entries to remove from the bounds
                    for k in 1:length(oedms_def["equalStep"]) # Loop over each entry of equalStep
                        tmpk = zeros(length(oedms_def["equalStep"][k])); # Temporal variable that will contain the same values as the equalStep entry except the ones that have been fixed for that step (this will happen if the user set an Any for the current inducer)
                        for h in 1:length(oedms_def["equalStep"][k]) # Loop over each value in the current entry
                            if oedms_def["equalStep"][k][h] in fs # Check if the current step considered is also in fixedStep
                                for w in 1:length(fs) # Account for ht euser not introducng steps sorted
                                    if oedms_def["equalStep"][k][h] == fs[w]
                                        if us[w][i-fixin] == Any # If the current step considered was defined as Any for the current inducer add the index, if it was fixed to a value ignore (this has been taken into account before)
                                            tmpk[h] = oedms_def["equalStep"][k][h];
                                        end
                                    end
                                end
                            else
                                tmpk[h] = oedms_def["equalStep"][k][h] # If the current step was not in fixedStep just add it
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

    sut = Symbol(string(oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"Utility"));
    opt = BOpt((@eval $sut),
          model,
          ExpectedImprovement(),                   # type of acquisition
          modeloptimizer,
          lowe,
          uppe,                                     # lowerbounds, upperbounds
          repetitions = 1,                          # evaluate the function for each input 1 times
          maxiterations = oedms_def["maxiter"],     # evaluate at 50 input positions
          sense = Max,                              # maximise the function
          verbosity = Progress);

    #################################### OPTIMISER SETINGS THAT CAN BE MODIFIED STOP HERE ####################################

    return opt

end

## Main function that runs the optimisation

function mainOEDMS(oedms_def)

    oedms_def = checkStructOEDMS(oedms_def);

    # Generate results directory
    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"_",today()))
    end

    try
        include(string(cudi, "\\Results\\", oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"_",
                today(), "\\OEDModelSelectionScripts\\",oedms_def["Model_1"]["NameF"], "_VS_",
                oedms_def["Model_2"]["NameF"],"_OEDMS.jl"))
    catch
        nothing
    end


    include(string(oedms_def["Model_1"]["modelpath"]));
    include(string(oedms_def["Model_2"]["modelpath"]));

    global tsMS = collect(0.0:round(oedms_def["finalTime"]))';
    global pD1MS = oedms_def["Theta_M1"];
    global pD2MS = oedms_def["Theta_M2"];
    global spMS = [convert(Int,v) for v in (round.(oedms_def["switchT"])')];
    global ivss1MS = oedms_def["y0_M1"];
    global ivss2MS = oedms_def["y0_M2"];
    global sampsMS = oedms_def["tsamps"];
    global pre1MS = oedms_def["preInd_M1"];
    global pre2MS = oedms_def["preInd_M2"];

    # Manage inputs order
    inpnam = oedms_def["Model_1"]["inpName"];
    if oedms_def["Model_1"]["nInp"] < oedms_def["Model_2"]["nInp"]
        tmp1 = [];
        try
            for i in 1:oedms_def["Model_2"]["nInp"]
                if oedms_def["Model_2"]["inpName"][i] in oedms_def["Model_1"]["inpName"]
                    nothing
                else
                    tmp1 = vcat(tmp1, oedms_def["Model_2"]["inpName"][i])
                end
            end
            inpnam = vcat(inpnam, tmp1)
        catch
            for i in 1:oedms_def["Model_2"]["nInp"]
                if oedms_def["Model_2"]["inpName"][i] in oedms_def["Model_1"]["inpName"]
                    nothing
                else
                    tmp1 = hcat(tmp1, oedms_def["Model_2"]["inpName"][i])
                end
            end
            inpnam = hcat(inpnam, tmp1)
        end
    end;

    # Optimisation section
    steps = length(spMS)-1;
    induc = maximum([oedms_def["Model_1"]["nInp"], oedms_def["Model_2"]["nInp"]]);
    fiii = length(oedms_def["fixedInp"]);

    doof = (steps*(induc-fiii))+fiii;


    opt = settingsBayesOpt(oedms_def);

    result = boptimize!(opt);

    # results

    # Extract results
    oedms_res = Dict()
    oedms_res["BestResOptim"] = result;
    oedms_res["ResOptim"] = Dict("acquisitionoptions"=>opt.acquisitionoptions,
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
    oedms_res["BestUtil"] = maximum(opt.model.y);

    convc = zeros(length(opt.model.y));
    convc[1] = opt.model.y[1];
    for i in 2:length(opt.model.y)
        if opt.model.y[i] > convc[i-1];
            convc[i] = opt.model.y[i];
        else
            convc[i] = convc[i-1];
        end
    end
    oedms_res["ConvCurv"] = hcat(collect(1:length(convc)), convc);

    # THE FOURTH MOUNTAINS!
    inducers = Dict(); # Library that will contain an entry for each of the inducers with its values ordered (including fixed values and repetitions)
    cnt2 = 1; # Index for the optimisation result (we move it by one once we have already extracted the current value for all the necessary inducers entries)
    fs = [oedms_def["fixedStep"][j][1] for j in 1:length(oedms_def["fixedStep"])]; # Fixed Steps (the actual step number)
    us = [oedms_def["fixedStep"][j][2] for j in 1:length(oedms_def["fixedStep"])]; # Input values (or Any entry) for each step specified in fixedStep
    fixin = 0; # Same use as in the other MOUNTAINS
    for i in 1:induc # Loop over each inducer
        if inpnam[i] in oedms_def["fixedInp"] # Check if the current inducer was fixed to a single value to be optimised
            inducers[inpnam[i]] = [result.observed_optimizer[cnt2]]; # Take the value
            fixin += 1; # Add this so there will not be bad indexing in fixedStep
            cnt2 += 1; # Move the index of the result vector to one since we already extracted this
        else
            allstp = Array{Any,1}(undef, steps); # Temporary variable that will have the inducer values for each step of the experiment.
            for j in 1:steps # Loop over the number of steps of the experiment
                if j in fs # Check if current step is in the list of fixedStep
                    for w in 1:length(fs) # Consider entries not sorted by step
                        if j==fs[w]
                            if oedms_def["fixedStep"][w][2][i-fixin] != Any # If inducer is in the fixedStep but set as Any (optimise), nothing. If it is fixed get the fixed value into the temporary variable allstp.
                                allstp[j] = oedms_def["fixedStep"][w][2][i-fixin];
                            end
                        end
                    end
                end
                if oedms_def["equalStep"] == [] # Case where equalStep is empty (no equal steps)
                    if typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{Any,1}},1} && # Case where there is no Any in the fixedStep case. Also check if the current step is nto assigned to a value. This can happen if there is also nothing in fixedStep. In this case just fill the allstp with the current value and move index.
                        typeof(oedms_def["fixedStep"]) != Array{Tuple{Int,Array{T,1} where T},1} &&
                        !isassigned(allstp, j)
                        allstp[j] = result.observed_optimizer[cnt2];
                        cnt2 += 1;
                    end
                else
                    for k in 1:length(oedms_def["equalStep"]) # Loop over each entry in equalStep
                        tmpk = zeros(length(oedms_def["equalStep"][k])); # Temporary variable to have same values as equalStep entry without the fixed steps
                        for h in 1:length(oedms_def["equalStep"][k]) # Loop over each step in the entry
                            if oedms_def["equalStep"][k][h] in fs # Check if current step is in fixedStep
                                for w in 1:length(fs) # Consider no sorted index in fixedStep
                                    if oedms_def["equalStep"][k][h] == fs[w]
                                        if us[w][i-fixin] == Any # If current inducer was not fixed (Any) then consider the step, if not no.
                                            tmpk[h] = oedms_def["equalStep"][k][h];
                                        end
                                    end
                                end
                            else
                                tmpk[h] = oedms_def["equalStep"][k][h] # Step not in fixedStep, so just consider it.
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

    oedms_res["uInpOpt"] = inducers;
    # Simulate!!!
    OEDsFunM1 = Symbol(string(join(oedms_def["Model_1"]["NameF"]),"_SolveAll"));
    OEDsFunM2 = Symbol(string(join(oedms_def["Model_2"]["NameF"]),"_SolveAll"));

    global inpuM1 = zeros(oedms_def["Model_1"]["nInp"]*(length(oedms_def["switchT"])-1));
    global inpuM2 = zeros(oedms_def["Model_2"]["nInp"]*(length(oedms_def["switchT"])-1));

    r1 = 1:(oedms_def["Model_1"]["nInp"]):length(inpuM1);
    r2 = 1:(oedms_def["Model_2"]["nInp"]):length(inpuM2);

    for i in 1:length(inpnam)
        if inpnam[i] in oedms_def["Model_1"]["inpName"]
            if inpnam[i] in oedms_def["fixedInp"]
                tmpin = oedms_res["uInpOpt"][inpnam[i]];
                inpuM1[r1.+(i-1)] .= tmpin;
            else
                tmpin = oedms_res["uInpOpt"][inpnam[i]];
                inpuM1[r1.+(i-1)] = tmpin;
            end
        end

        if inpnam[i] in oedms_def["Model_2"]["inpName"]
            if inpnam[i] in oedms_def["fixedInp"]
                tmpin = oedms_res["uInpOpt"][inpnam[i]];
                inpuM2[r2.+(i-1)] .= tmpin;
            else
                tmpin = oedms_res["uInpOpt"][inpnam[i]];
                inpuM2[r2.+(i-1)] = tmpin;
            end
        end
    end

    simulM1 = @eval $OEDsFunM1(tsMS, pD1MS, spMS, inpuM1, ivss1MS, sampsMS, pre1MS);
    simulM2 = @eval $OEDsFunM2(tsMS, pD2MS, spMS, inpuM2, ivss2MS, sampsMS, pre2MS);

    SimObsM1 = selectObsSim_te(simulM1, oedms_def["Obs"],oedms_def["Model_1"]["stName"]);
    SimObsM2 = selectObsSim_te(simulM2, oedms_def["Obs"],oedms_def["Model_2"]["stName"]);

    oedms_res["Simul_M1"] = simulM1;
    oedms_res["Simul_M2"] = simulM2;
    oedms_res["SimulObs_M1"] = SimObsM1;
    oedms_res["SimulObs_M2"] = SimObsM2;

    oedms_def["savepath"] = string(cudi, "\\Results\\", oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"_",today());
    oedms_def["savename"] = string("OEDModelSelectResults_",oedms_def["flag"],".jld");

    save(string(oedms_def["savepath"], "\\", oedms_def["savename"]), "oedms_res", oedms_res, "oedms_def", oedms_def);

    println("")
    println("----------------------------------------- RESULTS -----------------------------------------")
    println("OED for Model Selection results are saved in the directory: ")
    println(string("                 ", oedms_def["savepath"]))
    println(string("Under the name ",oedms_def["savename"]))
    println("--------------------------------------------------------------------------------------")
    println("")



    if oedms_def["plot"] == true
        plotOEDMSResults(oedms_res, oedms_def)
        println("")
        println("----------------------------------------- PLOTS -----------------------------------------")
        println("PLOTS are saved in the directory: ")
        println(string("                 ", oedms_def["savepath"]))
        println(string("Under the names PlotOEDMSResults_Exp1_",oedms_def["flag"],".png", " and Plot_OEDMSConvergence_", oedms_def["flag"], ".png"))
        println("--------------------------------------------------------------------------------------")
        println("")
    end

    return(oedms_res, oedms_def)

end
