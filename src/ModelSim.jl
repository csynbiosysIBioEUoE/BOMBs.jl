# using CSV
# using Dates
# using Plots
# using JLD

## Main simulation function

function simulateODEs(model_def, simul_def)

    # First make a folder where to save the results
    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
    end

    # Check which of the 2 types of simul_def have been given and make a choice on how to proceed
    if haskey(simul_def, "switchT")
        simul_def = checkStructSimul(model_def, simul_def)
    else
        simul_def = extractSimulCSV(model_def, simul_def)
    end

    # Make empty structure where the simulations will be
    simuls = Dict();

    # Get ODEs function
    ODEsFun = Symbol(string(join(model_def["NameF"]),"_SolveAll"));

    # loop over the different experiments to simulate
    for i in 1:simul_def["Nexp"]
        # Define simulation parameters
        global ts = collect(0.0:round(simul_def["finalTime"][i]))'; # 24 hour experiment
        global sp = [convert(Int,v) for v in (round.(simul_def["switchT"][i])')]; # 8 step experiment for now
        global ivss = simul_def["y0"][i]; # Experimental values for Y0 (IPTG, aTc, RFP, GFP)

        if model_def["Y0eqs"] != [] # ON inputs
            global pre = simul_def["preInd"][i];
        else
            global pre = [];
        end
        global samps = convert.(Int, round.(simul_def["tsamps"][i]));
        global theta = simul_def["theta"];

        # Inputs case
        if model_def["nInp"] != 0
            r = convert.(Int, 1:model_def["nInp"]:(length(simul_def["uInd"][i])));
            global inputs = zeros(convert.(Int,length(simul_def["uInd"][i])));
            for j in 1:convert.(Int,length(simul_def["uInd"][i])/model_def["nInp"])
                for k in 0:(model_def["nInp"]-1)
                    inputs[r[j]+k] = simul_def["uInd"][i][j,(k+1)];
                end
            end
        else
            r=0
        inputs=[];
        end

        simul = @eval $ODEsFun(ts, theta, sp, inputs, ivss, samps, pre);

        simuls[string("Exp_", i)] = simul;
    end

    simul_def["savepath"] = string(cudi, "\\Results\\", model_def["NameF"],"_",today());
    simul_def["savename"] = string(model_def["NameF"],"_",today(), "_SimulationResults_",simul_def["flag"],".jld");

    if simul_def["plot"] == true
        plotSimsODE(simuls,model_def,simul_def)
        println("")
        println("----------------------------------------- PLOTS -----------------------------------------")
        println("Simulation PLOTS are saved in the directory: ")
        println(string("                 ", simul_def["savepath"]))
        println(string("Under the name PlotSimulation_Exp(i)_",simul_def["flag"],".png"))
        println("--------------------------------------------------------------------------------------")
        println("")
    end

    save(string(simul_def["savepath"], "\\", simul_def["savename"]), "simuls", simuls, "model_def", model_def, "simul_def", simul_def)

    println("")
    println("----------------------------------------- RESULTS -----------------------------------------")
    println("Simulation results are saved in the directory: ")
    println(string("                 ", simul_def["savepath"]))
    println(string("Under the name ",simul_def["savename"]))
    println("--------------------------------------------------------------------------------------")
    println("")

    return simuls, model_def, simul_def

end

## Simulation Data structure Dict
function defSimulStruct()

    # The structure of things user should give
    simul_def = Dict()

    simul_def["Nexp"] = []; # Integer indicating the number of experiments to be simulated
    simul_def["finalTime"] = []; # -> Vector of final times for each simulation (initial time will allways be asumed as 0, so please consider that)
    simul_def["switchT"] = []; # Array with the switching times of the inducer in the simulation (time 0 and final time need to be considered)
    simul_def["y0"] = []; # Array of Y0s for the simulations for each experiment. If you are computing the steady state this vector might not be used, however you still need to introduce it
    simul_def["preInd"] = []; # Level of the inducer in the ON. This entry can be empty if no Steady State equations are given and only the Y0 given by the user is considered.
    simul_def["uInd"] = []; # Array with the values for the inducer for each experiment at each step
    simul_def["theta"] = []; # Vector/Matrix with the parameter samples or directory and file location of CSV file with them
    simul_def["tsamps"] = []; # Array of Sampling times vectors
    simul_def["plot"] = []; # Bollean or yes/no string to save the resulting simulations in the results directory (false will be considered as default)
    simul_def["flag"] = []; # String to attach a unique flag to the generated files so it is not overwritten. If empty, nothing will be added.

    return(simul_def)

end

# Check Simulations function
# Function to check the content of the fields the user has given
function checkStructSimul(model_def, simul_def)

    # Check taht all the dictionary entries are correct
    entries = ["Nexp", "finalTime", "switchT", "y0", "preInd", "uInd", "theta", "tsamps", "plot", "flag"]
    if symdiff(entries,keys(simul_def))!=[] && symdiff(entries,keys(simul_def)) != ["savepath", "savename"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is something wrong...")
        println(symdiff(entries,keys(simul_def)))
        return
    end

    # Check one to be sure there is no empty entry
    nee = ["Nexp", "finalTime", "switchT", "y0", "uInd", "theta", "tsamps"]; # No empty entries
    for i in 1:length(nee)
        if simul_def[nee[i]] == []
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Plese, check the field ", nee[i], ", this cannot be empty"))
            return
        end
    end

    # Checks to see if the content of the fields is the expected (also consier if single value entries have been given as a number or a vector)
    if ((typeof(simul_def["Nexp"])!=Array{Int,1}) || length(simul_def["Nexp"])!=1) && (typeof(simul_def["Nexp"])!=Int)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Nexp! This should be an integer. ")
        return
    elseif (typeof(simul_def["finalTime"]) != Array{Int,1}) && (typeof(simul_def["finalTime"]) != Array{Float64,1}) && (typeof(simul_def["finalTime"]) != Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field finalTime! This should be a vector of final times. ")
        return
    elseif (typeof(simul_def["switchT"]) != Array{Array{Int,1},1}) && (typeof(simul_def["switchT"]) != Array{Array{Float64,1},1}) && (typeof(simul_def["switchT"]) != Array{Array{Float32,1},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field switchT! This should be an array of arrays of numbers. ")
        return
    elseif (typeof(simul_def["y0"]) != Array{Array{Int,1},1}) && (typeof(simul_def["y0"]) != Array{Array{Float64,1},1}) && (typeof(simul_def["y0"]) != Array{Array{Float32,1},1}) &&
        (typeof(simul_def["y0"]) != Array{Array{Int,N} where N,1}) && (typeof(simul_def["y0"]) != Array{Array{Float64,N} where N,1}) && (typeof(simul_def["y0"]) != Array{Array{Float32,N} where N,1}) &&
        (typeof(simul_def["y0"]) != Array{Array{Int,2},1}) && (typeof(simul_def["y0"]) != Array{Array{Float64,2},1}) && (typeof(simul_def["y0"]) != Array{Array{Float32,2},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field y0! This should be an array of arrays of numbers. ")
        return
    elseif (typeof(simul_def["preInd"]) != Array{Array{Int,1},1}) && (typeof(simul_def["preInd"]) != Array{Array{Float64,1},1}) && (typeof(simul_def["preInd"]) != Array{Array{Float32,1},1}) && simul_def["preInd"] != []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field preInd! This should be an array of arrays of numbers. ")
        return
    # elseif (typeof(simul_def["uInd"]) != Array{Array{Int,1},1}) && (typeof(simul_def["uInd"]) != Array{Array{Float64,1},1}) && (typeof(simul_def["uInd"]) != Array{Array{Float32,1},1})
    #     println("-------------------------- Process STOPPED!!! --------------------------")
    #     println("Please, check the field uInd! This should be an array of arrays of numbers. ")
    #     return
    elseif (typeof(simul_def["theta"]) != Array{Float64,1}) && (typeof(simul_def["theta"]) != Array{Float32,1}) &&
        (typeof(simul_def["theta"]) != Array{Int,1}) && (typeof(simul_def["theta"]) != Array{Int,1}) &&
        (typeof(simul_def["theta"]) != Array{Float64,2}) && (typeof(simul_def["theta"]) != Array{Float32,2}) &&
        (typeof(simul_def["theta"]) != Array{String,1}) && (typeof(simul_def["theta"]) != String)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field theta! This should be a vector, matrix or string. ")
        return
    elseif (typeof(simul_def["tsamps"]) != Array{Array{Int,1},1}) && (typeof(simul_def["tsamps"]) != Array{Array{Float64,1},1}) && (typeof(simul_def["tsamps"]) != Array{Array{Float32,1},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field tsamps! This should be an array of arrays of numbers. ")
        return
    elseif (simul_def["plot"] != [])
        try
            if (typeof(simul_def["plot"]) != Bool) && (typeof(simul_def["plot"]) != Array{Bool,1}) &&
            (typeof(simul_def["plot"]) != String) && (typeof(simul_def["plot"]) != Array{String,1})
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check the field plot! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
                return
            end
        catch
             if (typeof(simul_def["plot"]) != Array{Bool,1}) &&
            (typeof(simul_def["plot"]) != Array{String,1})
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check the field plot! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
                return
            end
        end
    elseif (typeof(simul_def["flag"]) != Array{String,1}) && (typeof(simul_def["flag"]) != String) && ((simul_def["flag"]) != [])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field flag! This should be a vector, matrix or string. ")
    return

    end

    # Extract necessary elements to ease generalisation
    if typeof(simul_def["Nexp"]) == Array{Int,1}
        simul_def["Nexp"] = simul_def["Nexp"][1];
    end

    if simul_def["plot"]==[true] || simul_def["plot"]==["Yes"] || simul_def["plot"]==["yes"] || simul_def["plot"]=="Yes" || simul_def["plot"]=="yes"
        simul_def["plot"]=true
    elseif simul_def["plot"]==[false] || simul_def["plot"]==["No"] || simul_def["plot"]==["no"] || simul_def["plot"]=="No" || simul_def["plot"]=="no" || simul_def["plot"]==[]
        simul_def["plot"]=false
    end

    if (typeof(simul_def["theta"]) == Array{String,1})
        simul_def["theta"] = simul_def["theta"][1];
    end

    if typeof(simul_def["flag"]) == Array{String,1}
        simul_def["flag"] = simul_def["flag"][1];
    elseif typeof(simul_def["flag"]) == String
        simul_def["flag"] = simul_def["flag"];
    elseif ((simul_def["flag"]) == [])
        simul_def["flag"] = "";
    end


    # Check that all the contents make sense
    if (simul_def["Nexp"] != length(simul_def["finalTime"])) || (simul_def["Nexp"] != length(simul_def["switchT"])) ||
        (simul_def["Nexp"] != length(simul_def["y0"])) || (simul_def["Nexp"] != length(simul_def["uInd"])) ||
        (simul_def["Nexp"] != length(simul_def["tsamps"]))
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, some field does not match the number of experiments selected. Check that the number of")
        println("entries for finalTime, switchT, y0, uInd or tsamps matches the number of experiments in Nexp.")
        return
    end


    if (model_def["Y0eqs"] != []) && simul_def["preInd"]==[]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("You specified computation of Y0 as steady-state but you have not introduced an inducer value for it!")
        return
    end


    if (typeof(simul_def["theta"]) == Array{Float64,1}) || (typeof(simul_def["theta"]) == Array{Float32,1})
        if length(simul_def["theta"]) != model_def["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Number of parameters introduced does not match the specified")
            return
        end
    elseif (typeof(simul_def["theta"]) == Array{Float64,2}) || (typeof(simul_def["theta"]) == Array{Float32,2})
        if size(simul_def["theta"])[1] != model_def["nPar"] && size(simul_def["theta"])[2] != model_def["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Number of parameters introduced does not match the specified")
            return
        end
    elseif (typeof(simul_def["theta"]) == String)
        if simul_def["theta"][end-3:end] != ".csv"
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, add the .csv termination to the theta file! And no spaces after!")
            return
        end
        if isfile(simul_def["theta"])
            simul_def["theta"] = Matrix(CSV.read(simul_def["theta"]));
        else
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but the file path you introduced does not exist!")
            return
        end
    end


    for i in 1:simul_def["Nexp"]
        if simul_def["tsamps"][i][end] > simul_def["finalTime"][i]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check finalTime. You have selected a sampling point past it.")
            return
        end
        if model_def["nInp"] != 0
            if size(simul_def["uInd"][i])[1] != (length(simul_def["switchT"][i])-1)
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check uInd and switchT. Number of steps does not match the number of values for the inputs.")
                return
            end
        end

        if length(simul_def["y0"][i])/model_def["nStat"] == 1
            if length(simul_def["y0"][i]) != model_def["nStat"]
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but Y0 does not match the number or states in experiment ", i))
                return
            end
        else
            if size(simul_def["y0"][i])[1] != model_def["nStat"] && size(simul_def["y0"][i])[2] != model_def["nStat"]
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but Y0 does not match the number or states in experiment ", i))
                return
            end
            if length(simul_def["y0"][i])/model_def["nStat"] != length(simul_def["theta"])/model_def["nPar"]
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but Y0 does not match the dimensions of theta "))
                return
            end
            if size(simul_def["y0"][i])[1] == size(simul_def["y0"][i])[2]
                println("-------------------------- WARNING --------------------------")
                println(string("Sorry, but the number of rows and columns of the y0 matrix in experiment ",i," is the same, so the checks on "))
                println("correct orientation will not work. Please make sure that the dimensions follow: ")
                println("y0[samples, states]")
            end
        end
    end

    # Check warning in case the matrix introduced is symmetric
    if length(size(simul_def["theta"])) == 2
        if size(simul_def["theta"])[1] == size(simul_def["theta"])[2]
            println("-------------------------- WARNING --------------------------")
            println("Sorry, but the number of rows and columns of the theta matrix is the same, so the checks on ")
            println("correct orientation will not work. Please make sure that the dimensions follow: ")
            println("theta[samples, parameters]")
        end
    end

    # Check that no step is no smaller than 2 unit of time
    for i in 1:simul_def["Nexp"]
        for j in 1:size(simul_def["uInd"][i])[2]
            if (simul_def["switchT"][i][j+1]-simul_def["switchT"][i][j])<=4
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
            for i in 1:simul_def["Nexp"]
                if size(simul_def["uInd"][i])[2] != model_def["nInp"];
                    simul_def["uInd"][i] = simul_def["uInd"][i]';
                end
            end
        end
    catch
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but there is some issue with the contents of the field uInd."))
        return
    end

    return(simul_def)
end

## CSV introduction description
function fileStructInfo()

    println("ALL THE FILES MUST BE CSV FILES!!!")
    println("")
    println("---------------------------------------------------------------------------------------------------------")
    println("")
    println("The observables file entry should have the following structure: ")
    println("Column 1: Sampling Times for the simulations")
    println("Column 2 + Number of states: Y0 value for each state (in order) in each column. ")
    println("           Value might need to be repeated across all the column, but only the first row will be considered")
    println("")
    println("---------------------------------------------------------------------------------------------------------")
    println("")
    println("The event inputs file entry should have the following structure: ")
    println("Column 1: Each Switching time (final time not considered)")
    println("Column 2: Final time (can be repeated for as many entrances as column 1 has)")
    println("Column 3 + number of inducers: Repeated column for the value of the inducers in the ON. ")
    println("              If no pre-inducer will be used in the simulations, any number could be used")
    println("Column 4 + number of inducers: Value of each inducer at each switching time.")
    println("")
    println("This is an Example of an experiment with 3 steps and 2 inducers: ")
    println(" |Switching|  |FinalTime|  |IPTGpre|  |aTcPre|  |  IPTG  |  |  aTc  |")
    println("      0           1350         1          0          1          0    ")
    println("     400          1350         1          0         0.3         50   ")
    println("     700          1350         1          0         0.7         35   ")
    println("     950          1350         1          0         0.1        100   ")
    println("")
    println("---------------------------------------------------------------------------------------------------------")
    println("")
    println("For the main directory, remember that in Julia you have to type \\ twice!")

end

## CSV introduction structure
function defSimulStructFiles()

    fileStructInfo()

    # The structure of things user should give
    simul_def = Dict()

    simul_def["ObservablesFile"] = []; # File containing information about the states Y0 and sampling times
    simul_def["EventInputsFile"] = []; # -> File containing all the information about the stimuli
    simul_def["theta"] = []; # Vector/Matrix with the parameter samples or directory and file location of CSV file with them
    simul_def["MainDir"] = []; # Directory path where all files are located (to avoid repetition). This can be empty.
    simul_def["plot"] = [];
    simul_def["flag"] = [];

    return(simul_def)

end

## Extract data from CSV
# Function to extract contents of CSV files if user sellects that
function extractSimulCSV(model_def, simul_def)

    # Check taht all the dictionary entries are correct
    entries = ["ObservablesFile", "EventInputsFile", "theta", "MainDir", "plot", "flag"]
    if symdiff(entries,keys(simul_def))!=[]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is something wrong...")
        println(symdiff(entries,keys(simul_def)))
        return
    end

    # Checks on the different strings introduced
    if simul_def["MainDir"] == [] || simul_def["MainDir"] == ""
        simul_def["MainDir"] = "";
    elseif typeof(simul_def["MainDir"])==Array{String,1}
        simul_def["MainDir"] = simul_def["MainDir"][1]
    elseif typeof(simul_def["MainDir"]) == String
        nothing
    else
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, there is something wrong in the type of the field MainDir")
        return
    end

    if typeof(simul_def["ObservablesFile"]) != Array{String,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, there is something wrong in the type of the field ObservablesFile. This should be a Array{String,1}")
        return
    end

    if length(simul_def["ObservablesFile"]) != length(simul_def["EventInputsFile"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, there number of entries for ObservablesFile and EventInputsFile should be the same")
        return
    end

    for i in 1:length(simul_def["ObservablesFile"])
        if simul_def["ObservablesFile"][i][end-3:end] != ".csv"
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, add the .csv termination to the ObservablesFile file! And no spaces after!")
            return
        end
    end


    if typeof(simul_def["EventInputsFile"]) != Array{String,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, there is something wrong in the type of the field EventInputsFile. This should be a Array{String,1}")
        return
    end

    for i in 1:length(simul_def["EventInputsFile"])
        if simul_def["EventInputsFile"][i][end-3:end] != ".csv"
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, add the .csv termination to the EventInputsFile file! And no spaces after!")
            return
        end
    end



    simstu = defSimulStruct();

    simstu["theta"] = simul_def["theta"];

    tmpInp = Array{Any}(undef, 1,length(simul_def["EventInputsFile"]))
    for i in 1:length(simul_def["EventInputsFile"])
        if simul_def["MainDir"] != [] && simul_def["MainDir"] != ""
            tmpInp[i] = Matrix(CSV.read(string(simul_def["MainDir"], "\\\\",simul_def["EventInputsFile"][i])));
        else
            tmpInp[i] = Matrix(CSV.read(simul_def["EventInputsFile"][i]));
        end
    end


    tmpObs = Array{Any}(undef, 1,length(simul_def["ObservablesFile"]))
    for i in 1:length(simul_def["ObservablesFile"])
        if simul_def["MainDir"] != [] && simul_def["MainDir"] != ""
            tmpObs[i] = Matrix(CSV.read(string(simul_def["MainDir"], "\\\\",simul_def["ObservablesFile"][i])));
        else
            tmpObs[i] = Matrix(CSV.read(simul_def["ObservablesFile"][i]));
        end
    end

    try
        simstu["Nexp"] = length(simul_def["EventInputsFile"]);
        simstu["finalTime"] = [tmpInp[i][1,2] for i in 1:length(simul_def["EventInputsFile"])];
        simstu["switchT"] = [vcat(tmpInp[i][:,1], tmpInp[i][1,2]) for i in 1:length(simul_def["EventInputsFile"])];
        simstu["y0"] = [tmpObs[i][1,2:(1+model_def["nStat"])] for i in 1:length(simul_def["ObservablesFile"])];
        simstu["preInd"] = [tmpInp[i][1,3:(2+model_def["nInp"])] for i in 1:length(simul_def["EventInputsFile"])];
        simstu["uInd"] = [tmpInp[i][:,(3+model_def["nInp"]):(3+model_def["nInp"]+(model_def["nInp"]-1))] for i in 1:length(simul_def["EventInputsFile"])];
        simstu["theta"] = simul_def["theta"];
        simstu["tsamps"] = [tmpObs[i][:,1] for i in 1:length(simul_def["ObservablesFile"])];
        simstu["plot"] = simul_def["plot"];
        simstu["flag"] = simul_def["flag"];
    catch
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, but there seems to be some issue with the internal structure of the files you have given...")
        return
    end

    simstu = checkStructSimul(model_def, simstu);

    return(simstu)

end

## Plot simulations function
function plotSimsODE(simuls,model_def,simul_def)
    cudi = pwd(); # @__DIR__;
    for i in 1:simul_def["Nexp"]
        # Title for each subplot
        tit = "";
        yl1 = "";
        for k in 1:model_def["nStat"]
            tit = hcat(tit, string(model_def["stName"][k]));
            yl1 = hcat(yl1, "y");
        end
        tit = tit[:,2:end];
        yl1 = yl1[:,2:end];

        titu = "";
        yl2 = "";
        if model_def["nInp"] != 0
            for k in 1:model_def["nInp"]
                titu = hcat(titu, string(model_def["inpName"][k]));
                yl2 = hcat(yl2, "u");
            end
            titu = titu[:,2:end];
            yl2 = yl2[:,2:end];
        end

        tuu = hcat(tit, titu);
        yuu = hcat(yl1, yl2);

        # Elements for the inducers
        tu = Array{Array{Any,1}}(undef, model_def["nStat"]+model_def["nInp"])
        [tu[k] = [] for k in 1:model_def["nStat"]]
        su = Array{Array{Any,1}}(undef, model_def["nStat"]+model_def["nInp"])
        [su[k] = [] for k in 1:model_def["nStat"]]
        for k in 1:model_def["nInp"]
            tu[k+model_def["nStat"]] = round.(simul_def["switchT"][i]);
            su[k+model_def["nStat"]] = vcat(simul_def["uInd"][i][:,(k)], simul_def["uInd"][i][end,(k)])
        end

        t = round.(simul_def["tsamps"][i]);
        s = simuls[string("Exp_",i)];

        pl = plot(t, s[:,:,1], layout=model_def["nStat"]+model_def["nInp"], label = "", title = tit, xlabel = "time", ylabel = yuu,
        color = "blue", size = [2000, 1200])

        if model_def["nStat"]>1
            for p in 2:size(simuls[string("Exp_",i)])[3]
                pl = plot!(t,s[:,:,p], layout=model_def["nStat"], label = "", color = "blue");
            end
        end

        pl = plot!(tu, su, layout=model_def["nStat"]+model_def["nInp"], label = "", linetype=:step, title = tuu);

        savefig(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\PlotSimulation_Exp", i,"_", simul_def["flag"], ".png"));

    end

end
