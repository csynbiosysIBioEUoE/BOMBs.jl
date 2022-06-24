
## Main function to generate pseudo-Data

function GenPseudoDat(model_def, pseudo_def)

    # First make a folder where to save the results
    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
    end

    # Check which of the 2 types of simul_def have been given and make a choice on how to proceed
    if haskey(pseudo_def, "switchT")
        pseudo_def = checkStructPseudoDat(model_def, pseudo_def)
    else
        pseudo_def = extractPseudoDatCSV(model_def, pseudo_def)
    end

    simul_def = defSimulStruct();

    for (i,j) in simul_def
        simul_def[i] = pseudo_def[i];
    end
    # Make empty structure where the simulations will be
    simuls = Dict();

    # Get ODEs function
    OEDsFun = Symbol(string(join(model_def["NameF"]),"_SolveAll"));

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

        simul = @eval $OEDsFun(ts, theta, sp, inputs, ivss, samps, pre);

        simuls[string("Exp_", i)] = simul;
    end

    # Now extract the entries to which the pseudodata will be generated
    pseudoDme = Dict();
    pseudoDsd = Dict();
    simObs = Dict();
    for i in 1:simul_def["Nexp"]
        pseudoDme[string("PDExp_", i)] = zeros(size(simuls[string("Exp_", i)])[1], length(pseudo_def["Obs"]), size(simuls[string("Exp_", 1)])[3]);
        pseudoDsd[string("PDExp_", i)] = zeros(size(simuls[string("Exp_", i)])[1], length(pseudo_def["Obs"]), size(simuls[string("Exp_", 1)])[3]);
        simObs[string("PDExp_", i)] = zeros(size(simuls[string("Exp_", i)])[1], length(pseudo_def["Obs"]), size(simuls[string("Exp_", 1)])[3]);
        if typeof(pseudo_def["Obs"])==Array{Int,1}
            simObs[string("PDExp_", i)] = simuls[string("Exp_", i)][:,(pseudo_def["Obs"]),:];
        else
            coun = 1;
            for j in 1:length(pseudo_def["Obs"])
               if !(occursin.("+", pseudo_def["Obs"][j])) && !(occursin.("-", pseudo_def["Obs"][j])) &&
               !(occursin.("/", pseudo_def["Obs"][j])) && !(occursin.("*", pseudo_def["Obs"][j])) &&
                !(occursin.("^", pseudo_def["Obs"][j]))
                    tind = findall(x->x==pseudo_def["Obs"][j], model_def["stName"])[1];
                    simObs[string("PDExp_", i)][:,coun,:] = simuls[string("Exp_", i)][:,tind,:];
                    coun +=1;
                else #This is the complicated stuff where we need to evaluate the function written as a string
                    if (occursin(".+", pseudo_def["Obs"][j]))
                        pseudo_def["Obs"][j] = replace(pseudo_def["Obs"][j], ".+"=>" .+ ")
                    elseif (occursin("+", pseudo_def["Obs"][j]))
                        pseudo_def["Obs"][j] = replace(pseudo_def["Obs"][j], "+"=>" .+ ")
                    end

                    if (occursin(".-", pseudo_def["Obs"][j]))
                        pseudo_def["Obs"][j] = replace(pseudo_def["Obs"][j], ".-"=>" .- ")
                    elseif (occursin("-", pseudo_def["Obs"][j]))
                        pseudo_def["Obs"][j] = replace(pseudo_def["Obs"][j], "-"=>" .- ")
                    end

                    if (occursin("./", pseudo_def["Obs"][j]))
                        pseudo_def["Obs"][j] = replace(pseudo_def["Obs"][j], "./"=>" ./ ")
                    elseif (occursin("/", pseudo_def["Obs"][j]))
                        pseudo_def["Obs"][j] = replace(pseudo_def["Obs"][j], "/"=>" ./ ")
                    end

                    if (occursin(".*", pseudo_def["Obs"][j]))
                        pseudo_def["Obs"][j] = replace(pseudo_def["Obs"][j], ".*"=>" .* ")
                    elseif (occursin("*", pseudo_def["Obs"][j]))
                        pseudo_def["Obs"][j] = replace(pseudo_def["Obs"][j], "*"=>" .* ")
                    end

                    if (occursin(".^", pseudo_def["Obs"][j]))
                        pseudo_def["Obs"][j] = replace(pseudo_def["Obs"][j], ".^"=>" .^ ")
                    elseif (occursin("^", pseudo_def["Obs"][j]))
                        pseudo_def["Obs"][j] = replace(pseudo_def["Obs"][j], "^"=>" .^ ")
                    end

                    tmp1 = pseudo_def["Obs"][j];
                    for k in 1:model_def["nStat"]
                        if occursin(model_def["stName"][k], pseudo_def["Obs"][j])
                            tmp1 = replace(tmp1, model_def["stName"][k]=>string(simuls[string("Exp_",i)][:,k,:]))
                        end
                    end
                    simObs[string("PDExp_", i)][:,coun,:] = eval(Meta.parse(tmp1));
                    coun +=1;
                end
            end
        end
    end

    # Generate pseudo-data with noise

    for j in 1:simul_def["Nexp"]
        for stat in 1:length(pseudo_def["Obs"])
            if length(size(pseudo_def["theta"]))==2
                for p in 1:convert(Int, length(pseudo_def["theta"])/model_def["nPar"])
                    if pseudo_def["NoiseType"] == "hetero"
                        pseudoDsd[string("PDExp_", j)][:,stat,p] = simObs[string("PDExp_", j)][:,stat,p] .* pseudo_def["Noise"][stat];
                    elseif pseudo_def["NoiseType"] == "homo"
                        pseudoDsd[string("PDExp_", j)][:,stat,p] = ones(length(simObs[string("PDExp_", j)][:,stat,p])) .* pseudo_def["Noise"][stat];
                    end
                    pseudoDme[string("PDExp_", j)][:,stat,p] = simObs[string("PDExp_", j)][:,stat,p] .+ rand(MvNormal(zeros(size(simObs[string("PDExp_", j)][:,stat,p])[1]),pseudoDsd[string("PDExp_", j)][:,stat,p][:]))
                    if pseudoDsd[string("PDExp_", j)][1,stat,p] == 0 # This checks for the case wher ethe user gives y0 and this is 0. This will output a mean and standard deviation of 0 which will generate issues in lekelihood computation (MLE and Stan inference)
                        pseudoDsd[string("PDExp_", j)][1,stat,p] = 1e-20;
                    end
                    if pseudoDme[string("PDExp_", j)][1,stat,p] == 0
                        pseudoDme[string("PDExp_", j)][1,stat,p] = 1e-10;
                    end
                    for r in 1:size(pseudoDme[string("PDExp_", j)])[1]
                        if pseudoDme[string("PDExp_", j)][r,stat,p] == 0 # If you have homoscedastic noise and high, this can generate negative values which do not make much sense. A proxi to a truncation of the distribution I suppsose.
                            pseudoDme[string("PDExp_", j)][r,stat,p] = 1e-10;
                        end
                    end
                end
            else
                    p=1;
                    if pseudo_def["NoiseType"] == "hetero"
                        pseudoDsd[string("PDExp_", j)][:,stat,p] = simObs[string("PDExp_", j)][:,stat,p] .* pseudo_def["Noise"][stat];
                    elseif pseudo_def["NoiseType"] == "homo"
                        pseudoDsd[string("PDExp_", j)][:,stat,p] = ones(length(simObs[string("PDExp_", j)][:,stat,p])) .* pseudo_def["Noise"][stat];
                    end
                    pseudoDme[string("PDExp_", j)][:,stat,p] = simObs[string("PDExp_", j)][:,stat,p] .+ rand(MvNormal(zeros(size(simObs[string("PDExp_", j)][:,stat,p])[1]),pseudoDsd[string("PDExp_", j)][:,stat,p][:]))
                    if pseudoDsd[string("PDExp_", j)][1,stat,p] == 0
                        pseudoDsd[string("PDExp_", j)][1,stat,p] = 1e-20;
                    end
                    if pseudoDme[string("PDExp_", j)][1,stat,p] == 0
                        pseudoDme[string("PDExp_", j)][1,stat,p] = 1e-10;
                    end

                    for r in 1:size(pseudoDme[string("PDExp_", j)])[1]
                        if pseudoDme[string("PDExp_", j)][r,stat,p] == 0 # If you have homoscedastic noise and high, this can generate negative values which do not make much sense. A proxi to a truncation of the distribution I suppsose.
                            pseudoDme[string("PDExp_", j)][r,stat,p] = 1e-10;
                        end
                    end
            end
        end
    end
    # Save path
    pseudo_def["savepath"] = string(cudi, "\\Results\\", model_def["NameF"],"_",today());
    pseudo_def["savename"] = string(model_def["NameF"],"_",today(), "_PseudoDataResults_",pseudo_def["flag"],".jld");

    # Put all results into one same structure and save
    pseudo_res = Dict()
    pseudo_res["Sims"] = simuls;
    pseudo_res["SimsObs"] = simObs;
    pseudo_res["PData"] = pseudoDme;
    pseudo_res["PError"] = pseudoDsd;

    JLD.save(string(pseudo_def["savepath"], "\\", pseudo_def["savename"]), "PseudoData", pseudo_res, "model_def", model_def, "pseudo_def", pseudo_def)

    println("")
    println("----------------------------------------- RESULTS -----------------------------------------")
    println("Pseudo-Data results are saved in the directory: ")
    println(string("                 ", pseudo_def["savepath"]))
    println(string("Under the name ",pseudo_def["savename"]))
    println("--------------------------------------------------------------------------------------")
    println("")

    if simul_def["plot"] == true
        plotPseudoDatODE(pseudo_res,model_def,pseudo_def)
        println("")
        println("----------------------------------------- PLOTS -----------------------------------------")
        println("Pseudo-Data PLOTS are saved in the directory: ")
        println(string("                 ", pseudo_def["savepath"]))
        println(string("Under the name PlotPseudoDat_Exp(i)_",pseudo_def["flag"],".png"))
        println(string("If more than one sample from the parameters are used, the plot will be the average between all traces."))
        println("--------------------------------------------------------------------------------------")
        println("")
    end

    PDatCSVGen(pseudo_res,model_def,pseudo_def);
    println("")
    println("----------------------------------------- CSVs -----------------------------------------")
    println("Pseudo-Data CSVs with results have been generated in the directory: ")
    println(string("                 ", pseudo_def["savepath"], "\\PseudoDataFiles"))
    println(string("Under the names: "))
    println(string("        Simulations: ", model_def["NameF"],"_EXP(i)", "_",
                    pseudo_def["flag"],"_Simulations.csv"))
    println(string("        Observables: ", model_def["NameF"],"_EXP(i)", "_",
                    pseudo_def["flag"],"_Observables.csv"))
    println(string("        Event Inputs (only if the model has inputs): ", model_def["NameF"],"_EXP(i)", "_",
                    pseudo_def["flag"],"_Event_Inputs.csv"))
    println("--------------------------------------------------------------------------------------")
    println("")

    for i in 1:pseudo_def["Nexp"]
        pseudo_def["uInd"][i] = Array(pseudo_def["uInd"][i]');
    end

    return pseudo_res, model_def, pseudo_def

end

## Pseudodata general plot function

function plotPseudoDatODE(pseudo_res,model_def,pseudo_def)
    cudi = pwd(); # @__DIR__;
    for i in 1:pseudo_def["Nexp"]
        # Title for each subplot
        tit = "";
        yl1 = "";
        for k in 1:length(pseudo_def["Obs"])
            if typeof(pseudo_def["Obs"]) == Array{Int,1};
                tit = hcat(tit, string(model_def["stName"][pseudo_def["Obs"][k]]));
            else
                tit = hcat(tit, string(pseudo_def["Obs"][k]));
            end

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
        tu = Array{Array{Any,1}}(undef, length(pseudo_def["Obs"])+model_def["nInp"])
        [tu[k] = [] for k in 1:length(pseudo_def["Obs"])]
        su = Array{Array{Any,1}}(undef, length(pseudo_def["Obs"])+model_def["nInp"])
        [su[k] = [] for k in 1:length(pseudo_def["Obs"])]
        for k in 1:model_def["nInp"]
            tu[k+length(pseudo_def["Obs"])] = round.(pseudo_def["switchT"][i]);
            su[k+length(pseudo_def["Obs"])] = vcat(pseudo_def["uInd"][i][:,(k)], pseudo_def["uInd"][i][end,(k)])
        end

        t = round.(pseudo_def["tsamps"][i]);
        m = pseudo_res["PData"][string("PDExp_",i)];
        s = pseudo_res["PError"][string("PDExp_",i)];

        if length(size(pseudo_def["theta"])) == 1
            pl = plot(t, m[:,:,1],yerror = s[:,:,1], layout=length(pseudo_def["Obs"])+model_def["nInp"], label = "", title = tit, xlabel = "time", ylabel = yuu,
            color = "blue", size = [2000, 1200])

            pl = plot!(tu, su, layout=model_def["nStat"]+model_def["nInp"], label = "", linetype=:step, title = tuu);
        else
            m2 = mean(m, dims = 3);
            s2 = std(m, dims = 3);
            pl = plot(t, m2[:,:,1],yerror = s2[:,:,1], layout=length(pseudo_def["Obs"])+model_def["nInp"], label = "Average", title = tit, xlabel = "time", ylabel = yuu,
            color = "blue", size = [2000, 1200])

            pl = plot!(tu, su, layout=model_def["nStat"]+model_def["nInp"], label = "", linetype=:step, title = tuu);
        end

        savefig(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\PlotPseudoDat_Exp", i,"_", pseudo_def["flag"], ".png"));

    end

end


## Simulation Data structure Dict
function defPseudoDatStruct()

    # The structure of things user should give
    pseudo_def = Dict()

    pseudo_def["Nexp"] = []; # Integer indicating the number of experiments to be simulated
    pseudo_def["finalTime"] = []; # -> Vector of final times for each simulation (initial time will allways be asumed as 0, so please consider that)
    pseudo_def["switchT"] = []; # Array with the switching times of the inducer in the simulation (time 0 and final time need to be considered)
    pseudo_def["y0"] = []; # Array of Y0s for the simulations for each experiment. If you are computing the steady state this vector might not be used, however you still need to introduce it
    pseudo_def["preInd"] = []; # Level of the inducer in the ON. This entry can be empty if no Steady State equations are given and only the Y0 given by the user is considered.
    pseudo_def["uInd"] = []; # Array with the values for the inducer for each experiment at each step
    pseudo_def["theta"] = []; # Vector/Matrix with the parameter samples or directory and file location of CSV file with them. If a matrix is given, pseudo-data will be generated for each one of the theta samples
    pseudo_def["tsamps"] = []; # Array of Sampling times vectors
    pseudo_def["plot"] = []; # Bollean or yes/no string to save the resulting simulations in the results directory (false will be considered as default)
    pseudo_def["flag"] = []; # String to attach a unique flag to the generated files so it is not overwritten. If empty, nothing will be added.

    pseudo_def["Obs"] = []; # States that are observable. This is either a vecotr of strings (that has the same order asmodel_def["stName"] and model_def["eqns"]), a vector of integers insidating which states are observable or the string "All" if all states are observables. This can also be an experssion combining states (Only +,-,*,\ and ^ will be considered)
    pseudo_def["NoiseType"] = []; # String indicating if the desired noise is homoscedastic or heteroscedastic. The allowed strings are homo, homoscedastic or hetero, heteroscedastic (case independent). If empty, default will be heteroscedastic. Note that for high homoscedastic noise data could become negative (which to me does not make much sense), for now this will be set to a value close to 0.
    pseudo_def["Noise"] = []; # Percentage of heteroscedastic noise (introduced as value from 0 to 1). If empty 10% heteroscedastic noise will be assumed. In the homoscedastic case, if empty a standard deviation of 1 will be assumes. This has to be a vector of noise values for each observable.

    return(pseudo_def)

end

## Check pseudo-data structure function
# Function to check the content of the fields the user has given
function checkStructPseudoDat(model_def, pseudo_def)

    # Check taht all the dictionary entries are correct
    entries = ["Nexp", "finalTime", "switchT", "y0", "preInd", "uInd", "theta", "tsamps", "plot", "flag", "Obs", "NoiseType", "Noise"]
    if symdiff(entries,keys(pseudo_def))!=[] && symdiff(entries,keys(pseudo_def)) != ["savepath", "savename"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is something wrong...")
        println(symdiff(entries,keys(pseudo_def)))
        return
    end

    # Check one to be sure there is no empty entry
    nee = ["Nexp", "finalTime", "switchT", "y0", "uInd", "theta", "tsamps", "Obs"]; # No empty entries
    for i in 1:length(nee)
        if pseudo_def[nee[i]] == []
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Plese, check the field ", nee[i], ", this cannot be empty"))
            return
        end
    end

    # Checks to see if the content of the fields is the expected (also consier if single value entries have been given as a number or a vector)
    if ((typeof(pseudo_def["Nexp"][1])!=Int) || length(pseudo_def["Nexp"])!=1) && (typeof(pseudo_def["Nexp"])!=Int)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Nexp! This should be an integer. ")
        return
    elseif (typeof(pseudo_def["finalTime"]) != Array{Int,1}) && (typeof(pseudo_def["finalTime"]) != Array{Float64,1}) && (typeof(pseudo_def["finalTime"]) != Array{Int32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field finalTime! This should be a vector of final times. ")
        return
    elseif (typeof(pseudo_def["switchT"]) != Array{Array{Int,1},1}) && (typeof(pseudo_def["switchT"]) != Array{Array{Float64,1},1}) && (typeof(pseudo_def["switchT"]) != Array{Array{Float32,1},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field switchT! This should be an array of arrays of numbers. ")
        return
    elseif (typeof(pseudo_def["y0"]) != Array{Array{Int,1},1}) && (typeof(pseudo_def["y0"]) != Array{Array{Float64,1},1}) && (typeof(pseudo_def["y0"]) != Array{Array{Float32,1},1}) &&
        (typeof(pseudo_def["y0"]) != Array{Array{Int,N} where N,1}) && (typeof(pseudo_def["y0"]) != Array{Array{Float64,N} where N,1}) && (typeof(pseudo_def["y0"]) != Array{Array{Float32,N} where N,1}) &&
        (typeof(pseudo_def["y0"]) != Array{Array{Int,2},1}) && (typeof(pseudo_def["y0"]) != Array{Array{Float64,2},1}) && (typeof(pseudo_def["y0"]) != Array{Array{Float32,2},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field y0! This should be an array of arrays of numbers. ")
        return
    elseif (typeof(pseudo_def["preInd"]) != Array{Array{Int,1},1}) && (typeof(pseudo_def["preInd"]) != Array{Array{Float64,1},1}) && (typeof(pseudo_def["preInd"]) != Array{Array{Float32,1},1}) && pseudo_def["preInd"] != []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field preInd! This should be an array of arrays of numbers. ")
        return
    # elseif (typeof(pseudo_def["uInd"]) != Array{Array{Int,1},1}) && (typeof(pseudo_def["uInd"]) != Array{Array{Float64,1},1}) && (typeof(pseudo_def["uInd"]) != Array{Array{Float32,1},1})
    #     println("-------------------------- Process STOPPED!!! --------------------------")
    #     println("Please, check the field uInd! This should be an array of arrays of numbers. ")
    #     return
    elseif (typeof(pseudo_def["theta"]) != Array{Float64,1}) && (typeof(pseudo_def["theta"]) != Array{Float32,1}) &&
        (typeof(pseudo_def["theta"]) != Array{Float64,2}) && (typeof(pseudo_def["theta"]) != Array{Float32,2}) &&
        (typeof(pseudo_def["theta"]) != Array{String,1}) && (typeof(pseudo_def["theta"]) != String)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field theta! This should be a vector, matrix or string. ")
        return
    elseif (typeof(pseudo_def["tsamps"]) != Array{Array{Int,1},1}) && (typeof(pseudo_def["tsamps"]) != Array{Array{Float64,1},1}) && (typeof(pseudo_def["tsamps"]) != Array{Array{Float32,1},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field tsamps! This should be an array of arrays of numbers. ")
        return
    elseif (pseudo_def["plot"] != [])
        try
            if (typeof(pseudo_def["plot"]) != Bool) && (typeof(pseudo_def["plot"]) != Array{Bool,1}) &&
            (typeof(pseudo_def["plot"]) != String) && (typeof(pseudo_def["plot"]) != Array{String,1})
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check the field plot! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
                return
            end
        catch
             if (typeof(pseudo_def["plot"]) != Array{Bool,1}) &&
            (typeof(pseudo_def["plot"]) != Array{String,1})
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check the field plot! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
                return
            end
        end
    elseif (typeof(pseudo_def["flag"]) != Array{String,1}) && (typeof(pseudo_def["flag"]) != String) && ((pseudo_def["flag"]) != [])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field flag! This should be a vector, matrix or string. ")
        return
    elseif (typeof(pseudo_def["Obs"]) != Array{String,1}) && (typeof(pseudo_def["Obs"]) != Array{Int,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Obs! This should be a vector of strings or integers. ")
        return
    elseif (typeof(pseudo_def["NoiseType"]) != Array{String,1}) && (typeof(pseudo_def["NoiseType"]) != String) && ((pseudo_def["NoiseType"]) != [])

    elseif (typeof(pseudo_def["Noise"]) != Array{Float64,1}) && (typeof(pseudo_def["Noise"]) != Array{Float32,1})  &&
        ((pseudo_def["NoiseType"]) != []) && (typeof(pseudo_def["Noise"]) != Array{Int,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Noise! This should be a vector of error percentages from 0 to 1. ")
        return
    end

    # Extract necessary elements to ease generalisation

    if typeof(pseudo_def["Nexp"]) == Array{Int,1}
        pseudo_def["Nexp"] = pseudo_def["Nexp"][1];
    end

    if pseudo_def["plot"]==[true] || pseudo_def["plot"]==["Yes"] || pseudo_def["plot"]==["yes"] || pseudo_def["plot"]=="Yes" || pseudo_def["plot"]=="yes"
        pseudo_def["plot"]=true
    elseif pseudo_def["plot"]==[false] || pseudo_def["plot"]==["No"] || pseudo_def["plot"]==["no"] || pseudo_def["plot"]=="No" || pseudo_def["plot"]=="no" || pseudo_def["plot"]==[]
        pseudo_def["plot"]=false
    end

    if (typeof(pseudo_def["theta"]) == Array{String,1})
        pseudo_def["theta"] = pseudo_def["theta"][1];
    end

    if typeof(pseudo_def["flag"]) == Array{String,1}
        pseudo_def["flag"] = pseudo_def["flag"][1];
    elseif typeof(pseudo_def["flag"]) == String
        pseudo_def["flag"] = pseudo_def["flag"];
    elseif ((pseudo_def["flag"]) == [])
        pseudo_def["flag"] = "";
    end

    if (typeof(pseudo_def["NoiseType"]) == Array{String,1})
        pseudo_def["NoiseType"] = pseudo_def["NoiseType"][1];
    end

    # Check that all the contents make sense
    if (pseudo_def["Nexp"] != length(pseudo_def["finalTime"])) || (pseudo_def["Nexp"] != length(pseudo_def["switchT"])) ||
        (pseudo_def["Nexp"] != length(pseudo_def["y0"])) || (pseudo_def["Nexp"] != length(pseudo_def["uInd"])) ||
        (pseudo_def["Nexp"] != length(pseudo_def["tsamps"]))
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, some field does not match the number of experiments selected. Check that the number of")
        println("entries for finalTime, switchT, y0, uInd or tsamps matches the number of experiments in Nexp.")
        return
    end


    if (typeof(pseudo_def["theta"]) == Array{Float64,1}) || (typeof(pseudo_def["theta"]) == Array{Float32,1})
        if length(pseudo_def["theta"]) != model_def["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Number of parameters introduced does not match the specified")
            return
        end
    elseif (typeof(pseudo_def["theta"]) == Array{Float64,2}) || (typeof(pseudo_def["theta"]) == Array{Float32,2})
        if size(pseudo_def["theta"])[1] != model_def["nPar"] && size(pseudo_def["theta"])[2] != model_def["nPar"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Number of parameters introduced does not match the specified")
            return
        end
    elseif (typeof(pseudo_def["theta"]) == String)
        if pseudo_def["theta"][end-3:end] != ".csv"
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, add the .csv termination to the theta file! And no spaces after!")
            return
        end
        if isfile(pseudo_def["theta"])
            pseudo_def["theta"] = Matrix(CSV.read(pseudo_def["theta"]));
        else
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but the file path you introduced does not exist!")
            return
        end
    end

    for i in 1:pseudo_def["Nexp"]
        if pseudo_def["tsamps"][i][end] > pseudo_def["finalTime"][i]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check finalTime. You have selected a sampling point past it.")
            return
        end
        if model_def["nInp"] != 0
            if size(pseudo_def["uInd"][i])[2] != (length(pseudo_def["switchT"][i])-1)
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check uInd and switchT. Number of steps does not match the number of values for the inputs.")
                return
            end
        end

        if length(pseudo_def["y0"][i])/model_def["nStat"] == 1
            if length(pseudo_def["y0"][i]) != model_def["nStat"]
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but Y0 does not match the number or states in experiment ", i))
                return
            end
        else
            if size(pseudo_def["y0"][i])[1] != model_def["nStat"] && size(pseudo_def["y0"][i])[2] != model_def["nStat"]
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but Y0 does not match the number or states in experiment ", i))
                return
            end
            if length(pseudo_def["y0"][i])/model_def["nStat"] != length(pseudo_def["theta"])/model_def["nPar"]
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but Y0 does not match the dimensions of theta "))
                return
            end
            if size(pseudo_def["y0"][i])[1] == size(pseudo_def["y0"][i])[2]
                println("-------------------------- WARNING --------------------------")
                println(string("Sorry, but the number of rows and columns of the y0 matrix in experiment ",i," is the same, so the checks on "))
                println("correct orientation will not work. Please make sure that the dimensions follow: ")
                println("y0[samples, states]")
            end
        end
    end

    if (model_def["Y0eqs"] != []) && pseudo_def["preInd"]==[]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("You specified computation of Y0 as steady-state but you have not introduced an inducer value for it!")
        return
    end

    if pseudo_def["NoiseType"] == []
        pseudo_def["NoiseType"] = "hetero"
    end

    pseudo_def["NoiseType"] = lowercase(pseudo_def["NoiseType"]); # To avoid issues with uper cases and that

    if (pseudo_def["NoiseType"] == "homoscedastic")
        pseudo_def["NoiseType"] = "homo";
    elseif (pseudo_def["NoiseType"] == "heteroscedastic")
        pseudo_def["NoiseType"] = "hetero";
    end

    if pseudo_def["Noise"] == []
        pseudo_def["Noise"] = ones(length(pseudo_def["Obs"])) .* 0.1;
    end

    if (pseudo_def["NoiseType"] == "hetero") && (maximum(pseudo_def["Noise"]) > 1)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, but the noise percentage has to be between 0 and 1!")
        return
    end

    if (minimum(pseudo_def["Noise"]) < 0)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, but the noise percentage or standard deviation is negative and that CANNOT be!")
        return
    end

    if length(pseudo_def["Obs"]) != length(pseudo_def["Noise"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, but you need to assign an error value for each observable selected!")
        return
    end

    # Check warning in case the matrix introduced is symmetric
    if length(size(pseudo_def["theta"])) == 2
        if size(pseudo_def["theta"])[1] == size(pseudo_def["theta"])[2]
            println("-------------------------- WARNING --------------------------")
            println("Sorry, but the number of rows and columns of the theta matrix is the same, so the checks on ")
            println("correct orientation will not work. Please make sure that the dimensions follow: ")
            println("theta[samples, parameters]")
        end
    end

    # Check that no step is no smaller than 2 unit of time
    for i in 1:pseudo_def["Nexp"]
        for j in 1:size(pseudo_def["uInd"][i])[2]
            if (pseudo_def["switchT"][i][j+1]-pseudo_def["switchT"][i][j])<=4
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but 2 of the steps in experiment ", i, " are too close. This package cannot "))
                println("handle this for now. We are working on it!")
                return
            end
        end
    end

    # Check that the inputs for the experiment are given as collumns
    try
        for i in 1:pseudo_def["Nexp"]
            if size(pseudo_def["uInd"][i])[2] != model_def["nInp"];
                pseudo_def["uInd"][i] = pseudo_def["uInd"][i]';
            end
        end
    catch
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but there is some issue with the contents of the field uInd."))
        return
    end

    if pseudo_def["Obs"] == Array{Int,1}
        if (maximum(pseudo_def["Obs"]) > model_def["nStat"]) || (length(pseudo_def["Obs"]) > model_def["nStat"])
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but there is some issue with the contents of the field Obs. It seems that "))
            println("     you have selected more observables than states or an index higher than the number of observables")
            return
        end
    end

    if typeof(pseudo_def["Obs"]) == Array{String,1}
        if convert(Bool,sum(occursin.("+", pseudo_def["Obs"]))) || convert(Bool,sum(occursin.("-", pseudo_def["Obs"]))) ||
           convert(Bool,sum(occursin.("/", pseudo_def["Obs"]))) ||
           convert(Bool,sum(occursin.("*", pseudo_def["Obs"]))) || convert(Bool,sum(occursin.("^", pseudo_def["Obs"])))

            if sum([sum(occursin.(model_def["stName"], pseudo_def["Obs"][k])) for k in 1:length(pseudo_def["Obs"])]) == 0
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but there is some issue with the contents of the field Obs."))
                println("It seems that the observable(s) selected do not match any state")
                return
            end
        else
            if sum([sum(occursin.(model_def["stName"], pseudo_def["Obs"][k])) for k in 1:length(pseudo_def["Obs"])]) == 0
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but there is some issue with the contents of the field Obs."))
                println("It seems that the observable(s) selected do not match any state")
                return
            else
                pseudo_def["Obs"] = sort([findall(x->x==pseudo_def["Obs"][i], model_def["stName"])[1] for i in 1:length(pseudo_def["Obs"])])
            end
        end

    elseif typeof(pseudo_def["Obs"]) == Array{Int,1}
        pseudo_def["Obs"] = sort(pseudo_def["Obs"]);
    end

    return(pseudo_def)
end

## CSV introduction description
# function fileStructInfo()
#
#     println("ALL THE FILES MUST BE CSV FILES!!!")
#     println("")
#     println("---------------------------------------------------------------------------------------------------------")
#     println("")
#     println("The observables file entry should have the following structure: ")
#     println("Column 1: Sampling Times for the simulations")
#     println("Column 2 + Number of states: Y0 value for each each state (in order) in each column. ")
#     println("           Value might need to be repeated across all the column, but only the first row will be considered")
#     println("")
#     println("---------------------------------------------------------------------------------------------------------")
#     println("")
#     println("The event inputs file entry should have the following structure: ")
#     println("Column 1: Each Switching time (final time not considered)")
#     println("Column 2: Final time (can be repeated for as many entrances as column 1 has)")
#     println("Column 3 + number of inducers: Repeated column for the value of the inducers in the ON. ")
#     println("              If no pre-inducer will be used in the simulations, any number could be used")
#     println("Column 4 + number of inducers: Value of each inducer at each switching time.")
#     println("")
#     println("This is an Example of an experiment with 3 steps and 2 inducers: ")
#     println(" |Switching|  |FinalTime|  |IPTGpre|  |aTcPre|  |  IPTG  |  |  aTc  |")
#     println("      0           1350         1          0          1          0    ")
#     println("     400          1350         1          0         0.3         50   ")
#     println("     700          1350         1          0         0.7         35   ")
#     println("     950          1350         1          0         0.1        100   ")
#     println("")
#     println("---------------------------------------------------------------------------------------------------------")
#     println("")
#     println("For the main directory, remember that in Julia you have to type \\ twice!")
#
# end

## CSV introduction structure
function defPseudoDatStructFiles()

    fileStructInfo()

    # The structure of things user should give
    pseudo_def = Dict()

    pseudo_def["ObservablesFile"] = []; # File containing information about the states Y0 and sampling times
    pseudo_def["EventInputsFile"] = []; # -> File containing all the information about the stimuli
    pseudo_def["theta"] = []; # Vector/Matrix with the parameter samples or directory and file location of CSV file with them
    pseudo_def["MainDir"] = []; # Directory path where all files are located (to avoid repetition). This can be empty.
    pseudo_def["plot"] = [];
    pseudo_def["flag"] = [];

    pseudo_def["Obs"] = []; # States that are observable. This is either a vecotr of strings (that has the same order asmodel_def["stName"] and model_def["eqns"]), a vector of integers insidating which states are observable or the string "All" if all states are observables. This can also be an experssion combining states (Only +,-,*,\ and ^ will be considered)
    pseudo_def["NoiseType"] = []; # String indicating if the desired noise is homoscedastic or heteroscedastic. The allowed strings are homo, homoscedastic or hetero, heteroscedastic (case independent). If empty, default will be heteroscedastic. Note that for high homoscedastic noise data could become negative (which to me does not make much sense), for now this will be set to a value close to 0.
    pseudo_def["Noise"] = []; # Percentage of heteroscedastic noise (introduced as value from 0 to 1). If empty 10% heteroscedastic noise will be assumed. In the homoscedastic case, if empty a standard deviation of 1 will be assumes. This has to be a vector of noise values for each observable.

    return(pseudo_def)

end

## Extract data from CSV
# Function to extract contents of CSV files if user sellects that
function extractPseudoDatCSV(model_def, pseudo_def)

    # Check taht all the dictionary entries are correct
    entries = ["ObservablesFile", "EventInputsFile", "theta", "MainDir", "plot", "flag", "Obs", "Noise"]
    if symdiff(entries,keys(pseudo_def))!=[]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is something wrong...")
        println(symdiff(entries,keys(pseudo_def)))
        return
    end

    # Checks on the different strings introduced
    if pseudo_def["MainDir"] == [] || pseudo_def["MainDir"] == ""
        pseudo_def["MainDir"] = "";
    elseif typeof(pseudo_def["MainDir"])==Array{String,1}
        pseudo_def["MainDir"] = pseudo_def["MainDir"][1]
    elseif typeof(pseudo_def["MainDir"]) == String
        nothing
    else
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, there is something wrong in the type of the field MainDir")
        return
    end

    if typeof(pseudo_def["ObservablesFile"]) != Array{String,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, there is something wrong in the type of the field ObservablesFile. This should be a Array{String,1}")
        return
    end

    if length(pseudo_def["ObservablesFile"]) != length(pseudo_def["EventInputsFile"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, there number of entries for ObservablesFile and EventInputsFile should be the same")
        return
    end

    for i in 1:length(pseudo_def["ObservablesFile"])
        if pseudo_def["ObservablesFile"][i][end-3:end] != ".csv"
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, add the .csv termination to the ObservablesFile file! And no spaces after!")
            return
        end
    end


    if typeof(pseudo_def["EventInputsFile"]) != Array{String,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, there is something wrong in the type of the field EventInputsFile. This should be a Array{String,1}")
        return
    end

    for i in 1:length(pseudo_def["EventInputsFile"])
        if pseudo_def["EventInputsFile"][i][end-3:end] != ".csv"
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, add the .csv termination to the EventInputsFile file! And no spaces after!")
            return
        end
    end



    simstu = defPseudoDatStruct();

    simstu["theta"] = pseudo_def["theta"];

    tmpInp = Array{Any}(undef, 1,length(pseudo_def["EventInputsFile"]))
    for i in 1:length(pseudo_def["EventInputsFile"])
        if pseudo_def["MainDir"] != [] && pseudo_def["MainDir"] != ""
            tmpInp[i] = Matrix(CSV.read(string(pseudo_def["MainDir"], "\\\\",pseudo_def["EventInputsFile"][i])));
        else
            tmpInp[i] = Matrix(CSV.read(pseudo_def["EventInputsFile"][i]));
        end
    end


    tmpObs = Array{Any}(undef, 1,length(pseudo_def["ObservablesFile"]))
    for i in 1:length(pseudo_def["ObservablesFile"])
        if pseudo_def["MainDir"] != [] && pseudo_def["MainDir"] != ""
            tmpObs[i] = Matrix(CSV.read(string(pseudo_def["MainDir"], "\\\\",pseudo_def["ObservablesFile"][i])));
        else
            tmpObs[i] = Matrix(CSV.read(pseudo_def["ObservablesFile"][i]));
        end
    end

    try
        simstu["Nexp"] = length(pseudo_def["EventInputsFile"]);
        simstu["finalTime"] = [tmpInp[i][1,2] for i in 1:length(pseudo_def["EventInputsFile"])];
        simstu["switchT"] = [vcat(tmpInp[i][:,1], tmpInp[i][1,2]) for i in 1:length(pseudo_def["EventInputsFile"])];
        simstu["y0"] = [tmpObs[i][1,2:(1+model_def["nStat"])] for i in 1:length(pseudo_def["ObservablesFile"])];
        simstu["preInd"] = [tmpInp[i][1,3:(2+model_def["nInp"])] for i in 1:length(pseudo_def["EventInputsFile"])];
        simstu["uInd"] = [Array(tmpInp[i][:,(3+model_def["nInp"]):(3+model_def["nInp"]+(model_def["nInp"]-1))]') for i in 1:length(pseudo_def["EventInputsFile"])];
        simstu["theta"] = pseudo_def["theta"];
        simstu["tsamps"] = [tmpObs[i][:,1] for i in 1:length(pseudo_def["ObservablesFile"])];
        simstu["plot"] = pseudo_def["plot"];
        simstu["flag"] = pseudo_def["flag"];
        simstu["Obs"] = pseudo_def["Obs"];
        simstu["NoiseType"] = pseudo_def["NoiseType"];
        simstu["Noise"] = pseudo_def["Noise"];
    catch
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, but there seems to be some issue with the internal structure of the files you have given...")
        return
    end

    simstu = checkStructPseudoDat(model_def, simstu);

    return(simstu)

end


function PDatCSVGen(pseudo_res,model_def,pseudo_def)
    ## Generate CSV Files
    if !isdir(string(pseudo_def["savepath"], "\\PseudoDataFiles"))
        mkdir(string(pseudo_def["savepath"], "\\PseudoDataFiles"))
    end

    for expp in 1:pseudo_def["Nexp"]

        # Simulations Matrix
        a,b,c = size(pseudo_res["Sims"][string("Exp_", expp)]);
        Simulations = zeros(a, (b*c)+1);
        Simulations[:,1] = pseudo_def["tsamps"][expp];

        sims_head = ["time"];

        coun = 2;
        for i in 1:c
            Simulations[:,coun:(coun+b-1)] = pseudo_res["Sims"][string("Exp_", expp)][:,:,i];
            coun = coun+b;
            for j in 1:b
                sims_head = vcat(sims_head, string(model_def["stName"][j], "_Theta_", i))
            end
        end

        df1 = DataFrame[Simulations][1]
        rename!(df1, [sims_head[i] for i in 1:(b*c)+1])
        CSV.write(string(pseudo_def["savepath"], "\\PseudoDataFiles\\",model_def["NameF"],"_EXP",expp, "_",
                pseudo_def["flag"],"_Simulations.csv"), df1);

        # Observables matrix
        a,b,c = size(pseudo_res["PData"][string("PDExp_", expp)]);
        Observables = zeros(a, (3*b)*c);
        obs_head = Array{String,1}(undef,(3*b)*c);

        r = 1:3:(3*b)*c
        for i in 1:length(r)
            Observables[:,r[i]] = pseudo_def["tsamps"][expp];
            obs_head[r[i]] = string("time",i);
        end


        co = 2;
        for i in 1:c
            for j in 1:b
                Observables[:,co] = pseudo_res["PData"][string("PDExp_", expp)][:,j,i];
                Observables[:,co+1] = pseudo_res["PError"][string("PDExp_", expp)][:,j,i]
                obs_head[co] = string(pseudo_def["Obs"][j],"_theta", i, "_Mean");
                obs_head[co+1] = string(pseudo_def["Obs"][j],"_theta", i, "_Std");
                co +=3;
            end
        end

        df2 = DataFrame[Observables][1]
        rename!(df2, [obs_head[i] for i in 1:(3*b)*c])
        CSV.write(string(pseudo_def["savepath"], "\\PseudoDataFiles\\",model_def["NameF"],"_EXP",expp, "_",
                pseudo_def["flag"],"_Observables.csv"), df2);

        # Event Inputs matrix
        if model_def["nInp"] != 0
            EvnInputs = zeros(convert(Int, length(pseudo_def["uInd"][expp])/model_def["nInp"]),2+(model_def["nInp"]*2));
            evins_head = Array{String,1}(undef,2+(model_def["nInp"]*2));
            evins_head[1] = "Switchingtimes";
            evins_head[2] = "FinalTime";

            EvnInputs[:,1] = pseudo_def["switchT"][expp][1:end-1];
            EvnInputs[:,2] = repeat([pseudo_def["switchT"][expp][end]], outer = [convert(Int,length(pseudo_def["uInd"][expp])/model_def["nInp"])]);

            for i in 1:model_def["nInp"]
                if pseudo_def["preInd"] != []
                    EvnInputs[:,(2+i)] = repeat([pseudo_def["preInd"][expp][i]], outer = [convert(Int,length(pseudo_def["uInd"][expp])/model_def["nInp"])]);
                else
                    EvnInputs[:,(2+i)] = repeat([0], outer = [convert(Int,length(pseudo_def["uInd"][expp])/model_def["nInp"])]);
                end
                evins_head[(2+i)] = string(model_def["inpName"][i],"_Pre");
                if model_def["nInp"] == 1
                    EvnInputs[:,(2+model_def["nInp"])+i] = pseudo_def["uInd"][expp];
                else
                    EvnInputs[:,(2+model_def["nInp"])+i] = pseudo_def["uInd"][expp][:,i];
                end
                evins_head[(2+model_def["nInp"])+i] = string(model_def["inpName"][i]);
            end


            df3 = DataFrame[EvnInputs][1]
            rename!(df3, [evins_head[i] for i in 1:2+(model_def["nInp"]*2)])
            CSV.write(string(pseudo_def["savepath"], "\\PseudoDataFiles\\",model_def["NameF"],"_EXP",expp, "_",
                    pseudo_def["flag"],"_Event_Inputs.csv"), df3);
        end

    end

end
