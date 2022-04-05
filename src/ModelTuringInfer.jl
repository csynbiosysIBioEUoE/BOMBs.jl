
##
function defTurInfStruct()
    turinf_def = Dict();
    turinf_def["Priors"] = []; # No option to leave this empty
    turinf_def["Data"] = [];
    # turinf_def["TuringSettings"] = []; # Two options:
    #                         # 1) A dictionary with the basic fields from a turing run using NUTS. The structure of the dictionary can be
    #                             # extracted calling the function defBasicTuringSettingsStruct().
    #                         # 2) An empty array. If this is the case, no inference will be done after calling the main
    #                             # function. Instead, you will be given a TuringModel file and data structure and you will
    #                             # have to use this to run a call of Turing by yourself. An example will be provided.

    turinf_def["flag"] = []; # String to attach a unique flag to the generated files so it is not overwritten. If empty, nothing will be added.

    turinf_def["plot"] = []; # true, flase, "Yes", "yes", "No", "no" indicating if plots with the results will be generated. Default is false.

    # turinf_def["runInf"] = []; # true, flase, "Yes", "yes", "No", "no" indicating if inference wants to be performed or only stan structure is given. Default is true (if all the necessary elements are present)

    turinf_def["MultiNormFit"] = [];

    turinf_def["Trunc"] = []; # true, flase, "Yes", "yes", "No", "no" or empty vector indicating that if samples are used to fit the prior, if the user wants truncation in these or not. Truncations will be the minimum and maximum sample values for each parameter minus and plus 10 % of the value respectively. Highly recomended in Turing otherwise sampling issues arise at least with NUTS. Default will be true. Alternatively you can introduce a 2*nTheta martix with the desired truncation bounds.

    return(turinf_def)
end

##
# function defBasicTuringSettingsStruct()
#     turing_def = Dict();
#
#     turing_def["nchains"] = []; # Number of chains to be used. If none, 4 chains will be used.
#     turing_def["nsamples"] = []; # Number of samples for each chain. If none, 2000 will be used.
#     turing_def["nadaptations"] = []; # Number of samples for the adaptation step (warmup). If none, default from Turing will be used. Remember that these will be substracted from the total nsamples selected.
#     turing_def["acceptRatio"] = []; # Target acceptance ration for NUTS. Value between 0 and 1. If none, 0.45 will be assumed.
#     turing_def["printsummary"] = []; # Print inference progres in console. This can be either true or false. Default will be true.
#     turing_def["init"] = []; # Vector for initial theta. In this case this will be the same one for each chain. If none, just leave empty.
#
#     return(turing_def)
# end


##
function fitPriorSampsTuring(priorsamps, model_def, turinf_def)
    priorfit = Dict();

    if length(size(priorsamps)) == 2
        if size(priorsamps)[1] == size(priorsamps)[2]
            println("-------------------------- WARNING --------------------------")
            println("Sorry, but the number of rows and columns of the theta matrix is the same, so the checks on ")
            println("correct orientation will not work. Please make sure that the dimensions follow: ")
            println("theta[samples, parameters]")
        elseif size(priorsamps)[1] == model_def["nPar"]
            priorsamps = priorsamps';

        end

    elseif length(size(priorsamps)) == 1 && model_def["nPar"]>1
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, but the CSV file given only contains a vector, not a matrix.")
        return
    elseif length(size(priorsamps)) == 1 && model_def["nPar"]==1
        priorsamps = reshape(priorsamps, length(priorsamps),model_def["nPar"]);
    end


    dists = [Beta, Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform];
    dists2 = [Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform]; # For the parameters that are not between 0 and 1 (to avoid error)

    fitts = Dict(); # Where all the distribution fits will be stored
    bestfit = Dict(); # Where the best distribution fit will be stored
    bestfitInd = zeros(1,model_def["nPar"]); # Index for the best distribution fit for each parameter

    names = model_def["parName"];

    for i in 1:length(names)
        try
            fitts[names[i]] = fit.(dists, Ref(priorsamps[:,i]));
            bestfitInd[i] = findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2];
            bestfit[names[i]] = fitts[names[i]][findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2]];
        catch
            fitts[names[i]] = fit.(dists2, Ref(priorsamps[:,i]));
            bestfitInd[i] = findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2]+1;
            bestfit[names[i]] = fitts[names[i]][findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2]];
        end
    end

    if turinf_def["Trunc"] == true

        boundsup = maximum(priorsamps, dims=1)+(maximum(priorsamps, dims=1)*0.1);
        boundsdw = minimum(priorsamps, dims=1)-(minimum(priorsamps, dims=1)*0.1);
        for k in 1:length(boundsdw)
            if boundsdw[k] < 0
                boundsdw[k] = 0;
            end
        end

        newPri = [Truncated(bestfit[names[i]], boundsdw[i], boundsup[i]) for i in 1:model_def["nPar"]]
    elseif turinf_def["Trunc"] == false
        newPri = [bestfit[names[i]] for i in 1:model_def["nPar"]]
    else
        boundss = reshape(turinf_def["Trunc"], 2, model_def["nPar"]);
        boundsup = maximum(boundss, dims=1);
        boundsdw = minimum(boundss, dims=1);

        newPri = [Truncated(bestfit[names[i]], boundsdw[i], boundsup[i]) for i in 1:model_def["nPar"]]
    end


    return(newPri)
end

##
function fitPriorSampsMultiNormTuring(priorsamps, model_def, turinf_def)

    priorfit = Dict();

    if length(size(priorsamps)) == 2
        if size(priorsamps)[1] == size(priorsamps)[2]
            println("-------------------------- WARNING --------------------------")
            println("Sorry, but the number of rows and columns of the theta matrix is the same, so the checks on ")
            println("correct orientation will not work. Please make sure that the dimensions follow: ")
            println("theta[parameters,samples]")
        elseif size(priorsamps)[1] == model_def["nPar"]
            priorsamps = priorsamps';

        end

    elseif length(size(priorsamps)) == 1 && model_def["nPar"]>1
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, but the CSV file given only contains a vector, not a matrix.")
        return
    elseif length(size(priorsamps)) == 1 && model_def["nPar"]==1
        priorsamps = reshape(priorsamps, length(priorsamps),model_def["nPar"]);
    end


#     if size(priorsamps)[1]==model_def["nPar"]
#         newPri=fit(MvNormal, priorsamps)
#     else
#         newPri=fit(MvNormal, priorsamps')
#     end


    if turinf_def["Trunc"] == true

        boundsup = maximum(priorsamps, dims=1)+(maximum(priorsamps, dims=1)*0.1);
        boundsdw = minimum(priorsamps, dims=1)-(minimum(priorsamps, dims=1)*0.1);

        for k in 1:length(boundsdw)
            if boundsdw[k] < 0
                boundsdw[k] = 0;
            end
        end

        newPri = [Truncated(fit(Normal, priorsamps[:,i]), boundsdw[i], boundsup[i]) for i in 1:model_def["nPar"]]
    elseif turinf_def["Trunc"] == false
        newPri = [fit(Normal, priorsamps[:,i]) for i in 1:model_def["nPar"]]
    else
        boundss = reshape(turinf_def["Trunc"], 2, model_def["nPar"]);
        boundsup = maximum(boundss, dims=1);
        boundsdw = minimum(boundss, dims=1);

        newPri = [Truncated(fit(Normal, priorsamps[:,i]), boundsdw[i], boundsup[i]) for i in 1:model_def["nPar"]]
    end

    return(newPri)

end


##
function checkStructTurInf(model_def, turinf_def)

    # Check taht all the dictionary entries are correct
    # entries = ["Priors", "Data", "TuringSettings", "flag", "plot", "runInf", "Trunc", "MultiNormFit"]
    entries = ["Priors", "Data", "flag", "plot", "Trunc", "MultiNormFit"]
    if symdiff(entries,keys(turinf_def))!=[] && symdiff(entries,keys(turinf_def)) != ["ModelPath"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is something wrong...")
        println(symdiff(entries,keys(turinf_def)))
        return
    end

    # Check one to be sure there is no empty entry
    nee = ["Data"]; # No empty entries
    for i in 1:length(nee)
        if turinf_def[nee[i]] == []
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Plese, check the field ", nee[i], ", this cannot be empty"))
            return
        end
    end

    if turinf_def["flag"]!=[] && typeof(turinf_def["flag"])!=String && typeof(turinf_def["flag"])!=Array{String,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field flag!")
        return
    end

    if typeof(turinf_def["flag"]) ==Array{String,1}
        turinf_def["flag"] = turinf_def["flag"][1];
    end

    if (typeof(turinf_def["plot"]) != Bool) && (typeof(turinf_def["plot"]) != Array{Bool,1}) &&
    (typeof(turinf_def["plot"]) != String) && (typeof(turinf_def["plot"]) != Array{String,1}) && turinf_def["plot"] != []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field plot! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
        return
    end

    if (typeof(turinf_def["Trunc"]) != Bool) && (typeof(turinf_def["Trunc"]) != Array{Bool,1}) &&
    (typeof(turinf_def["Trunc"]) != String) && (typeof(turinf_def["Trunc"]) != Array{String,1}) && turinf_def["Trunc"] != [] &&
        typeof(turinf_def["Trunc"]) != Array{Int,2} && typeof(turinf_def["Trunc"]) != Array{Float32,2} &&
        typeof(turinf_def["Trunc"]) != Array{Float64,2} && typeof(turinf_def["Trunc"]) != Array{Any,2}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Trunc! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\", an empty array or a 2D matrix with teh desired truncation bounds")
        return
    end

    # if (typeof(turinf_def["runInf"]) != Bool) && (typeof(turinf_def["runInf"]) != Array{Bool,1}) &&
    # (typeof(turinf_def["runInf"]) != String) && (typeof(turinf_def["runInf"]) != Array{String,1}) && turinf_def["runInf"] != []
    #     println("-------------------------- Process STOPPED!!! --------------------------")
    #     println("Please, check the field runInf! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
    #     return
    # end

    if (typeof(turinf_def["MultiNormFit"]) != Bool) && (typeof(turinf_def["MultiNormFit"]) != Array{Bool,1}) &&
    (typeof(turinf_def["MultiNormFit"]) != String) && (typeof(turinf_def["MultiNormFit"]) != Array{String,1}) && turinf_def["MultiNormFit"] != []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field MultiNormFit! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
        return
    end

    if turinf_def["plot"]==[true] || turinf_def["plot"]==["Yes"] || turinf_def["plot"]==["yes"] || turinf_def["plot"]=="Yes" || turinf_def["plot"]=="yes"
        turinf_def["plot"]=true
    elseif turinf_def["plot"]==[false] || turinf_def["plot"]==["No"] || turinf_def["plot"]==["no"] || turinf_def["plot"]=="No" || turinf_def["plot"]=="no" || turinf_def["plot"]==[]
        turinf_def["plot"]=false
    end

    # if turinf_def["runInf"]==[true] || turinf_def["runInf"]==["Yes"] || turinf_def["runInf"]==["yes"] || turinf_def["runInf"]=="Yes" || turinf_def["runInf"]=="yes" || turinf_def["runInf"]==[]
    #     turinf_def["runInf"]=true
    # elseif turinf_def["runInf"]==[false] || turinf_def["runInf"]==["No"] || turinf_def["runInf"]==["no"] || turinf_def["runInf"]=="No" || turinf_def["runInf"]=="no"
    #     turinf_def["runInf"]=false
    # end

    if turinf_def["MultiNormFit"]==[true] || turinf_def["MultiNormFit"]==["Yes"] || turinf_def["MultiNormFit"]==["yes"] || turinf_def["MultiNormFit"]=="Yes" || turinf_def["MultiNormFit"]=="yes"
        turinf_def["MultiNormFit"]=true
    elseif turinf_def["MultiNormFit"]==[false] || turinf_def["MultiNormFit"]==["No"] || turinf_def["MultiNormFit"]==["no"] || turinf_def["MultiNormFit"]=="No" || turinf_def["MultiNormFit"]=="no" || turinf_def["MultiNormFit"]==[]
        turinf_def["MultiNormFit"]=false
    end

    if turinf_def["Trunc"]==[true] || turinf_def["Trunc"]==true || turinf_def["Trunc"]==["Yes"] || turinf_def["Trunc"]==["yes"] || turinf_def["Trunc"]=="Yes" || turinf_def["Trunc"]=="yes" || turinf_def["Trunc"]==[]
        turinf_def["Trunc"]=true
    elseif turinf_def["Trunc"]==[false] || turinf_def["Trunc"]==false || turinf_def["Trunc"]==["No"] || turinf_def["Trunc"]==["no"] || turinf_def["Trunc"]=="No" || turinf_def["Trunc"]=="no"
        turinf_def["Trunc"]=false
    else
        if size(turinf_def["Trunc"])!=(2,model_def["nPar"]) && size(turinf_def["Trunc"])!=(model_def["nPar"],2)

            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check the field Trunc! You have introduced a matrix with bounds but the size of this seems wrong.")
            return
        end
    end


    # Extract contents of Priors field
    if (typeof(turinf_def["Priors"]) <: Array) || (typeof(turinf_def["Priors"]) <: LinearAlgebra.Adjoint) ||
            typeof(turinf_def["Priors"]) == String# Object is an array. This could be an array of values or array of distributions

        if typeof(turinf_def["Priors"]) == Array{String,1}
            turinf_def["Priors"] = turinf_def["Priors"][1];
        end

        if typeof(turinf_def["Priors"]) == String
            if turinf_def["Priors"][end-3:end] == ".csv"
                if isfile(turinf_def["Priors"])
                    priorsamps = Matrix(CSV.read(turinf_def["Priors"]));
                    if turinf_def["MultiNormFit"] == false
                        turinf_def["Priors"] = fitPriorSampsTuring(priorsamps, model_def, turinf_def);
                    elseif turinf_def["MultiNormFit"] == true
                        turinf_def["Priors"] = fitPriorSampsMultiNormTuring(priorsamps, model_def, turinf_def);
                    end
                else
                    println("-------------------------- Process STOPPED!!! --------------------------")
                    println("Sorry, but the CSV file introduced for the prior fit does not exist.")
                    return
                end
            else
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check the entries of the dictionary in Priors. This has to be the path to a .csv")
                println("file. Also check that there are no spaces after the termination.")
                return
            end
        end


        if length(turinf_def["Priors"]) == 1
            if typeof(turinf_def["Priors"])<:Distribution
                nothing
            else
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("You have introduced an array with only one element but this is not a Distribution. ")
                println("Introduce a distribution, samples or the path to a CSV file with samples. ")
                return
            end
        else
            if (size(turinf_def["Priors"]) == (2,model_def["nPar"])) || (size(turinf_def["Priors"]) == (model_def["nPar"],2))
                if (size(turinf_def["Priors"]) == (model_def["nPar"],2))
                    turinf_def["Priors"]=turinf_def["Priors"]';
                end

                thmax = maximum(turinf_def["Priors"], dims = 1)
                thmin = minimum(turinf_def["Priors"], dims = 1)
                mess = (thmax+thmin)/2;
                errs = (thmax-thmin)/4;

                priors = Array{Any,1}(undef, model_def["nPar"]);

                for thet in 1:model_def["nPar"]
                    priors[thet] = Truncated(Normal(mess[thet],errs[thet]), thmin[thet],thmax[thet])
                end

                turinf_def["Priors"] = priors;


            elseif typeof(turinf_def["Priors"]) == Array{Float64,2} || typeof(turinf_def["Priors"]) == Array{Float32,2}
                if turinf_def["MultiNormFit"] == false
                    turinf_def["Priors"] = fitPriorSampsTuring(turinf_def["Priors"], model_def, turinf_def);
                elseif turinf_def["MultiNormFit"] == true
                    turinf_def["Priors"] = fitPriorSampsMultiNormTuring(turinf_def["Priors"], model_def, turinf_def);
                end
            end
        end
    elseif turinf_def["Priors"] == [];
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but in an inference using Turing you allways need to specify the priors in the main structure.")
            return

    end

    #     Check Data

    if haskey(turinf_def["Data"], "switchT")
        turinf_def["Data"] = checkStructBayInfData(model_def, turinf_def["Data"]);
    else
        turinf_def["Data"] = checkStructBayInfDataFiles(model_def, turinf_def["Data"])
    end
    #     Check Stan Setting
    # if turinf_def["TuringSettings"] != []
    #     turinf_def["TuringSettings"] = checkStructBayInfTuringSettings(model_def, turinf_def["TuringSettings"]);
    # end

    for i in 1:turinf_def["Data"]["Nexp"]
        if size(turinf_def["Data"]["uInd"][i])[1] != model_def["nInp"];
            turinf_def["Data"]["uInd"][i] = turinf_def["Data"]["uInd"][i]';
        end
    end

    return(turinf_def)
end



##
# function checkStructBayInfTuringSettings(model_def, turing_def)
#
#     # Check taht all the dictionary entries are correct
#     entries = ["acceptRatio", "nsamples", "init", "printsummary", "nchains", "nadaptations"]
#     if symdiff(entries,keys(turing_def))!=[] && symdiff(entries,keys(turing_def)) != ["savepath", "savename"]
#         println("-------------------------- Process STOPPED!!! --------------------------")
#         println("Please, check the entries of the dictionary, there is something wrong...")
#         println(symdiff(entries,keys(turing_def)))
#         return
#     end
#
#     # Checks to see if the content of the fields is the expected (also consier if single value entries have been given as a number or a vector)
#     if typeof(turing_def["nchains"]) != Int && typeof(turing_def["nchains"]) != Array{Int,1}
#         println("-------------------------- Process STOPPED!!! --------------------------")
#         println("Please, check the field nchains! This should be an string. ")
#         return
#     elseif typeof(turing_def["nsamples"]) != Int && typeof(turing_def["nsamples"]) != Array{Int,1}
#         println("-------------------------- Process STOPPED!!! --------------------------")
#         println("Please, check the field nsamples! This should be an string. ")
#         return
#     elseif typeof(turing_def["nadaptations"]) != Int && typeof(turing_def["nadaptations"]) != Array{Int,1}
#         println("-------------------------- Process STOPPED!!! --------------------------")
#         println("Please, check the field nadaptations! This should be an string. ")
#         return
#     elseif turing_def["printsummary"] != true && turing_def["printsummary"] != false && turing_def["printsummary"] != [true] &&
#         turing_def["printsummary"] != [false] && turing_def["printsummary"] != []
#         println("-------------------------- Process STOPPED!!! --------------------------")
#         println("Please, check the field printsummary! This should be an string. ")
#         return
#     elseif turing_def["init"] != [] && !(typeof(turing_def["init"]) <: Array)
#         println("-------------------------- Process STOPPED!!! --------------------------")
#         println("Please, check the field init! This should be an array of dictionaries or an empty array. ")
#         return
#     elseif !isempty(turing_def["init"])
#         if length(turing_def["init"]) != model_def["nPar"]
#             println("-------------------------- Process STOPPED!!! --------------------------")
#             println("Please, check the field init! This shoudl be a vector with a single sample for each parameter. ")
#             return
#         end
#         for i in 1:model_def["nPar"]
#             if typeof(turing_def["init"][i]) != Int && typeof(turing_def["init"][i]) != Float32 &&
#                 typeof(turing_def["init"][i]) != Float64
#                 println("-------------------------- Process STOPPED!!! --------------------------")
#                 println("Please, check the field init! This shoudl be a vector with a single sample for each parameter and this sample has to be a number!")
#                 return
#             end
#         end
#     end
#
#     # Extract necessary elements to ease generalisation
#
#     if typeof(turing_def["acceptRatio"]) == Array{Int,1} || typeof(turing_def["acceptRatio"]) == Array{Float32,1} ||
#         typeof(turing_def["acceptRatio"]) == Array{Float64,1}
#         turing_def["acceptRatio"] = turing_def["acceptRatio"][1];
#     end
#
#     if typeof(turing_def["nchains"]) == Array{Int,1}
#         turing_def["nchains"] = turing_def["nchains"][1];
#     end
#
#     if typeof(turing_def["nsamples"]) == Array{Int,1}
#         turing_def["nsamples"] = turing_def["nsamples"][1];
#     end
#
#     if typeof(turing_def["nadaptations"]) == Array{Int,1}
#         turing_def["nadaptations"] = turing_def["nadaptations"][1];
#     end
#
#     if typeof(turing_def["printsummary"]) == Array{Bool,1}
#         turing_def["printsummary"] = turing_def["printsummary"][1];
#     end
#
#     if turing_def["printsummary"] == []
#         turing_def["printsummary"] = true;
#     end
#
#
#     return(turing_def)
# end


##
function genTuringModel(model_def, turinf_def)

    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\ModelsFunctions"))
        mkdir(string(cudi,"\\ModelsFunctions"))
    end


    model_def = checkStruct(model_def);

    # Packages needed for the simulations
    Head = string("
    using DifferentialEquations
    using OrdinaryDiffEq
    using DiffEqBase
    using Sundials
    using ODEInterfaceDiffEq
    using Turing
    using MCMCChains
        ");

    # Generate ODEs function
    # String containing the equtions, ODEs and others
    tet = Array{Any,1}(undef, length(model_def["eqns"])); # Vector of strings that will contain the eauqtions for the function
    # odec = 1; # Counter for the ODEs to take into account equations that are not ODEs
    for i in 1:length(model_def["eqns"]) # Loop over number of equations
        if replace(model_def["eqns"][i], " "=>"")[1]=='d' # If the equation starts with a d (differential)
            for j in 1:model_def["nStat"]
                eqsg = findfirst("=", replace(model_def["eqns"][i], " "=>""))[1];
                if replace(model_def["eqns"][i], " "=>"")[1:eqsg-1] == string("d", model_def["stName"][j])
                    tet[i] = string("        du[", j, "] = ", model_def["eqns"][i],";\n");
                    break;
                else
                    tet[i] = string("        ", model_def["eqns"][i],";\n");
                end
            end
        else # Other equations
            tet[i] = string("            ", model_def["eqns"][i],";\n");
        end
    end;

    fun1 = string("
        # du -> Derivative
        # u -> State at time t
        # p-> paramter vector
        # t-> time as tuple (init, end)

        function ",join(model_def["NameF"]),"ODE!(du,u,p,t)
            ",join(model_def["stName"], ", ")," = u;

            ",join(model_def["parName"], ", "),", ",join(model_def["inpName"], ", ")," = p;

",join(tet),"

        end
        ");

    # Generate Steady State function
    y0eq = ""; # Vector that will contain the steady state equations (SSE)
    if model_def["Y0eqs"]==[] # If there is no equations, make the function return the same Y0 that the user gives
        y0eq = (string(join(model_def["stName"], ", ")," = I; \n        ",


                "       return(I) \n "));
    else
        tet = Array{Any,1}(undef, length(model_def["Y0eqs"])); # Vector of strings that will contain the eauqtions for the function
        odec = 1 # Counter for the SSEs to take into account equations that are not SSEs
        df = Array{Any,1}(undef, length(model_def["Y0eqs"]));
        for i in 1:length(model_def["Y0eqs"]) # Check for the equations that start with the name of the states (to tell apart from other equations)
            if !isnothing(findfirst.("=", model_def["Y0eqs"][i]))
                df[i] = model_def["Y0eqs"][i][1:findfirst.("=", model_def["Y0eqs"][i])[1]]
            else
                df[i] = "";
            end
        end;
        for i in 1:length(model_def["Y0eqs"]) # Loop over equations
            if sum(occursin.(model_def["stName"], df[i]))!=0 # If the equation looked begins with the name of a state
                tet[i] = string("        alp[", odec, "] = ", model_def["Y0eqs"][i],";\n");
                odec +=1;
            else
                tet[i] = string("        ", model_def["Y0eqs"][i],";\n");
            end
        end

        y0eq = (string(join([string("exp",model_def["stName"][k]) for k in 1:model_def["nStat"]], ", ")," = I; \n        ",

        join(model_def["parName"], ", "),", ",join(model_def["inpName"], ", ")," = p;\n \n \n ",

                "        alp = zeros(",length(model_def["stName"]),"); \n \n",

                join(tet),
                "        return(alp) \n "));

    end

    fun2 = string("
        # p -> parameter vector (plus inducer vector at the end)
        # I -> initial Y0 vector

        function ",join(model_def["NameF"]),"SteadyState(p,I)

            ",y0eq,"

        end
        ");

    # Generate Step Wise simulation function
    ONSim = "";
    if model_def["Y0Sim"]==false
        ONSim = "y0 = y_al;";
    elseif model_def["Y0Sim"]==true
        if model_def["solver"] == "CVODE_BDF"
            ONSim = string("
                prob = ODEProblem(",join(model_def["NameF"]),"ODE!,y_al,(0.0,Float64(24*60-1)),pSS);
                ssv = Sundials.solve(prob, CVODE_BDF(),reltol=",model_def["tols"][1],",abstol=",model_def["tols"][2],");
                y0 = ssv[:,end]; ");
        else
            ONSim = string("
                prob = ODEProblem(",join(model_def["NameF"]),"ODE!,y_al,(0.0,Float64(24*60-1)),pSS);
                ssv = DifferentialEquations.solve(prob, ",join(model_def["solver"]),"(),reltol=",model_def["tols"][1],",abstol=",model_def["tols"][2],");
                y0 = ssv[:,end]; ");
        end
    elseif typeof(model_def["Y0Sim"]) == Float64 || typeof(model_def["Y0Sim"]) == Float32 || typeof(model_def["Y0Sim"]) == Int
        if model_def["solver"] == "CVODE_BDF"
            ONSim = string("
                prob = ODEProblem(",join(model_def["NameF"]),"ODE!,y_al,(0.0,Float64(",model_def["Y0Sim"],"*24*60-1)),pSS);
                ssv = DifferentialEquations.solve(prob, CVODE_BDF(),reltol=",model_def["tols"][1],",abstol=",model_def["tols"][2],");
                y0 = ssv[:,end]; ");
        else
            ONSim = string("
                prob = ODEProblem(",join(model_def["NameF"]),"ODE!,y_al,(0.0,Float64(",model_def["Y0Sim"],"*24*60-1)),pSS);
                ssv = DifferentialEquations.solve(prob, ",join(model_def["solver"]),"(),reltol=",model_def["tols"][1],",abstol=",model_def["tols"][2],");
                y0 = ssv[:,end]; ");
        end
    end

    if model_def["solver"] == "CVOED_BDF"
        soso = "Sundials.";
    else
        soso = "DifferentialEquations.";
    end



    if typeof(turinf_def["Priors"])<:Distribution{Multivariate}

        tp1 = string("            p1 ~ ",string(typeof(turinf_def["Priors"]))[1:findfirst("{",
                    string(typeof(turinf_def["Priors"])))[1]-1],"(",string(turinf_def["Priors"].μ),",",
                    string(turinf_def["Priors"].Σ.mat),"); \n");
        tp2 = string("            p = p1;");
    else

        tp1 = [string("            p",k," ~ ", turinf_def["Priors"][k], "\n") for k in 1:length(turinf_def["Priors"])]
        for j in 1:length(tp1)
            if typeof(tp1[j])<:Distribution{Multivariate}
                tp1[j] = string("            p1 ~ ",string(typeof(tp1[j]))[1:findfirst("{",
                    string(typeof(tp1[j])))[1]-1],"(",string(tp1[j].μ),",",
                    string(tp1[j].Σ.mat),"); \n");
            else
                tp1[j]=replace(tp1[j], "{Float64}"=>"")
                if !isempty(findall("range=", tp1[j]))
                    tp1[j]=replace(tp1[j], "range=("=>"")
                    tp1[j] = replace(tp1[j], ")\n"=>"\n");
                end

                pst = findall("=", tp1[j])

                for i in 1:length(pst)
                    ps = findfirst("=", tp1[j]);
                    try
                        tp1[j]=replace(tp1[j], tp1[j][ps[1]-1:ps[1]]=>"")
                    catch
                        tp1[j]=replace(tp1[j], tp1[j][ps[1]-2:ps[1]]=>"")
                    end
                end
            end
        end
        tp2 = string("            p = [",join([string("p",k,",") for k in 1:length(turinf_def["Priors"])]),"];");
    end

    simsB = Array{Any,1}(undef, length(turinf_def["Data"]["Obs"]));
    if typeof(turinf_def["Data"]["Obs"]) == Array{Int,1}

        for i in 1:length(turinf_def["Data"]["Obs"])
            simsB[i] = string("                        finalBay2[",i,",d] = part1[",turinf_def["Data"]["Obs"][i],", convert.(Int, d-sp[q])]; \n");
        end

    else

        for i in 1:length(turinf_def["Data"]["Obs"])
            if (occursin(" .+ ", turinf_def["Data"]["Obs"][i]))
                nothing
            elseif (occursin(".+", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], ".+"=>" .+ ")
            elseif (occursin("+", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], "+"=>" .+ ")
            end

            if (occursin(" .- ", turinf_def["Data"]["Obs"][i]))
                nothing
            elseif (occursin(".-", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], ".-"=>" .- ")
            elseif (occursin("-", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], "-"=>" .- ")
            end

            if (occursin(" ./ ", turinf_def["Data"]["Obs"][i]))
                nothing
            elseif (occursin("./", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], "./"=>" ./ ")
            elseif (occursin("/", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], "/"=>" ./ ")
            end

            if (occursin(" .* ", turinf_def["Data"]["Obs"][i]))
                nothing
            elseif (occursin(".*", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], ".*"=>" .* ")
            elseif (occursin("*", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], "*"=>" .* ")
            end

            if (occursin(" .^ ", turinf_def["Data"]["Obs"][i]))
                nothing
            elseif (occursin(".^", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], ".^"=>" .^ ")
            elseif (occursin("^", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], "^"=>" .^ ")
            end

            tmp2 = string("                        finalBay2[",i,",d] = ",turinf_def["Data"]["Obs"][i],"; \n");
            for j in 1:model_def["nStat"]
                if occursin(model_def["stName"][j], tmp2)
                    simsB[i] = replace(tmp2, model_def["stName"][j]=> string("part1[",j,", convert.(Int, d-sp[q])]"))
                end
            end

        end

    end

    if typeof(turinf_def["Data"]["Obs"]) == Array{Int,1}
        if length(size(turinf_def["Data"]["DataError"][1][1])) == 1
            lik = string("

                r = ",turinf_def["Data"]["Obs"],";
                for ob in 1:",length(turinf_def["Data"]["Obs"]),"
                    if typeof(p) == Array{Float64, 1}
                        obs = final[convert.(Int,round.(data[","\"","tsamps","\"","][exp])).+1, [r[ob]]];
                    else
                        obs = finalBay2[ob, convert.(Int,round.(data[","\"","tsamps","\"","][exp])).+1];
                    end

                    for dat in 1:length(obs)
                        data[","\"","DataMean","\"","][exp][:,ob][dat] ~ Normal(obs[dat], data[","\"","DataError","\"","][1][ob][dat])
                    end
                end
                ");
        else
#             lik = string("


#                 for ob in 1:",length(turinf_def["Data"]["Obs"]),"
#                     if typeof(p) == Array{Float64, 1}
#                         obs = final[convert.(Int,round.(data[","\"","tsamps","\"","][exp])).+1, ",turinf_def["Data"]["Obs"],"];
#                     else
#                         obs = finalBay2[ob, convert.(Int,round.(data[","\"","tsamps","\"","][exp])).+1];
#                     end

#                     data[","\"","DataMean","\"","][ob] ~ MvNormal(reshape(obs, length(obs)), data[","\"","DataError","\"","][ob][1].+Diagonal(ones(length(obs)).*1e-10))
#                 end
#                 ");

            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but multivariate likelohoods are still not implemented. Please, provide a vector of errors instead of a covariance matrix")
            return

        end
    else
        Obsers1 = Array{Any,1}(undef, length(turinf_def["Data"]["Obs"]))
        Obsers2 = Array{Any,1}(undef, length(turinf_def["Data"]["Obs"]))
        for i in 1:length(turinf_def["Data"]["Obs"])
            if (occursin(" .+ ", turinf_def["Data"]["Obs"][i]))
                nothing
            elseif (occursin(".+", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], ".+"=>" .+ ")
            elseif (occursin("+", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], "+"=>" .+ ")
            end

            if (occursin(" .- ", turinf_def["Data"]["Obs"][i]))
                nothing
            elseif (occursin(".-", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], ".-"=>" .- ")
            elseif (occursin("-", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], "-"=>" .- ")
            end

            if (occursin(" ./ ", turinf_def["Data"]["Obs"][i]))
                nothing
            elseif (occursin("./", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], "./"=>" ./ ")
            elseif (occursin("/", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], "/"=>" ./ ")
            end

            if (occursin(" .* ", turinf_def["Data"]["Obs"][i]))
                nothing
            elseif (occursin(".*", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], ".*"=>" .* ")
            elseif (occursin("*", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], "*"=>" .* ")
            end

            if (occursin(" .^ ", turinf_def["Data"]["Obs"][i]))
                nothing
            elseif (occursin(".^", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], ".^"=>" .^ ")
            elseif (occursin("^", turinf_def["Data"]["Obs"][i]))
                turinf_def["Data"]["Obs"][i] = replace(turinf_def["Data"]["Obs"][i], "^"=>" .^ ")
            end

            tmp1 = string("obs[",i,"] = ",turinf_def["Data"]["Obs"][i],"; \n");
            for j in 1:model_def["nStat"]
                if occursin(model_def["stName"][j], tmp1)
                    Obsers1[i] = replace(tmp1, model_def["stName"][j]=> string("final[convert.(Int,round(data[","\"","tsamps","\"","][exp])).+1, ",j,"]"))
                    Obsers2[i] = string("obs[",i,"] = finalBay2[",i,", convert.(Int,round(data[","\"","tsamps","\"","][exp])).+1];")
                end
            end

        end
        if length(size(turinf_def["Data"]["DataError"][1][1])) == 1
            lik = string("

                obs = Array{Any,1}(undef, ",length(turinf_def["Data"]["Obs"]),")
                if typeof(p) == Array{Float64, 1}
                    ",join(Obsers1),"
                else
                    ",join(Obsers2),"
                end

                for ob in 1:",length(turinf_def["Data"]["Obs"]),"
                    for dat in 1:length(obs)
                        data[","\"","DataMean","\"","][exp][:,ob][dat] ~ Normal(obs[dat], data[","\"","DataError","\"","][1][ob][dat])
                    end
                end
                ")
        else
#             lik = string("

#                 obs = Array{Any,1}(undef, ",length(turinf_def["Data"]["Obs"]),")
#                 if typeof(p) == Array{Float64, 1}
#                     ",join(Obsers1),"
#                 else
#                     ",join(Obsers2),"
#                 end

#                 for ob in 1:",length(turinf_def["Data"]["Obs"]),"
#                     data[","\"","DataMean","\"","][ob] ~ MvNormal(obs[ob], data[","\"","DataError","\"","][ob][1].+Diagonal(ones(length(obs)).*1e-10))
#                 end
#                 ")
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but multivariate likelohoods are still not implemented. Please, provide a vector of errors instead of a covariance matrix")
            return
        end

    end

    fun3 = string("
        # ts -> time vector going from t=0 to t=end every 1
        # p -> parameter vector
        # sp -> vector with switching times (t=0 and t=end needs to be included). Length of vector will be length(inducer steps)+1
        # inputs -> Vecor with the inducers (if more than 1 inducer in the system, the shape will be a1,b1,c1,..., a2,b2,c2,... aend,bend,cend,... where each leter represents a different inducer)
        # ivss -> Initial Y0 vector
        # pre -> Vector with the inducers in the ON. It can be empty

        @model function ",join(model_def["NameF"]),"_TuringModel(data,prob1)

",join(tp1),"
",join(tp2),"

            Neq = ",model_def["nStat"],";

            for exp in 1:",turinf_def["Data"]["Nexp"],"

                ts=collect(0.:1:data[","\"","finalTime","\"","][exp])
                sp = data[","\"","switchT","\"","][exp];
                inputs = data[","\"","uInd","\"","][exp];
                ivss = data[","\"","y0","\"","][exp];
                pre=data[","\"","preInd","\"","][exp];

                maxtime = length(ts);
                Nsp = length(sp);
                Nevents = length(sp)-1;


                if typeof(p) == Array{Float64, 1}
                    if pre != []
                        pSS = vcat(p, pre);
                    else
                        pSS = p;
                    end
                else
                    if pre != []
                        pSS = vcat([p[i].value for i in 1:",model_def["nPar"],"], pre);
                    else
                        pSS = [p[i].value for i in 1:",model_def["nPar"],"];
                    end
                end


                final = zeros(maxtime, Neq);
                finalBay2 = Array{Any,2}(undef, ",length(turinf_def["Data"]["Obs"]),", maxtime)

                y_al = ",join(model_def["NameF"]),"SteadyState(pSS,ivss) # Calculation of initial guesses for steady state
                y0 = y_al;
                initialV = y0;
                i = 1;

                for q in collect(1:Nevents)

                    lts = length(ts[convert.(Int, (sp[q]+1):sp[q+1]+1)]);  # General way to define the number of elements in each event series

                    Tevent = ts[convert.(Int,(sp[q]+1):sp[q+1]+1)];  # General way to extract the times of each event
                    I = inputs[i:(i+(",model_def["nInp"]-1,"))];
                    pSte = vcat(p, I);

                    if q == 1
                        if typeof(p) == Array{Float64, 1}
                            prob = ODEProblem(",join(model_def["NameF"]),"ODE!,initialV,(ts[convert.(Int, (sp[q]+1))],ts[convert.(Int, sp[q+1]+1)]),pSte);
                            part1 = ",soso,"solve(prob, ",soso,model_def["solver"],"(),reltol=",model_def["tols"][1],",abstol=",model_def["tols"][2],",saveat=1);
                        else
                            prob = remake(prob1, p=pSte)
                            prob = remake(prob, u0=initialV)
                            prob = remake(prob, tspan = convert.(Int, (ts[convert.(Int,(sp[q]+1))],ts[convert.(Int,sp[q+1]+1)])))
                            part1 = ",soso,"solve(prob,",soso,model_def["solver"],"(),reltol=",model_def["tols"][1],",abstol=",model_def["tols"][2],",saveat=1)
                        end
                    else
                        if typeof(p) == Array{Float64, 1}
                            prob = ODEProblem(",join(model_def["NameF"]),"ODE!,initialV,(ts[convert.(Int, (sp[q]+1))],ts[convert.(Int,sp[q+1]+1)]),pSte);
                            part1 = ",soso,"solve(prob, ",soso,model_def["solver"],"(),reltol=",model_def["tols"][1],",abstol=",model_def["tols"][2],",saveat=1);
                        else
                            prob = remake(prob1, p=pSte)
                            prob = remake(prob, u0=initialV)
                            prob = remake(prob, tspan = (ts[convert.(Int,(sp[q]+1))],ts[convert.(Int,sp[q+1]+1)]))
                            part1 = ",soso,"solve(prob,",soso,model_def["solver"],"(),reltol=",model_def["tols"][1],",abstol=",model_def["tols"][2],",saveat=1)
                        end
                    end

                    if typeof(p) == Array{Float64, 1}
                        initialV = part1[lts];
                    else
                        initialV = [part1[:,lts][h].value for h in 1:",model_def["nStat"],"];
                    end

                    i+=",model_def["nInp"],";


                    if typeof(p) == Array{Float64, 1}
                        for d in collect((sp[q]+1):(sp[q]+lts))
                            final[convert.(Int,d),:] = part1[convert.(Int, d-sp[q])];
                        end
                    else
                        for d in collect(convert.(Int,(sp[q]+1):(sp[q]+lts)))
    ",join(simsB),"
                        end
                    end

                end

",join(lik),"

            end

        end
    ")

    turingmodel = string(Head, fun1, fun2, fun3);

    open(string(cudi, "\\ModelsFunctions\\", model_def["NameF"], "_TuringModel.jl"), "w") do io
       write(io, turingmodel);
    end;

    include(string(cudi,"\\ModelsFunctions\\", model_def["NameF"], "_TuringModel.jl"))

    turinf_def["ModelPath"] = string(cudi, "\\ModelsFunctions\\", model_def["NameF"], "_TuringModel.jl");

    println("")
    println("----------------------------------------- TURING MODEL GENERATION -----------------------------------------")
    println("The model has been generated in the directory: ")
    println(string("                 ", turinf_def["ModelPath"]))
    println("--------------------------------------------------------------------------------------")
    println("")



    return(turinf_def)
end


##
function getTuringInferenceElements(model_def, turinf_def)


    turinf_def = checkStructTurInf(model_def, turinf_def);
    turinf_def = genTuringModel(model_def, turinf_def);


    modelpath = turinf_def["ModelPath"];
    modelname = string(join(model_def["NameF"]),"_TuringModel(data,prob1)")
    inferdata = turinf_def["Data"];


    println("")
    println("----------------------------------------- TURING INFERENCE INFO -----------------------------------------")
    println("Remember, Turing inference will NOT be initialised. Once you are ready, please run Turing with your desired settings ")
    println("To include your generated functions into your environment just run: ")
    println(replace(string("         include(","\"",turinf_def["ModelPath"],"\"",")"), "\\"=> "\\\\"))
    println("The Turing model function is called: ")
    println("         ",join(model_def["NameF"]),"_TuringModel(data,prob1)")
    println("To automatically save your results, you can use: ")
    println("         saveTuringResults(chains, model_def, turinf_def, turing_res)")
    println(" ")
    println("An example on how to run the turing inference can be: ")
    println("         chainss = sample(model, HMC(0.1, 5), 4000)")
    println("Or if you want to run multiple chains in parallel using NUTS:")
    println("         chainss = sample(model, NUTS(300, .65), MCMCThreads(), 1000, 4, progress=true, init_theta = [0.09,1.4,0.02,5.7])")
    println("For more information visit the Turing.jl webpage https://turing.ml/dev/")
    println("--------------------------------------------------------------------------------------")
    println("")

    return modelpath, modelname, inferdata, model_def, turinf_def

end


##
function saveTuringResults(chainss, model_def, turinf_def)

    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\TuringResults"))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\TuringResults"))
    end

    JLD.save(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\TuringResults\\TuringResults_", model_def["NameF"],"_", turinf_def["flag"], ".jld"),"chainss",chainss)

    poster = chainss.value.data[:,11:11+model_def["nPar"]-1,1];
    for c in 2:size(chainss.value.data)[3]
        poster = vcat(poster, chainss.value.data[:,11:11+model_def["nPar"]-1,c])
    end

    dfP = DataFrame(poster);

    CSV.write(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\TuringResults\\TuringResults_", model_def["NameF"],"_", turinf_def["flag"], "_Posterior.csv"), dfP);

    println("")
    println("----------------------------------------- TURING RESULTS -----------------------------------------")
    println("Turing Inference results are saved in the directory: ")
    println(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\TuringResults"))
    println("--------------------------------------------------------------------------------------")
    println("")

    if turinf_def["plot"] == true
        plotTuringResults(poster,model_def,turinf_def)
        println("")
        println("----------------------------------------- PLOTS -----------------------------------------")
        println("Stan Inference Results PLOTS are saved in the directory: ")
        println(string("                 ", cudi, "\\Results\\", model_def["NameF"],"_",today() ))
        println(string("Under the names PlotTuringInfResults_Exp(i)_",turinf_def["flag"],".png for simulations"))
        println(string("And ", "PosteriorPlotsTur_", turinf_def["flag"], ".png for the posteriors."))
        println("--------------------------------------------------------------------------------------")
        println("")
    end

    return poster

end


##
# function runTuringInference(model_def, turinf_def)
#
#     turinf_def = checkStructTurInf(model_def, turinf_def);
#     turinf_def = genTuringModel(model_def, turinf_def);
#
#
#     modelpath = turinf_def["ModelPath"];
#     modelname = string(join(model_def["NameF"]),"_TuringModel(data,prob1)")
#
#     cudi = pwd(); #@__DIR__;
#     include(string(cudi,"\\ModelsFunctions\\", model_def["NameF"], "_TuringModel.jl"))
#
#     global inferdata = turinf_def["Data"];
#
#     mod2 = Symbol(string(join(model_def["NameF"]),"ODE!"))
# #     @eval $mod2
#     global prob1=ODEProblem((@eval $mod2),inferdata["y0"][1],(0,inferdata["finalTime"][1]),rand(turinf_def["Priors"],1))
#
#     mod = Symbol(string(join(model_def["NameF"]),"_TuringModel"))
#
#     model = @eval $mod(inferdata,prob1)
#
#     if turinf_def["TuringSettings"] == []
#         chainss = sample(model, NUTS(.45), MCMCThreads(), 2000, 4, progress=true);
#     else
#         if turinf_def["TuringSettings"]["acceptRatio"] == []
#             turinf_def["TuringSettings"]["acceptRatio"] = 0.45;
#         end
#
#         if turinf_def["TuringSettings"]["nsamples"] == []
#             turinf_def["TuringSettings"]["nsamples"] = 2000;
#         end
#
#         if turinf_def["TuringSettings"]["printsummary"] == []
#             turinf_def["TuringSettings"]["printsummary"] = true;
#         end
#
#         if turinf_def["TuringSettings"]["nchains"] == []
#             turinf_def["TuringSettings"]["nchains"] = 4;
#         end
#
#         if turinf_def["TuringSettings"]["init"] == []
#             if turinf_def["TuringSettings"]["nadaptations"] == []
#                 chainss = sample(model, NUTS(turinf_def["TuringSettings"]["acceptRatio"]),
#                         MCMCThreads(), turinf_def["TuringSettings"]["nsamples"],
#                         turinf_def["TuringSettings"]["nchains"], progress=turinf_def["TuringSettings"]["printsummary"]);
#             else
#                 chainss = sample(model, NUTS(turinf_def["TuringSettings"]["nadaptations"], turinf_def["TuringSettings"]["acceptRatio"]),
#                         MCMCThreads(), turinf_def["TuringSettings"]["nsamples"],
#                         turinf_def["TuringSettings"]["nchains"], progress=turinf_def["TuringSettings"]["printsummary"]);
#             end
#         else
#             if turinf_def["TuringSettings"]["nadaptations"] == []
#                 chainss = sample(model, NUTS(turinf_def["TuringSettings"]["acceptRatio"]),
#                         MCMCThreads(), turinf_def["TuringSettings"]["nsamples"],
#                         turinf_def["TuringSettings"]["nchains"], progress=turinf_def["TuringSettings"]["printsummary"],
#                         init_theta = turinf_def["TuringSettings"]["init"]);
#             else
#                 chainss = sample(model, NUTS(turinf_def["TuringSettings"]["nadaptations"], turinf_def["TuringSettings"]["acceptRatio"]),
#                         MCMCThreads(), turinf_def["TuringSettings"]["nsamples"],
#                         turinf_def["TuringSettings"]["nchains"], progress=turinf_def["TuringSettings"]["printsummary"],
#                         init_theta = turinf_def["TuringSettings"]["init"]);
#             end
#         end
#     end
#
#     Poster = saveTuringResults(chainss, model_def, turinf_def);
#
#     turinf_res = Dict();
#     turinf_res["chains"] = chainss;
#     turinf_res["posterior"] = Poster;
#     turinf_res["Model"] = modelname;
#     turinf_res["TuringModel"] = modelpath;
#     turinf_res["InferData"] = turinf_def["Data"];
#     turinf_res["prior"] = turinf_def["Priors"];
#     turinf_res["init"] = turinf_def["TuringSettings"]["init"];
#
#     cudi = pwd(); #@__DIR__;
#     if !isdir(string(cudi,"\\Results"))
#         mkdir(string(cudi,"\\Results"))
#     end
#     if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
#         mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
#     end
#
#     JLD.save(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\TuringInferenceResults_", model_def["NameF"],"_", turinf_def["flag"], ".jld"),"turinf_res",turinf_res, "model_def",model_def, "turinf_def", turinf_def)
#
#     println("")
#     println("----------------------------------------- RESULTS -----------------------------------------")
#     println("Turing Inference process results are saved in the directory: ")
#     println(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), ""))
#     println("Under the name: ")
#     println(string("TuringInferenceResults_", model_def["NameF"],"_", turinf_def["flag"], ".jld"))
#     println("--------------------------------------------------------------------------------------")
#     println("")
#
#     return turinf_res, model_def, turinf_def
#
# end


##
function plotTuringResults(poster, model_def, turinf_def)

    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
    end

    # Plot 1, the bivariate distributions
    pc = corrplot(poster, label = model_def["parName"]);
    savefig(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\PosteriorPlotsTur_", turinf_def["flag"], ".png"));

    # Plot 2, the simulations
    # Simulation Results
    simul_def = defSimulStruct()
    simul_def = SimToMle(simul_def, turinf_def["Data"])
    simul_def["plot"] = false
    simul_def["theta"] = convert(Array, poster);
    simul_def["flag"] = "TuringInferResults";

    simuls, model_def, simul_def = simulateODEs(model_def, simul_def);

    for i in 1:turinf_def["Data"]["Nexp"]
        tit = "";
        yl1 = "";
        for k in 1:length(turinf_def["Data"]["Obs"])
            if typeof(turinf_def["Data"]["Obs"]) == Array{Int,1};
                tit = hcat(tit, string(model_def["stName"][turinf_def["Data"]["Obs"][k]]));
            else
                tit = hcat(tit, string(turinf_def["Data"]["Obs"][k]));
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

        SimObs = selectObsSim_te(simuls[string("Exp_",i)], turinf_def["Data"]["Obs"],model_def["stName"]);

        errorss = zeros(size(SimObs)[1],length(turinf_def["Data"]["Obs"]));
        if length(size(turinf_def["Data"]["DataError"][i][1])) == 1
            for k in 1:length(turinf_def["Data"]["Obs"])
                errorss[:,k] = turinf_def["Data"]["DataError"][i][k];
            end
        else
            for k in 1:length(turinf_def["Data"]["Obs"])
                errorss[:,k] = [turinf_def["Data"]["DataError"][i][k][f,f] for f in 1:size(turinf_def["Data"]["DataError"][i][k])[1]];
            end
        end

        # Elements for the inducers
        tu = Array{Array{Any,1}}(undef, length(turinf_def["Data"]["Obs"])+model_def["nInp"])
        [tu[k] = [] for k in 1:length(turinf_def["Data"]["Obs"])]
        su = Array{Array{Any,1}}(undef, length(turinf_def["Data"]["Obs"])+model_def["nInp"])
        [su[k] = [] for k in 1:length(turinf_def["Data"]["Obs"])]
        for k in 1:model_def["nInp"]
            tu[k+length(turinf_def["Data"]["Obs"])] = round.(turinf_def["Data"]["switchT"][i]);
            su[k+length(turinf_def["Data"]["Obs"])] = vcat(turinf_def["Data"]["uInd"][i]'[:,(k)], turinf_def["Data"]["uInd"][i]'[end,(k)])
        end

        pl = plot(round.(turinf_def["Data"]["tsamps"][i]), turinf_def["Data"]["DataMean"][i],yerror = errorss,
        layout=length(turinf_def["Data"]["Obs"])+model_def["nInp"], label = "", title = tit, xlabel = "time", ylabel = yuu,
                color = "gray", size = [2000, 1200]);

        for p in 1:size(SimObs)[3]
            pl = plot!(round.(turinf_def["Data"]["tsamps"][i]),SimObs[:,:,p], layout=length(turinf_def["Data"]["Obs"]), label = "");
        end

        pl = plot!(tu, su, layout=length(turinf_def["Data"]["Obs"])+model_def["nInp"], label = "", linetype=:step, title = tuu);

        pl = plot!(round.(turinf_def["Data"]["tsamps"][i]), turinf_def["Data"]["DataMean"][i],yerror = errorss,
        layout=length(turinf_def["Data"]["Obs"])+model_def["nInp"], label = "", title = tit, xlabel = "time", ylabel = yuu,
                color = "gray", size = [2000, 1200]);

        savefig(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\PlotTuringInfResults_Exp", i,"_", turinf_def["flag"], ".png"));
    end


end


##
function TuringInfer(model_def, turinf_def)

    cudi = pwd(); #@__DIR__;
    turinf_def = checkStructTurInf(model_def, turinf_def);

    # if turinf_def["runInf"] == false

        # println("----------------------------------------- TURING INFERENCE INFO -----------------------------------------")
        # println("No inference will be performed since this is what has been selected. Structure elements")
        # println("to run the Turing inference will be given instead.")
        # println("--------------------------------------------------------------------------------------")

    println("----------------------------------------- TURING INFERENCE INFO -----------------------------------------")
    println("Remember, NO inference will be performed within BOMBs. Structure elements")
    println("to run the Turing inference will be given instead.")
    println("--------------------------------------------------------------------------------------")

    modelpath, modelname, inferdata, model_def, turinf_def = getTuringInferenceElements(model_def, turinf_def);
    turing_struct = Dict();
    turing_struct["modelpath"] = modelpath;
    turing_struct["modelname"] = modelname;
    turing_struct["inferdata"] = inferdata;
    # if isempty(turinf_def["TuringSettings"])
    #     turing_struct["init"] = [];
    # else
    #     turing_struct["init"] = turinf_def["TuringSettings"]["init"];
    # end
    turing_struct["model_def"] = model_def;
    turing_struct["turinf_def"] = turinf_def;
    mod2 = Symbol(string(join(model_def["NameF"]),"ODE!"));
    turing_struct["probExamp"]=ODEProblem((@eval $mod2),inferdata["y0"][1],(0,inferdata["finalTime"][1]),rand(turinf_def["Priors"],1))

    # elseif turinf_def["runInf"] == true
    #     println("----------------------------------------- TURING INFERENCE INFO -----------------------------------------")
    #     println("Inference will be performed since there is enough information from the user:")
    #     println("Turing temporary files will be stored in your current working directory in the folder tmp.")
    #     println("--------------------------------------------------------------------------------------")
    #
    #     turing_struct, model_def, turinf_def = runTuringInference(model_def, turinf_def);
    #
    #     if turinf_def["plot"] == true
    #         plotTuringResults(turing_struct,model_def,turinf_def)
    #         println("")
    #         println("----------------------------------------- PLOTS -----------------------------------------")
    #         println("Stan Inference Results PLOTS are saved in the directory: ")
    #         println(string("                 ", cudi, "\\Results\\", model_def["NameF"],"_",today() ))
    #         println(string("Under the names PlotTuringInfResults_Exp(i)_",turinf_def["flag"],".png for simulations"))
    #         println(string("And ", "PosteriorPlotsTur_", turinf_def["flag"], ".png for the posteriors."))
    #         println("--------------------------------------------------------------------------------------")
    #         println("")
    #     end
    #
    # end


    return turing_struct, model_def, turinf_def
end
