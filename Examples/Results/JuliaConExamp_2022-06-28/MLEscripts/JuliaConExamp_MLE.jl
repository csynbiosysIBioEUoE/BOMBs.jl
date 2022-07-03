
    using Distributed
    

    using DifferentialEquations
    using OrdinaryDiffEq
    using DiffEqBase
    using Sundials
    using Plots
    # using DiffEqGPU
    using ODEInterfaceDiffEq
    using DataFrames
    using CSV
    using Statistics
    using LinearAlgebra
    using JLD
    using StatsBase
    using Random
    using SharedArrays
    using BlackBoxOptim
    using BlackBoxOptim: num_func_evals
        
    include("e:\\UNI\\D_Drive\\PhD\\Year_1\\2022_06_23_JuliaCon2022\\ModelsFunctions\\JuliaConExamp_Model.jl")
        

    function selectObsSim(simul)

        simObs = zeros(size(simul)[1], length(Obs));
        if typeof(Obs)==Array{Int,1}
            simObs = simul[:,(Obs),1];
        else
            coun = 1;
            for j in 1:length(Obs)
               if !(occursin.("+", Obs[j])) && !(occursin.("-", Obs[j])) &&
               !(occursin.("/", Obs[j])) && !(occursin.("*", Obs[j])) &&
                !(occursin.("^", Obs[j]))
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


    
    function parOptim(opts)
        ress = convert(SharedArray{Float64,2}, zeros(4,length(opts)));
        for i in 1:length(opts)
            resT = bboptimize(opts[string("Opt_", i)]);
            ress[:,i] = resT.archive_output.best_candidate;
        end;
        return(ress);
    end
    
    function loglike(dats, mes, errs)

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

    

    function objectiveMLEJuliaConExamp(theta)

        LLK_Exps = zeros(1, Nexp);
        simobs = Array{Any,1}(undef,Nexp);
        for i in 1:Nexp
            simul = JuliaConExamp_SolveAll(tim[i], theta, esp[i], inp[i], ini[i], smp[i], prr[i])
            simobs[i] = selectObsSim(simul);
            llkobs = zeros(1,length(Obs));
            for j in 1:length(Obs)
                llkobs[1,j] = loglike(dataM[i][:,j], simobs[i][:,j], dataE[i][j]);
            end
            LLK_Exps[1,i] = sum(llkobs);
        end
        return(-sum(LLK_Exps))
    end

    

    function restructInputs(model_def, mle_def, expp)
        if model_def["nInp"] != 0
            r = convert.(Int, 1:model_def["nInp"]:(length(mle_def["uInd"][expp])));
            inputs = zeros(convert.(Int,length(mle_def["uInd"][expp])));
            for j in 1:convert.(Int,length(mle_def["uInd"][expp])/model_def["nInp"])
                for k in 0:(model_def["nInp"]-1)
                    inputs[r[j]+k] = mle_def["uInd"][expp][j,(k+1)];
                end
            end
        else
            inputs = [];
        end
        return(inputs)
    end

    

    function RunMLEJuliaConExamp(model_def2, mle_def2)




        global Nexp = mle_def2["Nexp"];
        global Obs = mle_def2["Obs"];
        global stName = model_def2["stName"];

        global dataM = mle_def2["DataMean"];
        global dataE = mle_def2["DataError"];

        global ts = Array{Any,1}(undef,mle_def2["Nexp"]);
        global sp = Array{Any,1}(undef,mle_def2["Nexp"]);
        global ivss = Array{Any,1}(undef,mle_def2["Nexp"]);
        global pre = Array{Any,1}(undef,mle_def2["Nexp"]);
        global samps = Array{Any,1}(undef,mle_def2["Nexp"]);
        global inputs = Array{Any,1}(undef,mle_def2["Nexp"]);

        for i in 1:mle_def2["Nexp"]
            ts[i] = collect(0.0:round(mle_def2["finalTime"][i]))';
            sp[i] = [convert(Int,v) for v in (round.(mle_def2["switchT"][i])')];
            ivss[i] = mle_def2["y0"][i];
            if model_def2["Y0eqs"] != [] # ON inputs
                pre[i] = mle_def2["preInd"][i];
            else
                pre[i] = [];
            end
            samps[i] = convert.(Int, round.(mle_def2["tsamps"][i]));
            inputs[i] = restructInputs(model_def2, mle_def2, i);
        end

        global tim = ts;
        global esp = sp;
        global inp = inputs;
        global ini = ivss;
        global prr = pre;
        global smp = samps;


        lb = mle_def2["thetaMIN"];
        up = mle_def2["thetaMAX"];

        rang = [Array{Tuple{Int, Int}}(undef, length(lb))];
        rang = [(lb[i], up[i]) for i in 1:length(lb)];

        opts = Dict()


        convcur = Array{Any,1}(undef, 1);
        for i in 1:1
            convcur[i] = Array{Tuple{Int, Float64},1}();
            callback = oc -> push!(convcur[i], (num_func_evals(oc), best_fitness(oc)))
            opts[string("Opt_", i)] = bbsetup(objectiveMLEJuliaConExamp; SearchRange = rang, NumDimensions = length(lb),
                Method = :adaptive_de_rand_1_bin_radiuslimited, MaxTime = 120, TraceMode = :silent, CallbackFunction = callback, CallbackInterval = 0.0);
        end


        println("----------------------------------------- OPTIMISATION STARTS! -----------------------------------------")

        @time begin
            tet= parOptim(opts);
        end

        println("----------------------------------------- OPTIMISATION ENDED -----------------------------------------")


        names = model_def2["parName"];

        alltog = Array{Dict{String,Any},1}(undef,1)
        for i in 1:(1)
            global tmp = Dict{String, Any}()
            for j in 1:length(names)
                tmp[names[j]] = tet[j,i];
            end
            alltog[i] = tmp
        end

        mle_res = Dict();
        mle_res["Theta"] = convert(Array,tet); # Best thetas results
        mle_res["convCurv"] = convcur; # Convergence curves
        mle_res["StanDict"] = alltog; # Best thetas in the structure Stan needs as initial guess.

        bcfv = zeros(length(convcur))
        for k in 1:length(convcur)
            bcfv[k] = convcur[k][end][2]
        end
        mle_res["BestCFV"] = bcfv; # Best cost function values for each run

        tpin = findfirst(x->x==minimum(bcfv), bcfv);
        mle_res["BestTheta"] = tet[:,tpin]; # Best theta value from all runs

        return(mle_res)

    end
    