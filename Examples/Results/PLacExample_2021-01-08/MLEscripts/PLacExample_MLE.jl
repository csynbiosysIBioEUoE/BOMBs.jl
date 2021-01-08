
    using Distributed
    

    @everywhere using DifferentialEquations
    @everywhere using OrdinaryDiffEq
    using DiffEqBase
    @everywhere using Sundials
    using Plots
    # using DiffEqGPU
    @everywhere using ODEInterfaceDiffEq
    using DataFrames
    using CSV
    @everywhere using Statistics
    @everywhere using LinearAlgebra
    using JLD
    @everywhere using StatsBase
    @everywhere using Random
    using SharedArrays
    @everywhere using BlackBoxOptim
    @everywhere using BlackBoxOptim: num_func_evals
        
    @everywhere include("E:\\UNI\\D_Drive\\PhD\\JULIAlang\\Generalisation_PLacExample\\Examples\\ModelsFunctions\\PLacExample_Model.jl")
        

    @everywhere function selectObsSim(simul)

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


    
    @everywhere function parOptim(opts)
        ress = convert(SharedArray{Float64,2}, zeros(9,length(opts)));
        arr = @sync @distributed (vcat) for i in 1:length(opts)
            resT = bboptimize(opts[string("Opt_", i)]);
            ress[:,i] = resT.archive_output.best_candidate;

            chr = (0,i);
            tmpcv = vcat(chr, convcur[i]);
            tmpcv
        end;
        return ress, arr;
    end

    
    @everywhere function loglike(dats, mes, errs)

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

    

    @everywhere function objectiveMLEPLacExample(theta)

        LLK_Exps = zeros(1, Nexp);
        simobs = Array{Any,1}(undef,Nexp);
        for i in 1:Nexp
            simul = PLacExample_SolveAll(tim[i], theta, esp[i], inp[i], ini[i], smp[i], prr[i])
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
        r = convert.(Int, 1:model_def["nInp"]:(length(mle_def["uInd"][expp])));
        inputs = zeros(convert.(Int,length(mle_def["uInd"][expp])));
        for j in 1:convert.(Int,length(mle_def["uInd"][expp])/model_def["nInp"])
            for k in 0:(model_def["nInp"]-1)
                inputs[r[j]+k] = mle_def["uInd"][expp][j,(k+1)];
            end
        end
        return(inputs)
    end

    

    function RunMLEPLacExample(model_def2, mle_def2)

        @everywhere @spawnat 2 mle_def2 
        @everywhere @spawnat 3 mle_def2 
        @everywhere @spawnat 4 mle_def2 
        @everywhere @spawnat 5 mle_def2 
        @everywhere @spawnat 6 mle_def2 
        @everywhere @spawnat 7 mle_def2 
        @everywhere @spawnat 8 mle_def2 
        @everywhere @spawnat 9 mle_def2 
        @everywhere @spawnat 10 mle_def2 
        @everywhere @spawnat 11 mle_def2 
        @everywhere @spawnat 12 mle_def2 

        @everywhere @spawnat 2 model_def2 
        @everywhere @spawnat 3 model_def2 
        @everywhere @spawnat 4 model_def2 
        @everywhere @spawnat 5 model_def2 
        @everywhere @spawnat 6 model_def2 
        @everywhere @spawnat 7 model_def2 
        @everywhere @spawnat 8 model_def2 
        @everywhere @spawnat 9 model_def2 
        @everywhere @spawnat 10 model_def2 
        @everywhere @spawnat 11 model_def2 
        @everywhere @spawnat 12 model_def2 


        @everywhere global Nexp = mle_def2["Nexp"];
        @everywhere global Obs = mle_def2["Obs"];
        @everywhere global stName = model_def2["stName"];

        @everywhere global dataM = mle_def2["DataMean"];
        @everywhere global dataE = mle_def2["DataError"];

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
        @spawnat 2 tim 
        @spawnat 2 esp 
        @spawnat 2 ini 
        @spawnat 2 prr 
        @spawnat 2 smp 
        @spawnat 2 inp 
        @spawnat 3 tim 
        @spawnat 3 esp 
        @spawnat 3 ini 
        @spawnat 3 prr 
        @spawnat 3 smp 
        @spawnat 3 inp 
        @spawnat 4 tim 
        @spawnat 4 esp 
        @spawnat 4 ini 
        @spawnat 4 prr 
        @spawnat 4 smp 
        @spawnat 4 inp 
        @spawnat 5 tim 
        @spawnat 5 esp 
        @spawnat 5 ini 
        @spawnat 5 prr 
        @spawnat 5 smp 
        @spawnat 5 inp 
        @spawnat 6 tim 
        @spawnat 6 esp 
        @spawnat 6 ini 
        @spawnat 6 prr 
        @spawnat 6 smp 
        @spawnat 6 inp 
        @spawnat 7 tim 
        @spawnat 7 esp 
        @spawnat 7 ini 
        @spawnat 7 prr 
        @spawnat 7 smp 
        @spawnat 7 inp 
        @spawnat 8 tim 
        @spawnat 8 esp 
        @spawnat 8 ini 
        @spawnat 8 prr 
        @spawnat 8 smp 
        @spawnat 8 inp 
        @spawnat 9 tim 
        @spawnat 9 esp 
        @spawnat 9 ini 
        @spawnat 9 prr 
        @spawnat 9 smp 
        @spawnat 9 inp 
        @spawnat 10 tim 
        @spawnat 10 esp 
        @spawnat 10 ini 
        @spawnat 10 prr 
        @spawnat 10 smp 
        @spawnat 10 inp 
        @spawnat 11 tim 
        @spawnat 11 esp 
        @spawnat 11 ini 
        @spawnat 11 prr 
        @spawnat 11 smp 
        @spawnat 11 inp 
        @spawnat 12 tim 
        @spawnat 12 esp 
        @spawnat 12 ini 
        @spawnat 12 prr 
        @spawnat 12 smp 
        @spawnat 12 inp 


        @everywhere lb = mle_def2["thetaMIN"];
        @everywhere up = mle_def2["thetaMAX"];

        @everywhere rang = [Array{Tuple{Int, Int}}(undef, length(lb))];
        @everywhere rang = [(lb[i], up[i]) for i in 1:length(lb)];

        @everywhere opts = Dict()


        @everywhere convcur = Array{Any,1}(undef, 3);
        for i in 1:3
            convcur[i] = Array{Tuple{Int, Float64},1}();
            callback = oc -> push!(convcur[i], (num_func_evals(oc), best_fitness(oc)))
            opts[string("Opt_", i)] = bbsetup(objectiveMLEPLacExample; SearchRange = rang, NumDimensions = length(lb),
                Method = :adaptive_de_rand_1_bin_radiuslimited, MaxTime = 10, TraceMode = :silent, CallbackFunction = callback, CallbackInterval = 0.0);
        end


        println("----------------------------------------- OPTIMISATION STARTS! -----------------------------------------")

        @time begin
            tet, convcur2 = parOptim(opts);

        convcur3 = Array{Any}(undef,3)
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
        @everywhere convcur = Array{Any}(undef,3);
        for k in 1:3
            tmpch = convcur3[k][1][2];
            convcur[tmpch] = convcur3[k][2:end];
        end


        
        end

        println("----------------------------------------- OPTIMISATION ENDED -----------------------------------------")


        names = model_def2["parName"];

        alltog = Array{Dict{String,Any},1}(undef,3)
        for i in 1:(3)
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
    