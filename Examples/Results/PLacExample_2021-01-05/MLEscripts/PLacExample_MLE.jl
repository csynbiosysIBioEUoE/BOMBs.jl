
    using Distributed
    addprocs(length(Sys.cpu_info())-1)

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



    @everywhere function restructInputs(model_def, mle_def, expp)
        r = convert.(Int, 1:model_def["nInp"]:(length(mle_def["uInd"][expp])));
        inputs = zeros(convert.(Int,length(mle_def["uInd"][expp])));
        for j in 1:convert.(Int,length(mle_def["uInd"][expp])/model_def["nInp"])
            for k in 0:(model_def["nInp"]-1)
                inputs[r[j]+k] = mle_def["uInd"][expp][j,(k+1)];
            end
        end
        return(inputs)
    end



    function RunMLEPLacExample(model_def, mle_def)

        @spawnat 1 mle_def
        @spawnat 2 mle_def
        @spawnat 3 mle_def
        @spawnat 4 mle_def
        @spawnat 5 mle_def
        @spawnat 6 mle_def
        @spawnat 7 mle_def
        @spawnat 8 mle_def
        @spawnat 9 mle_def
        @spawnat 10 mle_def
        @spawnat 11 mle_def
        @spawnat 12 mle_def

        @spawnat 1 model_def
        @spawnat 2 model_def
        @spawnat 3 model_def
        @spawnat 4 model_def
        @spawnat 5 model_def
        @spawnat 6 model_def
        @spawnat 7 model_def
        @spawnat 8 model_def
        @spawnat 9 model_def
        @spawnat 10 model_def
        @spawnat 11 model_def
        @spawnat 12 model_def


        @everywhere global Nexp = mle_def["Nexp"];
        @everywhere global Obs = mle_def["Obs"];
        @everywhere global stName = model_def["stName"];

        @everywhere global dataM = mle_def["DataMean"];
        @everywhere global dataE = mle_def["DataError"];

        @everywhere ts = Array{Any,1}(undef,mle_def["Nexp"]);
        @everywhere sp = Array{Any,1}(undef,mle_def["Nexp"]);
        @everywhere ivss = Array{Any,1}(undef,mle_def["Nexp"]);
        @everywhere pre = Array{Any,1}(undef,mle_def["Nexp"]);
        @everywhere samps = Array{Any,1}(undef,mle_def["Nexp"]);
        @everywhere inputs = Array{Any,1}(undef,mle_def["Nexp"]);

        for i in 1:mle_def["Nexp"]
            global ts[i] = collect(0.0:round(mle_def["finalTime"][i]))';
            global sp[i] = [convert(Int,v) for v in (round.(mle_def["switchT"][i])')];
            global ivss[i] = mle_def["y0"][i];
            if model_def["Y0eqs"] != [] # ON inputs
                global pre[i] = mle_def["preInd"][i];
            else
                global pre[i] = [];
            end
            global samps[i] = convert.(Int, round.(mle_def["tsamps"][i]));
            global inputs[i] = restructInputs(model_def, mle_def, i);
        end

        @spawnat 1 ts
        @spawnat 1 sp
        @spawnat 1 ivss
        @spawnat 1 pre
        @spawnat 1 samps
        @spawnat 1 inputs
        @spawnat 2 ts
        @spawnat 2 sp
        @spawnat 2 ivss
        @spawnat 2 pre
        @spawnat 2 samps
        @spawnat 2 inputs
        @spawnat 3 ts
        @spawnat 3 sp
        @spawnat 3 ivss
        @spawnat 3 pre
        @spawnat 3 samps
        @spawnat 3 inputs
        @spawnat 4 ts
        @spawnat 4 sp
        @spawnat 4 ivss
        @spawnat 4 pre
        @spawnat 4 samps
        @spawnat 4 inputs
        @spawnat 5 ts
        @spawnat 5 sp
        @spawnat 5 ivss
        @spawnat 5 pre
        @spawnat 5 samps
        @spawnat 5 inputs
        @spawnat 6 ts
        @spawnat 6 sp
        @spawnat 6 ivss
        @spawnat 6 pre
        @spawnat 6 samps
        @spawnat 6 inputs
        @spawnat 7 ts
        @spawnat 7 sp
        @spawnat 7 ivss
        @spawnat 7 pre
        @spawnat 7 samps
        @spawnat 7 inputs
        @spawnat 8 ts
        @spawnat 8 sp
        @spawnat 8 ivss
        @spawnat 8 pre
        @spawnat 8 samps
        @spawnat 8 inputs
        @spawnat 9 ts
        @spawnat 9 sp
        @spawnat 9 ivss
        @spawnat 9 pre
        @spawnat 9 samps
        @spawnat 9 inputs
        @spawnat 10 ts
        @spawnat 10 sp
        @spawnat 10 ivss
        @spawnat 10 pre
        @spawnat 10 samps
        @spawnat 10 inputs
        @spawnat 11 ts
        @spawnat 11 sp
        @spawnat 11 ivss
        @spawnat 11 pre
        @spawnat 11 samps
        @spawnat 11 inputs
        @spawnat 12 ts
        @spawnat 12 sp
        @spawnat 12 ivss
        @spawnat 12 pre
        @spawnat 12 samps
        @spawnat 12 inputs


        @everywhere global tim = ts;
        @everywhere global esp = sp;
        @everywhere global inp = inputs;
        @everywhere global ini = ivss;
        @everywhere global prr = pre;
        @everywhere global smp = samps;


        @everywhere lb = mle_def["thetaMIN"];
        @everywhere up = mle_def["thetaMAX"];

        @everywhere rang = [Array{Tuple{Int, Int}}(undef, length(lb))];
        @everywhere rang = [(lb[i], up[i]) for i in 1:length(lb)];

        @everywhere chains = mle_def["runs"];
        @everywhere opts = Dict()


        @everywhere convcur = Array{Any,1}(undef, chains);
        for i in 1:chains
            convcur[i] = Array{Tuple{Int, Float64},1}();
            callback = oc -> push!(convcur[i], (num_func_evals(oc), best_fitness(oc)))
            opts[string("Opt_", i)] = bbsetup(objectiveMLEPLacExample; SearchRange = rang, NumDimensions = length(lb),
                Method = :adaptive_de_rand_1_bin_radiuslimited, MaxTime = 10, TraceMode = :silent, CallbackFunction = callback, CallbackInterval = 0.0);
        end


        println("----------------------------------------- OPTIMISATION STARTS! -----------------------------------------")

        @time begin
            tet, convcur2 = parOptim(opts);

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
        @everywhere convcur = Array{Any}(undef,chains);
        for k in 1:chains
            tmpch = convcur3[k][1][2];
            convcur[tmpch] = convcur3[k][2:end];
        end



        end

        println("----------------------------------------- OPTIMISATION ENDED -----------------------------------------")


        names = model_def["parName"];

        alltog = Array{Dict{String,Any},1}(undef,chains)
        for i in 1:(chains)
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
