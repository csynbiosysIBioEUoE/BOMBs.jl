##
function defBayInfStruct()
    bayinf_def = Dict();
    bayinf_def["Priors"] = []; # Five options for this:
                            # 1) A 2*N array containing the bounds for the parameters. Order of parameters will be assumed as
                                    # the one introduced in model_def["parName"]. As a prior a truncated Normal covering 2
                                    # standard deviations in the bounds given will be generated as prior.
                            # 2) Path to a CSV file containing samples for the parameters. Fitting of the samples to
                                    # different type of distributions will be done to generate the priors. Order of parameters
                                    # will be assumed as the one introduced in model_def["parName"].
                            # 3) Dictionary with fields pars, transpars and pridis defining the
                                    # parameters, subsequent desired transformations and prior distributions.
                                    # In the field pars, parameters
                                    # have to be defined with the same names and order as in model_def["parName"],
                                    # otherwise the script will not proceed.
                            # 4) An empty array. If this is the case, the stan model will be generated with nothing in
                                    # the parameters and transformed parameters section. The path to the stan file will
                                    # be given so the user can fill these sections.
                            # 5) Path to a Stan model file if you already have one. For example, if first you introduce an
                                    # empty array, then you can fill the parameters section and run from that script instead
                                    # of having to copy the things in here.

    bayinf_def["Data"] = []; # Two options:
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

    bayinf_def["StanSettings"] = []; # Two options:
                            # 1) A dictionary with the basic fields from a stan run. The structure of the dictionary can be
                                # extracted calling the function defBasicStanSettingsStruct().
                            # 2) An empty array. If this is the case, no inference will be done after calling the main
                                # function. Instead, you will be given a StanModel file and data structure and you will
                                # have to use this to run a call of Stan by yourself. An example will be provided.

    bayinf_def["flag"] = []; # String to attach a unique flag to the generated files so it is not overwritten. If empty, nothing will be added.

    bayinf_def["plot"] = []; # true, flase, "Yes", "yes", "No", "no" indicating if plots with the results will be generated. Default is false.

    bayinf_def["runInf"] = []; # true, flase, "Yes", "yes", "No", "no" indicating if inference wants to be performed or only stan structure is given. Default is true (if all the necessary elements are present)

    bayinf_def["MultiNormFit"] = [];

    return(bayinf_def)
end

##
function defBayInfDataStruct()
    data_def = Dict();

    data_def["Nexp"] = []; # Integer indicating the number of experiments to be simulated
    data_def["finalTime"] = []; # -> Vector of final times for each simulation (initial time will allways be asumed as 0, so please consider that)
    data_def["switchT"] = []; # Array with the switching times of the inducer in the simulation (time 0 and final time need to be considered)
    data_def["y0"] = []; # Array of Y0s for the simulations for each experiment. If you are computing the steady state this vector might not be used, however you still need to introduce it
    data_def["preInd"] = []; # Level of the inducer in the ON. This entry can be empty if no Steady State equations are given and only the Y0 given by the user is considered.
    data_def["uInd"] = []; # Array with the values for the inducer for each experiment at each step
    data_def["tsamps"] = []; # Array of Sampling times vectors

    data_def["DataMean"] = []; # Array containin ght evector of means for each experiment.
    data_def["DataError"] = []; # Array containing the vector or matrices (covariance included) of errors for the data for each experiment
    data_def["Obs"] = []; # Vector containing the observables for the data.

    return(data_def)
end
##
function defBayInfDataFromFilesStruct()
    data_def = Dict();

    data_def["Observables"] = []; # Array containin ght evector of means for each experiment.
    data_def["Inputs"] = []; # Array containing the vector or matrices (covariance included) of errors for the data for each experiment
    data_def["Obs"] = []; # Vector containing the observables for the data.
    data_def["y0"] = []; # Array of Y0s for the simulations for each experiment. If you are computing the steady state this vector might not be used, however you still need to introduce it

    return(data_def)
end
##
function defBasicStanSettingsStruct()
    stan_def = Dict();

    stan_def["cmdstan_home"] = []; # Full path to the installation directory for cmdstan
    stan_def["nchains"] = []; # Number of chains for the inference
    stan_def["nsamples"] = []; # Number of post-warmup samples for each chain
    stan_def["nwarmup"] = []; # Number of warm-up samples for each chain.
    stan_def["printsummary"] = []; # Print sumary of inference at the end. This can be either true or false. Default will be true.
    stan_def["init"] = []; # Initial point for the parameters in each chain. See output of MLE results. This field can be empty. Please introduce the parameter values in their true parameter range
    stan_def["maxdepth"] = []; # Maximum tree-depth. This field can be empty. Check Stan documentation for more information.
    stan_def["adaptdelta"] = []; # Delta value between 0 and 1. This field can be empty.  Check Stan documentation for more information.
    stan_def["jitter"] = []; # Jitter value between 0 and 1. This field can be empty.  Check Stan documentation for more information.

    return(stan_def)
end

##
# ma -> maximum true value
# mi -> minimum true value
# bo -> bottom transformed range
# up -> upper transformed range
function convertBoundTo2(x, bo, up)

    if bo>=up
        println("Please, give correct bounds")
        return
    end
    a1, a2, b1, b2 = minimum(x), maximum(x), bo, up

    t = b1.+ ((x.-a1).*(b2-b1))./(a2-a1)
    return(t)

end

    # function convertBoundToReal(x, bo, up)
    #
    #     if bo>=up
    #         println("Please, give correct bounds")
    #         return
    #     end
    #     a1, a2, b1, b2 = bo, up, minimum(x), maximum(x)
    #
    #     t = b1.+ ((x.-a1).*(b2-b1))./(a2-a1)
    #     return(t)
    #
    # end

##
function fitPriorSamps(priorsamps, model_def)
    priorfit = Dict();

    if length(size(priorsamps)) == 2
        if size(priorsamps)[1] == size(priorsamps)[2]
            println("-------------------------- WARNING --------------------------")
            println("Sorry, but the number of rows and columns of the theta matrix is the same, so the checks on correct ")
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

    Posterior2 = zeros(size(priorsamps)); # Temporary mapped samples to assess beta distribution with the rest
    for j in 1:model_def["nPar"]
        Posterior2[:,j] = convertBoundTo2(priorsamps[:,j], 0.1, 1)
    end

    dists = [Beta, Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform];
    dists2 = [Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform]; # For the parameters that are not between 0 and 1 (to avoid error)

    fitts = Dict(); # Where all the distribution fits will be stored
    bestfit = Dict(); # Where the best distribution fit will be stored
    bestfitInd = zeros(1,model_def["nPar"]); # Index for the best distribution fit for each parameter

    fittsNorm = Dict();
    bestfitNorm = Dict();
    bestfitIndNorm = zeros(1,model_def["nPar"]);

    names = model_def["parName"];

    for i in 1:length(names)
        fittsNorm[names[i]] = fit.(dists, Ref(Posterior2[:,i]));
        bestfitIndNorm[i] = findmax(loglikelihood.(fittsNorm[names[i]], Ref(Posterior2[:,i])))[2];
        bestfitNorm[names[i]] = fittsNorm[names[i]][findmax(loglikelihood.(fittsNorm[names[i]], Ref(Posterior2[:,i])))[2]];
    end

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

    newPri = Vector{String}(undef,model_def["nPar"]);
    parDef = Vector{String}(undef,model_def["nPar"]);
    transPar = Vector{String}(undef,model_def["nPar"]);

    for i in 1:model_def["nPar"]
        parDef[i] = string("real ", names[i], "; \n");
        transPar[i] = string("theta[", i, "] = ", names[i], "; \n");

        if bestfitInd[i] == 1 # Beta
            conSamp = convertBoundTo2(priorsamps[:,i], 0, 1) # Map samples to appropiate range
            fitSamp = fit(dists[1], conSamp) # Fit distribution to mapped samples
            newPri[i] = string(names[i], " ~ beta(", fitSamp.α, ", ", fitSamp.β, "); \n"); # Define string with pdf for prior
            transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",1," - (",0,")); \n"); # Define remapping to true range of values

        elseif bestfitInd[i] == 2 # Exponential
            conSamp = convertBoundTo2(priorsamps[:,i], 0, 2)
            fitSamp = fit(dists[2], conSamp)
            newPri[i] = string(names[i], " ~ exponential(", fitSamp.θ, "); \n");
            transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0,")); \n");

        elseif bestfitInd[i] == 3 # LogNormal
            parDef[i] = string("real ", names[i], "; \n");
            transPar[i] = string("theta[",i,"] = exp(((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n"); # Map logNormal to univariate Normal. This is how to transform the samples back
            newPri[i] = string(names[i], " ~ normal(", 0, ", ", 1, "); \n"); # Univariate normal

        elseif bestfitInd[i] == 4 # Normal
            if ((0-bestfit[names[i]].μ)/bestfit[names[i]].σ) >=-2 # Truncation only if the actual 0 for the parameter is inside the defined distribution.
                parDef[i] = string("real<lower=",(0-bestfit[names[i]].μ)/bestfit[names[i]].σ,"> ", names[i], "; \n");
            else
                parDef[i] = string("real ", names[i], "; \n");
            end
            transPar[i] = string("theta[",i,"] = (((",names[i],")*(",bestfit[names[i]].σ,"))+(",bestfit[names[i]].μ,")); \n") # Map Normal to univariate Normal. This is how to transform the samples back
            newPri[i] = string(names[i], " ~ normal(", 0, ", ", 1, "); \n"); # Univariate normal

        elseif bestfitInd[i] == 5 # Gamma
            conSamp = convertBoundTo2(priorsamps[:,i], 0.001, 2)
            fitSamp = fit(dists[5], conSamp)
            newPri[i] = string(names[i], " ~ gamma(", fitSamp.α, ", ", fitSamp.θ, "); \n");
            transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0.001,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0.001,")); \n");

        elseif bestfitInd[i] == 6 # Laplace
            conSamp = convertBoundTo2(priorsamps[:,i], -2, 2)
            fitSamp = fit(dists[6], conSamp)
            newPri[i] = string(names[i], " ~ double_exponential(", fitSamp.μ, ", ", fitSamp.θ, "); \n");
            transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");

        elseif bestfitInd[i] == 7 # Pareto
            conSamp = convertBoundTo2(priorsamps[:,i], 0.001, 2)
            fitSamp = fit(dists[7], conSamp)
            newPri[i] = string(names[i], " ~ pareto(", fitSamp.θ, ", ", fitSamp.α, "); \n");
            transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0.001,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0.001,")); \n");

        elseif bestfitInd[i] == 8 # Rayleigh
            conSamp = convertBoundTo2(priorsamps[:,i], 0, 2)
            fitSamp = fit(dists[8], conSamp)
            newPri[i] = string(names[i], " ~ rayleigh(", fitSamp.σ, "); \n");
            transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",0,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",0,")); \n");

        elseif bestfitInd[i] == 9 # Cauchy
            conSamp = convertBoundTo2(priorsamps[:,i], -2, 2)
            fitSamp = fit(dists[9], conSamp)
            newPri[i] = string(names[i], " ~ cauchy(", fitSamp.μ, ", ", fitSamp.σ, "); \n");
            transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");

        elseif bestfitInd[i] == 10
            conSamp = convertBoundTo2(priorsamps[:,i], -2, 2)
            fitSamp = fit(dists[10], conSamp)
            newPri[i] = string(names[i], " ~ uniform(", -2, ", ", 2, "); \n");
            transPar[i] = string("theta[",i,"] = ",minimum(priorsamps[:,i]),"+ ((",names[i]," - (",-2,"))*(",maximum(priorsamps[:,i])," - ",minimum(priorsamps[:,i]),"))/(",2," - (",-2,")); \n");

        end
    end

    pushfirst!(transPar, "real theta[nParms]; \n");

    priorfit["pars"] = parDef;
    priorfit["transpars"] = transPar;
    priorfit["pridis"] = newPri;

    return(priorfit)
end

##
function fitPriorSampsMultiNorm(priorsamps, model_def)

    priorfit = Dict();

    if length(size(priorsamps)) == 2
        if size(priorsamps)[1] == size(priorsamps)[2]
            println("-------------------------- WARNING --------------------------")
            println("Sorry, but the number of rows and columns of the theta matrix is the same, so the checks on correct ")
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

    dists = [Uniform, Normal, LogNormal];

    fitts = Dict(); # Where all the distribution fits will be stored
    bestfit = Dict(); # Where the best distribution fit will be stored
    bestfitInd = zeros(1,model_def["nPar"]); # Index for the best distribution fit for each parameter

    names = model_def["parName"];

    for i in 1:length(names)
        fitts[names[i]] = fit.(dists, Ref(priorsamps[:,i]));
        bestfitInd[i] = findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2];
        bestfit[names[i]] = fitts[names[i]][findmax(loglikelihood.(fitts[names[i]], Ref(priorsamps[:,i])))[2]];
    end

    s = 1:model_def["nPar"]
    lognorms = s[findall(x->x==3, reshape(convert.(Int, bestfitInd), model_def["nPar"]))]
    uninorms = s[findall(x->x==1, reshape(convert.(Int, bestfitInd), model_def["nPar"]))]

    # Transform samples to log values
    Posterior2 = zeros(size(priorsamps));
    Posterior2[:,:] = priorsamps[:,:];
    for i in 1:length(lognorms)
        Posterior2[:,lognorms[i]] = log.(Posterior2[:,lognorms[i]]);
    end

    me = mean(Matrix(Posterior2[:,:]),dims=1);
    sd = cov(Matrix(Posterior2[:,:]));

    # Remove uniform parameters
    Pos = Posterior2;
    Pos = Pos[:,setdiff(1:end, uninorms)];

    unst = Vector{String}(undef,length(uninorms)+1)
    for i in 1:length(uninorms)
        bou = convertBoundTo2(priorsamps[:,uninorms[i]], -2, 2);
        unst[i] = string("real theta_Un",i,"; \n");
    end
    unst[end] = string("vector[nParms-",length(uninorms),"] theta2; \n");

    Posterior3 = zeros(size(Posterior2));
    for i in 1:size(Posterior2)[2]
       Posterior3[:,i] = (Posterior2[:,i].-me[i])./sqrt(sd[i,i]);
    end

    Pos3 = Posterior3;
    Pos3 = Pos3[:,setdiff(1:end, uninorms)];

    if length(uninorms) != model_def["nPar"]
        me2 = mean(Matrix(Pos3[:,:]),dims=1);
        sd2 = cov(Matrix(Pos3[:,:]));
    end

    stst = Vector{String}(undef,model_def["nPar"])
    count = 1;
    count2 = 1;
    for i in 1:model_def["nPar"]
        if i in lognorms
            stst[i] = string("theta[",i,"] = exp(((theta2[",count,"])*(",sqrt(sd[i,i]),"))+(",me[i],")); \n");
            count +=1;
        elseif i in uninorms
            par = fit(Uniform, priorsamps[:,i]);
            a1 = -2
            a2 = 2
            b1 = par.a
            b2 = par.b
            stst[i] = string("theta[",i,"] = ",b1,"+ ((theta_Un",count2," - (",a1,"))*(",b2," - ",b1,"))/(",a2," - (",a1,")); \n");
            count2 +=1;
        else
            stst[i] = string("theta[",i,"] = (((theta2[",count,"])*(",sqrt(sd[i,i]),"))+(",me[i],")); \n");
            count +=1;
        end
    end

    unst2 = Vector{String}(undef,length(uninorms))
    for i in 1:length(uninorms)
        par = fit(Uniform, priorsamps[:,uninorms[i]]);
        unst2[i] = string("theta_Un",i," ~ uniform(",-2,", ",2,"); \n");
    end

    sts2 = Array{String,1}(undef, length(stst)+1);
    sts2[1] = "real theta[nParms]; \n";
    sts2[2:end] = stst;

    unst3 = Array{String,1}(undef, length(unst2)+1);
    unst3[1] = "theta2 ~ multi_normal(mera,cora); \n";
    unst3[2:end] = unst2;

    priorfit["pars"] = unst;
    priorfit["transpars"] = sts2;
    priorfit["pridis"] = unst3;
    priorfit["numN"] = model_def["nPar"]-length(uninorms);
    if length(uninorms) != model_def["nPar"]
        priorfit["mera"] = vec(me2);
        priorfit["cora"] = sd2;
    else
        priorfit["mera"] = [];
        priorfit["cora"] = [];
    end
    return(priorfit)

end


## Function to check the contents of bayinf_def. Mostly the contents of the field priors plus call the functions checkStructBayInfData and checkStructBayInfStanSettings
function checkStructBayInf(model_def, bayinf_def)

    # Check taht all the dictionary entries are correct
    entries = ["Priors", "Data", "StanSettings", "flag", "plot", "runInf", "MultiNormFit"]
    if symdiff(entries,keys(bayinf_def))!=[] && symdiff(entries,keys(bayinf_def)) != ["ModelPath"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is soemthign wrong...")
        println(symdiff(entries,keys(bayinf_def)))
        return
    end

    # Check one to be sure there is no empty entry
    nee = ["Data"]; # No empty entries
    for i in 1:length(nee)
        if bayinf_def[nee[i]] == []
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Plese, check the field ", nee[i], ", this cannot be empty"))
            return
        end
    end

    if bayinf_def["flag"]!=[] && typeof(bayinf_def["flag"])!=String && typeof(bayinf_def["flag"])!=Array{String,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field flag!")
        return
    end

    if typeof(bayinf_def["flag"]) ==Array{String,1}
        bayinf_def["flag"] = bayinf_def["flag"][1];
    end

    if (typeof(bayinf_def["plot"]) != Bool) && (typeof(bayinf_def["plot"]) != Array{Bool,1}) &&
    (typeof(bayinf_def["plot"]) != String) && (typeof(bayinf_def["plot"]) != Array{String,1}) && bayinf_def["plot"] != []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field plot! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
        return
    end

    if (typeof(bayinf_def["runInf"]) != Bool) && (typeof(bayinf_def["runInf"]) != Array{Bool,1}) &&
    (typeof(bayinf_def["runInf"]) != String) && (typeof(bayinf_def["runInf"]) != Array{String,1}) && bayinf_def["runInf"] != []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field runInf! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
        return
    end

    if (typeof(bayinf_def["MultiNormFit"]) != Bool) && (typeof(bayinf_def["MultiNormFit"]) != Array{Bool,1}) &&
    (typeof(bayinf_def["MultiNormFit"]) != String) && (typeof(bayinf_def["MultiNormFit"]) != Array{String,1}) && bayinf_def["MultiNormFit"] != []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field MultiNormFit! This should be true, fals, \"Yes\", \"yes\", \"No\", \"no\" or an empty array ")
        return
    end

    if bayinf_def["plot"]==[true] || bayinf_def["plot"]==["Yes"] || bayinf_def["plot"]==["yes"] || bayinf_def["plot"]=="Yes" || bayinf_def["plot"]=="yes"
        bayinf_def["plot"]=true
    elseif bayinf_def["plot"]==[false] || bayinf_def["plot"]==["No"] || bayinf_def["plot"]==["no"] || bayinf_def["plot"]=="No" || bayinf_def["plot"]=="no" || bayinf_def["plot"]==[]
        bayinf_def["plot"]=false
    end

    if bayinf_def["runInf"]==[true] || bayinf_def["runInf"]==["Yes"] || bayinf_def["runInf"]==["yes"] || bayinf_def["runInf"]=="Yes" || bayinf_def["runInf"]=="yes" || bayinf_def["runInf"]==[]
        bayinf_def["runInf"]=true
    elseif bayinf_def["runInf"]==[false] || bayinf_def["runInf"]==["No"] || bayinf_def["runInf"]==["no"] || bayinf_def["runInf"]=="No" || bayinf_def["runInf"]=="no"
        bayinf_def["runInf"]=false
    end

    if bayinf_def["MultiNormFit"]==[true] || bayinf_def["MultiNormFit"]==["Yes"] || bayinf_def["MultiNormFit"]==["yes"] || bayinf_def["MultiNormFit"]=="Yes" || bayinf_def["MultiNormFit"]=="yes"
        bayinf_def["MultiNormFit"]=true
    elseif bayinf_def["MultiNormFit"]==[false] || bayinf_def["MultiNormFit"]==["No"] || bayinf_def["MultiNormFit"]==["no"] || bayinf_def["MultiNormFit"]=="No" || bayinf_def["MultiNormFit"]=="no" || bayinf_def["MultiNormFit"]==[]
        bayinf_def["MultiNormFit"]=false
    end

    # Extract contents of Priors field
    if typeof(bayinf_def["Priors"]) <: Dict # Object is a dictionary
        entries = ["pars", "transpars", "pridis"]
        if symdiff(entries,keys(bayinf_def["Priors"]))!=[] && symdiff(entries,keys(bayinf_def["Priors"])) != ["mera", "cora", "numN"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check the entries of the dictionary in Priors. This have to be pars, transpars and pridis.")
            println(symdiff(entries,keys(bayinf_def["Priors"])))
            return
        elseif bayinf_def["Priors"]["pars"] == [] && bayinf_def["Priors"]["pars"] != Array{String,1}
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check the entries of the dictionary in Priors. More presicely the field pars.")
            return
        elseif bayinf_def["Priors"]["pridis"] == [] && bayinf_def["Priors"]["pridis"] != Array{String,1}
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check the entries of the dictionary in Priors. More presicely the field pridis.")
            return
        else
            println("-------------------------- WARNING!!! --------------------------")
            println("Please, be carefull with the definition of parameters, transformations and prior distributions.")
            println("No check will be done in this section.")
        end
    elseif bayinf_def["Priors"] == [];
        nothing
    elseif typeof(bayinf_def["Priors"]) <: LinearAlgebra.Adjoint || typeof(bayinf_def["Priors"]) <: Array || typeof(bayinf_def["Priors"]) == String
        if typeof(bayinf_def["Priors"]) == Array{String,1}
            bayinf_def["Priors"] = bayinf_def["Priors"][1];
        end

        if typeof(bayinf_def["Priors"]) == String
            if bayinf_def["Priors"][end-3:end] == ".csv"
                if isfile(bayinf_def["Priors"])
                    priorsamps = Matrix(CSV.read(bayinf_def["Priors"]));
                    if bayinf_def["MultiNormFit"] == false
                        bayinf_def["Priors"] = fitPriorSamps(priorsamps, model_def);
                    elseif bayinf_def["MultiNormFit"] == true
                        bayinf_def["Priors"] = fitPriorSampsMultiNorm(priorsamps, model_def);
                    end
                else
                    println("-------------------------- Process STOPPED!!! --------------------------")
                    println("Sorry, but the CSV file introduced for the prior fit does not exist.")
                    return
                end
            elseif bayinf_def["Priors"][end-4:end] == ".stan"
                if isfile(bayinf_def["Priors"])
                    nothing
                else
                    println("-------------------------- Process STOPPED!!! --------------------------")
                    println("Sorry, but the STAN file introduced for the model does not exist.")
                    return
                end
            else
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Please, check the entries of the dictionary in Priors. This has to be the path to a .csv or a .stan")
                println("file. Also check that there is no spaces after the termination.")
                return
            end

        elseif size(bayinf_def["Priors"]) == (2,model_def["nPar"])
            thmax = maximum(bayinf_def["Priors"], dims = 1)
            thmin = minimum(bayinf_def["Priors"], dims = 1)
            mess = (thmax+thmin)/2;
            errs = (thmax-thmin)/4;
            priors = Dict();
            priors["pars"] = Array{String,1}(undef, model_def["nPar"]);
            priors["transpars"] = Array{String,1}(undef, model_def["nPar"]+1);
            priors["pridis"] = Array{String,1}(undef, model_def["nPar"]);

            priors["transpars"][1] = "real theta[nParms]; \n";

            for i in 1:model_def["nPar"]
                priors["pars"][i] = string("real<lower=-2,upper=2> ",model_def["parName"][i],"; \n")
                priors["transpars"][i+1] = string("theta[",i,"] = ((",model_def["parName"][i],")*(",errs[i],"))+",mess[i],"; \n");
                priors["pridis"][i] = string(model_def["parName"][i], " ~ normal(0,1); \n");
            end

            bayinf_def["Priors"] = priors
        elseif typeof(bayinf_def["Priors"]) == Array{Float64,2} || typeof(bayinf_def["Priors"]) == Array{Float32,2}
            if bayinf_def["MultiNormFit"] == false
                bayinf_def["Priors"] = fitPriorSamps(bayinf_def["Priors"], model_def);
            elseif bayinf_def["MultiNormFit"] == true
                bayinf_def["Priors"] = fitPriorSampsMultiNorm(bayinf_def["Priors"], model_def);
            end
        end
    end

    #     Check Data

    if haskey(bayinf_def["Data"], "switchT")
        bayinf_def["Data"] = checkStructBayInfData(model_def, bayinf_def["Data"]);
    else
        bayinf_def["Data"] = checkStructBayInfDataFiles(model_def, bayinf_def["Data"])
    end
    #     Check Stan Setting
    if bayinf_def["StanSettings"] != []
        bayinf_def["StanSettings"] = checkStructBayInfStanSettings(model_def, bayinf_def["StanSettings"]);
    end

    return(bayinf_def)
end

## Check the contents of data_def containing the data files or data info
function checkStructBayInfData(model_def, data_def)
    # Check taht all the dictionary entries are correct
    entries = ["Nexp", "finalTime", "switchT", "y0", "preInd", "uInd", "tsamps", "Obs", "DataMean", "DataError"]
    if symdiff(entries,keys(data_def))!=[] && symdiff(entries,keys(data_def)) != ["savepath", "savename"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is soemthign wrong...")
        println(symdiff(entries,keys(data_def)))
        return
    end

    # Checks to see if the content of the fields is the expected (also consier if single value entries have been given as a number or a vector)
    if ((typeof(data_def["Nexp"][1])!=Int) || length(data_def["Nexp"])!=1) && (typeof(data_def["Nexp"])!=Int)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Nexp! This should be an integer. ")
        return
    elseif (typeof(data_def["finalTime"]) != Array{Int,1}) && (typeof(data_def["finalTime"]) != Array{Float64,1}) && (typeof(data_def["finalTime"]) != Array{Int32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field finalTime! This should be a vector of final times. ")
        return
    elseif (typeof(data_def["switchT"]) != Array{Array{Int,1},1}) && (typeof(data_def["switchT"]) != Array{Array{Float64,1},1}) && (typeof(data_def["switchT"]) != Array{Array{Float32,1},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field switchT! This should be an array of arrays of numbers. ")
        return
    elseif (typeof(data_def["y0"]) != Array{Array{Int,1},1}) && (typeof(data_def["y0"]) != Array{Array{Float64,1},1}) && (typeof(data_def["y0"]) != Array{Array{Float32,1},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field y0! This should be an array of arrays of numbers. ")
        return
    elseif (typeof(data_def["preInd"]) != Array{Array{Int,1},1}) && (typeof(data_def["preInd"]) != Array{Array{Float64,1},1}) && (typeof(data_def["preInd"]) != Array{Array{Float32,1},1}) && data_def["preInd"] != []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field preInd! This should be an array of arrays of numbers. ")
        return
    # elseif (typeof(data_def["uInd"]) != Array{Array{Int,1},1}) && (typeof(data_def["uInd"]) != Array{Array{Float64,1},1}) && (typeof(data_def["uInd"]) != Array{Array{Float32,1},1})
    #     println("-------------------------- Process STOPPED!!! --------------------------")
    #     println("Please, check the field uInd! This should be an array of arrays of numbers. ")
    #     return
    elseif (typeof(data_def["tsamps"]) != Array{Array{Int,1},1}) && (typeof(data_def["tsamps"]) != Array{Array{Float64,1},1}) && (typeof(data_def["tsamps"]) != Array{Array{Float32,1},1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field tsamps! This should be an array of arrays of numbers. ")
        return
    elseif (typeof(data_def["Obs"]) != Array{String,1}) && (typeof(data_def["Obs"]) != Array{Int,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Obs! This should be a vector of strings or integers. ")
        return
    end

    # Extract necessary elements to ease generalisation
    if typeof(data_def["Nexp"]) == Array{Int,1}
        data_def["Nexp"] = data_def["Nexp"][1];
    end

    # Check that all the contents make sense
    if (data_def["Nexp"] != length(data_def["finalTime"])) || (data_def["Nexp"] != length(data_def["switchT"])) ||
        (data_def["Nexp"] != length(data_def["y0"])) || (data_def["Nexp"] != length(data_def["uInd"])) ||
        (data_def["Nexp"] != length(data_def["tsamps"]))
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, some field does not match the number of experiments selected. Check that the number of")
        println("entries for finalTime, switchT, y0, uInd or tsamps matches the number of experiments in Nexp.")
        return
    end

    for i in 1:data_def["Nexp"]
        if data_def["tsamps"][i][end] > data_def["finalTime"][i]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check finalTime. You have selected a sampling point past it.")
            return
        end

        if length(data_def["uInd"][i]) != (length(data_def["switchT"][i])-1)
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check uInd and switchT. Number of steps does not match the number of values for the inputs.")
            return
        end

        if length(data_def["y0"][i]) != model_def["nStat"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but Y0 does not match the number or states in experiment ", i))
            return
        end
    end

    if (model_def["Y0eqs"] != []) && data_def["preInd"]==[]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("You specified computation of Y0 as steady state but you have not introduced an inducer value for it!")
        return
    end

    if length(data_def["DataMean"]) != length(data_def["DataError"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, the dimensions the data means and errors do not match!")
        return
    end

    if length(data_def["DataMean"]) != data_def["Nexp"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, the dimensions the data and number of experiments indicated do not match!")
        return
    end

    if (typeof(data_def["DataMean"]) == Array{String,1} && typeof(data_def["DataError"]) != Array{String,1} ) ||
        (typeof(data_def["DataMean"]) != Array{String,1} && typeof(data_def["DataError"]) == Array{String,1} )
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, if files are introduced for the data please introduce the same in DataMean and DataError!")
        return
    end

    for i in 1:data_def["Nexp"]
        if size(data_def["DataMean"][i])[1] != length(data_def["tsamps"][i])
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, the dimensions the data and number of samples indicated do not match!")
            return
        end
    end

    for i in 1:data_def["Nexp"]

        if (size(data_def["DataMean"][i])[2] != length(data_def["Obs"]))
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but the data dimensions do not match the number of observables!")
            return
        end

        if (size(data_def["DataMean"][i])[1] != length(data_def["tsamps"][i]))
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but the number of data points do not match the number of sampling points!")
            return
        end

        for j in 1:length(data_def["Obs"])
            if size(data_def["DataError"][i][j])[1] != length(data_def["tsamps"][i])
                println("-------------------------- Process STOPPED!!! --------------------------")
                println("Sorry, but the number of error points do not match the number of sampling points!")
                return
            end

            if length(size(data_def["DataError"][i][j])) == 2
                if (size(data_def["DataError"][i][j]))[1] != (size(data_def["DataError"][i][j]))[2]
                    println("-------------------------- Process STOPPED!!! --------------------------")
                    println("Sorry, but it seems that you have introduced a covariance matrix for the data but this is not simetric!")
                    return
                end
            end
        end

    end

    # Check that no step is no smaller than 2 unit of time
    for i in 1:data_def["Nexp"]
        for j in 1:length(data_def["uInd"][i])
            if (data_def["switchT"][i][j+1]-data_def["switchT"][i][j])<=4
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
            for i in 1:data_def["Nexp"]
                if size(data_def["uInd"][i])[2] != model_def["nInp"];
                    data_def["uInd"][i] = data_def["uInd"][i]';
                end
            end
        end
    catch
        println("-------------------------- Process STOPPED!!! --------------------------")
        println(string("Sorry, but there is some issue with the contents of the field uInd."))
        return
    end

    if data_def["Obs"] == Array{Int,1}
        if (maximum(data_def["Obs"]) > model_def["nStat"]) || (length(data_def["Obs"]) > model_def["nStat"])
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Sorry, but there is some issue with the contents of the field Obs. It seems that "))
            println("     you have selected more observables than states or an index higher than the number of observables")
            return
        end
    end

    if typeof(data_def["Obs"]) == Array{String,1}
        if convert(Bool,sum(occursin.("+", data_def["Obs"]))) || convert(Bool,sum(occursin.("-", data_def["Obs"]))) ||
           convert(Bool,sum(occursin.("/", data_def["Obs"]))) ||
           convert(Bool,sum(occursin.("*", data_def["Obs"]))) || convert(Bool,sum(occursin.("^", data_def["Obs"])))
            nothing
        else
            if sum(occursin.(model_def["stName"], data_def["Obs"])) == 0
                println("-------------------------- Process STOPPED!!! --------------------------")
                println(string("Sorry, but there is some issue with the contents of the field Obs."))
                println("It seems that the observable(s) selected do not match any state")
                return
            else
                data_def["Obs"] = sort([findall(x->x==data_def["Obs"][i], model_def["stName"])[1] for i in 1:length(data_def["Obs"])])
            end
        end

    elseif typeof(data_def["Obs"]) == Array{Int,1}
        data_def["Obs"] = sort(data_def["Obs"]);
    end


    return(data_def)
end


## Check the contents of data_def containing the data files or data info
function checkStructBayInfDataFiles(model_def, data_def)
    entries = ["Obs", "Observables", "Inputs", "y0"]
    if symdiff(entries,keys(data_def))!=[]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is soemthign wrong...")
        println(symdiff(entries,keys(data_def)))
        return
    end

    # Check one to be sure there is no empty entry
    nee = ["Obs", "Observables", "Inputs", "y0"]; # No empty entries
    for i in 1:length(nee)
        if data_def[nee[i]] == []
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Plese, check the field ", nee[i], ", this cannot be empty"))
            return
        end
    end

    if (typeof(data_def["Obs"]) != Array{String,1}) && (typeof(data_def["Obs"]) != Array{Int,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Obs! This should be a vector of strings or integers. ")
        return
    end

    if typeof(data_def["Observables"]) != Array{String,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, there is something wrong in the type of the field Observables. This should be a Array{String,1}")
        return
    end

    if typeof(data_def["Inputs"]) != Array{String,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, there is something wrong in the type of the field Inputs. This should be a Array{String,1}")
        return
    end

    if length(data_def["Observables"]) != length(data_def["Inputs"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, there number of entries for ObservablesFile and EventInputsFile should be the same")
        return
    end

    for i in 1:length(data_def["Observables"])
        if data_def["Observables"][i][end-3:end] != ".csv"
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, add the .csv termination to the ObservablesFile file! And no spaces after!")
            return
        end
    end

    for i in 1:length(data_def["Inputs"])
        if data_def["Inputs"][i][end-3:end] != ".csv"
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, add the .csv termination to the ObservablesFile file! And no spaces after!")
            return
        end
    end

    datatu = defBayInfDataStruct();
    datatu["Obs"] = data_def["Obs"];
    datatu["y0"] = data_def["y0"];

    tmpInp = Array{Any}(undef, 1,length(data_def["Inputs"]))
    for i in 1:length(data_def["Inputs"])
        tmpInp[i] = Matrix(CSV.read(data_def["Inputs"][i]));
    end

    tmpObs = Array{Any}(undef, 1,length(data_def["Observables"]))
    for i in 1:length(data_def["Observables"])
        tmpObs[i] = Matrix(CSV.read(data_def["Observables"][i]));
    end

    try
        datatu["Nexp"] = length(data_def["Observables"]);
        datatu["finalTime"] = [tmpInp[i][1,2] for i in 1:length(data_def["Inputs"])];
        datatu["switchT"] = [vcat(tmpInp[i][:,1], tmpInp[i][1,2]) for i in 1:length(data_def["Inputs"])];
        datatu["preInd"] = [tmpInp[i][1,3:(2+model_def["nInp"])] for i in 1:length(data_def["Inputs"])];
        datatu["uInd"] = [tmpInp[i][:,(3+model_def["nInp"]):(3+model_def["nInp"]+(model_def["nInp"]-1))] for i in 1:length(data_def["Inputs"])];
        datatu["tsamps"] = [tmpObs[i][:,1] for i in 1:length(data_def["Inputs"])];

        tmp1 = Array{Any,1}(undef,datatu["Nexp"])
        tmp2 = Array{Any,1}(undef,datatu["Nexp"])
        for i in 1:datatu["Nexp"]

            tmpD = Matrix(CSV.read(data_def["Observables"][i]));
            if size(tmpD)[2] != 3*length(data_def["Obs"])

                tmpdat = zeros(size(tmpD)[1],convert(Int,size(tmpD)[2]/(3*length(data_def["Obs"]))), length(data_def["Obs"]));
                co = 2;
                for k in 1:convert(Int,size(tmpD)[2]/(3*length(data_def["Obs"])))
                    for j in 1:length(data_def["Obs"])
                        tmpdat[:,k] = tmpD[:,co,j];
                        co +=3;
                    end
                end

                tmco = Array{Any}(undef,length(data_def["Obs"]))
                for j in 1:length(data_def["Obs"])
                    tmco[j] = cov(tmpdat[:,:,j]');
                end

                tmp1[i] = mean(tmpdat, dims=2);
                tmp2[i] = tmco;

            else
                nobs = length(datatu["Obs"])

                mees = zeros(size(tmpObs[i])[1], nobs);
                errs =  Array{Any,1}(undef,nobs);#zeros(size(tmpObs[i])[1], nobs);
                r = 1:3:nobs*3
                for j in 1:nobs
                    mees[:,j] = tmpObs[i][:,r[j]+1];
                    errs[j] = tmpObs[i][:,r[j]+2];
                end
                tmp1[i] = mees;
                tmp2[i] = errs;
            end
        end

        datatu["DataMean"] = tmp1; # Array containin ght evector of means for each experiment.
        datatu["DataError"] = tmp2;

    catch
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, but there seems to be some issue with the internal structure of the files you have given...")
        return
    end

    datatu = checkStructBayInfData(model_def, datatu);

    return(datatu)
end


## Check the contents of stan_def containing the settings for stan.
function checkStructBayInfStanSettings(model_def, stan_def)

    # Check taht all the dictionary entries are correct
    entries = ["cmdstan_home", "nchains", "nsamples", "nwarmup", "printsummary", "init", "maxdepth", "adaptdelta", "jitter"]
    if symdiff(entries,keys(stan_def))!=[] && symdiff(entries,keys(stan_def)) != ["savepath", "savename"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is soemthign wrong...")
        println(symdiff(entries,keys(stan_def)))
        return
    end

    # Check one to be sure there is no empty entry
    nee = ["cmdstan_home", "nchains", "nsamples", "nwarmup"]; # No empty entries
    for i in 1:length(nee)
        if stan_def[nee[i]] == []
            println("-------------------------- Process STOPPED!!! --------------------------")
            println(string("Plese, check the field ", nee[i], ", this cannot be empty"))
            return
        end
    end

    # Checks to see if the content of the fields is the expected (also consier if single value entries have been given as a number or a vector)
    if typeof(stan_def["cmdstan_home"]) != String && typeof(stan_def["cmdstan_home"]) != Array{String,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field cmdstan_home! This should be an string. ")
        return
    elseif typeof(stan_def["nchains"]) != Int && typeof(stan_def["nchains"]) != Array{Int,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field nchains! This should be an string. ")
        return
    elseif typeof(stan_def["nsamples"]) != Int && typeof(stan_def["nsamples"]) != Array{Int,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field nsamples! This should be an string. ")
        return
    elseif typeof(stan_def["nwarmup"]) != Int && typeof(stan_def["nwarmup"]) != Array{Int,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field nwarmup! This should be an string. ")
        return
    elseif stan_def["printsummary"] != true && stan_def["printsummary"] != false && stan_def["printsummary"] != [true] &&
        stan_def["printsummary"] != [false] && stan_def["printsummary"] != []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field printsummary! This should be an string. ")
        return
    elseif stan_def["init"] != [] && !(typeof(stan_def["init"]) <: Array)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field init! This should be an array of dictionaries or an empty array. ")
        return
    elseif !isempty(stan_def["init"])
        if (typeof(stan_def["init"]) <: Array) && !(typeof(stan_def["init"][1]) <: Dict)
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check the field init! This should be an array of dictionaries or an empty array. ")
            return
        end
    elseif typeof(stan_def["maxdepth"]) != Int && typeof(stan_def["maxdepth"]) != Array{Int,1} && stan_def["maxdepth"] != [] &&
    typeof(stan_def["maxdepth"]) != Float64 && typeof(stan_def["maxdepth"]) != Array{Float64,1} &&
    typeof(stan_def["maxdepth"]) != Float32 && typeof(stan_def["maxdepth"]) != Array{Float32,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field maxdepth! This should be an string. ")
        return
    elseif typeof(stan_def["adaptdelta"]) != Int && typeof(stan_def["adaptdelta"]) != Array{Int,1} && stan_def["adaptdelta"] != [] &&
    typeof(stan_def["adaptdelta"]) != Float64 && typeof(stan_def["adaptdelta"]) != Array{Float64,1} &&
    typeof(stan_def["adaptdelta"]) != Float32 && typeof(stan_def["adaptdelta"]) != Array{Float32,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field adaptdelta! This should be an string. ")
        return
    elseif typeof(stan_def["jitter"]) != Int && typeof(stan_def["jitter"]) != Array{Int,1} && stan_def["jitter"] != [] &&
    typeof(stan_def["jitter"]) != Float64 && typeof(stan_def["jitter"]) != Array{Float64,1} &&
    typeof(stan_def["jitter"]) != Float32 && typeof(stan_def["jitter"]) != Array{Float32,1}
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field jitter! This should be an string. ")
        return

    end

    # Extract necessary elements to ease generalisation

    if typeof(stan_def["cmdstan_home"]) == Array{String,1}
        stan_def["cmdstan_home"] = stan_def["cmdstan_home"][1];
    end

    if typeof(stan_def["nchains"]) == Array{Int,1}
        stan_def["nchains"] = stan_def["nchains"][1];
    end

    if typeof(stan_def["nsamples"]) == Array{Int,1}
        stan_def["nsamples"] = stan_def["nsamples"][1];
    end

    if typeof(stan_def["nwarmup"]) == Array{Int,1}
        stan_def["nwarmup"] = stan_def["nwarmup"][1];
    end

    if typeof(stan_def["printsummary"]) == Array{Bool,1}
        stan_def["printsummary"] = stan_def["printsummary"][1];
    end

    if typeof(stan_def["maxdepth"]) == Array{Int,1} || typeof(stan_def["maxdepth"]) == Array{Float64,1} || typeof(stan_def["maxdepth"]) == Array{Float32,1}
        stan_def["maxdepth"] = stan_def["maxdepth"][1];
    end

    if typeof(stan_def["adaptdelta"]) == Array{Int,1} || typeof(stan_def["adaptdelta"]) == Array{Float64,1} || typeof(stan_def["adaptdelta"]) == Array{Float32,1}
        stan_def["adaptdelta"] = stan_def["adaptdelta"][1];
    end

    if typeof(stan_def["jitter"]) == Array{Int,1} || typeof(stan_def["jitter"]) == Array{Float64,1} || typeof(stan_def["jitter"]) == Array{Float32,1}
        stan_def["jitter"] = stan_def["jitter"][1];
    end

    if stan_def["printsummary"] == []
        stan_def["printsummary"] = true;
    end

    return(stan_def)
end


## Generate Dictionary for stan initial samples
function genStanInitDict(samps, names, chains)
    alltog = Array{Dict{String,Any},1}(undef,chains);
    for i in 1:(chains)
        global tmp = Dict{String, Any}()
        for j in 1:length(names)
            tmp[names[j]] = samps[j,i];
        end
        alltog[i] = tmp
    end

    return(alltog)
end

## Transform dictionary entries to correct transformed range defined in transformed parameters
function reparamDictStan(standict, bayinf_def)

    fuus = Array{String,2}(undef, length(standict), length(bayinf_def["Priors"]["transpars"]));
    standict2 = Array{Dict{String,Any},1}(undef, length(standict));

    for dis in 1:length(standict)
        standict2[dis] = Dict();
        for pa in 1:length(bayinf_def["Priors"]["transpars"])
            fuus[dis,pa] = "";
            for k in keys(standict[dis])
                if occursin(k, bayinf_def["Priors"]["transpars"][pa])
                    t1 = replace(bayinf_def["Priors"]["transpars"][pa], bayinf_def["Priors"]["transpars"][pa][1:findfirst("=", bayinf_def["Priors"]["transpars"][pa])[1]] => "")
                    t1 = replace(t1, "\n" => "")
                    t1 = replace(t1, ";" => "")
                    t1 = replace(t1, k => "x")
                    tt = (string("f",dis,pa,"(x) = ", string("(", t1, ")-", standict[dis][k]), "; val = find_zero(f",dis,pa,", 0);"))
                    fuus[dis,pa] = tt;
                    s1 = Meta.parse(fuus[dis,pa]);
                    standict2[dis][k] = @eval $s1;
                end
            end
        end
    end
    return(standict2)

end



## Generation of the stan model
# Function that generates the Stan Model
function genStanModel(model_def, bayinf_def)

    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\ModelsFunctions"))
        mkdir(string(cudi,"\\ModelsFunctions"))
    end


    if typeof(bayinf_def["Priors"]) == String
        bayinf_def["ModelPath"] = bayinf_def["Priors"];
    else
        bayinf_def["ModelPath"] = string(cudi, "\\ModelsFunctions\\", model_def["NameF"], "_StanModel.stan");

        inps = Array{String,1}(undef, model_def["nInp"]);
        [inps[i] = string("      real ", model_def["inpName"][i], " = x_r[",i,"]; \n") for i in 1:model_def["nInp"]];

        parsss = Array{String,1}(undef, model_def["nPar"]);
        [parsss[i] = string("      real ", model_def["parName"][i], " = p[",i,"]; \n") for i in 1:model_def["nPar"]];

        stats = Array{String,1}(undef, model_def["nStat"]);
        [stats[i] = string("      real ", model_def["stName"][i], " = y[",i,"]; \n") for i in 1:model_def["nStat"]];

        eqs = Array{String,1}(undef, length(model_def["eqns"]));
        for i in 1:length(model_def["eqns"])
            for j in 1:model_def["nStat"]
                if occursin(string("d", model_def["stName"][j]), model_def["eqns"][i][1:findfirst("=", model_def["eqns"][i])[1]])
                    eqs[i] = string("      dInd_dt[",j,"] = ", model_def["eqns"][i][(findfirst("=", model_def["eqns"][i])[1]+1):end], "; \n");
                end
            end
            if !isassigned(eqs, i)
                eqs[i] = string("      ", model_def["eqns"][i], "; \n")
            end
            for j in 1:model_def["nStat"]
                if occursin(string("d", model_def["stName"][j]), eqs[i][(findfirst("=", eqs[i])[1]+1):end])
                    sta = (string("d", model_def["stName"][j]));
                    tmp1 =  replace((eqs[i][(findfirst("=", eqs[i])[1]+1):end]), sta => string("dInd_dt[",j,"]"))
                    eqs[i] = replace(eqs[i], (eqs[i][(findfirst("=", eqs[i])[1]+1):end]) => tmp1)

                end
            end
        end

        if "alp" in model_def["parName"] || "alp" in model_def["inpName"] || "alp" in model_def["stName"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Sorry, but the name alp is reserved for something else in Stan...")
            return
        end

        parinp1 = Array{String,1}(undef, model_def["nPar"]);
        for i in 1:model_def["nPar"]
            parinp1[i] = string("      real ", model_def["parName"][i], " = p[",i,"]; \n")
        end

        parinp2 = Array{String,1}(undef, model_def["nInp"]);
        for i in 1:model_def["nInp"]
            parinp2[i] = string("      real ", model_def["inpName"][i], " = x_r[",i,"]; \n")
        end

        # ssstates = Array{String,1}(undef, model_def["nStat"]);
        # [ssstates[i] = string("      real ", model_def["stName"][i], " = alp[",i,"]; \n") for i in 1: model_def["nStat"]];

        sseqs = Array{String,1}(undef, length(model_def["Y0eqs"]));
        for i in 1:length(model_def["Y0eqs"])
            for j in 1:model_def["nStat"]
                if !isnothing(findfirst("=", model_def["Y0eqs"][i]))
                    if occursin(string(model_def["stName"][j]), model_def["Y0eqs"][i][1:findfirst("=", model_def["Y0eqs"][i])[1]])
                        sseqs[i] = string("      alp[",j,"] = ", model_def["Y0eqs"][i][(findfirst("=", model_def["Y0eqs"][i])[1]+1):end], "; \n");
                    end
                end
            end
            if !isassigned(sseqs, i)
                sseqs[i] = string("      ", model_def["Y0eqs"][i], "; \n")
            end
        end

        for i in 1:length(sseqs)
            for j in 1:model_def["nStat"]
                sseqs[i] = replace(sseqs[i], string("exp",model_def["stName"][j])=>string("init[",j,"]"));
                sseqs[i] = replace(sseqs[i], string(model_def["stName"][j])=>string("alp[",j,"]"));
            end
        end

        # ssinps = Array{String,1}(undef, (model_def["nStat"]));
        # [ssinps[i] = string("      real exp", model_def["stName"][i], " = init[",i,"]; \n") for i in 1:model_def["nStat"]];

        # Steady State Equations
        sseqsdef = "";
        if model_def["Y0eqs"] == []
            sseqsdef = "      return init; \n";
        else
            sseqsdef = string("
      vector[",model_def["nStat"],"] alp;

      // Parameters and Inputs
",join(parinp1),"
",join(parinp2),"

      // Equations

",join(sseqs),"

      // Results
      return alp;
                ");
        end;

        functs = string("

functions{

    real[] ",model_def["NameF"],"_ODEs(real t, real[] y, real[] p, real[] x_r, int[] x_i){

      // Inputs
",join(inps),"

      // Parameters
",join(parsss),"

      // ODEs
      real dInd_dt[",model_def["nStat"],"];

",join(stats),"

",join(eqs),"

      // Results
      return dInd_dt;

    }


    vector SteadyState(vector init, vector p, real[] x_r, int[] x_i){

",join(sseqsdef),"

    }

}
    ");

        if model_def["Y0ON"] == true
            onvec = string("
    // 24h incubation times for steady state calculation
    int tonil; //-------------------------------------------------------------------> Introduce a check in case Y0 needs steady state simulation
    real toni[tonil];
                ")
        else
            onvec = "";
        end

        if length(size(bayinf_def["Data"]["DataError"][1][1])) == 1
            erros = "      real Erros[stslm,m,obser];";
        else
    #             erros = Array{Any,1}(undef, bayinf_def["Data"]["Nexp"])
    #                 st = length(bayinf_def["Data"]["tsamps"][j]);
            erros = string("    real ErrosMul[stslm,stslm,m,obser]; \n")
        end

        if bayinf_def["MultiNormFit"]
            mnp = string("
    // Prior parameters
    int numN;
    vector[numN] mera; // Mean
    matrix[numN,numN] cora; // Covariance matrix
    ");
        else
            mnp="";
        end

        dats = string("

data {

    // Observables
    int m; // Total number of data series
    int stslm; // Maximum number of rows for all the observable matrices
    int stsl[1,m]; // Number of elements at each time series for each series m
    int sts[stslm,m]; // Sampling times for each series m
    int obser;//-> Introduce this so we have all the data in one same array (easier generalisation). Work on generalisation in case different experiments have different obsevables?
    int obSta[1,obser]; // -> This variable will be to know which are the observable states

    real Means[stslm,m,obser]; // ---> General arrays of means and errors
",join(erros),"

    // Inputs
    int elm; // Number of rows in the matrices for IPTG and aTc, half for the inputs and -1 for the total number of events
    int tml; // Maximum length of the rows for the sampling times matrix
    int Nsp[1,m]; // Number of event switching points (including final time) times for each series m
    real ts[tml, m]; // Time series for each serie m
    int tsl[1,m]; // Length of sampling time series per serie m

    int nindu; // -> Number of inducers/stimuly
    real preInd[nindu,m]; // Values of inputs for each serie m for the ON incubation

    real inputs[(elm*nindu),m]; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...

    int evnT[(elm+1),m]; // Event change time points for each serie m
    real Y0us[",model_def["nStat"],",m]; // Y0 vectors

",join(onvec),"
",join(mnp),"
}


    ");

        transdat = string("

transformed data {

    int nParms = ",model_def["nPar"],"; // Number of parameters of the model //--------> Introduce number in generation of script
    int Neq = ",model_def["nStat"],"; // Total number of equations of the model //-----> Introduce number in generation of script
    int x_i[0]; // Empty x_i object (needs to be defined)
    real x_r[(elm*nindu),m]=inputs; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
    real ivss[Neq,m] = Y0us; // Initial experimental values for the calculation of the steady state ordered as LacI+RFP, TetR+GFP -------------------------> Carefull to how I define this
    real pre[nindu,m]; // Input values during the 24h incubation ordered as IPTG, aTc

    for(i in 1:m){
        for(k in 1:nindu){
          pre[k,i] = preInd[k,i]; //-> If people only give a Y0 value, this will be NO needed
        };
    };

}
    ");


        if typeof(bayinf_def["Priors"]) <: Dict
            for i in 1:length(bayinf_def["Priors"]["transpars"])
                if bayinf_def["Priors"]["transpars"][i][end] !='\n'
                    bayinf_def["Priors"]["transpars"][i] = string(bayinf_def["Priors"]["transpars"][1], "\n");
                end
            end

            for i in 1:length(bayinf_def["Priors"]["pars"])
                if bayinf_def["Priors"]["pars"][i][end] !='\n'
                    bayinf_def["Priors"]["pars"][i] = string(bayinf_def["Priors"]["pars"][1], "\n");
                end
            end

            for i in 1:length(bayinf_def["Priors"]["pridis"])
                if bayinf_def["Priors"]["pridis"][i][end] !='\n'
                    bayinf_def["Priors"]["pridis"][i] = string(bayinf_def["Priors"]["pridis"][1], "\n");
                end
            end
        elseif bayinf_def["Priors"] == [];
            bayinf_def["Priors"]["transpars"] = "";
            bayinf_def["Priors"]["pars"] = "";
            bayinf_def["Priors"]["pridis"] = "";
        end

        pars = string("

parameters {

",join(bayinf_def["Priors"]["pars"]),"

}

    ");

        transpars = string("

transformed parameters {

",join(bayinf_def["Priors"]["transpars"]),"

}

    ");

        if model_def["Y0ON"] == true
            ssv = "real ssv[tonil,Neq]; // Real that will include the solution of the ODE for the 24h incubation ";
            ing = string("
            // 24h incubation calculation for the steady state
            ssv = integrate_ode_bdf(",model_def["NameF"],"_ODEs, Y0[,j],0,toni,theta,pre[,j], x_i, ",model_def["tols"][1],", ",model_def["tols"][2],", 1e7);
            Y0[,j] = ssv[tonil];
                ")
        else
            ssv = "";
            ing = "";
        end;


        if length(size(bayinf_def["Data"]["DataError"][1][1])) == 1
            if typeof(bayinf_def["Data"]["Obs"]) == Array{Int,1}
                lik = string("
            for(ob in 1:obser){
              yhat[t,j,ob] = y_hat[(sts[t,j]+1),obSta[ob]]; //---> Will need to double check if this is enough or there will be soem confusion on selecting the observable states and aching it with data
              Means[t,j,ob] ~ normal(yhat[t,j,ob],Erros[t,j,ob]); //---> Normal as default, but should be able to include other optios. Also have a check in case the user whants to define different ones for different observables
            }

                    ");
                meees = "";
                errrs = "";
            else
                meees = "";
                errrs = "";
                lik1 = Array{String,1}(undef,length(bayinf_def["Data"]["Obs"]));
                lik2 = Array{String,1}(undef,length(bayinf_def["Data"]["Obs"]));
                for i in 1:length(bayinf_def["Data"]["Obs"])
                    bayinf_def["Data"]["Obs"][i] = replace(bayinf_def["Data"]["Obs"][i], ".+"=>"+")
                    bayinf_def["Data"]["Obs"][i] = replace(bayinf_def["Data"]["Obs"][i], ".-"=>"-")
                    bayinf_def["Data"]["Obs"][i] = replace(bayinf_def["Data"]["Obs"][i], ".*"=>"*")
                    bayinf_def["Data"]["Obs"][i] = replace(bayinf_def["Data"]["Obs"][i], "./"=>"/")
                    bayinf_def["Data"]["Obs"][i] = replace(bayinf_def["Data"]["Obs"][i], ".^"=>"^")
                    tmp2 = string("        yhat[t,j,",i,"] = ", bayinf_def["Data"]["Obs"][i], "; \n");
                    for j in 1:model_def["nStat"]
                        if occursin(model_def["stName"][j], tmp2)
                            tmp2 = replace(tmp2, model_def["stName"][j]=> string("y_hat[(sts[t,j]+1),",j,"]"))
                        end
                    end
                    lik1[i] = tmp2;

                    lik2[i] = string("        Means[t,j,",i,"] ~ normal(yhat[t,j,",i,"],Erros[t,j,",i,"]);")
                    lik = join([join(lik1),join(lik2)]);
                end

            end
        else
            if typeof(bayinf_def["Data"]["Obs"]) == Array{Int,1}

                meees = Array{Any,1}(undef, length(bayinf_def["Data"]["Obs"]));
                errrs = Array{Any,1}(undef, length(bayinf_def["Data"]["Obs"]));
                for j in 1:length(bayinf_def["Data"]["Obs"])
                    meees[j] = string("  to_vector(Means[1:stsl[1,j],j,",j,"]) ~ multi_normal(to_vector(yhat[1:stsl[1,j],j,ob]), (Erros_",j,")); \n")
                    errrs[j] = string("  matrix[stsl[1,j],stsl[1,j]] Erros_",j," = to_matrix(ErrosMul[1:stsl[1,j],1:stsl[1,j],j,",j,"]); \n");
                end

                lik = string("
            for(ob in 1:obser){
              yhat[t,j,ob] = y_hat[(sts[t,j]+1),obSta[ob]]; //---> Will need to double check if this is enough or there will be soem confusion on selecting the observable states and aching it with data
            }
                    ");
            else
                errrs = Array{String,1}(undef,length(bayinf_def["Data"]["Obs"]));
                lik = Array{String,1}(undef,length(bayinf_def["Data"]["Obs"]));
                meees = Array{String,1}(undef,length(bayinf_def["Data"]["Obs"]));
                for i in 1:length(bayinf_def["Data"]["Obs"])
                    errrs[i] = string("   matrix[stsl[1,j],stsl[1,j]] Erros_",i," = to_matrix(ErrosMul[1:stsl[1,j],1:stsl[1,j],j,",i,"]); \n");

                    bayinf_def["Data"]["Obs"][i] = replace(bayinf_def["Data"]["Obs"][i], "+."=>"+")
                    bayinf_def["Data"]["Obs"][i] = replace(bayinf_def["Data"]["Obs"][i], "-."=>"-")
                    bayinf_def["Data"]["Obs"][i] = replace(bayinf_def["Data"]["Obs"][i], "*."=>"*")
                    bayinf_def["Data"]["Obs"][i] = replace(bayinf_def["Data"]["Obs"][i], "/."=>"/")
                    bayinf_def["Data"]["Obs"][i] = replace(bayinf_def["Data"]["Obs"][i], "^."=>"^")
                    tmp2 = string("        yhat[t,j,",i,"] = ", bayinf_def["Data"]["Obs"][i], "; \n");
                    for j in 1:model_def["nStat"]
                        if occursin(model_def["stName"][j], tmp2)
                            tmp2 = replace(tmp2, model_def["stName"][j]=> string("y_hat[(sts[t,j]+1),",j,"]"))
                        end
                    end
                    lik[i] = tmp2;

                    meees[i] = string("  to_vector(Means[1:stsl[1,j],j,",i,"]) ~ multi_normal(to_vector(yhat[1:stsl[1,j],j,",i,"]),(Erros_",i,"));")
                end


            end
        end

        mode = string("

model {

  // Intermediate parameters
  int i; // Increasing index for the inputs
  vector[Neq] ing; // Vector that will include the solution of the algebraic solution for the steady state of the model
  ",join(ssv),"
  real Y0[Neq,m]; // Initial values for the ODEs variables at the first event
  real yhat[stslm,m,obser]; // ---> Generall array to include all the observables (easier generalisation)

  // Reparameterised priors definition
",join(bayinf_def["Priors"]["pridis"]),"

  // Likelihood
  for (j in 1:m){
",join(errrs),"
    real ivst[Neq]; // Initial value of the states
    real y_hat[(tsl[1,j]),Neq]; // Object to include the ODEs solutions for each state

    // Calculation of Y0

    ing = SteadyState(to_vector(ivss[,j]), to_vector(theta), pre[,j], x_i);
    for(g in 1:Neq){
      Y0[g,j] = ing[g];
    };
",join(ing),"

    i = 1;

    // Loop (over the number of events/inputs) to solve the ODE for each event stopping the solver and add them to the final object y_hat

    for (q in 1:Nsp[1,j]-1){
      int itp = evnT[q,j];  // Initial time points of each event
      int lts = num_elements(ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j]);  // Length of the time series for each event
      real part1[lts,Neq]; // Temporary object that will include the solution of the ODE for each event at each loop

      // Calculation of the solution for the ODEs where for events that are not the first one. The time series starts one minute before the original point of the time serie overlaping with the last point of the previous event with same state values at the time
      if (q == 1){
        ivst = Y0[,j];
        part1 = integrate_ode_bdf(",model_def["NameF"],"_ODEs,ivst,itp,ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j],theta,to_array_1d(inputs[i:(i+(nindu-1)),j]), x_i, ",model_def["tols"][1],", ",model_def["tols"][2],", 1e7);
      }
      else{
        part1 = integrate_ode_bdf(",model_def["NameF"],"_ODEs, ivst,(itp-1e-7),ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j],theta,to_array_1d(inputs[i:(i+(nindu-1)),j]), x_i, ",model_def["tols"][1],", ",model_def["tols"][2],", 1e7);
      }

      // Modification of the initial state values for the next event
      ivst = part1[lts];
      // Increase index for inputs
      i=i+nindu;

      // Introduction of the result of part1 into the object y_hat
      for (y in (itp+1):(itp+lts)){
        y_hat[(y),]=(part1)[(y-itp),];
      };
    };

    // Likelihood definition (residuals) at each sampling time

    for (t in 1:stsl[1,j]){ //----> General form
",join(lik),"
    }
",join(meees),"
  }
}

    ");

        stanmodel = join([functs, dats, transdat, pars, transpars, mode]);
        open(string(cudi, "\\ModelsFunctions\\", model_def["NameF"], "_StanModel.stan"), "w") do io
           write(io, stanmodel);
        end;

        println("")
        println("----------------------------------------- STAN MODEL GENERATION -----------------------------------------")
        println("The model has been generated in the directory: ")
        println(string("                 ", bayinf_def["ModelPath"]))
        println("--------------------------------------------------------------------------------------")
        println("")

    end


    return(bayinf_def)
end


## Function to get the data info and re-structure it in a dictionary for Stan.
function restructureDataInference(model_def, bayinf_def)

    ######################## OBSERVABLES ########################
    mdp = zeros(bayinf_def["Data"]["Nexp"]); # Maximum data point for each file (multiexperimental case, which should be a vector as long as number of experiments)
    for i in collect(1:bayinf_def["Data"]["Nexp"]) # Loop that opens a file for each iteration to extract the file with maximum data points
        mdp[i] = length(bayinf_def["Data"]["tsamps"][i]);
    end
    mdp = convert.(Int, mdp);

    mrow = convert(Int,maximum(mdp)); # Maximum number of rows
    mcol = bayinf_def["Data"]["Nexp"]; # Maximum number of columns

    # Definition of the matrices with 0s on them
    samplingT = zeros(mrow,mcol);
    for i in collect(1:bayinf_def["Data"]["Nexp"]) # Loop that opens a file for each iteration to extract the data
        tss = bayinf_def["Data"]["tsamps"][i];
        for j in 1:length(tss[:,1])
          samplingT[j,i] = round(tss[j,1]);

        end
    end

    obser = length(bayinf_def["Data"]["Obs"]);
    MEANS = zeros(mrow,mcol,obser);
    if length(size(bayinf_def["Data"]["DataError"][1][1])) == 1
        STDS = zeros(mrow,mcol,obser);
        for i in 1:mcol
            for j in 1:obser
                MEANS[1:length(bayinf_def["Data"]["tsamps"][i]),i,j] = bayinf_def["Data"]["DataMean"][i][:,j];
                STDS[1:length(bayinf_def["Data"]["tsamps"][i]),i,j] = bayinf_def["Data"]["DataError"][i][j];
            end
        end
    else
        STDS = zeros(mrow,mrow,mcol,obser);
        for i in 1:mcol
            for j in 1:obser
                MEANS[1:length(bayinf_def["Data"]["tsamps"][i]),i,j] = bayinf_def["Data"]["DataMean"][i][:,j];
                corrmat = Diagonal(ones(length(bayinf_def["Data"]["tsamps"][i]))*0.1); # Necessary to ensure the matrix be positive definite
                STDS[1:length(bayinf_def["Data"]["tsamps"][i]),1:length(bayinf_def["Data"]["tsamps"][i]),i,j] = bayinf_def["Data"]["DataError"][i][j].+corrmat;
            end
        end
    end




    mdp2 = zeros(1,mcol);
    mdp2[1,:] = mdp;

    ######################## INPUTS ########################
    mdpI  = zeros(bayinf_def["Data"]["Nexp"]); # Maximum data point for each file
    mtimes  = zeros(bayinf_def["Data"]["Nexp"]); # Maximum time point for each file

    for i in collect(1:bayinf_def["Data"]["Nexp"]) # Loop that opens a file for each iteration to extract the file with maximum data points
        mdpI[i] = length(bayinf_def["Data"]["uInd"][i]);
        mtimes[i] = convert(Int,round(bayinf_def["Data"]["finalTime"][i]));
    end

    mrowI = convert(Int,maximum(mdpI)); # Maximum number of rows
    mcolI = bayinf_def["Data"]["Nexp"]; # Maximum number of columns
    mt = convert(Int,maximum(mtimes)); # Maximum number of time points

    mdpI2 = zeros(1,mcol);
    mdpI2[1,:] = mdpI;

    evnT = zeros(mrowI+1,mcolI); # Event chance time points
    preI = zeros(model_def["nInp"],mcolI); # Inducer initial values
    time = zeros(mt+1,mcolI); # Time points vector
    inps = zeros(mrowI*model_def["nInp"],mcolI); # Inputs vector

    ltimes = zeros(bayinf_def["Data"]["Nexp"]); # Length of time series

    for i in 1:bayinf_def["Data"]["Nexp"] # Loop that opens a file for each iteration to extract the data

        tempT = collect(1e-9:round(bayinf_def["Data"]["finalTime"][i])+1)';
        ltimes[i] = length(tempT);
        evnT[1:length(bayinf_def["Data"]["switchT"][i]),i] = bayinf_def["Data"]["switchT"][i];

        r = convert.(Int, 1:model_def["nInp"]:(length(bayinf_def["Data"]["uInd"][i])));
        inputs = zeros(convert.(Int,length(bayinf_def["Data"]["uInd"][i])));
        for j in 1:convert.(Int,length(bayinf_def["Data"]["uInd"][i])/model_def["nInp"])
            for k in 0:(model_def["nInp"]-1)
                inputs[r[j]+k] = bayinf_def["Data"]["uInd"][i][j,(k+1)];
            end
        end
        inps[1:length(inputs),i] = inputs;

        for l in 1:length(tempT)
            time[l,i] = tempT[l];
        end

        try
            if bayinf_def["Data"]["preInd"] == [];
                 bayinf_def["Data"]["preInd"] = zeros(model_def["nInp"]);
            end
        catch
            if bayinf_def["Data"]["preInd"][i] == [];
                bayinf_def["Data"]["preInd"][i] = zeros(model_def["nInp"]);
            end
        end
        preI[:,i] = bayinf_def["Data"]["preInd"][i];

    end

    ltimes2 = zeros(1,length(ltimes));
    ltimes2[1,:] = ltimes;

    toni = collect(1e-9:24*60);

    obsta = zeros(1,obser);
    if typeof(bayinf_def["Data"]["Obs"]) != Array{String,1}
        obsta[1,:] = bayinf_def["Data"]["Obs"];
    end

    y0us = zeros(model_def["nStat"],mcol)
    for i in 1:mcol
        y0us[:,i] = bayinf_def["Data"]["y0"][i];
    end

    ######################## DICTIONARY ########################

    data_multi = Dict(
        "elm" => mrowI, # Maximum length of the rows of the matrices except for time and evnT and pres
        "tml" => trunc(mt+1), # Maximum length of the rows for the time matrix
        "ts" => time.+(1e-9), # Experimental time vector
        "tsl" => convert.(Int,ltimes2), # length of time series per event
        "tsmax" => convert.(Int,round.(mtimes)), # maximum time per event
        "Nsp" => convert.(Int, mdpI2.+1), # length(evnT) or number of event switching times
        "inputs" => inps, # Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
        "evnT" => convert.(Int,round.(evnT)), # event switching times
        "m" => mcol, # Number of time series
        "stsl" => convert.(Int,mdp2), # Number of elements at each time series
        "stslm" => mrow, # Maximum number of time points
        "sts" => convert.(Int,round.(samplingT)), # sampling time points
        "obser" => obser,
        "obSta" => convert.(Int, obsta),
        "nindu" => model_def["nInp"],
        "preInd" => preI,
        "Means" => MEANS,
        #         "Erros" => STDS,
        "Y0us" => y0us
    );
    if length(size(bayinf_def["Data"]["DataError"][1][1])) == 1
        data_multi["Erros"] = STDS;
    else
        data_multi["ErrosMul"] = STDS;
    end

    if model_def["Y0ON"]==true
        data_multi["toni"] = toni;
        data_multi["tonil"] = length(toni);
    end

    if bayinf_def["MultiNormFit"] == true && haskey(bayinf_def["Priors"], "cora") && haskey(bayinf_def["Priors"], "mera") && haskey(bayinf_def["Priors"], "numN")
        data_multi["mera"] = bayinf_def["Priors"]["mera"];
        data_multi["cora"] = bayinf_def["Priors"]["cora"];
        data_multi["numN"] = bayinf_def["Priors"]["numN"];
    end

    return(data_multi)
end

## function that returns all the necessary elements for the stan inference if more control over is desired.
function getStanInferenceElements(model_def, bayinf_def)


    bayinf_def = checkStructBayInf(model_def, bayinf_def);
    bayinf_def = genStanModel(model_def, bayinf_def);


    modelpath = bayinf_def["ModelPath"];

    inferdata = restructureDataInference(model_def, bayinf_def);

    stream = open(modelpath,"r");
    Model = read(stream,String);
    close(stream);

    if bayinf_def["StanSettings"] == []
        StanModel = Stanmodel(name=string(model_def["NameF"], "_Stan_",bayinf_def["flag"]), model=Model, nchains =4,
            num_samples  = 2000, num_warmup = 1000, printsummary=true);
    else
        set_cmdstan_home!(bayinf_def["StanSettings"]["cmdstan_home"]);

        StanModel = Stanmodel(name=string(model_def["NameF"], "_Stan_",bayinf_def["flag"]),
            model=Model, nchains =bayinf_def["StanSettings"]["nchains"],
            num_samples  = bayinf_def["StanSettings"]["nsamples"],
            num_warmup = bayinf_def["StanSettings"]["nwarmup"], printsummary=bayinf_def["StanSettings"]["printsummary"]);
        if bayinf_def["StanSettings"]["maxdepth"] != []
            StanModel.method.algorithm.engine.max_depth = bayinf_def["StanSettings"]["maxdepth"];
        end
        if bayinf_def["StanSettings"]["adaptdelta"] != []
            StanModel.method.adapt.delta = bayinf_def["StanSettings"]["adaptdelta"];
        end
        if bayinf_def["StanSettings"]["jitter"] != []
            StanModel.method.algorithm.stepsize_jitter = bayinf_def["StanSettings"]["jitter"];
        end
    end

    if bayinf_def["StanSettings"]["init"] != [] && bayinf_def["Priors"] <: Dict
        if bayinf_def["Priors"]["transpars"] == []
            init = bayinf_def["StanSettings"]["init"];
        else
            init = reparamDictStan(bayinf_def["StanSettings"], bayinf_def);
        end
    else
        init = [];
    end

    println("")
    println("----------------------------------------- STAN INFERENCE INFO -----------------------------------------")
    println("Stan inference will NOT be initialised due to your input. Once you are ready, please run Stan by using the command: ")
    println("         rc, chns, cnames = stan(StanModel, inferdata, init = init)")
    println("To automatically save your results, you can use: ")
    println("         saveStanResults(rc, chns, cnames, model_def, bayinf_def)")
    println("--------------------------------------------------------------------------------------")
    println("")

    return modelpath, Model, StanModel, inferdata, init, model_def, bayinf_def

end

## Function to save stan RESULTS
function saveStanResults(rc, chns, cnames, model_def, bayinf_def)
    # Generate results directory
    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\StanResults"))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\StanResults"))
    end

    JLD.save(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\StanResults\\StanResults_", model_def["NameF"],"_", bayinf_def["flag"], ".jld"),"cnames",cnames, "chns",chns, "rc", rc)


    poster = chns[:,end-(model_def["nPar"]-1):end,1];
    for c in 2:size(chns)[3]
        poster = vcat(poster, chns[:,end-(model_def["nPar"]-1):end,c])
    end

    dfP = DataFrame(poster);
    CSV.write(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\StanResults\\StanResults_", model_def["NameF"],"_", bayinf_def["flag"], "_Posterior.csv"), dfP);

    println("")
    println("----------------------------------------- STAN RESULTS -----------------------------------------")
    println("Stan Inference results are saved in the directory: ")
    println(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\StanResults"))
    println("--------------------------------------------------------------------------------------")
    println("")

    return poster

end


# Function that runs a stan inference if the settings are given
function runStanInference(model_def, bayinf_def)

    bayinf_def = checkStructBayInf(model_def, bayinf_def);
    bayinf_def = genStanModel(model_def, bayinf_def);

    modelpath = bayinf_def["ModelPath"];

    inferdata = restructureDataInference(model_def, bayinf_def);

    stream = open(modelpath,"r");
    Model = read(stream,String);
    close(stream);

    if bayinf_def["StanSettings"] == []
        StanModel = Stanmodel(name=string(model_def["NameF"], "_Stan_",bayinf_def["flag"]), model=Model, nchains =4,
            num_samples  = 2000, num_warmup = 1000, printsummary=true);
    else
        set_cmdstan_home!(bayinf_def["StanSettings"]["cmdstan_home"]);

        StanModel = Stanmodel(name=string(model_def["NameF"], "_Stan_",bayinf_def["flag"]),
            model=Model, nchains =bayinf_def["StanSettings"]["nchains"],
            num_samples  = bayinf_def["StanSettings"]["nsamples"],
            num_warmup = bayinf_def["StanSettings"]["nwarmup"], printsummary=bayinf_def["StanSettings"]["printsummary"]);
        if bayinf_def["StanSettings"]["maxdepth"] != []
            StanModel.method.algorithm.engine.max_depth = bayinf_def["StanSettings"]["maxdepth"];
        end
        if bayinf_def["StanSettings"]["adaptdelta"] != []
            StanModel.method.adapt.delta = bayinf_def["StanSettings"]["adaptdelta"];
        end
        if bayinf_def["StanSettings"]["jitter"] != []
            StanModel.method.algorithm.stepsize_jitter = bayinf_def["StanSettings"]["jitter"];
        end
    end

    if bayinf_def["StanSettings"]["init"] != [] && typeof(bayinf_def["Priors"]) <: Dict
        if bayinf_def["Priors"]["transpars"] == []
            init = bayinf_def["StanSettings"]["init"];
        else
            init = reparamDictStan(bayinf_def["StanSettings"]["init"], bayinf_def);
        end
    else
        init = [];
    end

    # Run stan inference
    if init == []
        rc, chns, cnames = stan(StanModel, inferdata);
    else
        rc, chns, cnames = stan(StanModel, inferdata, init = init);
    end

    # Save stan results

    Poster = saveStanResults(rc, chns, cnames, model_def, bayinf_def);

    staninf_res = Dict();
    staninf_res["chains"] = chns;
    staninf_res["cnames"] = cnames;
    staninf_res["posterior"] = Poster;
    staninf_res["Model"] = Model;
    staninf_res["StanModel"] = StanModel;
    staninf_res["InferData"] = inferdata;
    staninf_res["Init"] = init;


    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
    end

    JLD.save(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\StanInferenceResults_", model_def["NameF"],"_", bayinf_def["flag"], ".jld"),"staninf_res",staninf_res, "model_def",model_def, "bayinf_def", bayinf_def)

    println("")
    println("----------------------------------------- RESULTS -----------------------------------------")
    println("Stan Inference process results are saved in the directory: ")
    println(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), ""))
    println("Under the name: ")
    println(string("StanInferenceResults_", model_def["NameF"],"_", bayinf_def["flag"], ".jld"))
    println("--------------------------------------------------------------------------------------")
    println("")

    return staninf_res, model_def, bayinf_def

end


function plotStanResults(staninf_res, model_def, bayinf_def)

    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
    end

    # Plot 1, the bivariate distributions
    pc = corrplot(staninf_res["posterior"], label = model_def["parName"], size = [4000, 4000]);
    savefig(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\PosteriorPlots_", bayinf_def["flag"], ".png"));

    # Plot 2, the simulations
    # Simulation Results
    simul_def = defSimulStruct()
    simul_def = SimToMle(simul_def, bayinf_def["Data"])
    simul_def["plot"] = false
    simul_def["theta"] = convert(Array, staninf_res["posterior"]);
    simul_def["flag"] = "StanInferResults";

    simuls, model_def, simul_def = simulateODEs(model_def, simul_def);

    for i in 1:bayinf_def["Data"]["Nexp"]
        tit = "";
        yl1 = "";
        for k in 1:length(bayinf_def["Data"]["Obs"])
            if typeof(bayinf_def["Data"]["Obs"]) == Array{Int,1};
                tit = hcat(tit, string(model_def["stName"][bayinf_def["Data"]["Obs"][k]]));
            else
                tit = hcat(tit, string(bayinf_def["Data"]["Obs"][k]));
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

        SimObs = selectObsSim_te(simuls[string("Exp_",i)], bayinf_def["Data"]["Obs"],model_def["stName"]);

        errorss = zeros(size(SimObs)[1],length(bayinf_def["Data"]["Obs"]));
        if length(size(bayinf_def["Data"]["DataError"][i][1])) == 1
            for k in 1:length(bayinf_def["Data"]["Obs"])
                errorss[:,k] = bayinf_def["Data"]["DataError"][i][k];
            end
        else
            for k in 1:length(bayinf_def["Data"]["Obs"])
                errorss[:,k] = [bayinf_def["Data"]["DataError"][i][k][f,f] for f in 1:size(bayinf_def["Data"]["DataError"][i][k])[1]];
            end
        end

        # Elements for the inducers
        tu = Array{Array{Any,1}}(undef, length(bayinf_def["Data"]["Obs"])+model_def["nInp"])
        [tu[k] = [] for k in 1:length(bayinf_def["Data"]["Obs"])]
        su = Array{Array{Any,1}}(undef, length(bayinf_def["Data"]["Obs"])+model_def["nInp"])
        [su[k] = [] for k in 1:length(bayinf_def["Data"]["Obs"])]
        for k in 1:model_def["nInp"]
            tu[k+length(bayinf_def["Data"]["Obs"])] = round.(bayinf_def["Data"]["switchT"][i]);
            su[k+length(bayinf_def["Data"]["Obs"])] = vcat(bayinf_def["Data"]["uInd"][i][:,(k)], bayinf_def["Data"]["uInd"][i][end,(k)])
        end

        pl = plot(round.(bayinf_def["Data"]["tsamps"][i]), bayinf_def["Data"]["DataMean"][i],yerror = errorss,
        layout=length(bayinf_def["Data"]["Obs"])+model_def["nInp"], label = "", title = tit, xlabel = "time", ylabel = yuu,
                color = "gray", size = [2000, 1200]);

        for p in 1:size(SimObs)[3]
            pl = plot!(round.(bayinf_def["Data"]["tsamps"][i]),SimObs[:,:,p], layout=length(bayinf_def["Data"]["Obs"]), label = "");
        end

        pl = plot!(tu, su, layout=length(bayinf_def["Data"]["Obs"])+model_def["nInp"], label = "", linetype=:step, title = tuu);

        savefig(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\PlotStanInfResults_Exp", i,"_", bayinf_def["flag"], ".png"));
    end

end

# Main Stan inference section function
function StanInfer(model_def, bayinf_def)
    cudi = pwd(); #@__DIR__;
    bayinf_def = checkStructBayInf(model_def, bayinf_def);

    if bayinf_def["runInf"] == false

        println("----------------------------------------- STAN INFERENCE INFO -----------------------------------------")
        println("No inference will be performed since this is what has been selected. Structure elements")
        println("to run a Stan inference will be given instead.")
        println("--------------------------------------------------------------------------------------")

        modelpath, Model, StanModel, inferdata, init, model_def, bayinf_def = getStanInferenceElements(model_def, bayinf_def);
        stan_struct = Dict();
        stan_struct["modelpath"] = modelpath;
        stan_struct["Model"] = Model;
        stan_struct["StanModel"] = StanModel;
        stan_struct["inferdata"] = inferdata;
        stan_struct["init"] = init;

    elseif bayinf_def["runInf"] == true
        if typeof(bayinf_def["Priors"]) <: Dict
            if bayinf_def["Priors"]["pridis"] == []

                println("----------------------------------------- STAN INFERENCE INFO -----------------------------------------")
                println("No inference will be performed since there is no priors definition. Structure elements")
                println("to run a Stan inference will be given instead.")
                println("--------------------------------------------------------------------------------------")

                modelpath, Model, StanModel, inferdata, init, model_def, bayinf_def = getStanInferenceElements(model_def, bayinf_def);
                stan_struct = Dict();
                stan_struct["modelpath"] = modelpath;
                stan_struct["Model"] = Model;
                stan_struct["StanModel"] = StanModel;
                stan_struct["inferdata"] = inferdata;
                stan_struct["init"] = init;

            else

                println("----------------------------------------- STAN INFERENCE INFO -----------------------------------------")
                println("Inference will be performed since there is enough information from the user:")
                println("Stan temporary files will be stored in your current working directory in the folder tmp.")
                println("--------------------------------------------------------------------------------------")

                stan_struct, model_def, bayinf_def = runStanInference(model_def, bayinf_def);

                if bayinf_def["plot"] == true
                    plotStanResults(stan_struct,model_def,bayinf_def)
                    println("")
                    println("----------------------------------------- PLOTS -----------------------------------------")
                    println("Stan Inference Results PLOTS are saved in the directory: ")
                    println(string("                 ", cudi, "\\Results\\", model_def["NameF"],"_",today() ))
                    println(string("Under the names PlotStanInfResults_Exp(i)_",bayinf_def["flag"],".png for simulations"))
                    println(string("And ", "PosteriorPlots_", bayinf_def["flag"], ".png for the posteriors."))
                    println("--------------------------------------------------------------------------------------")
                    println("")
                end

            end
        elseif bayinf_def["Priors"] == []
            println("----------------------------------------- STAN INFERENCE INFO -----------------------------------------")
            println("No inference will be performed since there is no priors definition. Structure elements")
            println("to run a Stan inference will be given instead.")
            println("--------------------------------------------------------------------------------------")

            modelpath, Model, StanModel, inferdata, init, model_def, bayinf_def = getStanInferenceElements(model_def, bayinf_def);
            stan_struct = Dict();
            stan_struct["modelpath"] = modelpath;
            stan_struct["Model"] = Model;
            stan_struct["StanModel"] = StanModel;
            stan_struct["inferdata"] = inferdata;
            stan_struct["init"] = init;

        elseif bayinf_def["Priors"] == String
            println("----------------------------------------- STAN INFERENCE INFO -----------------------------------------")
            println("Inference will be performed since there is enough information from the user:")
            println("Stan temporary files will be stored in your current working directory in the folder tmp.")
            println("--------------------------------------------------------------------------------------")

            stan_struct, model_def, bayinf_def = runStanInference(model_def, bayinf_def);

            if bayinf_def["plot"] == true
                plotStanResults(stan_struct,model_def,bayinf_def)
                println("")
                println("----------------------------------------- PLOTS -----------------------------------------")
                println("Stan Inference Results PLOTS are saved in the directory: ")
                println(string("                 ", cudi, "\\Results\\", model_def["NameF"],"_",today() ))
                println(string("Under the names PlotStanInfResults_Exp(i)_",bayinf_def["flag"],".png for simulations"))
                println(string("And ", "PosteriorPlots_", bayinf_def["flag"], ".png for the posteriors."))
                println("--------------------------------------------------------------------------------------")
                println("")
            end

        end
    end


    return stan_struct, model_def, bayinf_def
end
