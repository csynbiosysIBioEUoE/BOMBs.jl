
## Generate samples from the prior distribution
function genSamplesPrior(model_def, bayinf_def, nsamps, mu=[],coo=[])

    sam = ["lognormal", "normal", "uniform", "bernoulli", "binomial", "neg_binomial", "poisson", "student_t", "cauchy",
    "double_exponential", "gumbel", "chi_square", "exponential", "gamma", "weibull", "frechet",  "rayleigh", "pareto", "beta", "von_mises", "dirichlet"];

    println("-------------------------------------- WARNING --------------------------------------")
    println("This function will only work if the contents from bayinf_def[","\"","Priors","\"","][","\"","pridis","\"","] are")
    println("distributions wiht numbers defined in the parameters. If the parameters given to the ")
    println("distributions are variables instead of numbers an error will be prompted since sampling will")
    println("not work. So no hierarchical priors. Sorry for the inconvenience. ")
    println(" ")
    println("The only exception is the multivariate normal case where you can introduce the mean vector")
    println("and the covariance matrix as 4th and 5th arguments in the function in this order.")
    println("However, truncations will not be possible and only 1 MultiVariate normal can be assessed. ")
    println("")
    println("Also, sampling from all types of stan distributions is not implemnted yet. For now the ")
    println("available options are: ")
    println(join(sam, ", "))
    println(" ")
    println("-------------------------------------------------------------------------------------")
    println(" ")

    origsamp = Array{String,1}(undef, length(bayinf_def["Priors"]["pridis"]));

    for i in 1:length(bayinf_def["Priors"]["pridis"])
        origsamp[i] = replace(bayinf_def["Priors"]["pridis"][i], bayinf_def["Priors"]["pridis"][i][1:findfirst("~", bayinf_def["Priors"]["pridis"][i])[1]] => "")
        origsamp[i] = replace(origsamp[i], ";" => "")
        origsamp[i] = replace(origsamp[i], "\n" => "")
        if occursin("multi_normal", bayinf_def["Priors"]["pridis"][i])
            if mu == [] || coo ==[]
                println("Process Stopped, please give means and covariances for the multinormal");
                return
            else
                origsamp[i] = replace(bayinf_def["Priors"]["pridis"][i], bayinf_def["Priors"]["pridis"][i][findfirst("~", bayinf_def["Priors"]["pridis"][i])[1]:end] => string("= MvNormal(",mu,", ",coo,")"))
            end
        elseif occursin("lognormal", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "lognormal" => "LogNormal");
        elseif occursin("normal", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "normal" => "Normal");
        elseif occursin("uniform", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "uniform" => "Uniform");
        elseif occursin("bernoulli", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "bernoulli" => "Bernoulli");
        elseif occursin("binomial", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "binomial" => "Binomial");
        elseif occursin("neg_binomial", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "neg_binomial" => "NegativeBinomial");
        elseif occursin("poisson", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "poisson" => "Poisson");
        elseif occursin("student_t", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "student_t" => "TDist");
        elseif occursin("cauchy", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "cauchy" => "Cauchy");
        elseif occursin("double_exponential", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "double_exponential" => "Laplace");
    #     elseif occursin("logistic", bayinf_def["Priors"]["pridis"][i])
    #         origsamp[i] = replace(origsamp[i], "" => "");
        elseif occursin("gumbel", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "gumbel" => "Gumbel");
        elseif occursin("chi_square", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "chi_square" => "Chisq");
        elseif occursin("exponential", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "exponential" => "Exponential");
        elseif occursin("gamma", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "gamma" => "Gamma");
        elseif occursin("weibull", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "weibull" => "Weibull");
        elseif occursin("frechet", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "frechet" => "Frechet");
        elseif occursin("rayleigh", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "rayleigh" => "Rayleigh");
        elseif occursin("pareto", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "pareto" => "Pareto");
        elseif occursin("beta", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "beta" => "Beta");
        elseif occursin("von_mises", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "von_mises" => "VonMises");
        elseif occursin("dirichlet", bayinf_def["Priors"]["pridis"][i])
            origsamp[i] = replace(origsamp[i], "dirichlet" => "Dirichlet");
        end
    end

    for i in 1:length(bayinf_def["Priors"]["pars"])
        if !occursin("MvNormal", origsamp[i])
            if occursin("lower=", bayinf_def["Priors"]["pars"][i])
                ind1 = findfirst("lower=", bayinf_def["Priors"]["pars"][i])[end]+1;
                ind2 = "";
                try
                    ind2 = findfirst(",", bayinf_def["Priors"]["pars"][i])[1]-1;
                catch
                    ind2 = findfirst(">", bayinf_def["Priors"]["pars"][i])[1]-1;
                end
                lob = bayinf_def["Priors"]["pars"][i][ind1:ind2]
            else
                lob = [];
            end

            if occursin("upper=", bayinf_def["Priors"]["pars"][i])
                ind3 = findfirst("upper=", bayinf_def["Priors"]["pars"][i])[end]+1;
                ind4 = findfirst(">", bayinf_def["Priors"]["pars"][i])[1]-1;
                upb = bayinf_def["Priors"]["pars"][i][ind3:ind4]
            else
                upb = [];
            end

            if lob == [] && upb !=[]
                lob = "convert(Float64,*Inf)"
            elseif lob != [] && upb ==[]
                upb = "convert(Float64,Inf)"
            end

            if lob != [] && upb != []
                origsamp[i] = string("truncated(",origsamp[i],", ",lob,", ",upb,")")
            else
                origsamp[i] = "";
            end
        end
    end

    origsampTR = Array{String,1}(undef, length(bayinf_def["Priors"]["pars"]));

    for i in 1:length(bayinf_def["Priors"]["pars"])
        for j in 1:length(bayinf_def["Priors"]["transpars"])
            if occursin(model_def["parName"][i], bayinf_def["Priors"]["transpars"][j])
                origsampTR[i] = replace(bayinf_def["Priors"]["transpars"][j], ";"=>"")
                origsampTR[i] = replace(origsampTR[i], "\n"=>"")
                origsampTR[i] = replace(origsampTR[i], origsampTR[i][1:findfirst("=", origsampTR[i])[1]]=>"")
            end
        end
        if !isassigned(origsampTR,i)
            origsampTR[i] = ""
        end
    end

    pretransP = zeros(model_def["nPar"],nsamps);
    for i in 1:length(origsamp)
        tpp = Meta.parse(string("rand(", origsamp[i], ",", nsamps, ")"));
        valp = @eval $tpp;
        pretransP[i,:] = valp;
    end

    priorSamps = zeros(nsamps, model_def["nPar"]);
    for i in 1:model_def["nPar"]
        if origsampTR != ""
            tmpex = origsampTR[i];
            if occursin("*", origsampTR[i])
                tmpex = replace(tmpex, "*"=>".*")
            end
            if occursin("/", origsampTR[i])
                tmpex = replace(tmpex, "/"=>"./")
            end
            if occursin("+", origsampTR[i])
                tmpex = replace(tmpex, "+"=>".+")
            end
            if occursin("-", origsampTR[i])
                tmpex = replace(tmpex, "-"=>".-")
            end
            if occursin("^", origsampTR[i])
                tmpex = replace(tmpex, "^"=>".^")
            end

            for j in 1:model_def["nPar"]
                if occursin(model_def["parName"][j], tmpex)
                    tmpex = replace(tmpex, model_def["parName"][j] => string(pretransP[j,:]));
                end
            end
            tppT = Meta.parse(tmpex);
            valpT = @eval $tppT;
            priorSamps[:,i] = valpT;
        else
            priorSamps[:,i] = pretransP[i,:];
        end
    end

    return priorSamps
end

##
# --------------------------------- RELATIVE ENTROPY FUNCTION --------------------------------- #

# Calculation of the entropy of prior and posterior results for an inference result using the Taylor series approximation
# for Gaussian Mixtures proposed in "On Entropy Approximation for Gaussian Mixture Random Vectors". For all entropy calculations,
# x represents a random vector, MU the vectors of means for the Gaussian mixtures, E the covariance matrices for the Gaussian
# mixtures and w the vector of weights for the Gaussian Mixtures. As inputs the function takes a character string with the
# name of the experimental profile to be analysed.


## Calculation of upper bound for the posterior Entropy

function H_Upper(w,E)

    comp = length(w);
    dims = size(E[1])[1];
    Hu = zeros(1,comp);

    for i in 1:comp
        hu = w[i]*(-log(w[i])+0.5*log(det(Matrix(E[i]))*(2*pi*exp(1))^dims))
        Hu[1,i] = hu;
    end

    return(sum(Hu))
end

## PDF of a multivariate Gaussian distribution

function mvGauss(x, MU, E)

    gx = (1/sqrt(det(Matrix(E))*(2*pi)^2))*exp((-1/2)*(x-MU)'*inv(Matrix(E))*(x-MU))

    return(gx)

end

## Lower bound for the posterior entropy
function H_Lower(w, E, MU)

    comp = length(w);
    Hl = zeros(1,comp);

    for i in 1:comp
        inter = zeros(1,comp);
        for j in 1:comp
            insu = w[j]*BOMBS.mvGauss(MU[i,:], MU[j,:], E[i]+E[j]);
            inter[1,j] = insu;
        end
        Hl[1,i] = w[i]*log(sum(inter))
    end

    return(-sum(Hl))

end

## Function for a multivariate PDF of a Gaussian Mixture
function GaussMix(x, MU, E, w)

    comp = length(w);
    FX = zeros(1,comp);

    for i in 1:comp
        fx = w[i]*(1/sqrt(det(Matrix(E[i]))*(2*pi)^2))*exp((-1/2)*(x-MU[i,:])'*inv(Matrix(E[i]))*(x-MU[i,:]));
        FX[i] = fx;
    end

    return(sum(FX))

end

## Zero-order Taylor Series Expansion

function ZOTSE(MU, E, w)

    comp = length(w);
    ZO = zeros(1,comp);

    for i in 1:comp
        zo = w[i]*log(GaussMix(MU[i,:], MU, E, w));
        ZO[i] = zo;
    end

    return(-sum(ZO))
end

## Function to calculate the most computationaly expensive part of the Taylor Expansion component
function GaussMix2(x)

    comp = length(mvPro);
    FX = zeros(1,comp);

    for i in 1:comp
        fx = mvPro[i]*(1/sqrt(det(Matrix(mvCovar[i]))*(2*pi)^2))*exp((-1/2)*(x-mvMean[i,:])'*inv(Matrix(mvCovar[i]))*(x-mvMean[i,:]));
        FX[i] = fx;
    end

    return(sum(FX))

end

function FMix(x, MU, E, w)

    comp = length(w);
    dims = size(E[1])[1];

    F_x = Array{Any,1}(undef,comp);

    for j in 1:comp
        inter = (1/GaussMix(x, MU, E, w))*(x-MU[j,:])*(Calculus.gradient(GaussMix2, x))'+(x-MU[j,:])*(inv(Matrix(E[j]))*(x-MU[j,:]))'-Diagonal(ones(dims));
        fx = w[j]*inv(Matrix(E[j]))*(inter)*(mvGauss(x , MU[j,:], E[j]));
        F_x[j] = fx;
    end
    enn = sum(F_x)/GaussMix(x, MU, E, w)

    return(enn)
end

## Calculation of the Second-Order Taylor Series expanssion term

function SOTSE(MU, E, w)

    comp = length(w);
    dims = size(E[1])[1];

    H2 = zeros(1,comp);

    for i in 1:comp
        println(string("              Iteration ", i, " of ", comp))
        inter2 = (w[i]/2)*sum(FMix(MU[i,:],MU,E,w).*Matrix(E[i]))
        H2[i] = inter2;
    end

    return(sum(H2))

end

## Function to estimate the ENTROPY
# Function to compute entropy from distribution samples
function computeH(sampl, model_def, tag = "")

    if length(size(sampl)) == 1
        sampl = reshape(sampl,length(sampl), 1);
    end

    if size(sampl)[1]==model_def["nPar"]
        sampl = convert(Array, sampl');
    elseif size(sampl)[2]==model_def["nPar"]
        nothing
    else
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the samples that you have introduced. There is something odd in there. ")
        return
    end

    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
    end

    # Estimate joint posterior as Gaussian Mixtures
    # Perform clustering in a range between 1 to 40 components and use LLK to asses the optimum value
    grid_search = ScikitLearnBase.fit!(GridSearchCV(GMM(n_components=1, kind=:full), Dict(:n_components=>collect(1:40))), sampl)

    # Get and store best results
    comps = length(grid_search.best_estimator_.w); # Number of components
    weis = grid_search.best_estimator_.w; # Weights for each component
    meaas = grid_search.best_estimator_.μ; # Means for each component
    covas = [inv(Matrix(grid_search.best_estimator_.Σ[i]))*inv(Matrix(grid_search.best_estimator_.Σ[i])') for i in 1:length(grid_search.best_estimator_.Σ)]; # Covariance Matrix for each component

    entro_res = Dict();
    entro_res["BestGC"] = Dict();
    entro_res["BestGC"]["c"] = comps;
    entro_res["BestGC"]["w"] = weis;
    entro_res["BestGC"]["mu"] = meaas;
    entro_res["BestGC"]["cov"] = covas;

    global mvMean = meaas;
    global mvCovar = covas;
    global mvPro = weis;

    # Upper bound for the posterior entropy
    H_U_Post = H_Upper(weis, covas);
    entro_res["HUpper"] = H_U_Post;

    # Lower bound for the posterior entropy
    H_L_Post = H_Lower(weis, covas, meaas);
    entro_res["HLower"] = H_L_Post;

    # Zero-order Taylor series result for the approximation of the posterior entropy
    ZO_Expan = ZOTSE(meaas, covas, weis);

    # Second-order Taylor series result term for the approximation of the posterior entropy
    SO_Expan = SOTSE(meaas, covas, weis);

    # Aproximation of the Entropy for the posterior
    H_posterior = ZO_Expan-SO_Expan;
    entro_res["HTheta"] = H_posterior;

    # Assess if entropy estimation is between bounds
    if H_posterior>=H_L_Post && H_posterior<=H_U_Post
        eebb = true;
    else
        eebb = false;
    end
    entro_res["BetweenBounds"] = eebb;

    entro_res["Cov"] = cov(sampl);

    JLD.save(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\EntropyResults_", tag, ".jld"), "entro_res", entro_res)

    println("")
    println("----------------------------------------- ENTROPY RESULTS -----------------------------------------")
    println("Entropy results are saved in the directory: ")
    println(string("    ",cudi, "\\Results\\", model_def["NameF"],"_",today(), ""))
    println("Under the name: ")
    println(string("            EntropyResults_", tag, ".jld"))
    println("--------------------------------------------------------------------------------------")
    println("")

    return(entro_res)

end

## Function to compute entropy gain between prior and posterior (based on samples)
# In this case bayinf_def can only have the entry Priors (nothing else will be checked)
function computeHgain(prior, posterior, model_def, tag = "")

    cudi = pwd(); #@__DIR__;
    if !isdir(string(cudi,"\\Results"))
        mkdir(string(cudi,"\\Results"))
    end
    if !isdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
        mkdir(string(cudi, "\\Results\\", model_def["NameF"],"_",today()))
    end

    if typeof(prior) <: Dict
        pris = genSamplesPrior(model_def, prior, 10000);
    else
        pris = prior;
    end;

    entro_prior = computeH(pris, model_def, "prior");
    entro_posterior = computeH(posterior, model_def, "posterior");

    Relative_H = entro_prior["HTheta"]-entro_posterior["HTheta"];

    Hdiff_res = Dict();
    Hdiff_res["HPrior"] = entro_prior;
    Hdiff_res["HPosterior"] = entro_posterior;
    Hdiff_res["HDiff"] = Relative_H;

    JLD.save(string(cudi, "\\Results\\", model_def["NameF"],"_",today(), "\\EntropyGainResults_", tag, ".jld"), "Hdiff_res", Hdiff_res)

    println("")
    println("----------------------------------------- ENTROPY RESULTS -----------------------------------------")
    println("Entropy results are saved in the directory: ")
    println(string("    ",cudi, "\\Results\\", model_def["NameF"],"_",today(), ""))
    println("Under the name: ")
    println(string("            EntropyGainResults_", tag, ".jld"))
    println("--------------------------------------------------------------------------------------")
    println("")

    return(Hdiff_res)

end
