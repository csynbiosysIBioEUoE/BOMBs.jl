
        

    function PLacExampleUtility(ins)

        # Definition of the inputs for the ODEs
        if length(ins) == 1
            IPTG2,IPTG3 = ins[1];
        else
            IPTG2,IPTG3 = ins;
        end
        inputsMC = [0,IPTG2,IPTG3,IPTG2,IPTG3,IPTG2];

        # Solve ODEs
        solMC = PLacExample_SolveAll(tsMC, pD1MC, spMC, inputsMC, ivss1MC, sampsMC, pre1MC);

        # Extracte wanted vectors (observables) with time reduction
        Obs1_MC = 3 .*solMC[:,4,:]; 


    
        EntObs1_MC = zeros(size(solMC)[1]); 


        dists = [Beta, Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform];
        dists2 = [Exponential, LogNormal, Normal, Gamma, Laplace, Pareto, Rayleigh, Cauchy, Uniform];

        names = ["timePoint"];

        for j in 1:size(solMC)[1]
            fitts1 = Dict(); 
        bestfit1 = Dict(); 
        bestfitInd1 = zeros(1,1); 

    
            try
                fitts1[names[1]] = fit.(dists, Ref(Obs1_MC[j,:]));
                bestfitInd1[1] = findmax(loglikelihood.(fitts1[names[1]], Ref(Obs1_MC[j,:])))[2];
                bestfit1[names[1]] = fitts1[names[1]][findmax(loglikelihood.(fitts1[names[1]], Ref(Obs1_MC[j,:])))[2]];
            catch
                fitts1[names[1]] = fit.(dists2, Ref(Obs1_MC[j,:]));
                bestfitInd1[1] = findmax(loglikelihood.(fitts1[names[1]], Ref(Obs1_MC[j,:])))[2]+1;
                bestfit1[names[1]] = fitts1[names[1]][findmax(loglikelihood.(fitts1[names[1]], Ref(Obs1_MC[j,:])))[2]];
            end

            EntObs1_MC[j] = entropy(bestfit1[names[1]]); 

        end

        HES = zeros(1,1);
        HES[1,1] = sum(EntObs1_MC); 

        util = sum(HES)
    

        return(util)
    end

        