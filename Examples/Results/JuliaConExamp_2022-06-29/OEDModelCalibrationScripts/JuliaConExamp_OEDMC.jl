
        

    function JuliaConExampUtility(ins)

        # Definition of the inputs for the ODEs
        if length(ins) == 1
            U1,U2,U3 = ins[1];
        else
            U1,U2,U3 = ins;
        end
        inputsMC = [U1,U2,U3];

        # Solve ODEs
        solMC = JuliaConExamp_SolveAll(tsMC, pD1MC, spMC, inputsMC, ivss1MC, sampsMC, pre1MC);

        # Extracte wanted vectors (observables) with time reduction
        Obs1_MC = solMC[:,2,:]; 


    
        LowqObs1_MC = zeros(size(solMC)[1]); 

        HighqObs1_MC = zeros(size(solMC)[1]); 

        for i in 1:size(solMC)[1]
            LowqObs1_MC[i] = percentile(Obs1_MC[i,:],0.5); 

            HighqObs1_MC[i] = percentile(Obs1_MC[i,:],99.5); 

        end

        # Compute Euclidean distances
        EuObs1_MC = sqrt(sum((LowqObs1_MC-HighqObs1_MC).^2)); 

        util = (EuObs1_MC)*(1/1);
    

        return(util)
    end

        