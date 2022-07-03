
        

function JuliaConExamp_VS_JuliaConExampModel2Utility(ins)

    # Definition of the inputs for the ODEs
    if length(ins) == 1
        U1,U2,U3 = ins[1];
    else
        U1,U2,U3 = ins;
    end
    inputsM1 = [U1,U2,U3];
    inputsM2 = [U1,U2,U3];

    # Solve ODEs

    solM1 = JuliaConExamp_SolveAll(tsMS, pD1MS, spMS, inputsM1, ivss1MS, sampsMS, pre1MS);
    solM2 = JuliaConExampModel2_SolveAll(tsMS, pD2MS, spMS, inputsM2, ivss2MS, sampsMS, pre2MS);

    # Extracte wanted vectors (observables) with time reduction
    Obs1_M1 = solM1[:,2,:]; 

    Obs1_M2 = solM2[:,2,:]; 

    
    mu1_M1 = mean(Obs1_M1, dims=2); 

    sd1_M1 = cov(Obs1_M1, dims=2).+(0.1*Matrix((Diagonal(ones(size(Obs1_M1)[1]))))); 


    mu1_M2 = mean(Obs1_M2, dims=2); 

    sd1_M2 = cov(Obs1_M2, dims=2).+(0.1*Matrix((Diagonal(ones(size(Obs1_M2)[1]))))); 


    BHD1 = BhattacharyyaDist(mu1_M1[:,1], mu1_M2[:,1], sd1_M1, sd1_M2); 


    util = mean([BHD1])[1];

        

    return(util)
end

    