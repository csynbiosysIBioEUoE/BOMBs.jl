
        

function PLacExample_VS_TestPLac_M2Utility(ins)

    # Definition of the inputs for the ODEs
    IPTG2,aTc1 = ins;
    inputsM1 = [0,IPTG2,0,IPTG2,0,IPTG2];
    inputsM2 = [0,aTc1,IPTG2,0,0,aTc1,IPTG2,0,0,aTc1,IPTG2,0];

    # Solve ODEs

    solM1 = PLacExample_SolveAll(tsMS, pD1MS, spMS, inputsM1, ivss1MS, sampsMS, pre1MS);
    solM2 = TestPLac_M2_SolveAll(tsMS, pD2MS, spMS, inputsM2, ivss2MS, sampsMS, pre2MS);

    # Extracte wanted vectors (observables) with time reduction
    Obs1_M1 = 3 .*solM1[:,4,:]; 

    Obs1_M2 = 3 .*solM2[:,4,:]; 

    
    mu1_M1 = mean(Obs1_M1, dims=2); 

    sd1_M1 = cov(Obs1_M1, dims=2).+(0.1*Matrix((Diagonal(ones(size(Obs1_M1)[1]))))); 


    mu1_M2 = mean(Obs1_M2, dims=2); 

    sd1_M2 = cov(Obs1_M2, dims=2).+(0.1*Matrix((Diagonal(ones(size(Obs1_M2)[1]))))); 


    BHD1 = BhattacharyyaDist(mu1_M1[:,1], mu1_M2[:,1], sd1_M1, sd1_M2); 


    util = mean(BHD1);

        

    return(util)
end

    