


function PLacExample_VS_TestPLac_M2Utility(ins)

    # Definition of the inputs for the ODEs
    if length(ins) == 1
        IPTG2,IPTG3,aTc = ins[1];
    else
        IPTG2,IPTG3,aTc = ins;
    end
    inputsM1 = [0,IPTG2,IPTG3];
    inputsM2 = [0,aTc,IPTG2,aTc,IPTG3,aTc];

    # Solve ODEs

    solM1 = PLacExample_SolveAll(tsMS, pD1MS, spMS, inputsM1, ivss1MS, sampsMS, pre1MS);
    solM2 = TestPLac_M2_SolveAll(tsMS, pD2MS, spMS, inputsM2, ivss2MS, sampsMS, pre2MS);

    # Extracte wanted vectors (observables) with time reduction
    Obs1_M1 = solM1[:,4,:]./maximum(solM1[:,4,:]);

    Obs1_M2 = solM2[:,4,:]./maximum(solM2[:,4,:]);


    regtmp_m1 = Array{Any,1}(undef, 1);
    reg_m1 = zeros(1);
    regtmp_m2 = Array{Any,1}(undef, 1);
    reg_m2 = zeros(1);

    regtmp_m1[1] = diag(cov(Obs1_M1, dims=2));

    regtmp_m2[1] = diag(cov(Obs1_M2, dims=2));

    for k in 1:length(regtmp_m1)
        regtmp_m1[k][regtmp_m1[k].==0] .= maximum(regtmp_m1[k]); # This is to avoid the case of having a variance of 0 which would cause issues (semi-definite matrix)!
        reg_m1[k] = (maximum(regtmp_m1[k])*0.1);
        if reg_m1[k] > 0.1
            reg_m1[k] = 0.1;
        end
        regtmp_m2[k][regtmp_m2[k].==0] .= maximum(regtmp_m2[k]);
        reg_m2[k] = (maximum(regtmp_m2[k])*0.1);
        if reg_m2[k] > 0.1
            reg_m2[k] = 0.1;
        end
    end

    println(1)

    mu1_M1 = mean(Obs1_M1, dims=2);

    sd1_M1 = cov(Obs1_M1, dims=2).+(0.1*Matrix((Diagonal(ones(size(Obs1_M1)[1])))));


    mu1_M2 = mean(Obs1_M2, dims=2);

    sd1_M2 = cov(Obs1_M2, dims=2).+(0.1*Matrix((Diagonal(ones(size(Obs1_M2)[1])))));


    BHD1 = BhattacharyyaDist(mu1_M1[:,1], mu1_M2[:,1], sd1_M1, sd1_M2);


    util = mean(BHD1);



    return(util)
end
