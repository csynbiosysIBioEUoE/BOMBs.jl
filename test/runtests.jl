using BOMBS
using Test
using CSV
using DataFrames
using LinearAlgebra

# -------------------------------------------------------- MODEL DEFINITION TESTS
model_def = Dict();
model_def["NameF"] = [];
model_def["nStat"] = [];
model_def["nPar"] = [];
model_def["nInp"] = [];
model_def["stName"] = [];
model_def["parName"] = [];
model_def["inpName"] = [];
model_def["eqns"] = [];
model_def["Y0eqs"] = [];
model_def["Y0ON"] = [];
model_def["tols"] = [];
model_def["solver"] = [];

model_def2 = defModStruct();
model_def2["NameF"] = ["Test"];
model_def2["nStat"] = 2;
model_def2["nPar"] = 4;
model_def2["nInp"] = 1;
model_def2["stName"] = ["A", "B"];
model_def2["parName"] = ["k1", "k2", "k3", "k4"];
model_def2["inpName"] = ["inp1"];
model_def2["eqns"] = ["dA = k1*A - k2*A*inp1", "dB = k3*(A^2) - k4*B"];
model_def2["Y0eqs"] = ["A = 1", "B = k3*(A^2)/k4"];
model_def2["Y0ON"] = false;
model_def2["tols"] = [1e-5,1e-5];
model_def2["solver"] = "CVODE_BDF";

# -------------------------------------------------------- SIMULATION TESTS
simul_def = Dict()
simul_def["Nexp"] = [];
simul_def["finalTime"] = [];
simul_def["switchT"] = [];
simul_def["y0"] = [];
simul_def["preInd"] = [];
simul_def["uInd"] = [];
simul_def["theta"] = [];
simul_def["tsamps"] = [];
simul_def["plot"] = [];
simul_def["flag"] = [];

simul_def2 = defSimulStruct();
simul_def2["Nexp"] = 2;
simul_def2["finalTime"] = [40, 40];
simul_def2["switchT"] = [[0,20,40], [0,20,40]];
simul_def2["y0"] = [[0,0], [0,0]];
simul_def2["preInd"] = [[0.1], [0.1]];
simul_def2["uInd"] = [[1,1], [1,1]];
simul_def2["theta"] = [0.1 0.2 0.2 0.02; 0.2 0.1 0.2 0.01; 1 1 1 0.1];
simul_def2["tsamps"] = [[0,5,10,15,20,25,30,35,40], [0,3,5,10,12,15,20,25,27, 30,33, 35,40]];
simul_def2["plot"] = false;
simul_def2["flag"] = "testsim";

simul_def3 = Dict()
simul_def3["ObservablesFile"] = [];
simul_def3["EventInputsFile"] = [];
simul_def3["theta"] = [];
simul_def3["MainDir"] = [];
simul_def3["plot"] = [];
simul_def3["flag"] = [];

simul_defCSV1 = Dict();
simul_defCSV1["ObservablesFile"] = ["TestSimulCSVObs.csv"];
simul_defCSV1["EventInputsFile"] = ["TestSimulCSVInps.csv"];
simul_defCSV1["theta"] = ["Theta1.csv"];
simul_defCSV1["MainDir"] = [];
simul_defCSV1["plot"] = false;
simul_defCSV1["flag"] = "testsimCSV";

simul_defCSV2 = Dict();
simul_defCSV2["Nexp"] = 1;
simul_defCSV2["finalTime"] = [40];
simul_defCSV2["switchT"] = [[0.0,20.0,40.0]];
simul_defCSV2["y0"] = [[0,0]];
simul_defCSV2["preInd"] = [[0.1]];
simul_defCSV2["uInd"] = [reshape([1, 0.5], 2,1)];
simul_defCSV2["theta"] = [0.1 0.2 0.2 0.002];
simul_defCSV2["tsamps"] = [[0,5,10,15,20,25,30,35,40]];
simul_defCSV2["plot"] = false;
simul_defCSV2["flag"] = "testsimCSV";


# -------------------------------------------------------- PSEUDO-DATA TESTS

pseudo_def = Dict();
pseudo_def["Nexp"] = [];
pseudo_def["finalTime"] = [];
pseudo_def["switchT"] = [];
pseudo_def["y0"] = [];
pseudo_def["preInd"] = [];
pseudo_def["uInd"] = [];
pseudo_def["theta"] = [];
pseudo_def["tsamps"] = [];
pseudo_def["plot"] = [];
pseudo_def["flag"] = [];
pseudo_def["Obs"] = [];
pseudo_def["Noise"] = [];

pseudo_def2 = defPseudoDatStruct();
pseudo_def2["Nexp"] = 2;
pseudo_def2["finalTime"] = [40, 40];
pseudo_def2["switchT"] = [[0,20,40], [0,20,40]];
pseudo_def2["y0"] = [[0,0], [0,0]];
pseudo_def2["preInd"] = [[0.1], [0.1]];
pseudo_def2["uInd"] = [[1,1], [1,1]];
pseudo_def2["theta"] = [0.1 0.2 0.2 0.02; 0.2 0.1 0.2 0.01; 1 1 1 0.1];
pseudo_def2["tsamps"] = [[0,5,10,15,20,25,30,35,40], [0,3,5,10,12,15,20,25,27, 30,33, 35,40]];
pseudo_def2["plot"] = false;
pseudo_def2["flag"] = "testPD";
pseudo_def2["Obs"] = ["B*2"];
pseudo_def2["Noise"] = [0.1];

pseudo_def3 = Dict()
pseudo_def3["ObservablesFile"] = [];
pseudo_def3["EventInputsFile"] = [];
pseudo_def3["theta"] = [];
pseudo_def3["MainDir"] = [];
pseudo_def3["plot"] = [];
pseudo_def3["flag"] = [];
pseudo_def3["Obs"] = [];
pseudo_def3["Noise"] = [];

pseudo_defCSV1 = Dict();
pseudo_defCSV1["ObservablesFile"] = ["TestSimulCSVObsPD.csv"];
pseudo_defCSV1["EventInputsFile"] = ["TestSimulCSVInpsPD.csv"];
pseudo_defCSV1["theta"] = ["Theta1PD.csv"];
pseudo_defCSV1["MainDir"] = [];
pseudo_defCSV1["plot"] = false;
pseudo_defCSV1["flag"] = "testPDCSV";
pseudo_defCSV1["Obs"] = ["B*2"];
pseudo_defCSV1["Noise"] = [0.1];

pseudo_defCSV2 = Dict();
pseudo_defCSV2["Nexp"] = 1;
pseudo_defCSV2["finalTime"] = [40];
pseudo_defCSV2["switchT"] = [[0.0,20.0,40.0]];
pseudo_defCSV2["y0"] = [[0,0]];
pseudo_defCSV2["preInd"] = [[0.1]];
pseudo_defCSV2["uInd"] = [reshape([1, 0.5], 2,1)];
pseudo_defCSV2["theta"] = [0.1 0.2 0.2 0.002];
pseudo_defCSV2["tsamps"] = [[0,5,10,15,20,25,30,35,40]];
pseudo_defCSV2["plot"] = false;
pseudo_defCSV2["flag"] = "testPDCSV";
pseudo_defCSV2["Obs"] = ["B*2"];
pseudo_defCSV2["Noise"] = [0.1];

# -------------------------------------------------------- MLE TESTS
mle_def = Dict();
mle_def["Nexp"] = [];
mle_def["finalTime"] = [];
mle_def["switchT"] = [];
mle_def["y0"] = [];
mle_def["preInd"] = [];
mle_def["uInd"] = [];
mle_def["tsamps"] = [];
mle_def["plot"] = [];
mle_def["flag"] = [];
mle_def["thetaMAX"] = [];
mle_def["thetaMIN"] = [];
mle_def["runs"] = [];
mle_def["parallel"] = [];
mle_def["DataMean"] = [];
mle_def["DataError"] = [];
mle_def["Obs"] = [];
mle_def["OPTsolver"] = [];
mle_def["MaxTime"] = [];
mle_def["MaxFuncEvals"] = [];

mle_def2 = Dict();
mle_def2["Nexp"] = [1];
mle_def2["finalTime"] = [40];
mle_def2["switchT"] = [[0,20,40]];
mle_def2["y0"] = [[0,0]];
mle_def2["preInd"] = [[0.1]];
mle_def2["uInd"] = [[1,1]];
mle_def2["tsamps"] = [[0,5,10,15,20,25,30,35,40]];
mle_def2["plot"] = false;
mle_def2["flag"] = "testmle";
mle_def2["thetaMAX"] = [0.15 0.25 0.25 0.025];
mle_def2["thetaMIN"] = [0.1 0.2 0.2 0.02];
mle_def2["runs"] = 1;
mle_def2["parallel"] = false;
mle_def2["DataMean"] = [[0 1 2 3 4 5 6 7 8; 0 1 2 3 4 5 6 7 8]'];
mle_def2["DataError"] = [[[0,1,2,3,4,5,6,7,8], [0,1,2,3,4,5,6,7,8]]];
mle_def2["Obs"] = ["A", "B"];
mle_def2["OPTsolver"] = "adaptive_de_rand_1_bin_radiuslimited";
mle_def2["MaxTime"] = 4;
mle_def2["MaxFuncEvals"] = [];

cvmle_def = Dict();
cvmle_def["Nexp"] = [];
cvmle_def["finalTime"] = [];
cvmle_def["switchT"] = [];
cvmle_def["y0"] = [];
cvmle_def["preInd"] = [];
cvmle_def["uInd"] = [];
cvmle_def["theta"] = [];
cvmle_def["tsamps"] = [];
cvmle_def["plot"] = [];
cvmle_def["flag"] = [];
cvmle_def["DataMean"] = [];
cvmle_def["DataError"] = [];
cvmle_def["Obs"] = [];

cvmle_def2 = Dict();
cvmle_def2["Nexp"] = [1];
cvmle_def2["finalTime"] = [40];
cvmle_def2["switchT"] = [[0,20,40]];
cvmle_def2["y0"] = [[0,0]];
cvmle_def2["preInd"] = [[0.1]];
cvmle_def2["uInd"] = [[1,1]];
cvmle_def2["theta"] = [0.1 0.2 0.2 0.02];
cvmle_def2["tsamps"] = [[0,5,10,15,20,25,30,35,40]];
cvmle_def2["plot"] = false;
cvmle_def2["flag"] = "cvmletest";
cvmle_def2["DataMean"] = [[0 1 2 3 4 5 6 7 8]'];
cvmle_def2["DataError"] = [[[0,1,2,3,4,5,6,7,8]]];;
cvmle_def2["Obs"] = ["B"];

testdict1 = Dict();
testdict2 = Dict();
testdict1["a"] = "a";
testdict1["b"] = "b";
testdict2["a"] = [];
testdict2["b"] = 1;


@testset "ModelGenTests" begin

    # In case any modification to the structure is made, this will be the test to see if I have forgot some entry
    @test model_def == defModStruct()

    # Check function does some modifications, so it is good to see that it does what it should
    @test model_def2 == checkStruct(model_def2)

    # Check that a file has been generated (the functions containing the ODEs)
    GenerateModel(model_def2)
    @test isfile(string(pwd(),"\\ModelsFunctions\\", model_def2["NameF"], "_Model.jl"))

    rm(string(pwd(),"\\ModelsFunctions\\", model_def2["NameF"], "_Model.jl"))
    # rm(string(pwd(),"\\ModelsFunctions"))

end

@testset "ModelSimulationTests" begin

    CSV.write("Theta1.csv", DataFrame([0.1 0.2 0.2 0.002]))
    CSV.write("TestSimulCSVInps.csv", DataFrame([0 40 0.1 1; 20 40 0.1 0.5]))
    CSV.write("TestSimulCSVObs.csv", DataFrame([0 0 0; 5 0 0; 10 0 0; 15 0 0; 20 0 0; 25 0 0; 30 0 0; 35 0 0; 40 0 0]))

    # In case any modification to the structure is made, this will be the test to see if I have forgot some entry
    @test simul_def == defSimulStruct()

    # Check function does some modifications, so it is good to see that it does what it should
    @test simul_def2 == checkStructSimul(model_def2, simul_def2)

    # In case any modification to the structure is made, this will be the test to see if I have forgot some entry
    @test simul_def3 == defSimulStructFiles()

    # Might need to add a test for this, but then I need some example CSV files
    @test simul_defCSV2 == extractSimulCSV(model_def2, simul_defCSV1)

    rm("Theta1.csv")
    rm("TestSimulCSVInps.csv")
    rm("TestSimulCSVObs.csv")

    GenerateModel(model_def2)
    simuls, ~, ~ = simulateODEs(model_def2, simul_def2);

    @test length(simuls) == simul_def2["Nexp"]
    for i in 1:2
        a,b,c = size(simuls[string("Exp_", i)])
        @test a == length(simul_def2["tsamps"][i])
        @test b == (model_def2["nStat"])
        @test c == length(simul_def2["theta"])/model_def2["nPar"]
        @test sum(simuls[string("Exp_", i)]) != 0
    end

    plotSimsODE(simuls,model_def2,simul_def2)
    for i in 1:2
        @test isfile(string(simul_def2["savepath"], "\\PlotSimulation_Exp", i,"_", simul_def2["flag"], ".png"))
        rm(string(simul_def2["savepath"], "\\PlotSimulation_Exp", i,"_", simul_def2["flag"], ".png"))
    end

    rm(string(pwd(),"\\ModelsFunctions\\", model_def2["NameF"], "_Model.jl"))
    rm(string(simul_def2["savepath"], "\\", simul_def2["savename"]))

end


@testset "PseudoDataGenerationTests" begin

    CSV.write("Theta1PD.csv", DataFrame([0.1 0.2 0.2 0.002]))
    CSV.write("TestSimulCSVInpsPD.csv", DataFrame([0 40 0.1 1; 20 40 0.1 0.5]))
    CSV.write("TestSimulCSVObsPD.csv", DataFrame([0 0 0; 5 0 0; 10 0 0; 15 0 0; 20 0 0; 25 0 0; 30 0 0; 35 0 0; 40 0 0]))


    # In case any modification to the structure is made, this will be the test to see if I have forgot some entry
    @test pseudo_def == defPseudoDatStruct()

    # Check function does some modifications, so it is good to see that it does what it should
    @test pseudo_def2 == checkStructPseudoDat(model_def2, pseudo_def2)

    # In case any modification to the structure is made, this will be the test to see if I have forgot some entry
    @test pseudo_def3 == defPseudoDatStructFiles()

    # Might need to add a test for this, but then I need some example CSV files
    @test pseudo_defCSV2 == extractPseudoDatCSV(model_def2, pseudo_defCSV1)

    rm("Theta1PD.csv")
    rm("TestSimulCSVInpsPD.csv")
    rm("TestSimulCSVObsPD.csv")

    GenerateModel(model_def2)
    pdat, ~, ~ = GenPseudoDat(model_def2, pseudo_def2);

    @test length(pdat["SimsObs"]) == pseudo_def2["Nexp"]
    @test length(pdat["Sims"]) == pseudo_def2["Nexp"]
    @test length(pdat["PData"]) == pseudo_def2["Nexp"]
    @test length(pdat["PError"]) == pseudo_def2["Nexp"]

    for i in 1:2
        a,b,c = size(pdat["Sims"][string("Exp_", i)])
        @test a == length(pseudo_def2["tsamps"][i])
        @test b == (model_def2["nStat"])
        @test c == length(pseudo_def2["theta"])/model_def2["nPar"]
        @test sum(pdat["Sims"][string("Exp_", i)]) != 0

        a,b,c = size(pdat["SimsObs"][string("PDExp_", i)])
        @test b == length(pseudo_def2["Obs"]);
        @test pdat["SimsObs"][string("PDExp_", i)][:,1,:] == pdat["Sims"][string("Exp_", i)][:,2,:].*2

        a,b,c = size(pdat["PData"][string("PDExp_", i)])
        @test b == length(pseudo_def2["Obs"]);

        @test length(pdat["PData"][string("PDExp_", i)]) == length(pdat["PError"][string("PDExp_", i)])
    end

    plotPseudoDatODE(pdat,model_def2,pseudo_def2)
    for i in 1:2
        @test isfile(string(pseudo_def2["savepath"], "\\PlotPseudoDat_Exp", i,"_", pseudo_def2["flag"], ".png"))
        rm(string(pseudo_def2["savepath"], "\\PlotPseudoDat_Exp", i,"_", pseudo_def2["flag"], ".png"))
    end

    rm(string(pwd(),"\\ModelsFunctions\\", model_def2["NameF"], "_Model.jl"))
    rm(string(pseudo_def2["savepath"], "\\", pseudo_def2["savename"]))

end

@testset "MLESeriesTests" begin

    # First check that both structure calls give what it is supposed to
    @test mle_def == defMLEStruct();
    @test cvmle_def == defCrossValMLEStruct();

    # Check functionality of helpping functions
    testdict3 = SimToMle(testdict1, testdict2);
    @test testdict3 == testdict2;

    simst = reshape([1 2 3 4; 1 2 3 4; 1 2 3 4], 3,4,1);
    obst = ["c*2"];
    stnamest = ["a","b","c","d"];
    @test sum(selectObsSim_te(simst, obst, stnamest)) == sum(simst[:,3,:]*2);
    @test (selectObsSim_te(simst, [2], stnamest)) == reshape([2,2,2], 3,1,1);

    ttess1 = Dict();
    ttess1["nInp"] = 2;
    ttess2 = Dict();
    ttess2["uInd"] = [[1 1 1; 2 2 2]'];
    @test restructInputs_te(ttess1, ttess2, 1) == [1,2,1,2,1,2];
    ttess1["nInp"] = 3;
    ttess2["uInd"] = [[1 1 1; 2 2 2; 3 3 3]'];
    @test restructInputs_te(ttess1, ttess2, 1) == [1,2,3,1,2,3,1,2,3];

    # More tests might be needed?
    m=10;
    s=1;
    d=9;
    @test (-0.5)*(log(2*pi) + 0 + (1)) == UVloglike(d, m, s);
    m=[10, 10, 10, 10];
    s=[1,1,1,1];
    d=[9, 11, 9, 11];
    @test ((-0.5)*(log(2*pi) + 0 + (1)))*4 == UVloglike(d, m, s);
    m=[10, 15, 10, 15];
    s=[0.5,0.5,0.5,0.5];
    d=[12, 12, 12, 12];
    @test ((-0.5)*(log(2*pi) + log(0.25) + (2^2)/0.25)) + ((-0.5)*(log(2*pi) + log(0.25) + (2^2)/0.25)) +
    ((-0.5)*(log(2*pi) + log(0.25) + (3^2)/0.25)) + ((-0.5)*(log(2*pi) + log(0.25) + (3^2)/0.25)) == UVloglike(d, m, s);


    m=[10, 10];
    s=[1 1; 1 1]; # Diagonal element of 0.1 is added to avoid Infs due to non positive definite covariances.
    d=[9, 11];
    correcmat = Diagonal(ones(length(d))).*0.1;
    @test [MVloglike(d, m, s)] == (-0.5) * ((length(d)*log(2*pi)) .+ log(det(s.+correcmat)) .+ ([-1 1]*inv(s.+correcmat)*[-1,1]));

    # Check structure check functions (with and without CSVs)
    @test mle_def2 == checkStructMLE(model_def2, mle_def2)
    @test cvmle_def2 == checkStructCrossValMLE(model_def2, cvmle_def2)

    # Check MLE
    # MLEtheta




    # Check CV for MLE results
    # CrossValMLE
    # finishMLEres




    # Test Plotting
    # plotMLEResults
    # plotCrossValMLEResults




end
