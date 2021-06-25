using BOMBs
using Test
using CSV
using DataFrames
using LinearAlgebra
using Dates
using Distributions
using Random
using CmdStan

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
model_def["Y0Sim"] = [];
model_def["tols"] = [];
model_def["solver"] = [];

model_def2 = defModStruct();
model_def2["NameF"] = "Test";
model_def2["nStat"] = 2;
model_def2["nPar"] = 4;
model_def2["nInp"] = 1;
model_def2["stName"] = ["A", "B"];
model_def2["parName"] = ["k1", "k2", "k3", "k4"];
model_def2["inpName"] = ["inp1"];
model_def2["eqns"] = ["dA = k1*A - k2*A*inp1", "dB = k3*(A^2) - k4*B"];
model_def2["Y0eqs"] = ["A = 1", "B = k3*(A^2)/k4"];
model_def2["Y0Sim"] = false;
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
mle_def2["Nexp"] = 1;
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
mle_def2["runs"] = 2;
mle_def2["parallel"] = false;
mle_def2["DataMean"] = [[0 1 2 3 4 5 6 7 8; 0 1 2 3 4 5 6 7 8]'];
mle_def2["DataError"] = [[[0,1,2,3,4,5,6,7,8], [0,1,2,3,4,5,6,7,8]]];
mle_def2["Obs"] = ["A", "B"];
mle_def2["OPTsolver"] = "adaptive_de_rand_1_bin_radiuslimited";
mle_def2["MaxTime"] = [];
mle_def2["MaxFuncEvals"] = 10;

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
cvmle_def2["Nexp"] = 1;
cvmle_def2["finalTime"] = [40];
cvmle_def2["switchT"] = [[0,20,40]];
cvmle_def2["y0"] = [[0,0]];
cvmle_def2["preInd"] = [[0.1]];
cvmle_def2["uInd"] = [[1,1]];
cvmle_def2["theta"] = convert(Array, [0.1 0.2 0.2 0.02;0.12 0.22 0.22 0.022]');
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

    GenerateModel(model_def2);
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
    rm(string(simul_def2["savepath"]))

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

    for i in 1:2
        rm(string(pseudo_def2["savepath"], "\\PseudoDataFiles\\", model_def2["NameF"], "_EXP", i, "_", pseudo_def2["flag"], "_Events_Inputs.csv"))
        rm(string(pseudo_def2["savepath"], "\\PseudoDataFiles\\", model_def2["NameF"], "_EXP", i, "_", pseudo_def2["flag"], "_Observables.csv"))
        rm(string(pseudo_def2["savepath"], "\\PseudoDataFiles\\", model_def2["NameF"], "_EXP", i, "_", pseudo_def2["flag"], "_Simulations.csv"))
    end
    rm(string(pseudo_def2["savepath"], "\\PseudoDataFiles"))
    rm(string(pseudo_def2["savepath"]))
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
    GenerateModel(model_def2);
    mle_res, model_def3, mle_def3 = MLEtheta(model_def2, mle_def2);

    @test !isempty(mle_res);
    @test typeof(mle_res["StanDict"]) <: Array;
    @test typeof(mle_res["StanDict"][1]) <: Dict;
    @test typeof(mle_res["StanDict"][2]) <: Dict;
    @test symdiff(["k1","k2","k3","k4"],keys(mle_res["StanDict"][1])) == [];
    @test symdiff(["k1","k2","k3","k4"],keys(mle_res["StanDict"][2])) == [];

    @test length(mle_res["Theta"]) == 4*2;
    @test size(mle_res["Theta"]) == (4,2);

    @test length(mle_res["convCurv"]) == 2;

    @test (typeof(mle_res["convCurv"][1][1]) == Tuple{Int,Float64} || typeof(mle_res["convCurv"][1][1]) == Tuple{Int,Float32});
    @test (typeof(mle_res["convCurv"][1][2]) == Tuple{Int,Float64} || typeof(mle_res["convCurv"][1][2]) == Tuple{Int,Float32});
    @test (typeof(mle_res["convCurv"][2][1]) == Tuple{Int,Float64} || typeof(mle_res["convCurv"][2][1]) == Tuple{Int,Float32});
    @test (typeof(mle_res["convCurv"][2][2]) == Tuple{Int,Float64} || typeof(mle_res["convCurv"][2][2]) == Tuple{Int,Float32});

    @test length(mle_res["BestTheta"]) == 4;
    @test length(mle_res["BestCFV"]) == 2;

    #----------------------------- the parallel stuff!!!!
    mle_def2["parallel"] = true;
    mle_res2, model_def3, mle_def3 = MLEtheta(model_def2, mle_def2);
    @test isempty(mle_res2);
    @test isfile(string(mle_def3["savepath"], "\\MLEScripts\\", model_def3["NameF"], "_MLE.jl"));
    mle_res3, ~, ~ = finishMLEres(mle_res, model_def2, mle_def2);
    @test !isempty(mle_res3);
    mle_def2["parallel"] = false;


    # Check CV for MLE results
    cvmle_res, ~, ~ = CrossValMLE(model_def2, cvmle_def2);
    @test !isempty(cvmle_res);
    @test !isempty(cvmle_res["BestSimulations"][1]);
    @test size(cvmle_res["BestSimulations"][1])[2] == 2;
    @test length(cvmle_res["Costs"]) == 2;
    @test !isempty(cvmle_res["BestSimObservables"][1]);
    @test size(cvmle_res["BestSimObservables"][1])[2] == 1;
    @test length(cvmle_res["BestTheta"]) == 4;
    @test typeof(cvmle_res["SimObservables"]) <: Dict;
    @test size(cvmle_res["SimObservables"]["ExpObs_1"])[2] == 1;
    @test typeof(cvmle_res["Simulations"]) <: Dict;
    @test size(cvmle_res["Simulations"]["Exp_1"])[2] == 2;


    # Test Plotting
    plotMLEResults(mle_res,model_def2,mle_def2)
    @test isfile(string(mle_def2["savepath"], "\\PlotMLEResults_Exp", 1,"_", mle_def2["flag"], ".png"))
    @test isfile(string(mle_def2["savepath"], "\\Plot_MLEConvergence", "_", mle_def2["flag"], ".png"))
    rm(string(mle_def2["savepath"], "\\PlotMLEResults_Exp", 1,"_", mle_def2["flag"], ".png"))
    rm(string(mle_def2["savepath"], "\\Plot_MLEConvergence", "_", mle_def2["flag"], ".png"))

    rm(string(pwd(),"\\ModelsFunctions\\", model_def2["NameF"], "_Model.jl"))
    rm(string(mle_def2["savepath"], "\\", mle_def2["savename"]))

    rm(string(mle_def2["savepath"], "\\MLEScripts\\", model_def2["NameF"], "_MLE.jl"))
    rm(string(mle_def2["savepath"], "\\MLEScripts"))

    rm(string(cvmle_def2["savepath"], "\\", cvmle_def2["savename"]))
    rm(string(mle_def2["savepath"], "\\",model_def2["NameF"], "_", today(), "_SimulationResults_MLEsimulations1.jld"))
    rm(string(mle_def2["savepath"], "\\",model_def2["NameF"], "_", today(), "_SimulationResults_MLEsimulations2.jld"))
    rm(string(mle_def2["savepath"]))

end

@testset "StanInferenceTests" begin

    # Tests on structure definitions
    bayinf_def = defBayInfStruct();
    @test typeof(bayinf_def) <: Dict;
    entries1 = ["Priors", "Data", "StanSettings", "flag", "plot", "runInf", "MultiNormFit"];
    @test isempty(symdiff(entries1,keys(bayinf_def)));
    datinf_def = defBayInfDataStruct();
    @test typeof(datinf_def) <: Dict;
    entries2 = ["Nexp", "finalTime", "switchT", "y0", "preInd", "uInd", "tsamps", "Obs", "DataMean", "DataError"];
    @test isempty(symdiff(entries2,keys(datinf_def)));
    csvinf_def = defBayInfDataFromFilesStruct();
    @test typeof(csvinf_def) <: Dict;
    entries3 = ["Obs", "Observables", "Inputs", "y0"];
    @test isempty(symdiff(entries3,keys(csvinf_def)));
    stainf_def = defBasicStanSettingsStruct();
    @test typeof(stainf_def) <: Dict;
    entries4 = ["cmdstan_home", "nchains", "nsamples", "nwarmup", "printsummary", "init", "maxdepth", "adaptdelta", "jitter"];
    @test isempty(symdiff(entries4,keys(stainf_def)));

    # Check structures
    stainf_def["cmdstan_home"] = "C:/Users/David/.cmdstanpy/cmdstan-2.20.0";
    stainf_def["nchains"] = 2;
    stainf_def["nsamples"] = 20;
    stainf_def["nwarmup"] = 20;
    stainf_def["printsummary"] = false;
    stainf_def["init"] = [];
    stainf_def["maxdepth"] = 13;
    stainf_def["adaptdelta"] = 0.95;
    stainf_def["jitter"] = 0.5;
    @test stainf_def == checkStructBayInfStanSettings(model_def2, stainf_def);

    datinf_def["Nexp"] = 1;
    datinf_def["finalTime"] = [40];
    datinf_def["switchT"] = [[0,20,40]];
    datinf_def["y0"] = [[0,0]];
    datinf_def["preInd"] = [[0.1]];
    datinf_def["uInd"] = [[1,1]];
    datinf_def["tsamps"] = [[0,5,10,15,20,25,30,35,40]];
    datinf_def["DataMean"] = [[0 1 2 3 4 5 6 7 8; 0 1 2 3 4 5 6 7 8]'];
    datinf_def["DataError"] = [[[0,1,2,3,4,5,6,7,8], [0,1,2,3,4,5,6,7,8]]];
    datinf_def["Obs"] = ["A", "B"];
    @test datinf_def == checkStructBayInfData(model_def2, datinf_def);

    bayinf_def["Priors"] = [0.15 0.25 0.25 0.025; 0.1 0.2 0.2 0.02]'';
    bayinf_def["Data"] = datinf_def;
    bayinf_def["StanSettings"] = stainf_def;
    bayinf_def["flag"] = "stantest";
    bayinf_def["plot"] = false;
    bayinf_def["runInf"] = false;
    bayinf_def["MultiNormFit"] = false;
    bayinf_def2 =  checkStructBayInf(model_def2, bayinf_def);
    @test bayinf_def["Data"] == bayinf_def2["Data"]
    @test bayinf_def["StanSettings"] == bayinf_def2["StanSettings"]
    @test bayinf_def["flag"] == bayinf_def2["flag"]
    @test bayinf_def["plot"] == bayinf_def2["plot"]
    @test bayinf_def["runInf"] == bayinf_def2["runInf"]
    @test bayinf_def["MultiNormFit"] == bayinf_def2["MultiNormFit"]
    @test length(bayinf_def2["Priors"]) == 3;
    @test typeof(bayinf_def2["Priors"]) <: Dict
    # checkStructBayInfDataFiles(model_def, data_def)

    # Helpper functions tests
    s1 = [0,1,3,9,10]
    @test convertBoundTo2(s1, 10, 20) == s1.+10;
    @test convertBoundTo2(s1.+10, 0, 10) == s1;

    samp1 = rand(Normal(10,1), 80000,4);
    fipri1 = fitPriorSamps(samp1, model_def2);
    @test length(fipri1["pars"]) == 4;
    @test length(fipri1["transpars"]) == 5;
    @test length(fipri1["pridis"]) == 4;
    fipri2 = fitPriorSampsMultiNorm(samp1, model_def2);
    @test length(fipri2["pars"]) < 4;
    @test length(fipri2["transpars"]) == 5;
    @test length(fipri2["pridis"]) < 4;
    kee = ["transpars","mera","cora","pars","pridis","numN"]
    @test isempty(symdiff(kee,keys(fipri2)));

    dis1 = genStanInitDict(samp1', model_def2["parName"], 10);
    @test typeof(dis1[1]) <: Dict;
    @test isempty(symdiff(model_def2["parName"],keys(dis1[1])));
    @test length(dis1) == 10;

    # reparamDictStan(dis1, bayinf_def)

    # Stan
    genStanModel(model_def2, bayinf_def)
    @test isfile(string(pwd(),"\\ModelsFunctions\\", model_def2["NameF"], "_StanModel.stan"))
    sdat = restructureDataInference(model_def2, bayinf_def);
    datks = ["elm","tml","ts","tsl","tsmax","Nsp","inputs","evnT","m","stsl","stslm","sts","obser","obSta","nindu","preInd","Means","Y0us","Erros"]
    for i in 1:length(datks)
        @test datks[i] in keys(sdat)
        @test !isempty(sdat[datks[i]])
    end
    modelpath, Model, StanModel, inferdata, init, ~, ~ = getStanInferenceElements(model_def2, bayinf_def);
    @test isfile(modelpath);
    rm(modelpath)
    @test typeof(Model) == String;
    @test typeof(StanModel) == CmdStan.Stanmodel;
    @test inferdata == sdat;
    @test isempty(init);

    # saveStanResults(rc, chns, cnames, model_def, bayinf_def)
    # runStanInference(model_def, bayinf_def)

    stin, ~, ~ = StanInfer(model_def2, bayinf_def);
    kssm = ["init", "StanModel", "inferdata", "Model", "modelpath"]
    @test isempty(symdiff(kssm,keys(stin)));

    # Entropy tests
    ps = genSamplesPrior(model_def2, bayinf_def, 100);
    @test size(ps) == (100, 4);

    w = [1];
    E = [[1 0; 0 1]];
    MU = reshape([10, 10], 2, 1);
    x = reshape([10, 10], 2, 1);
    @test round(BOMBs.H_Upper(w,E)) == 3;
    @test round(BOMBs.mvGauss(x, MU, E[1])[1], digits=2) == 0.16;
    @test round(BOMBs.H_Lower(w, E, MU'), digits = 1) == 2.5;
    @test round(BOMBs.GaussMix(x[:,1], convert(Array, MU'), E, w), digits=2) == 0.16;
    @test round(BOMBs.ZOTSE(MU', E, w)) == 2;
    # Due to the use of global variables I cannot test the other functions... But I can test that the main script runs.
    # Do not know why, but I cannot make ScikitLearn work in the test file (works in scripts and jupyter)??? Need to have a look at it.
    # computeH(ps, model_def2, "test")




    # Plots
    # plotStanResults(staninf_res, model_def, bayinf_def)

    rm(string(pwd(),"\\ModelsFunctions\\", model_def2["NameF"], "_StanModel.stan"));
    rm(string(pwd(),"\\tmp\\", model_def2["NameF"], "_Stan_",bayinf_def["flag"],".stan"));
end

@testset "OEDModelSelectionTests" begin

    # structure definition
    oedms_def = defODEModelSelectStruct();
    @test typeof(oedms_def) <: Dict;
    entries1 = ["Model_1", "Model_2", "Obs", "Theta_M1", "Theta_M2", "y0_M1", "y0_M2", "preInd_M1", "preInd_M2",
                "finalTime", "switchT", "tsamps", "equalStep",
                "fixedInp", "fixedStep", "plot", "flag", "uUpper", "uLower", "maxiter"];
    @test isempty(symdiff(entries1,keys(oedms_def)));

    #structure check
    model_def3 = defModStruct();
    model_def3["NameF"] = "Test2";
    model_def3["nStat"] = 2;
    model_def3["nPar"] = 4;
    model_def3["nInp"] = 1;
    model_def3["stName"] = ["A", "B"];
    model_def3["parName"] = ["k1", "k2", "k3", "k4"];
    model_def3["inpName"] = ["inp1"];
    model_def3["eqns"] = ["dA = k1*A - k2*A*inp1", "dB = k3*(A^2) - k4*B"];
    model_def3["Y0eqs"] = ["A = 1", "B = k3*(A^2)/k4"];
    model_def3["Y0Sim"] = false;
    model_def3["tols"] = [1e-5,1e-5];
    model_def3["solver"] = "CVODE_BDF";

    oedms_def = Dict()
    oedms_def["Model_1"] = model_def2;
    oedms_def["Model_2"] = model_def3;
    oedms_def["Obs"] = ["B"];
    oedms_def["Theta_M1"] = [0.1 0.2 0.2 0.02; 0.11 0.21 0.21 0.021; 0.12 0.22 0.22 0.022];
    oedms_def["Theta_M2"] = [0.13 0.23 0.23 0.023; 0.14 0.24 0.24 0.024; 0.15 0.25 0.25 0.025];
    oedms_def["y0_M1"] = [0,0];
    oedms_def["y0_M2"] = [0,0];
    oedms_def["preInd_M1"] = [0.1];
    oedms_def["preInd_M2"] = [0.1];
    oedms_def["finalTime"] = [40];
    oedms_def["switchT"] = [0,20,30,40];
    oedms_def["tsamps"] = [0,5,10,15,20,25,30,35,40];
    oedms_def["fixedInp"] = [];
    oedms_def["fixedStep"] = [(1,[0])];
    oedms_def["equalStep"] = [[2,3]];
    oedms_def["plot"] = false;
    oedms_def["flag"] = "testoedms";
    oedms_def["uUpper"] = [1];
    oedms_def["uLower"] = [0];
    oedms_def["maxiter"] = 5;
    @test oedms_def == checkStructOEDMS(oedms_def);

    # Distance functions
    m1 = [10, 10];
    s1 = [1 0; 0 1];
    m2 = [11, 11];
    s2 = [1 0; 0 1];
    m3 = [15, 15];
    s3 = [1 0; 0 1];
    @test BhattacharyyaDist(m1, m1, s1, s1) == 0;
    @test BhattacharyyaDist(m1, m2, s1, s2) == 0.25;
    @test BhattacharyyaDist(m1, m3, s1, s3) == 6.25;

    @test EuclideanDist(m1, m1) == 0;
    @test EuclideanDist(m1, m2) == sqrt(2);
    @test EuclideanDist(m1, m3) == sqrt(50);

    # Generate utility script
    oedms_def = genOptimMSFuncts(oedms_def);
    @test isfile(string(pwd(), "\\ModelsFunctions\\", model_def2["NameF"], "_Model.jl"));
    @test isfile(string(pwd(), "\\ModelsFunctions\\", model_def3["NameF"], "_Model.jl"));
    @test isfile(string(pwd(), "\\Results\\", model_def2["NameF"], "_VS_", model_def3["NameF"], "_", today(),
                        "\\OEDModelSelectionScripts\\", model_def2["NameF"], "_VS_", model_def3["NameF"], "_OEDMS.jl"));

    # Bayes settings
    # Need to figure out if there is any test I can do here
    # opt = settingsBayesOpt(oedms_def);

    # Main function
    oedms_res, oedms_def = mainOEDMS(oedms_def);
    @test !isempty(oedms_res);
    @test typeof(oedms_res) <: Dict;
    @test isfile(string(pwd(), "\\Results\\", model_def2["NameF"], "_VS_", model_def3["NameF"], "_", today(),"\\OEDModelSelectResults_", oedms_def["flag"], ".jld"));
    @test oedms_res["uInpOpt"]["inp1"][1] == 0;
    @test oedms_res["uInpOpt"]["inp1"][2] == oedms_res["uInpOpt"]["inp1"][3];
    @test length(oedms_res["BestUtil"]) == 1;
    @test size(oedms_res["ConvCurv"])[1] == 5;
    @test size(oedms_res["Simul_M1"]) == size(oedms_res["Simul_M2"]);
    @test size(oedms_res["SimulObs_M1"])[2] == size(oedms_res["Simul_M1"])[2]-1;
    @test size(oedms_res["SimulObs_M2"])[2] == size(oedms_res["Simul_M2"])[2]-1;


    # Plot results
    plotOEDMSResults(oedms_res, oedms_def);
    @test isfile(string(oedms_def["savepath"], "\\Plot_OEDMSConvergence_", oedms_def["flag"], ".png"));
    rm(string(oedms_def["savepath"], "\\Plot_OEDMSConvergence_", oedms_def["flag"], ".png"));
    @test isfile(string(oedms_def["savepath"], "\\PlotOEDMSResults_Exp1_", oedms_def["flag"], ".png"));
    rm(string(oedms_def["savepath"], "\\PlotOEDMSResults_Exp1_", oedms_def["flag"], ".png"));

    rm(string(pwd(), "\\Results\\", model_def2["NameF"], "_VS_", model_def3["NameF"], "_", today(),"\\OEDModelSelectResults_", oedms_def["flag"], ".jld"));
    rm(string(pwd(), "\\ModelsFunctions\\", model_def2["NameF"], "_Model.jl"))
    rm(string(pwd(), "\\ModelsFunctions\\", model_def3["NameF"], "_Model.jl"))
    rm(string(pwd(), "\\Results\\", model_def2["NameF"], "_VS_", model_def3["NameF"], "_", today(),
                        "\\OEDModelSelectionScripts\\", model_def2["NameF"], "_VS_", model_def3["NameF"], "_OEDMS.jl"));
    rm(string(pwd(), "\\Results\\", model_def2["NameF"], "_VS_", model_def3["NameF"], "_", today(),
                        "\\OEDModelSelectionScripts"));
    rm(string(pwd(), "\\Results\\", model_def2["NameF"], "_VS_", model_def3["NameF"], "_", today()));

end

@testset "OEDModelCalibrationTests" begin

    # structure definition
    oedmc_def = defODEModelCalibrStruct();
    @test typeof(oedmc_def) <: Dict;
    entries1 = ["Model", "Obs", "Theta", "y0", "preInd", "finalTime", "switchT", "tsamps", "equalStep",
                "fixedInp", "fixedStep", "plot", "flag", "uUpper", "uLower", "maxiter", "util"];
    @test isempty(symdiff(entries1,keys(oedmc_def)));

    oedmc_def = Dict()
    oedmc_def["Model"] = model_def2;
    oedmc_def["Obs"] = ["B"];
    oedmc_def["Theta"] = [0.1 0.2 0.2 0.02; 0.11 0.21 0.21 0.021; 0.12 0.22 0.22 0.022];
    oedmc_def["y0"] = [0,0];
    oedmc_def["preInd"] = [0.1];
    oedmc_def["finalTime"] = [40];
    oedmc_def["switchT"] = [0,20,30,40];
    oedmc_def["tsamps"] = [0,5,10,15,20,25,30,35,40];
    oedmc_def["fixedInp"] = [];
    oedmc_def["fixedStep"] = [(1,[0])];
    oedmc_def["equalStep"] = [[2,3]];
    oedmc_def["plot"] = false;
    oedmc_def["flag"] = "testoedmc";
    oedmc_def["uUpper"] = [1];
    oedmc_def["uLower"] = [0];
    oedmc_def["maxiter"] = 5;
    oedmc_def["util"] = "perc";
    @test oedmc_def == checkStructOEDMC(oedmc_def);

    # Generate utility script
    oedmc_def = genOptimMCFuncts(oedmc_def);
    @test isfile(string(pwd(), "\\ModelsFunctions\\", model_def2["NameF"], "_Model.jl"));
    @test isfile(string(pwd(), "\\Results\\", model_def2["NameF"], "_", today(),
                        "\\OEDModelCalibrationScripts\\", model_def2["NameF"], "_OEDMC.jl"));

    # Main function
    oedmc_res, oedmc_def = mainOEDMC(oedmc_def);
    @test !isempty(oedmc_res);
    @test typeof(oedmc_res) <: Dict;
    @test isfile(string(pwd(), "\\Results\\", model_def2["NameF"], "_", today(),"\\OEDModelCalibrationResults_", oedmc_def["flag"], ".jld"));
    @test oedmc_res["uInpOpt"]["inp1"][1] == 0;
    @test oedmc_res["uInpOpt"]["inp1"][2] == oedmc_res["uInpOpt"]["inp1"][3];
    @test length(oedmc_res["BestUtil"]) == 1;
    @test size(oedmc_res["ConvCurv"])[1] == 5;
    @test size(oedmc_res["SimulObs_MC"])[2] == size(oedmc_res["Simul_MC"])[2]-1;

    # Plot results
    plotOEDMCResults(oedmc_res, oedmc_def);
    @test isfile(string(oedmc_def["savepath"], "\\Plot_OEDMCConvergence_", oedmc_def["flag"], ".png"));
    rm(string(oedmc_def["savepath"], "\\Plot_OEDMCConvergence_", oedmc_def["flag"], ".png"));
    @test isfile(string(oedmc_def["savepath"], "\\PlotOEDMCResults_Exp1_", oedmc_def["flag"], ".png"));
    rm(string(oedmc_def["savepath"], "\\PlotOEDMCResults_Exp1_", oedmc_def["flag"], ".png"));


    rm(string(pwd(), "\\ModelsFunctions\\", model_def2["NameF"], "_Model.jl"));
    rm(string(pwd(), "\\Results\\", model_def2["NameF"], "_", today(),
                        "\\OEDModelCalibrationScripts\\", model_def2["NameF"], "_OEDMC.jl"));
    rm(string(pwd(), "\\Results\\", model_def2["NameF"], "_", today(),"\\OEDModelCalibrationResults_", oedmc_def["flag"], ".jld"));

    rm(string(pwd(), "\\Results\\", model_def2["NameF"],  "_", today(),
                        "\\OEDModelCalibrationScripts"));
    rm(string(pwd(), "\\Results\\", model_def2["NameF"], "_", today()));
end



# rm(string(pwd(), "\\ModelsFunctions"))
# rm(string(pwd(), "\\Results"))
# rm(string(pwd(), "\\tmp"))
