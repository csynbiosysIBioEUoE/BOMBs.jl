using BOMBS
using Test
using CSV
using DataFrames

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
