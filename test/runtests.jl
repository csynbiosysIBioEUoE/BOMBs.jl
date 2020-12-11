using BOMBS
using Test

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
model_def2["nStat"] = 1;
model_def2["nPar"] = 1;
model_def2["nInp"] = 1;
model_def2["stName"] = ["Tst"];
model_def2["parName"] = ["thta1"];
model_def2["inpName"] = ["inp1"];
model_def2["eqns"] = ["dTst = thta1*Tst-inp1"];
model_def2["Y0eqs"] = ["Tst = inp1/(thta1)"];
model_def2["Y0ON"] = false;
model_def2["tols"] = [1e-5,1e-5];
model_def2["solver"] = "CVODE_BDF";

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
simul_def2["Nexp"] = 1;
simul_def2["finalTime"] = [40];
simul_def2["switchT"] = [[0,20,40]];
simul_def2["y0"] = [[0]];
simul_def2["preInd"] = [[0]];
simul_def2["uInd"] = [[5,10]];
simul_def2["theta"] = [3.0];
simul_def2["tsamps"] = [[0,5,10,15,20,25,30,35,40]];
simul_def2["plot"] = false;
simul_def2["flag"] = "testsim";

simul_def3 = Dict()
simul_def3["ObservablesFile"] = [];
simul_def3["EventInputsFile"] = [];
simul_def3["theta"] = [];
simul_def3["MainDir"] = [];
simul_def3["plot"] = [];
simul_def3["flag"] = [];


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

    # In case any modification to the structure is made, this will be the test to see if I have forgot some entry
    @test simul_def == defSimulStruct()

    # Check function does some modifications, so it is good to see that it does what it should
    @test simul_def2 == checkStructSimul(model_def2, simul_def2)

    # In case any modification to the structure is made, this will be the test to see if I have forgot some entry
    @test simul_def3 == defSimulStructFiles()

    # Might need to add a test for this, but then I need some example CSV files
    # extractSimulCSV(model_def, simul_def)

    # Test generating a plot?
    # plotSimsODE(simuls,model_def,simul_def)

    # Add a test to check simulation results structure
    # simulateODEs(model_def, simul_def)

end
