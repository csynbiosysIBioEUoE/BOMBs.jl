using BOMBS
using Test



@testset "ModelGenTests" begin

    # In case any modification to the structure is made, this will be the test to see if I have forgot some entry
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

    @test model_def == defModStruct()

    # Check function does some modifications, so it is good to see that it does what it should
    model_def = defModStruct();
    model_def["NameF"] = ["Test"];
    model_def["nStat"] = 1;
    model_def["nPar"] = 1;
    model_def["nInp"] = 1;
    model_def["stName"] = ["Tst"];
    model_def["parName"] = ["thta1"];
    model_def["inpName"] = ["inp1"];
    model_def["eqns"] = ["dTst = thta1*Tst-inp1"];
    model_def["Y0eqs"] = ["Tst = inp1/(thta1)"];
    model_def["Y0ON"] = false;
    model_def["tols"] = [1e-5,1e-5];
    model_def["solver"] = "CVODE_BDF";

    @test model_def == checkStruct(model_def)

    # Check that a file has been generated (the functions containing the ODEs)
    GenerateModel(model_def)
    @test isfile(string(pwd(),"\\ModelsFunctions\\", model_def["NameF"], "_Model.jl"))

    rm(string(pwd(),"\\ModelsFunctions\\", model_def["NameF"], "_Model.jl"))
    # rm(string(pwd(),"\\ModelsFunctions"))

end

# @testset "ModelGenTests2" begin
#
#
#
# end
