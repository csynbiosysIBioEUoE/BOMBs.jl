

## Function so the user does not need to type it every time
# --------------------------------------> For now, resolution of simulation is every 1. Need to add modification in case user wants to go lower
function defModStruct()
    model_def = Dict();

    model_def["NameF"] = []; # Name for the model (scripts will be stored with this name) --> String
    model_def["nStat"] = []; # Nuber of states in the model --> Integer
    model_def["nPar"] = []; # Number of parameters of the mode --> Integer
    model_def["nInp"] = []; # Nuber of stimuly (inducer) in the model --> Integer
    model_def["stName"] = []; # Name for all the states --> Vector of strings
    model_def["parName"] = []; # Name for all the parameters --> Vector of strings
    model_def["inpName"] = []; # Name for all the inputs --> Vector of strings
    model_def["eqns"] = []; # Equations defining the model (left and right hand side) --> Vector of strings
    model_def["Y0eqs"] = []; # Steady state equations --> Vector of strings or empty vector if not used. If some element of the equation requires an experimental value for the calculation, please add exp at the beginning of the state(exemple: Cmrna -> expCmrna). Please do not name anything alp.
    model_def["Y0ON"] = []; # If the ON simulation is required for the steady state --> true, false or a string saying Yes/yes, No/no. If empty the default would be No.
    model_def["tols"] = []; # Relative and absolute tolerances (in that order) --> vector of 2 floats (can be left empty, where 1e-6 will be taken as default for both)
    model_def["solver"] = []; # IVP solver to solve the ODEs. If nothing specified, the default will be Tsit5(). For more info check https://diffeq.sciml.ai/v2.0/tutorials/ode_example.html

    return(model_def)
end


## Function with all the checks
function checkStruct(model_def)

    # Check taht all the dictionary entries are correct
    entries = ["stName","inpName","eqns","nPar","NameF","nInp","parName","nStat","Y0eqs","Y0ON","tols","solver"]
    if symdiff(entries,keys(model_def))!=[] && symdiff(entries,keys(model_def))!=["modelpath"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the entries of the dictionary, there is soemthign wrong...")
        println(symdiff(entries,keys(model_def)))
        return
    end

    # Check one to be sure there is no empty entry
    if model_def["NameF"] == []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, introduce a name for the model!")
        return
    elseif model_def["nStat"] == []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, introduce the number of states!")
        return
    elseif model_def["nPar"] == []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, introduce the number of paramters!")
        return
    elseif model_def["nInp"] == []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, introduce the number of stimuli!")
        return
    elseif model_def["stName"] == []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, introduce the state names!")
        return
    elseif model_def["parName"] == []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, introduce the parameter names!")
        return
    elseif model_def["inpName"] == []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, introduce the stimuli names!")
        return
    elseif model_def["eqns"] == []
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, introduce the equations!")
        return
    end

    if "alp" in model_def["parName"] || "alp" in model_def["inpName"] || "alp" in model_def["stName"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, but the name alp is reserved for something else in the model...")
        return
    end

    # Checks to see if the content of the fields is the expected (also consier if single value entries have been given as a number or a vector)
    if ((typeof(model_def["NameF"])!=Array{String,1}) && (typeof(model_def["NameF"])!=String))
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field NameF!")
        return
    elseif (typeof(model_def["nStat"][1])!=Int) && (typeof(model_def["nStat"])!=Int)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field nStat!")
        return
    elseif (typeof(model_def["nPar"][1])!=Int) && (typeof(model_def["nPar"])!=Int)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field nPar!")
        return
    elseif (typeof(model_def["nInp"][1])!=Int) && (typeof(model_def["nInp"])!=Int)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field nInp!")
        return
    elseif (typeof(model_def["stName"])!=Array{String,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field stName!")
        return
    elseif (typeof(model_def["parName"])!=Array{String,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field parName!")
        return
    elseif (typeof(model_def["inpName"])!=Array{String,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field inpName!")
        return
    elseif (typeof(model_def["eqns"])!=Array{String,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field eqns!")
        return
    elseif (typeof(model_def["Y0eqs"])!=Array{Any,1}) && (typeof(model_def["Y0eqs"])!=Array{String,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Y0eqs!")
        return
    elseif ((typeof(model_def["solver"])!=Array{String,1}) && (typeof(model_def["solver"])!=String)) && ((model_def["solver"])!=[])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field solver!")
        return
    elseif model_def["Y0ON"]!=true && model_def["Y0ON"]!="Yes" && model_def["Y0ON"]!="yes" &&
    model_def["Y0ON"]!= false && model_def["Y0ON"] != "No" && model_def["Y0ON"] != "no"&&
    model_def["Y0ON"]!=[true] && model_def["Y0ON"]!=["Yes"] && model_def["Y0ON"]!=["yes"] &&
    model_def["Y0ON"]!= [false] && model_def["Y0ON"] != ["No"] && model_def["Y0ON"] != ["no"]&&
    model_def["Y0ON"]!= [] && model_def["Y0ON"] != ""

        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field Y0ON!")
        return
    elseif (typeof(model_def["tols"])!=Array{Any,1}) && (typeof(model_def["tols"])!=Array{Float64,1}) && (typeof(model_def["tols"])!=Array{Float32,1})
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field tols!")
        return
    elseif ((typeof(model_def["tols"])==Array{Float64,1}) || (typeof(model_def["tols"])==Array{Float32,1})) && (length(model_def["tols"]) != 2)
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the field tols!")
        return
    end


    # Take single elements outside vector (if so) to generalise things after
    if typeof(model_def["NameF"]) == Array{String,1}
        model_def["NameF"] = model_def["NameF"][1];
    end

    if typeof(model_def["nStat"]) == Array{Int,1}
        model_def["nStat"] = model_def["nStat"][1];
    end

    if typeof(model_def["nPar"]) == Array{Int,1}
        model_def["nPar"] = model_def["nPar"][1];
    end

    if typeof(model_def["nInp"]) == Array{Int,1}
        model_def["nInp"] = model_def["nInp"][1];
    end

    if model_def["Y0ON"]==true || model_def["Y0ON"]=="Yes" || model_def["Y0ON"]=="yes" ||
    model_def["Y0ON"]==[true] || model_def["Y0ON"]==["Yes"] || model_def["Y0ON"]==["yes"]
        model_def["Y0ON"]=true;
    elseif model_def["Y0ON"]== false || model_def["Y0ON"] == "No" || model_def["Y0ON"] == "no"||
    model_def["Y0ON"]== [false] || model_def["Y0ON"] == ["No"] || model_def["Y0ON"] == ["no"]||
    model_def["Y0ON"]== [] || model_def["Y0ON"] == ""
        model_def["Y0ON"]=false;
    end

    if model_def["tols"] == []
        model_def["tols"] = [1e-6, 1e-6];
    end

    if model_def["solver"] == []
        model_def["solver"] = "Tsit5";
    elseif model_def["solver"] == Array{String,1}
        model_def["solver"] = model_def["solver"][1];
    end

    # Check that the user has introduced a correct solver
    if (model_def["solver"] != "AutoTsit5(Rosenbrock23())") && (model_def["solver"] != "AutoVern7(Rodas5())") &&
        (model_def["solver"] != "Tsit5") && (model_def["solver"] != "BS3") && (model_def["solver"] != "Vern7") &&
        (model_def["solver"] != "Rodas4") && (model_def["solver"] != "Rodas5") && (model_def["solver"] != "KenCarp4") &&
        (model_def["solver"] != "TRBDF2") && (model_def["solver"] != "RadauIIA") && (model_def["solver"] != "CVODE_BDF")
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Sorry, but you have introduced a wrong solver!")
        println("The options are AutoTsit5(Rosenbrock23()), AutoVern7(Rodas5()), Tsit5, BS3, Vern7, Rodas4, Rodas5, ")
        println("          KenCarp4, TRBDF2, RadauIIA, CVODE_BDF.")
        println("For more info, please check https://diffeq.sciml.ai/stable/tutorials/ode_example/#ode_example-1")
        return
    end


    # Check tha there are no spaces in determined fields
    if findfirst(" ", join(model_def["stName"])) != nothing
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, do NOT include spaces in the entries of the field stName!")
        return
    end

    if findfirst(" ", join(model_def["inpName"])) != nothing
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, do NOT include spaces in the entries of the field inpName!")
        return
    end

    if findfirst(" ", join(model_def["NameF"])) != nothing
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, do NOT include spaces in the entries of the field NameF!")
        return
    end

    if findfirst(" ", join(model_def["parName"])) != nothing
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, do NOT include spaces in the entries of the field parName!")
        return
    end


    # Check that contents make sense
    if length(model_def["stName"]) != model_def["nStat"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, chec stName and nStat, they do not make sense")
        return
    end
    if length(model_def["inpName"]) != model_def["nInp"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, chec inpName and nInp, they do not make sense")
        return
    end

    if length(model_def["parName"]) != model_def["nPar"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, chec parName and nPar, they do not make sense")
        return
    end

#     if length(model_def["eqns"]) != model_def["nStat"]
#         println("-------------------------- Process STOPPED!!! --------------------------")
#         println("Please, chec eqns and nStat, the number of equations is not the same as the number of states")
#         return
#     end

    # ----------------------------------------> Might need to double check this part
    # Check if in the ODEs all the states are in the equations and with a d in front
    # First check that the equations contain an = sign and only 1

    # Check that in the ODEs there is only 1 = sign
    if (findfirst.("=", model_def["eqns"])!=findlast.("=", model_def["eqns"]))
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the number of = signs in your OED equations! You might have put 2 equations in one same string!")
        return
    end
    if model_def["Y0eqs"] != [] # Avoind the check if no Steady state equations are given
        # Check that in the Steady State Equations there is only 1 equal sign
        if (findfirst.("=", model_def["Y0eqs"])!=findlast.("=", model_def["Y0eqs"]))
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check the number of = signs in your Steady State equations! You might have put 2 equations in one same string!")
            return
        end
        # Check that in the Steady state equations all entries have left and right hand side
        if sum(occursin.("=", model_def["Y0eqs"])) != length(model_def["Y0eqs"])
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check the equations! You might have forgoten a left hand side or add 2 equations in one same string in the Steady State Equation.")
            return
        end
        # Check that in the Steady State equations all the states are given as equations
        ds2 = [model_def["Y0eqs"][i][1:findfirst.("=", model_def["Y0eqs"][i])[1]] for i in 1:model_def["nStat"]]
        if sum(occursin.(model_def["stName"], join(ds2))) !=  model_def["nStat"]
            println("-------------------------- Process STOPPED!!! --------------------------")
            println("Please, check the equations! You migh have forgoten a Steady State Equation!.")
            return
        end
    end

    # Check that in the ODEs all entries have left and right hand side
    if sum(occursin.("=", model_def["eqns"])) != length(model_def["eqns"])
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the equations! You might have forgoten a left hand side of the ODEs or add 2 equations in one same string.")
        return
    end

    # Check that in the ODEs all the states are given as equations
    ds = [model_def["eqns"][i][1:findfirst.("=", model_def["eqns"][i])[1]] for i in 1:length(model_def["eqns"])];
    if sum(occursin.(string.("d", model_def["stName"]), join(ds))) != model_def["nStat"]
        println("-------------------------- Process STOPPED!!! --------------------------")
        println("Please, check the equations! You migh have forgoten an ODE!.")
        println("Remember, ODEs should start with a d and no other equation can.")
        return
    end


    return(model_def)
end


## Function to generate all the necesary Julia scripts for the model simulations

function GenerateModel(model_def)

    model_def = checkStruct(model_def);

    # Packages needed for the simulations
    Head = string("
    using DifferentialEquations
    using OrdinaryDiffEq
    using DiffEqBase
    using Sundials
    using ODEInterfaceDiffEq
        ");

    # Generate ODEs function
    # String containing the equtions, ODEs and others
    tet = Array{Any,1}(undef, length(model_def["eqns"])); # Vector of strings that will contain the eauqtions for the function
    odec = 1; # Counter for the ODEs to take into account equations that are not ODEs
    for i in 1:length(model_def["eqns"]) # Loop over number of equations
        if replace(model_def["eqns"][i], " "=>"")[1]=='d' # If the equation starts with a d (differential)
            tet[i] = string("        du[", odec, "] = ", model_def["eqns"][i],";\n");
            odec +=1;
            else # Other equations
            tet[i] = string("        ", model_def["eqns"][i],";\n");
        end
    end;

    fun1 = string("
        # du -> Derivative
        # u -> State at time t
        # p-> paramter vector
        # t-> time as tuple (init, end)

        function ",join(model_def["NameF"]),"ODE!(du,u,p,t)
            ",join(model_def["stName"], ", ")," = u;

            ",join(model_def["parName"], ", "),", ",join(model_def["inpName"], ", ")," = p;

    ",join(tet),"

        end
        ");

    # Generate Steady State function
    y0eq = ""; # Vector that will contain the steady state equations (SSE)
    if model_def["Y0eqs"]==[] # If there is no equations, make the function return the same Y0 that the user gives
        y0eq = (string(join(model_def["stName"], ", ")," = I; \n        ",


                "       return(I) \n "));
    else
        tet = Array{Any,1}(undef, length(model_def["Y0eqs"])); # Vector of strings that will contain the eauqtions for the function
        odec = 1 # Counter for the SSEs to take into account equations that are not SSEs
        df = [model_def["Y0eqs"][i][1:findfirst.("=", model_def["Y0eqs"][i])[1]] for i in 1:length(model_def["Y0eqs"])]; # Check for the equations that start with the name of the states (to tell apart from other equations)
        for i in 1:length(model_def["Y0eqs"]) # Loop over equations
            if sum(occursin.(model_def["stName"], df[i]))!=0 # If the equation looked begins with the name of a state
                tet[i] = string("        alp[", odec, "] = ", model_def["Y0eqs"][i],";\n");
                odec +=1;
            else
                tet[i] = string("        ", model_def["Y0eqs"][i],";\n");
            end
        end

        y0eq = (string(join([string("exp",model_def["stName"][k]) for k in 1:model_def["nStat"]], ", ")," = I; \n        ",

        join(model_def["parName"], ", "),", ",join(model_def["inpName"], ", ")," = p;\n \n \n ",

                "        alp = zeros(4); \n \n",

                join(tet),
                "        return(alp) \n "));

    end

    fun2 = string("
        # p -> parameter vector (plus inducer vector at the end)
        # I -> initial Y0 vector

        function ",join(model_def["NameF"]),"SteadyState(p,I)

            ",y0eq,"

        end
        ");

    # Generate Step Wise simulation function
    ONSim = "";
    if model_def["Y0ON"]==false
        ONSim = "y0 = y_al;";
    elseif model_def["Y0ON"]==true
        if model_def["solver"] == "CVODE_BDF"
            ONSim = string("
                prob = ODEProblem(",join(model_def["NameF"]),"ODE!,y_al,(0.0,Float64(24*60-1)),pSS);
                ssv = Sundials.solve(prob, CVODE_BDF(),reltol=",model_def["tols"][1],",abstol=",model_def["tols"][2],");
                y0 = ssv[:,end]; ");
        else
            ONSim = string("
                prob = ODEProblem(",join(model_def["NameF"]),"ODE!,y_al,(0.0,Float64(24*60-1)),pSS);
                ssv = solve(prob, ",join(model_def["solver"]),"(),reltol=",model_def["tols"][1],",abstol=",model_def["tols"][2],");
                y0 = ssv[:,end]; ");
        end


    end

    if model_def["solver"] == "CVOED_BDF"
        soso = "Sundials.";
    else
        soso = "";
    end

    fun3 = string("
        # ts -> time vector going from t=0 to t=end every 1
        # p -> parameter vector
        # sp -> vector with switching times (t=0 and t=end needs to be included). Length of vector will be length(inducer steps)+1
        # inputs -> Vecor with the inducers (if more than 1 inducer in the system, the shape will be a1,b1,c1,..., a2,b2,c2,... aend,bend,cend,... where each leter represents a different inducer)
        # ivss -> Initial Y0 vector
        # pre -> Vector with the inducers in the ON. It can be empty

        function ",join(model_def["NameF"]),"_solvecoupledODE(ts, p, sp, inputs, ivss, pre=[])

            maxtime = length(ts);
            Nsp = length(sp);
            Nevents = length(sp)-1;
            Neq = ",join(model_def["nStat"]),";
            if pre != []
                pSS = vcat(p, pre);
            else
                pSS = p;
            end

            final = zeros(maxtime, Neq);

            y_al = ",join(model_def["NameF"]),"SteadyState(pSS,ivss) # Calculation of initial guesses for steady state

            ",join(ONSim),"
            initialV = y0;
            i = 1;

            for q in collect(1:Nevents)
                lts = length(ts[(sp[q]+1):sp[q+1]+1]);  # General way to define the number of elements in each event series
                Tevent = ts[(sp[q]+1):sp[q+1]+1];  # General way to extract the times of each event
                I = inputs[i:(i+(",model_def["nInp"]-1,"))];
                pSte = vcat(p, I);

                if q == 1
                    prob = ODEProblem(",join(model_def["NameF"]),"ODE!,initialV,(ts[(sp[q]+1)],ts[sp[q+1]+1]),pSte);
                    part1 = ",join(soso),"solve(prob, ",join(model_def["solver"]),"(),reltol=",model_def["tols"][1],",abstol=",model_def["tols"][2],",saveat=1);
                else
                    prob = ODEProblem(",join(model_def["NameF"]),"ODE!,initialV,(ts[(sp[q]+1)],ts[sp[q+1]+1]),pSte);
                    part1 = ",join(soso),"solve(prob, ",join(model_def["solver"]),"(),reltol=",model_def["tols"][1],",abstol=",model_def["tols"][2],",saveat=1);
                end

                initialV = part1[lts];
                i+=",model_def["nInp"],";

                for d in collect((sp[q]+1):(sp[q]+lts))
                    final[d,:] = part1[d-sp[q]];
                end
            end
            return(final)

        end
        ");

    # Function to simulate the system as many times as parameter draws are given
    fun4 = string("
        # ts -> time vector going from t=0 to t=end every 1
        # pD -> Matrix containing all the draws for the parameters (can be a single vector if nly 1 draw of parameters is considered)
        # sp -> vector with switching times (t=0 and t=end needs to be included). Length of vector will be length(inducer steps)+1
        # inputs -> Vecor with the inducers (if more than 1 inducer in the system, the shape will be a1,b1,c1,..., a2,b2,c2,... aend,bend,cend,... where each leter represents a different inducer)
        # ivss -> Initial Y0 vector
        # samps -> Vector of Time points that want to be extracted from the simulations (it considers the initial time point as 0)
        # pre -> Vector with the inducers in the ON. It can be empty

        function ",join(model_def["NameF"]),"_SolveAll(ts, pD, sp, inputs, ivss, samps, pre=[])

            if length(size(pD)) == 1
                pD = reshape(pD,size(pD)[1],1);
            end

            if size(pD)[2] != ",model_def["nPar"],"
                pD = pD';
            end

            if length(ivss)/",model_def["nStat"]," > 1
                if size(ivss)[2] != ",model_def["nStat"],"
                    ivss = ivss';
                end
            end

            AllSolTest = zeros(length(samps), ",model_def["nStat"],", length(pD[:,1]));

            for drawInd in collect(1:length(pD[:,1]))
                p = pD[drawInd,:];
                if length(ivss)/",model_def["nStat"]," > 1
                    ivss2 = ivss[drawInd,:];
                else
                    ivss2 = ivss;
                end
                temp = ",join(model_def["NameF"]),"_solvecoupledODE(ts, p, sp, inputs, ivss2, pre);
                AllSolTest[:,:,drawInd] = temp[convert.(Int,samps.+1),:];
            end

            return(AllSolTest);
        end
        ");

    finscr = string(Head, fun1, fun2, fun3, fun4);

    fidi = pwd(); #@__DIR__;

    if !isdir(string(fidi,"\\ModelsFunctions"))
        mkdir(string(fidi,"\\ModelsFunctions"));
    end

    open(string(fidi,"\\ModelsFunctions\\", model_def["NameF"], "_Model.jl"), "w") do io
       write(io, finscr);
    end;

    include(string(fidi,"\\ModelsFunctions\\", model_def["NameF"], "_Model.jl"))

    model_def["modelpath"] = string(fidi,"\\ModelsFunctions\\", model_def["NameF"], "_Model.jl");

    println("")
    println("----------------------------------------- MODEL GENERATION -----------------------------------------")
    println("The model has been generated in the directory: ")
    println(string("                 ", model_def["modelpath"]))
    println("--------------------------------------------------------------------------------------")
    println("")

    return(model_def)

end
