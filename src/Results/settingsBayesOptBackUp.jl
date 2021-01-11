## Main function that runs the optimisation

function settingsBayesOpt(oedms_def)

    # Manage inputs order
    inpnam = oedms_def["Model_1"]["inpName"];
    if oedms_def["Model_1"]["nInp"] < oedms_def["Model_2"]["nInp"]
        tmp1 = [];
        try
            for i in 1:oedms_def["Model_2"]["nInp"]
                if oedms_def["Model_2"]["inpName"][i] in oedms_def["Model_1"]["inpName"]
                    nothing
                else
                    tmp1 = vcat(tmp1, oedms_def["Model_2"]["inpName"][i])
                end
            end
            inpnam = vcat(inpnam, tmp1)
        catch
            for i in 1:oedms_def["Model_2"]["nInp"]
                if oedms_def["Model_2"]["inpName"][i] in oedms_def["Model_1"]["inpName"]
                    nothing
                else
                    tmp1 = hcat(tmp1, oedms_def["Model_2"]["inpName"][i])
                end
            end
            inpnam = hcat(inpnam, tmp1)
        end
    end;

    # Optimisation section
    steps = length(oedms_def["switchT"])-1;
    induc = maximum([oedms_def["Model_1"]["nInp"], oedms_def["Model_2"]["nInp"]]);
    fiii = length(oedms_def["fixedInp"]);

    doof = (steps*(induc-fiii))+fiii;

    # THE THIRD MOUNTAINS!
    uppe = Array{Any,1}(undef, doof); # Vector of upper bounds. Length set with doof but deppending on the user settings from fixedStep and equalStep this can change.
    lowe = Array{Any,1}(undef, doof); # Vector of lower bounds. Same ^
    coun = 1; # Index for the entry of ustr that will advance only if certain conditions are meet (that that position has been already checked and a variable or number has been introduced)
                # Tecnichally this indicated the first entry of uppe and lowwe for each inducer variables.
    fs = [oedms_def["fixedStep"][j][1] for j in 1:length(oedms_def["fixedStep"])]; # Fixed Steps (the actual step number)
    us = [oedms_def["fixedStep"][j][2] for j in 1:length(oedms_def["fixedStep"])]; # Input values (or Any entry) for each step specified in fixedStep
    fixin = 0; # Index to account for fixedInp. If the fixedInp is the first one in the inpam list this is the best way (that I have found) to account for it. This is because of the use of i (input index). This is because first I consider the case where there is a fisedInp and then check the rest.

    for i in 1:length(inpnam) # Loop over the number of inducers
        if inpnam[i] in oedms_def["fixedInp"] # If the current inducer is fixed then add only one bound and move index (coun). Fixin works as described in the other MOUNTAINS
            uppe[coun] = oedms_def["uUpper"][i];
            lowe[coun] = oedms_def["uLower"][i];
            coun += 1;
            fixin += 1;
        else
            try # Using coun you take all the indexes corresponding to the current inducer. The you fill it with the specified uper and lower bounds. Added the try just in case someoen runns it in a machien that does not support Float 64? (I am not sure)
                uppe[coun:coun+(steps-1)] = repeat([convert(Float64, oedms_def["uUpper"][i])], outer = steps);
                lowe[coun:coun+(steps-1)] = repeat([convert(Float64, oedms_def["uLower"][i])], outer = steps);
            catch
                uppe[coun:coun+(steps-1)] = repeat([convert(Float32, oedms_def["uUpper"][i])], outer = steps);
                lowe[coun:coun+(steps-1)] = repeat([convert(Float32, oedms_def["uLower"][i])], outer = steps);
            end

            for j in 1:length(uppe[coun:coun+(steps-1)]) # Now we loop over each entry that we filled to check if we need to delete it according to the users fixedStep and equalStep
                if j in fs # Since this is not a fixed input j also indicates the step index, so check if this is in the list of fixedStep
                    for w in 1:length(fs) # As commentetd in other loops, this for and if is to consider the case the user does not introduce fixedStep entries sorted by step.
                        if j == fs[w]
                            if us[w][i-fixin] != Any # Check if current input contained in a fixedStep is defined as Any (to be optimised) or not (has been fixed. )
                                uppe[coun-1+j] = "";# Since it was fixed (different from Any) ad an empty string there to remove later in both bounds.
                                lowe[coun-1+j] = "";
                            end
                        end
                    end
                end
                if oedms_def["equalStep"] != [] # If equalStep is empty then we do not need to do anything else, if it has something then there are entries to remove from the bounds
                    for k in 1:length(oedms_def["equalStep"]) # Loop over each entry of equalStep
                        tmpk = zeros(length(oedms_def["equalStep"][k])); # Temporal variable that will contain the same values as the equalStep entry except the ones that have been fixed for that step (this will happen if the user set an Any for the current inducer)
                        for h in 1:length(oedms_def["equalStep"][k]) # Loop over each value in the current entry
                            if oedms_def["equalStep"][k][h] in fs # Check if the current step considered is also in fixedStep
                                for w in 1:length(fs) # Account for ht euser not introducng steps sorted
                                    if oedms_def["equalStep"][k][h] == fs[w]
                                        if us[w][i-fixin] == Any # If the current step considered was defined as Any for the current inducer add the index, if it was fixed to a value ignore (this has been taken into account before)
                                            tmpk[h] = oedms_def["equalStep"][k][h];
                                        end
                                    end
                                end
                            else
                                tmpk[h] = oedms_def["equalStep"][k][h] # If the current step was not in fixedStep just add it
                            end
                        end
                        tmpk = convert.(Int, filter!(e->e∉[0],tmpk)); # Remove 0 entries (fixed values)
                        if j in tmpk && j != tmpk[1] # Now, check if the current step is in our temporary index. If it is the first entry in the vector then do nothing, if not remove (it will be the same as the first, so no bound for it)
                            uppe[coun-1+j] = "";
                            lowe[coun-1+j] = "";
                        end
                    end
                end
            end
            coun += (steps); # Add the number of steps to the varable coun so we can move to the next inducer
        end
    end

    uppe = [uppe[i] for i in 1:length(uppe) if uppe[i] != ""]; # Remove empty strings so we are left with a vector of numbers with the bounds for the inducers.
    lowe = [lowe[i] for i in 1:length(lowe) if lowe[i] != ""];
    #################################### OPTIMISER SETINGS THAT CAN BE MODIFIED START HERE ####################################

    model = ElasticGPE(convert(Int, length(uppe)),                            #
           mean = MeanConst(0.),
           kernel = SEArd(zeros(length(uppe)), 5.),
           logNoise = 0.,
           capacity = 3000);

    modeloptimizer = MAPGPOptimizer(every = 50, noisebounds = [-4, 3],       # bounds of the logNoise
                   kernbounds = [[-1*ones(length(uppe)); 0], [4*ones(length(uppe)); 10]],  # bounds of the 3 parameters GaussianProcesses.get_param_names(model.kernel)
                   maxeval = 40);

    sut = Symbol(string(oedms_def["Model_1"]["NameF"], "_VS_", oedms_def["Model_2"]["NameF"],"Utility"));
    opt = BOpt((@eval $sut),
          model,
          ExpectedImprovement(),                   # type of acquisition
          modeloptimizer,
          lowe,
          uppe,                                     # lowerbounds, upperbounds
          repetitions = 1,                          # evaluate the function for each input 1 times
          maxiterations = oedms_def["maxiter"],     # evaluate at 50 input positions
          sense = Max,                              # maximise the function
          verbosity = Progress);

    #################################### OPTIMISER SETINGS THAT CAN BE MODIFIED STOP HERE ####################################

    return opt

end
