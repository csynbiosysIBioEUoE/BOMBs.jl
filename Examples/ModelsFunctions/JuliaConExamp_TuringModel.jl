
    using DifferentialEquations
    using OrdinaryDiffEq
    using DiffEqBase
    using Sundials
    using ODEInterfaceDiffEq
    using Turing
    using MCMCChains
        
        # du -> Derivative
        # u -> State at time t
        # p-> paramter vector
        # t-> time as tuple (init, end)

        function JuliaConExampODE!(du,u,p,t)
            mRNA, Prot = u;

            n1, KM2, d2, km, U = p;

        du[1] = dmRNA = 1 + (200 / (1 + ((Prot/20)*((U^2)/((KM2^2)+(U^2))))^n1)) - 0.5*mRNA;
        du[2] = dProt = km*mRNA - d2*Prot;


        end
        
        # p -> parameter vector (plus inducer vector at the end)
        # I -> initial Y0 vector

        function JuliaConExampSteadyState(p,I)

            mRNA, Prot = I; 
               return(I) 
 

        end
        
        # ts -> time vector going from t=0 to t=end every 1
        # p -> parameter vector
        # sp -> vector with switching times (t=0 and t=end needs to be included). Length of vector will be length(inducer steps)+1
        # inputs -> Vecor with the inducers (if more than 1 inducer in the system, the shape will be a1,b1,c1,..., a2,b2,c2,... aend,bend,cend,... where each leter represents a different inducer)
        # ivss -> Initial Y0 vector
        # pre -> Vector with the inducers in the ON. It can be empty

        @model function JuliaConExamp_TuringModel(data,prob1)

            p1 ~ Truncated(Normal(3.0, 1.0), 1.0, 5.0)
            p2 ~ Truncated(Normal(15.5, 7.25), 1.0, 30.0)
            p3 ~ Truncated(Normal(0.505, 0.2475), 0.01, 1.0)
            p4 ~ Truncated(Normal(3.0, 1.0), 1.0, 5.0)

            p = [p1,p2,p3,p4,];

            Neq = 2;

            finalExps = Array{Any,1}(undef, 2);
            finalBay2EXPS = Array{Any,1}(undef, 2);

            for exp in 1:2

                ts=collect(0.:1:data["finalTime"][exp])
                sp = data["switchT"][exp];
                inputs = data["uInd"][exp];
                ivss = data["y0"][exp];
                pre=data["preInd"][exp];

                maxtime = length(ts);
                Nsp = length(sp);
                Nevents = length(sp)-1;


                if typeof(p) == Array{Float64, 1}
                    if pre != []
                        pSS = vcat(p, pre);
                    else
                        pSS = p;
                    end
                else
                    if pre != []
                        pSS = vcat([p[i].value for i in 1:4], pre);
                    else
                        pSS = [p[i].value for i in 1:4];
                    end
                end


                final = zeros(maxtime, Neq);
                finalBay2 = Array{Any,2}(undef, 1, maxtime)

                y_al = JuliaConExampSteadyState(pSS,ivss) # Calculation of initial guesses for steady state
                y0 = y_al;
                initialV = y0;
                i = 1;

                for q in collect(1:Nevents)

                    lts = length(ts[convert.(Int, (sp[q]+1):sp[q+1]+1)]);  # General way to define the number of elements in each event series

                    Tevent = ts[convert.(Int,(sp[q]+1):sp[q+1]+1)];  # General way to extract the times of each event
                    I = inputs[i:(i+(0))];
                    pSte = vcat(p, I);

                    if q == 1
                        if typeof(p) == Array{Float64, 1}
                            prob = ODEProblem(JuliaConExampODE!,initialV,(ts[convert.(Int, (sp[q]+1))],ts[convert.(Int, sp[q+1]+1)]),pSte);
                            part1 = DifferentialEquations.solve(prob, DifferentialEquations.Tsit5(),reltol=1.0e-9,abstol=1.0e-9,saveat=1);
                        else
                            prob = remake(prob1, p=pSte)
                            prob = remake(prob, u0=initialV)
                            prob = remake(prob, tspan = convert.(Int, (ts[convert.(Int,(sp[q]+1))],ts[convert.(Int,sp[q+1]+1)])))
                            part1 = DifferentialEquations.solve(prob,DifferentialEquations.Tsit5(),reltol=1.0e-9,abstol=1.0e-9,saveat=1)
                        end
                    else
                        if typeof(p) == Array{Float64, 1}
                            prob = ODEProblem(JuliaConExampODE!,initialV,(ts[convert.(Int, (sp[q]+1))],ts[convert.(Int,sp[q+1]+1)]),pSte);
                            part1 = DifferentialEquations.solve(prob, DifferentialEquations.Tsit5(),reltol=1.0e-9,abstol=1.0e-9,saveat=1);
                        else
                            prob = remake(prob1, p=pSte)
                            prob = remake(prob, u0=initialV)
                            prob = remake(prob, tspan = (ts[convert.(Int,(sp[q]+1))],ts[convert.(Int,sp[q+1]+1)]))
                            part1 = DifferentialEquations.solve(prob,DifferentialEquations.Tsit5(),reltol=1.0e-9,abstol=1.0e-9,saveat=1)
                        end
                    end

                    if typeof(p) == Array{Float64, 1}
                        initialV = part1[lts];
                    else
                        initialV = [part1[:,lts][h].value for h in 1:2];
                    end

                    i+=1;


                    if typeof(p) == Array{Float64, 1}
                        for d in collect((sp[q]+1):(sp[q]+lts))
                            final[convert.(Int,d),:] = part1[convert.(Int, d-sp[q])];
                        end
                    else
                        for d in collect(convert.(Int,(sp[q]+1):(sp[q]+lts)))
                            finalBay2[1,d] = part1[2, convert.(Int, d-sp[q])]; 

                        end
                    end

                end

                if typeof(p) == Array{Float64, 1}
                    finalExps[exp] = final;
                else
                    finalBay2EXPS[exp] = finalBay2;
                end

            end



                r = [2];
                obs = Array{Any,1}(undef, 2)
                for ob in 1:1
                    if typeof(p) == Array{Float64, 1}
                        for exp in 1:2
                            obs[exp] = finalExps[exp][convert.(Int,round.(data["tsamps"][exp])).+1, [r[ob]]];
                        end
                    else
                        for exp in 1:2
                            obs[exp] = finalBay2EXPS[exp][ob, convert.(Int,round.(data["tsamps"][exp])).+1];
                        end
                    end

                    for exp in 1:2
                        for dat in 1:length(obs[exp])
                            data["DataMean"][exp][:,ob][dat] ~ Normal(obs[exp][dat], data["DataError"][exp][ob][dat])
                        end
                    end
                end
                

        end
    