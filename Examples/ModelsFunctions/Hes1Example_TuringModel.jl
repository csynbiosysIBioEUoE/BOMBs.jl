
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

        function Hes1ExampleODE!(du,u,p,t)
            m, p1, p2 = u;

            k1, P0, v, h, inps = p;

        du[1] = dm = -0.03*m + (1/(1+(p2/P0)^h))*inps;
        du[2] = dp1 = -0.03*p1 + v*m - k1*p1;
        du[3] = dp2 = -0.03*p2 + k1*p1;


        end
        
        # p -> parameter vector (plus inducer vector at the end)
        # I -> initial Y0 vector

        function Hes1ExampleSteadyState(p,I)

            m, p1, p2 = I; 
               return(I) 
 

        end
        
        # ts -> time vector going from t=0 to t=end every 1
        # p -> parameter vector
        # sp -> vector with switching times (t=0 and t=end needs to be included). Length of vector will be length(inducer steps)+1
        # inputs -> Vecor with the inducers (if more than 1 inducer in the system, the shape will be a1,b1,c1,..., a2,b2,c2,... aend,bend,cend,... where each leter represents a different inducer)
        # ivss -> Initial Y0 vector
        # pre -> Vector with the inducers in the ON. It can be empty

        @model function Hes1Example_TuringModel(data,prob1)

            p1 ~ Truncated(Distributions.Normal(0.05, 0.025), 0.0, 0.1)
            p2 ~ Truncated(Distributions.Normal(1.0, 0.5), 0.0, 2.0)
            p3 ~ Truncated(Distributions.Normal(0.05, 0.025), 0.0, 0.1)
            p4 ~ Truncated(Distributions.Normal(5.5, 2.25), 1.0, 10.0)

            p = [p1,p2,p3,p4,];

            Neq = 3;

            for exp in 1:1

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

                y_al = Hes1ExampleSteadyState(pSS,ivss) # Calculation of initial guesses for steady state
                y0 = y_al;
                initialV = y0;
                i = 1;

                for q in collect(1:Nevents)

                    lts = length(ts[(sp[q]+1):sp[q+1]+1]);  # General way to define the number of elements in each event series

                    Tevent = ts[(sp[q]+1):sp[q+1]+1];  # General way to extract the times of each event
                    I = inputs[i:(i+(0))];
                    pSte = vcat(p, I);

                    if q == 1
                        if typeof(p) == Array{Float64, 1}
                            prob = ODEProblem(Hes1ExampleODE!,initialV,(ts[(sp[q]+1)],ts[sp[q+1]+1]),pSte);
                            part1 = DifferentialEquations.solve(prob, DifferentialEquations.Tsit5(),reltol=1.0e-9,abstol=1.0e-9,saveat=1);
                        else
                            prob = remake(prob1, p=pSte)
                            prob = remake(prob, u0=initialV)
                            prob = remake(prob, tspan = convert.(Int, (ts[(sp[q]+1)],ts[sp[q+1]+1])))
                            part1 = DifferentialEquations.solve(prob,DifferentialEquations.Tsit5(),reltol=1.0e-9,abstol=1.0e-9,saveat=1)
                        end
                    else
                        if typeof(p) == Array{Float64, 1}
                            prob = ODEProblem(Hes1ExampleODE!,initialV,(ts[(sp[q]+1)],ts[sp[q+1]+1]),pSte);
                            part1 = DifferentialEquations.solve(prob, DifferentialEquations.Tsit5(),reltol=1.0e-9,abstol=1.0e-9,saveat=1);
                        else
                            prob = remake(prob1, p=pSte)
                            prob = remake(prob, u0=initialV)
                            prob = remake(prob, tspan = (ts[(sp[q]+1)],ts[sp[q+1]+1]))
                            part1 = DifferentialEquations.solve(prob,DifferentialEquations.Tsit5(),reltol=1.0e-9,abstol=1.0e-9,saveat=1)
                        end
                    end

                    if typeof(p) == Array{Float64, 1}
                        initialV = part1[lts];
                    else
                        initialV = [part1[:,lts][h].value for h in 1:3];
                    end

                    i+=1;


                    if typeof(p) == Array{Float64, 1}
                        for d in collect((sp[q]+1):(sp[q]+lts))
                            final[d,:] = part1[d-sp[q]];
                        end
                    else
                        for d in collect((sp[q]+1):(sp[q]+lts))
                            finalBay2[1,d] = part1[1, d-sp[q]]; 

                        end
                    end

                end




                for ob in 1:1
                    if typeof(p) == Array{Float64, 1}
                        obs = final[convert.(Int,round.(data["tsamps"][exp])).+1, [1]];
                    else
                        obs = finalBay2[ob, convert.(Int,round.(data["tsamps"][exp])).+1];
                    end

                    for dat in 1:length(obs)
                        data["DataMean"][ob][dat] ~ Normal(obs[dat], data["DataError"][ob][1][dat])
                    end
                end
                

            end

        end
    