
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

        function PLacExampleODE!(du,u,p,t)
            Cit_mrna, Cit_foldedP, Cit_fluo, Cit_AU = u;

            alpha1, Vm1, h1, Km1, d1, alpha2, d2, Kf, sc_molec, IPTG = p;

        du[1] = dCit_mrna=alpha1+Vm1*(IPTG^h1/(Km1^h1+IPTG^h1))-d1*Cit_mrna;
        du[2] = dCit_foldedP=alpha2*Cit_mrna-(d2+Kf)*Cit_foldedP;
        du[3] = dCit_fluo=Kf*Cit_foldedP-d2*Cit_fluo;
        du[4] = dCit_AU = sc_molec*dCit_fluo;


        end
        
        # p -> parameter vector (plus inducer vector at the end)
        # I -> initial Y0 vector

        function PLacExampleSteadyState(p,I)

            expCit_mrna, expCit_foldedP, expCit_fluo, expCit_AU = I; 
        alpha1, Vm1, h1, Km1, d1, alpha2, d2, Kf, sc_molec, IPTG = p;
 
 
         alp = zeros(4); 
 
        alp[1] = Cit_mrna = (alpha1 + Vm1*(IPTG^h1/(Km1^h1+IPTG^h1)))/d1;;
        alp[2] = Cit_foldedP = (alpha2*Cit_mrna)/(Kf+d2);;
        alp[3] = Cit_fluo = (Kf*Cit_foldedP)/d2;;
        alp[4] = Cit_AU= sc_molec*Cit_fluo;;
        return(alp) 
 

        end
        
        # ts -> time vector going from t=0 to t=end every 1
        # p -> parameter vector
        # sp -> vector with switching times (t=0 and t=end needs to be included). Length of vector will be length(inducer steps)+1
        # inputs -> Vecor with the inducers (if more than 1 inducer in the system, the shape will be a1,b1,c1,..., a2,b2,c2,... aend,bend,cend,... where each leter represents a different inducer)
        # ivss -> Initial Y0 vector
        # pre -> Vector with the inducers in the ON. It can be empty

        @model function PLacExample_TuringModel(data,prob1)

            p1 ~ Truncated(Distributions.Normal(0.2475194, 0.1237403), 3.88e-5, 0.495)
            p2 ~ Truncated(Distributions.Normal(0.2669, 0.11405), 0.0388, 0.495)
            p3 ~ Truncated(Distributions.Normal(2.7, 1.1), 0.5, 4.9)
            p4 ~ Truncated(Distributions.Normal(6.0, 2.0), 2.0, 10.0)
            p5 ~ Truncated(Distributions.Normal(0.11885000000000001, 0.055575), 0.0077, 0.23)
            p6 ~ Truncated(Distributions.Normal(3.525, 1.6408500000000001), 0.2433, 6.8067)
            p7 ~ Truncated(Distributions.Normal(0.1224799, 0.06121005), 5.98e-5, 0.2449)
            p8 ~ Truncated(Distributions.Normal(0.01685, 0.002425), 0.012, 0.0217)
            p9 ~ Truncated(Distributions.Normal(5.0005, 2.49975), 0.001, 10.0)

            p = [p1,p2,p3,p4,p5,p6,p7,p8,p9,];

            Neq = 4;

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
                        pSS = vcat([p[i].value for i in 1:9], pre);
                    else
                        pSS = [p[i].value for i in 1:9];
                    end
                end


                final = zeros(maxtime, Neq);
                finalBay2 = Array{Any,2}(undef, 1, maxtime)

                y_al = PLacExampleSteadyState(pSS,ivss) # Calculation of initial guesses for steady state
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
                            prob = ODEProblem(PLacExampleODE!,initialV,(ts[(sp[q]+1)],ts[sp[q+1]+1]),pSte);
                            part1 = DifferentialEquations.solve(prob, DifferentialEquations.Tsit5(),reltol=1.0e-9,abstol=1.0e-9,saveat=1);
                        else
                            prob = remake(prob1, p=pSte)
                            prob = remake(prob, u0=initialV)
                            prob = remake(prob, tspan = convert.(Int, (ts[(sp[q]+1)],ts[sp[q+1]+1])))
                            part1 = DifferentialEquations.solve(prob,DifferentialEquations.Tsit5(),reltol=1.0e-9,abstol=1.0e-9,saveat=1)
                        end
                    else
                        if typeof(p) == Array{Float64, 1}
                            prob = ODEProblem(PLacExampleODE!,initialV,(ts[(sp[q]+1)],ts[sp[q+1]+1]),pSte);
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
                        initialV = [part1[:,lts][h].value for h in 1:4];
                    end

                    i+=1;


                    if typeof(p) == Array{Float64, 1}
                        for d in collect((sp[q]+1):(sp[q]+lts))
                            final[d,:] = part1[d-sp[q]];
                        end
                    else
                        for d in collect((sp[q]+1):(sp[q]+lts))
                            finalBay2[1,d] = 3 .* part1[4, d-sp[q]]; 

                        end
                    end

                end



                obs = Array{Any,1}(undef, 1)
                if typeof(p) == Array{Float64, 1}
                    obs[1] = 3 .* final[convert.(Int,round(data["tsamps"][exp])).+1, 4]; 

                else
                    obs[1] = finalBay2[1, convert.(Int,round(data["tsamps"][exp])).+1];
                end

                for ob in 1:1
                    for dat in 1:length(obs)
                        data["DataMean"][ob][dat] ~ Normal(obs[ob][dat], data["DataError"][ob][1][dat])
                    end
                end
                

            end

        end
    