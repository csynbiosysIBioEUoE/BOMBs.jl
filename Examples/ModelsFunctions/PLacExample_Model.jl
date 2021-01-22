
    using DifferentialEquations
    using OrdinaryDiffEq
    using DiffEqBase
    using Sundials
    using ODEInterfaceDiffEq
        
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

        function PLacExample_solvecoupledODE(ts, p, sp, inputs, ivss, pre=[])

            maxtime = length(ts);
            Nsp = length(sp);
            Nevents = length(sp)-1;
            Neq = 4;
            if pre != []
                pSS = vcat(p, pre);
            else
                pSS = p;
            end

            final = zeros(maxtime, Neq);

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
                    prob = ODEProblem(PLacExampleODE!,initialV,(ts[(sp[q]+1)],ts[sp[q+1]+1]),pSte);
                    part1 = solve(prob, Tsit5(),reltol=1.0e-9,abstol=1.0e-9,saveat=1);
                else
                    prob = ODEProblem(PLacExampleODE!,initialV,(ts[(sp[q]+1)],ts[sp[q+1]+1]),pSte);
                    part1 = solve(prob, Tsit5(),reltol=1.0e-9,abstol=1.0e-9,saveat=1);
                end

                initialV = part1[lts];
                i+=1;

                for d in collect((sp[q]+1):(sp[q]+lts))
                    final[d,:] = part1[d-sp[q]];
                end
            end
            return(final)

        end
        
        # ts -> time vector going from t=0 to t=end every 1
        # pD -> Matrix containing all the draws for the parameters (can be a single vector if nly 1 draw of parameters is considered)
        # sp -> vector with switching times (t=0 and t=end needs to be included). Length of vector will be length(inducer steps)+1
        # inputs -> Vecor with the inducers (if more than 1 inducer in the system, the shape will be a1,b1,c1,..., a2,b2,c2,... aend,bend,cend,... where each leter represents a different inducer)
        # ivss -> Initial Y0 vector
        # samps -> Vector of Time points that want to be extracted from the simulations (it considers the initial time point as 0)
        # pre -> Vector with the inducers in the ON. It can be empty

        function PLacExample_SolveAll(ts, pD, sp, inputs, ivss, samps, pre=[])

            if length(size(pD)) == 1
                pD = reshape(pD,size(pD)[1],1);
            end

            if size(pD)[2] != 9
                pD = pD';
            end

            if length(ivss)/4 > 1
                if size(ivss)[2] != 4
                    ivss = ivss';
                end
            end

            AllSolTest = zeros(length(samps), 4, length(pD[:,1]));

            for drawInd in collect(1:length(pD[:,1]))
                p = pD[drawInd,:];
                if length(ivss)/4 > 1
                    ivss2 = ivss[drawInd,:];
                else
                    ivss2 = ivss;
                end
                temp = PLacExample_solvecoupledODE(ts, p, sp, inputs, ivss2, pre);
                AllSolTest[:,:,drawInd] = temp[convert.(Int,samps.+1),:];
            end

            return(AllSolTest);
        end
        