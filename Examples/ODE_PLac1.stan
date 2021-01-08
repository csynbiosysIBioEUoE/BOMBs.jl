// ------------------------- TOGGLE SWITCH STAN MODEL 3 ------------------------- //


// Stan model script containing the inducer exchange model for the Toggle Switch developed by L. Bandiera presented 
// in the paper "Information content analysis reveals desirable aspects of in vivo experiments of a synthetic circuit".
// The script can be used for inference on a single experimental result or for a multiexperimental inference

functions{
  
  // Function containing the ODEs to be used for the inference
  
  real[] Toogle(real t, real[] y, real[] p, real[] x_r, int[] x_i){
    
    // Inputs (stimuly) definition
    // real u_IPTG = x_r[1];-------------------> x_r
    // real u_aTc = x_r[2];
    real IPTG = x_r[1];
    
    // Parameters definition
    // real k_IPTG = p[1]; -----------------------> Definition of parameter names with p
    // real k_aTc = p[2];
    // real k_L_pm0 = p[3];
    // real k_L_pm = p[4];
    // real theta_T = p[5];
    // real theta_aTc = p[6];
    // real n_aTc = p[7];
    // real n_T = p[8];
    // real k_T_pm0 = p[9];
    // real k_T_pm = p[10];
    // real theta_L = p[11];
    // real theta_IPTG = p[12];
    // real n_IPTG = p[13];
    // real n_L = p[14];
    
    real alpha1 = p[1];
    real Vm1 = p[2];
    real h1 = p[3];
    real Km1 = p[4];
    real d1 = p[5];
    real alpha2 = p[6];
    real d2 = p[7];
    real Kf = p[8];
    real sc_molec = p[9];
    
    

    // ODEs right-hand side
    // Order of equations(dInd_dt) follows as dIPTG/dt, daTc/dt, dLacI/dt and dTetR/dt
    real dInd_dt[4]; // -----------------------------------------------------------------------> Indicate number of states/equations amd generalise if someone whants to put an equation that is not a state
    
    real Cit_mrna = y[1];  // -------------------------------> Try if this works to make it more generable
    real Cit_foldedP = y[2];
    real Cit_fluo = y[3];
    real Cit_AU = y[4];
    
    dInd_dt[1] = alpha1+Vm1*(IPTG^h1/(Km1^h1+IPTG^h1))-d1*Cit_mrna;
    dInd_dt[2] = alpha2*Cit_mrna-(d2+Kf)*Cit_foldedP;
    dInd_dt[3] = Kf*Cit_foldedP-d2*Cit_fluo;
    dInd_dt[4] = sc_molec*dInd_dt[3]; // --------------------------------> Check if this actually works and find a way to generalise it (probably check if there is ad in fromnt of the state, or just tell people to write it as it is)
    
    
    
    // RESULTS
    return dInd_dt;
  }
  
  // Function type vector containing the equations where the root needs to be calculated for the ODEs steady states
  vector SteadyState(vector init, vector p, real[] x_r, int[] x_i){
    
    vector[4] alpha; // -----------------------------------------> Again, check this number when generating the file. 
    // Parameters definition
    real alpha1 = p[1];
    real Vm1 = p[2];
    real h1 = p[3];
    real Km1 = p[4];
    real d1 = p[5];
    real alpha2 = p[6];
    real d2 = p[7];
    real Kf = p[8];
    real sc_molec = p[9];
    
    // ODEs steady state equations. Order of initial guesses init is u_IPTG, u_aTc, LacI-RFP, TetR-GFP.
    // Order of equations (alpha) follows as dIPTG/dt, daTc/dt, dLacI/dt and dTetR/dt
    // ---------------------------------------------------------> init will just be IPTG here
    real IPTG = init[1];
    
    real Cit_mrna = alpha[1];  // -------------------------------> Try if this works to make it more generable
    real Cit_foldedP = alpha[2];
    real Cit_fluo = alpha[3];
    real Cit_AU = alpha[4];
    
    alpha[1] = (alpha1 + Vm1*(IPTG^h1/(Km1^h1+IPTG^h1)))/d1;
    alpha[2] = (alpha2*Cit_mrna)/(Kf+d2);
    alpha[3] = (Kf*Cit_foldedP)/d2;
    alpha[4] = sc_molec*Cit_fluo;
    
    // Results
    return alpha;
  }
  
}

data {
  
    // Observables
    int m; // Total number of data series
    int stslm; // Maximum number of rows for all the observable matrices
    int stsl[1,m]; // Number of elements at each time series for each series m
    int sts[stslm,m]; // Sampling times for each series m
    int obser;//-----------------------------------------------------------------------------------> Introduce this so we have all the data in one same array (easier generalisation). Work on generalisation in case different experiments have different obsevables?
    int obSta[obser]; // -----------------------------------------------------------------> This variable will be to know which are the observable states
    
    real Means[stslm,m,obser]; // -----------------> General arrays of means and errors
    real Erros[stslm,m,obser];
    
    // real GFPmean[stslm,m]; // estimated observables for TetR+GFP and LacI+RFP at each sampling time
    // real RFPmean[stslm,m]; 
    // real GFPstd[stslm,m]; // standard error for TetR+GFP and LacI+RFP at each sampling time
    // real RFPstd[stslm,m]; 
    
    // Inputs
    int elm; // Number of rows in the matrices for IPTG and aTc, half for the inputs and -1 for the total number of events
    int tml; // Maximum length of the rows for the sampling times matrix
    int Nsp[1,m]; // Number of event switching points (including final time) times for each series m
    real ts[tml, m]; // Time series for each serie m
    int tsl[1,m]; // Length of sampling time series per serie m
    
    int nindu; // -------------------------------------------------------------------> Number of inducers/stimuly
    real preInd[nindu,m];
    // real preIPTG[1,m]; // Values of inputs for each serie m for the ON incubation 
    // real preaTc[1,m];
    
    real Indu[elm,m,nindu];
    // real IPTG[elm,m]; // Values of inputs at each event for each serie m
    // real aTc[elm,m];
    
    real inputs[(elm*nindu),m]; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
    // real inputs[(elm*2),m]; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
    int evnT[(elm+1),m]; // Event change time points for each serie m
    
    // 24h incubation times for steady state calculation
    // int tonil; //-------------------------------------------------------------------> Introduce a check in case Y0 needs steady state simulation
    // real toni[tonil];
    // real Y0us[4,m]; //-----------------------------------------------------------------> In case user introduces a vector of Y0 instead of computing it (the number of states shoudl be introduced in the generation of the script when the user gives the model equations)
    

}

transformed data {
  int nParms = 9; // Number of parameters of the model //-------------------------------------------------------------------> Introduce number in generation of script
  int Neq = 4; // Total number of equations of the model //-------------------------------------------------------------------> Introduce number in generation of script
  int x_i[0]; // Empty x_i object (needs to be defined)
  real x_i2[0]; //--------------------------------------------------------------------------------------------------> Added to have a second empty vector
  real x_r[(elm*nindu),m]=inputs; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
  // real x_r[(elm*2),m]=inputs; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
  real ivss[Neq,m]; // Initial experimental values for the calculation of the steady state ordered as LacI+RFP, TetR+GFP -------------------------> Carefull to how I define this
  real pre[nindu,m]; // Input values during the 24h incubation ordered as IPTG, aTc
  // real pre[2,m]; // Input values during the 24h incubation ordered as IPTG, aTc
  
  
  // Introduction of the initial experimental values into the vector ivss and the inducer values for the 24h incubation
  // into the vector pre
  // -------------------------------------------> Here the options will deppend on if the user gives y0, steady state equation or ON needs to be simulated
  
  
  for(i in 1:m){
    for(j in 1:Neq){
      ivss[j,i] = 0; // ---------------------------------------> Here will need to add a check to give Y0 or stuff from preInd and data deppending on the scenario. (Perhaps in that scenario people should add Y0 vector too??) 
    };
    for(k in 1:nindu){
      pre[k,i] = preInd[k,i]; //--------------------------------------------> If people only give a Y0 value, this will be NO needed
    };
  };
  
  
  // for(i in 1:m){
  //   ivss[1,i] = preIPTG[1,i];
  //   ivss[2,i] = preaTc[1,i];
  //   ivss[3,i] = RFPmean[1,i];
  //   ivss[4,i] = GFPmean[1,i];
  //   pre[1,i] = preIPTG[1,i];
  //   pre[2,i] = preaTc[1,i];
  // };
}

parameters { //-------------------------------------------------------------------> This will change a lot depending on how the user wants to define the priors (Need to think about it)
    // ------------------------------------------------------------------> There should also be an option where the user gives a matrix of samples and we fit it
    // -----------------------------------------------------------------> If nothing, ask for at least some bounds and apply as default the wide normals (Almost uniform)
    // Parameters to be infered in the model
    real<lower=-2.1,upper=2.1> alpha1_raw;
    real<lower=-2.1,upper=2.1> Vm1_raw;
    real<lower=-2.1,upper=2.1> h1_raw;
    real<lower=-2.1,upper=2.1> Km1_raw;
    real<lower=-2.1,upper=2.1> d1_raw;
    real<lower=-2.1,upper=2.1> alpha2_raw;
    real<lower=-2.1,upper=2.1> d2_raw;
    real<lower=-2.1,upper=2.1> Kf_raw;
    real<lower=-2.1,upper=2.1> sc_molec_raw;
    
}



transformed parameters { // ---------------------------------------------------------------> This will need to be checked too. It could be defined by the user, or automatically if priors are normal, or apply the remapping script I have. Need to think a bit more.
  // Introduction of the parameters in an indexed object with the pertinent reparameterisation to obtain the parameter
  // values to be passed to the ODE and steady state functions
  real theta[nParms];
  theta[1] = exp(((alpha1_raw)*(1.15129254649702))+(-3.2188758248682));
  theta[2] = exp(((Vm1_raw)*(1.15129254649702))+(-2.30258509299405));
  theta[3] = exp(((h1_raw)*(1.15129254649702))+(-3.50655789731998));
  theta[4] = exp(((Km1_raw)*(1.15129254649702))+(2.30258509299405));
  theta[5] = exp(((d1_raw)*(1.15129254649702))+(3.40119738166216));
  theta[6] = exp(((alpha2_raw)*(1.15129254649702))+(2.30258509299405));
  theta[7] = (((d2_raw)*(1.25))+(2.5));
  theta[8] = (((Kf_raw)*(1.25))+(2.5));
  theta[9] = exp((((sc_molec_raw))*(1.15129254649702))-2.30258509299405);

}

model {
  
  // Intermediate parameters
  int i; // Increasing index for the inputs
  vector[4] ing; // Vector that will include the solution of the algebraic solution for the steady state of the model
  // real ssv[tonil,Neq]; // Real that will include the solution of the ODE for the 24h incubation ----------------------------> Not allways necessary
  real Y0[Neq,m]; // Initial values for the ODEs variables at the first event
  
  real yhat[stslm,m,obser]; // ---------------------------------------------------------------------------> Generall array to include all the observables (easier generalisation)
  // real yhat3[stslm,m]; // Reals that will include the values of RFP and GFP over time separately
  // real yhat4[stslm,m];
  
  // Reparameterised Priors definition 
  alpha1_raw ~ normal(0,1);
  Vm1_raw ~ normal(0,1);
  h1_raw ~ normal(0,1);
  Km1_raw ~ normal(0,1);
  d1_raw ~ normal(0,1);
  alpha2_raw ~ normal(0,1);
  d2_raw ~ normal(0,1);
  Kf_raw ~ normal(0,1);
  sc_molec_raw ~ normal(0,1);

  
  // Likelihood
  
  for (j in 1:m){
    real ivst[Neq]; // Initial value of the states 
    real y_hat[(tsl[1,j]),Neq]; // Object to include the ODEs solutions for each state
    
    // Calculation of initial guesses//----------------------------------------------------------------------------------> Alternative 1
    ing = SteadyState(to_vector(preInd[1:nindu,j]), to_vector(theta), x_i2, x_i);
    for(g in 1:Neq){
      Y0[g,j] = ing[g];
    };
    
    
    // // Calculation of initial guesses for steady state ------------------------------------------> This will change deppending on the scenario the user defines
    // ing = SteadyState(to_vector(ivss[1:4,j]), to_vector(theta), pre[1:2,j], x_i); 
    // Y0[1,j] = ing[1];
    // Y0[2,j] = ing[2];
    // Y0[3,j] = ing[3];
    // Y0[4,j] = ing[4];
    // // 24h incubation calculation for the steady state
    // ssv = integrate_ode_bdf(Toogle, Y0[,j],0,toni,theta,pre[1:2,j], x_i, 1e-9, 1e-9, 1e7); 
    // 
    // Y0[,j] = ssv[tonil];
    i = 1;
    
    // Loop (over the number of events/inputs) to solve the ODE for each event stopping the solver and add them to the final object y_hat
    for (q in 1:Nsp[1,j]-1){
      
      int itp = evnT[q,j];  // Initial time points of each event
      int lts = num_elements(ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j]);  // Length of the time series for each event
      real part1[lts,Neq]; // Temporary object that will include the solution of the ODE for each event at each loop
      
      // Calculation of the solution for the ODEs where for events that are not the first one. The time series starts one minute before the original point of the time serie overlaping with the last point of the previous event with same state values at the time
      if (q == 1){
        ivst = Y0[,j];
        part1 = integrate_ode_bdf(Toogle,ivst,itp,ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j],theta,to_array_1d(inputs[i:(i+(nindu-1)),j]), x_i, 1e-9, 1e-9, 1e7);
        // part1 = integrate_ode_bdf(Toogle,ivst,itp,ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j],theta,to_array_1d(inputs[i:(i+1),j]), x_i, 1e-9, 1e-9, 1e7);
      }
      else{
        part1 = integrate_ode_bdf(Toogle, ivst,(itp-1e-7),ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j],theta,to_array_1d(inputs[i:(i+(nindu-1)),j]), x_i, 1e-9, 1e-9, 1e7);
        // part1 = integrate_ode_bdf(Toogle, ivst,(itp-1e-7),ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j],theta,to_array_1d(inputs[i:(i+1),j]), x_i, 1e-9, 1e-9, 1e7);
      }

      // Modification of the initial state values for the next event
      ivst = part1[lts];
      // Increase index for inputs
      i=i+nindu;
      
      // Introduction of the result of part1 into the object y_hat
      for (y in (itp+1):(itp+lts)){
        y_hat[(y),]=(part1)[(y-itp),];
      };
    };

    // Likelihood definition (residuals) at each sampling time
    
    for (t in 1:stsl[1,j]){ //--------------------------------------------------------------------------------------> General form
      for(ob in 1:obser){
        yhat[t,j,ob] = y_hat[(sts[t,j]+1),obSta[ob]]; //----------------------------------------------------------------> Will need to double check if this is enough or there will be soem confusion on selecting the observable states and aching it with data
        Means[t,j,ob] ~ normal(yhat[t,j,ob],Erros[t,j,ob]); //------------------------------------------------------> Normal as default, but should be able to include other optios. Also have a check in case the user whants to define different ones for different observables 
      }
    }
    
    
    // for (t in 1:stsl[1,j]){
    //   yhat3[t,j] = y_hat[(sts[t,j]+1),3];
    //   yhat4[t,j] = y_hat[(sts[t,j]+1),4];
    //   RFPmean[t,j] ~ normal(yhat3[t,j],RFPstd[t,j]);
    //   GFPmean[t,j] ~ normal(yhat4[t,j],GFPstd[t,j]);
    // }

  };

}
