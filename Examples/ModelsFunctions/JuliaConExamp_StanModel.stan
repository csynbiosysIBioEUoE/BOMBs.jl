

functions{

    real[] JuliaConExamp_ODEs(real t, real[] y, real[] p, real[] x_r, int[] x_i){

      // Inputs
      real U = x_r[1]; 


      // Parameters
      real n1 = p[1]; 
      real KM2 = p[2]; 
      real d2 = p[3]; 
      real km = p[4]; 


      // ODEs
      real dInd_dt[2];

      real mRNA = y[1]; 
      real Prot = y[2]; 


      dInd_dt[1] =  1 + (200 / (1 + ((Prot/20)*((U^2)/((KM2^2)+(U^2))))^n1)) - 0.5*mRNA; 
      dInd_dt[2] =  km*mRNA - d2*Prot; 


      // Results
      return dInd_dt;

    }


    vector SteadyState(vector init, vector p, real[] x_r, int[] x_i){

      return init; 


    }

}
    

data {

    // Observables
    int m; // Total number of data series
    int stslm; // Maximum number of rows for all the observable matrices
    int stsl[1,m]; // Number of elements at each time series for each series m
    int sts[stslm,m]; // Sampling times for each series m
    int obser;//-> Introduce this so we have all the data in one same array (easier generalisation). Work on generalisation in case different experiments have different obsevables?
    int obSta[1,obser]; // -> This variable will be to know which are the observable states

    real Means[stslm,m,obser]; // ---> General arrays of means and errors
      real Erros[stslm,m,obser];

    // Inputs
    int elm; // Number of rows in the matrices for IPTG and aTc, half for the inputs and -1 for the total number of events
    int tml; // Maximum length of the rows for the sampling times matrix
    int Nsp[1,m]; // Number of event switching points (including final time) times for each series m
    real ts[tml, m]; // Time series for each serie m
    int tsl[1,m]; // Length of sampling time series per serie m

    int nindu; // -> Number of inducers/stimuly
    real preInd[nindu,m]; // Values of inputs for each serie m for the ON incubation

    real inputs[(elm*nindu),m]; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...

    int evnT[(elm+1),m]; // Event change time points for each serie m
    real Y0us[2,m]; // Y0 vectors



}


    

transformed data {

    int nParms = 4; // Number of parameters of the model //--------> Introduce number in generation of script
    int Neq = 2; // Total number of equations of the model //-----> Introduce number in generation of script
    int x_i[0]; // Empty x_i object (needs to be defined)
    real x_r[(elm*nindu),m]=inputs; // Input values for each event ordered as IPTG, aTc, IPTG, aTc, ...
    real ivss[Neq,m] = Y0us; // Initial experimental values for the calculation of the steady state ordered as LacI+RFP, TetR+GFP -------------------------> Careful to how I define this
    real pre[nindu,m]; // Input values during the 24h incubation ordered as IPTG, aTc

    for(i in 1:m){
        for(k in 1:nindu){
          pre[k,i] = preInd[k,i]; //-> If people only give a Y0 value, this will be NO needed
        };
    };

}
    

parameters {

real n1; 
real KM2; 
real d2; 
real km; 


}

    

transformed parameters {

real theta[nParms]; 
theta[1] = exp(((n1)*(0.04143871887549682))+(1.1197892379138183)); 
theta[2] = (((KM2)*(6.389203819134893))+(15.683296447500005)); 
theta[3] = exp(((d2)*(0.041890435442655666))+(-0.23346539422180088)); 
theta[4] = exp(((km)*(0.047232512131290254))+(1.0962967683697238)); 


}

    

model {

  // Intermediate parameters
  int i; // Increasing index for the inputs
  vector[Neq] ing; // Vector that will include the solution of the algebraic solution for the steady state of the model
  
  real Y0[Neq,m]; // Initial values for the ODEs variables at the first event
  real yhat[stslm,m,obser]; // ---> Generall array to include all the observables (easier generalisation)

  // Reparameterised priors definition
n1 ~ normal(0, 1); 
KM2 ~ normal(0, 1); 
d2 ~ normal(0, 1); 
km ~ normal(0, 1); 


  // Likelihood
  for (j in 1:m){

    real ivst[Neq]; // Initial value of the states
    real y_hat[(tsl[1,j]),Neq]; // Object to include the ODEs solutions for each state

    // Calculation of Y0

    ing = SteadyState(to_vector(ivss[,j]), to_vector(theta), pre[,j], x_i);
    for(g in 1:Neq){
      Y0[g,j] = ing[g];
    };


    i = 1;

    // Loop (over the number of events/inputs) to solve the ODE for each event stopping the solver and add them to the final object y_hat

    for (q in 1:Nsp[1,j]-1){
      int itp = evnT[q,j];  // Initial time points of each event
      int lts = num_elements(ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j]);  // Length of the time series for each event
      real part1[lts,Neq]; // Temporary object that will include the solution of the ODE for each event at each loop

      // Calculation of the solution for the ODEs where for events that are not the first one. The time series starts one minute before the original point of the time serie overlaping with the last point of the previous event with same state values at the time
      if (q == 1){
        ivst = Y0[,j];
        part1 = integrate_ode_bdf(JuliaConExamp_ODEs,ivst,itp,ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j],theta,to_array_1d(inputs[i:(i+(nindu-1)),j]), x_i, 1.0e-9, 1.0e-9, 1e7);
      }
      else{
        part1 = integrate_ode_bdf(JuliaConExamp_ODEs, ivst,(itp-1e-7),ts[(evnT[q,j]+1):(evnT[q+1,j]+1),j],theta,to_array_1d(inputs[i:(i+(nindu-1)),j]), x_i, 1.0e-9, 1.0e-9, 1e7);
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

    for (t in 1:stsl[1,j]){ //----> General form

            for(ob in 1:obser){
              yhat[t,j,ob] = y_hat[(sts[t,j]+1),obSta[1,ob]]; //---> Will need to double check if this is enough or there will be soem confusion on selecting the observable states and aching it with data
              Means[t,j,ob] ~ normal(yhat[t,j,ob],Erros[t,j,ob]); //---> Normal as default, but should be able to include other optios. Also have a check in case the user whants to define different ones for different observables
            }

                    
    }

  }
}

    