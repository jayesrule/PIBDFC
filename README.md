# Predictor Informed Bayesian Dynamic Functional Connectivity
This is a repository for the Matlab implementation of the Predictor-Informed Dynamic Functional Connectivity model of Lee et al, along with an implementation of the graph-based Bayesian Dynamic Functional Connectivity Model by Warnick et al. 

## Examples
There are 2 matlab scripts walking through how to fit the models using simulated data: [pibdfc_example.m](https://github.com/jayesrule/PIBDFC/blob/public/pibdfc_example.m) and [bdfc_example.m](https://github.com/jayesrule/PIBDFC/blob/public/bdfc_example.m)

Following will be a walkthrough of a simulated data analysis using pibdfc.

# Simulating the Data

The model expects data to be in a TxRxN array. T is the number of timepoints, R is the number of ROIs varying across time, and N is the number of subjects measured. We do not currently support data where the subjects are measured with differing amounts of timepoints. Additionally, the model expects the covariate data to be a TxKxN array as well. With K being the number of covariates measured along side the outcome time series.
We have a function in the main directory, [create_gaussian_simset.m](https://github.com/jayesrule/PIBDFC/blob/public/create_gaussian_simset.m), that will simulate data in the appropriate format for example purposes.

```
run('create_gaussian_simset.m');
```

Once this function is run, time series for 5 subjects, 16 Regions, over 300 timepoints will be generated and put into an array: Y. Additionally, the TxKxN design matrix xDat is generated. 

# Setting the Prior Hyperparameters

There are prior hyper parameters that can be tuned to fit the model.

* S: [int] the number of functional connectivity states to consider
* a0: [positive num] the scale parameters for the prior on $\rho$ and $xi$, the individual level state transition parameters. A bigger number represenet a prior belief in more heterogeneity of state transition behavior between subjects.
* b0: [positive num] the scale parameters for the prior on $\eta$ and $Z$, the group level state transition parameters. A bigger number represents a prior belief in more rapid state dynamics.
* Z_prior: [SxS numeric matrix] the $i$, $j$, element of this matrix is the log odds of transitioning from state i to state j, relative to transitioning from state i to state 1. It is recommended to have Z_prior = $z*I_S$ where z is some positive number to encourage self-transitions, i.e a more stable estimated state sequence.
* eta_prior: [SxK matrix] the mean parameter for the effect of covariates on the state transition probabilities. The default is a matrix of 0s.

```
%% Setting Up for the MCMC
nsims = 1000; %Number of simulations total
nburn = 1000; %Number of burn-in sims

S = 3;
a0 = 0.01; %Prior variance for rho/xi
b0 = 0.01; %prior variance for eta/Z
tau_0 = 1; %prior scale parameter for graphical horseshoe
Z_prior = 2*eye(S); %prior mean for Z
Z_prior(1,1) = 0;
eta_prior = zeros(S,size(xDat,2)); %prior mean for eta
n_report = 100; %report current state assignments and parameter values every n_report iterations
thin = 4; %keep only 1/thin iterations 
parallel = false; 
```

# Fitting the Model

Once the prior hyper parameters are set, we only need to run the model using the pibdfc function.
```
%Perform MCMC
[post_draws] = pibdfc(zscore(Y), xDat, S, Z_prior, eta_prior, a0, b0, tau_0, nsims, nburn, n_report, thin, parallel);
```

# Posterior Summaries
The posterior samples are stored in a structure, labeled post_draws in the example. We can generate many posterior summaries from this output.

### Estimated State sequence
We can output a NxT matrix of the most likely state sequence for each subject. We can do so by plugging the posterior samples into the [mapstates_mode](https://github.com/jayesrule/PIBDFC/blob/public/pibdfc_functions/mapstates_mode.m) function.
```
map_states = mapstates_mode(post_draws.states_post); %compute state sequence estimate
```

### Graph Selection
We can output the location of the 0s in the conditional indpendence graph for each state using our false discovery rate approach. The [fdr_est](https://github.com/jayesrule/PIBDFC/blob/public/pibdfc_functions/fdr_est.m) function will output these graphs at a prespecified fdr level (0.05 in our example).

```
%% Creating summaries
est_tot = fdr_est(post_draws.kappa_post, 0.05); %compute selection graph
```

### The State Specific Partial Correlation Matrices
Our model generates posterior samples of the state specific precision matrices $\Omega$. If, for interpretability reasons, we wish to transform these into partial correlation matrices, we can use the included [prec2parcor](https://github.com/jayesrule/PIBDFC/blob/public/pibdfc_functions/prec2parcor.m) function.
```
parcor = prec2parcor(mean(post_draws.Omega_post,4));
```


