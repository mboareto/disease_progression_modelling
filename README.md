# Model Implementation

Our modeling approach consists of three parts: 1) estimate the disease trajectory of each clinical score, 2) infer disease progression fro natural history data, and 3) infer placebo response from clinical trial data. 


The model is implemented using the probabilistic programming language Turing.jl (_Julia 1.10.5_). 
The `Project.toml` and `Manifest.toml` files contain information about dependencies, versions, package names, UUIDs, etc. 
Documentation on how to start a Julia environment can be found [here](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project).



### 1. Disease trajectories
We begin by estimating the disease trajectories for each specific clinical score throughout the progression of the disease. T
o achieve this, we employ a generalized logistic model, where the lower and upper bounds are parameters that are inferred from the available data for each score.
We found that this is a convenient approach that allows us to screen multiple potential trajectories such as linear, exponential, logistic or a combination of them, 
while keeping the model simple and dependent on only two parameters to be inferred. 

The files [`code/1_generating_simulated_data.jl`](https://github.com/mboareto/disease_progression_modelling/blob/main/code/1_generating_simulated_data.jl) is used to generate the simulated data
[`code/1_simulated_data.csv`](https://github.com/mboareto/disease_progression_modelling/blob/main/code/1_simulated_data.csv). Three scenarios are simulated (linear, exponential, and logistic). 
The file [`code/1_estimate_disease_trajectories.jl`](https://github.com/mboareto/disease_progression_modelling/blob/main/code/1_estimate_disease_trajectories.jl) is used to fit the model to the data. 



### 2. Natural history progression
Once the lower and upper bound values are estimated from data, the following equation can be used to transform the values of the clinical score (x), to a scaled value (y) so that the scaled value is now linear time. 
With a linear model, progression rate is described as a linear combination of predictive covariates, and the predictive value of each covariate ($\beta$) is estimated.  


The files [`code/2_generating_simulated_data.jl`](https://github.com/mboareto/disease_progression_modelling/blob/main/code/2_generating_simulated_data.jl) is used to generate the simulated data
[`code/2_simulated_data.csv`](https://github.com/mboareto/disease_progression_modelling/blob/main/code/2_simulated_data.csv) and simulated covariate matrix 
[`code/2_simulated_covariatematrix.csv`](https://github.com/mboareto/disease_progression_modelling/blob/main/code/2_simulated_covariatematrix.csv). 
The file [`code/2_natural_history_model.jl`](https://github.com/mboareto/disease_progression_modelling/blob/main/code/2_natural_history_model.jl) is used to fit the model to the data. 


### 3. Placebo response
Once the the parameters related to natural history progression are inferred, we used clinical trial data to infer placebo response. 

The files [`code/3_generating_simulated_data.jl`](https://github.com/mboareto/disease_progression_modelling/blob/main/code/3_generating_simulated_data.jl) is used to generate the simulated data
[`code/3_simulated_data.csv`](https://github.com/mboareto/disease_progression_modelling/blob/main/code/3_simulated_data.csv) and simulated covariate matrix 
[`code/3_simulated_covariatematrix.csv`](https://github.com/mboareto/disease_progression_modelling/blob/main/code/3_simulated_covariatematrix.csv). 
The file [`code/3_placebo_response.jl`](https://github.com/mboareto/disease_progression_modelling/blob/main/code/3_placebo_response.jl) is used to fit the model to the data. 

