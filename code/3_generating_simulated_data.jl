# Activate project environment
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# Loading packages
using DataFrames, CSV, Distributions

## -----------------------------------------------------------------------------------------------------------------------------
# Generate simulated data for testing the model

n_individuals = 1000
n_timepoints  = rand(1:5, n_individuals)  # Random number of time points (betwen 1 and 5) per individual
n_datapoints  = sum(n_timepoints)

# Initialize arrays to store the data
id = Vector{  Int64}(undef, n_datapoints)
x0 = Vector{Float64}(undef, n_datapoints)
t  = Vector{Float64}(undef, n_datapoints)

# Populate the arrays with simulated data
idx = 1
for i in 1:n_individuals
    x0_i = rand(0:20)           # Initial score for individual i ranging from 0 to 20
    for j in 1:n_timepoints[i]
        id[idx] = i
        x0[idx] = x0_i
        t[ idx] = j
        idx += 1
    end
end
## -----------------------------------------------------------------------------------------------------------------------------


# Adding noise
ϵ = rand(Normal(), n_datapoints) 

# weight of each covariate in predicting progression rate (β)
# 30 non relevant covariates, 20 relevant covariates with increased relevance (higher β)
β = vcat(zeros(30), collect(0.01:0.01:0.2))  

# Generating a covariate matrix 
n_covariates = length(β)
Z = zeros(n_individuals, n_covariates)
for j in 1:n_covariates
    Z[:, j] = rand(Normal(), n_individuals)
end

# progression rate for each individual
β0 = rand(n_individuals)          
r = β0 .+ Z * β           

# ----------- same for placebo response -----------------

# Generating a covariate matrix 
n_covariates = length(β)
Z_α = zeros(n_individuals, n_covariates)
for j in 1:n_covariates
    Z_α[:, j] = rand(Normal(), n_individuals)
end

# placebo response for each individual
α = Z_α * β 
# --------------------------------------------------------

# Generate the observed changes in scores from baseline
δy = α[id] .+ r[id] .* t .+ ϵ  

## -----------------------------------------------------------------------------------------------------------------------------

# data into a DataFrame
data = DataFrame(δy = δy, t = t, id = id)
# saving longitudinal data
CSV.write("./simulated_data/3_simulated_data.csv", data)
# saving covariate matrix data
CSV.write("./simulated_data/3_simulated_covariatematrix.csv", DataFrame(Z_α, :auto))

## -----------------------------------------------------------------------------------------------------------------------------
