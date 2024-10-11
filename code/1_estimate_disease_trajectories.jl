# Activate project environment
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# Loading packages
using Turing, ReverseDiff, Distributions, DataFrames, CSV

# ------------------------------ loading specific functions -------------------------------------------------------

# non-linear temporal progression transformed to linear temporal progression
logistic_to_linear(x::Float64, l::Float64, u::Float64) = -log((u - l) / (x - l) - 1.0)

# linear temporal progression transformed to non-linear temporal progression
linear_to_logistic(y::Float64, l::Float64, u::Float64) = l + (u - l) / (exp.(-y) + 1.0)

# Define the f_score function to model disease progression in a linear space
function f_score(x0::Float64, r_i::Float64, t::Float64, l::Float64, u::Float64)
    y0 = logistic_to_linear(x0, l, u)    # Covert from non-linear temporal progression transformed to linear temporal progression
    y = y0 + r_i * t                     # Linear progression over time
    return linear_to_logistic(y, l, u)   # Convert back to non-linear temporal progression
end


# Define the Turing model for estimating disease trajectories
@model function estimate_disease_trajectories(x::Vector{Float64}, x0::Vector{Float64}, t::Vector{Float64}, id::Vector{Int64},
    xmin::Float64, xmax::Float64)

    # Precompute parameters
    n = length(unique(id))  # Number of unique individuals

    # Define the model priors
    σ ~ truncated(Normal(0, 1), 0, Inf) # Standard deviation of the noise
    Ω_r ~ truncated(Normal(0, 1), 0, Inf) # Standard deviation of the random effects

    lower ~ truncated(Normal(0.5, 1.0), 0.01, 2.0) # Lower bound scaling factor
    upper ~ truncated(Normal(0.5, 1.0), 0.01, 2.0) # Upper bound scaling factor

    δ = xmax - xmin       # Range of the observed data
    l = xmin - δ * lower  # Lower bound for logistic transformation
    u = xmax + δ * upper  # Upper bound for logistic transformation

    r_dist ~ Normal(0.0, 5.0)  # Population-level rate of progression
    r_pop = fill(r_dist, n)    # Fill the population with the same rate distribution
    r_i ~ MvNormal(r_pop, Ω_r) # Individual-level random effects

    # Compute the predicted scores based on choice of l and u 
    μ = f_score.(x0, r_i[id], t, l, u)

    # Likelihood of the observed data
    x ~ MvNormal(μ, σ)
end
## -----------------------------------------------------------------------------------------------------------------------------

# loading simulated data
data = CSV.read("./simulated_data/1_simulated_data.csv", DataFrame)

# save model output as a DataFrame for different types of models (linear, exponential, logistic)
model_output = DataFrame(lower=Float64[], upper=Float64[], model=String[])
for model in ["linear", "exponential", "logistic"]
    # model input
    x  = data[data.model.==model, :x ]
    t  = data[data.model.==model, :t ]
    x0 = data[data.model.==model, :x0]
    id = data[data.model.==model, :id]

    # Estimate the lower and upper bound (parameters of disease trajectory model) using the Turing model
    xmin, xmax = extrema(x)

    # Establish the Turing model
    m = estimate_disease_trajectories(x, x0, t, id, xmin, xmax)         

    # Run sampling using the NUTS algorithm
    chn = sample(m, NUTS(500, 0.65, adtype=AutoReverseDiff(false)), 200) # run sampling 

    # save model output (inferred lower and upper parameters)
    model_output = vcat(model_output, DataFrame(lower=mean(chn[:lower]), upper=mean(chn[:upper]), model=model))
end
model_output
## -----------------------------------------------------------------------------------------------------------------------------

# Expected results
# lower, upper >> 1 :  linear model
# lower, upper << 1 :  logistic model
# lower >> upper    :  exponential model
# upper >> lower    :  exponential model

