# Activate project environment
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# Loading packages
using Turing, ReverseDiff, Distributions, DataFrames, CSV

# ------------------------------ loading specific functions -------------------------------------------------------

# compute covariate score (β) for R2D2 HorseShoe prior
function compute_betas(chn, n_cov)
    n_chn = size(chn, 1) # number of chains
    n_thr = size(chn, 3) # number of threads

    B = zeros(n_cov, n_chn, n_thr)  # Matrix to store the final covariate scores
    z = zeros(n_cov, n_chn, n_thr)  # Matrix to store z values
    ϕ = zeros(n_cov, n_chn, n_thr)  # Matrix to store ϕ values

    # Extract z and ϕ values from the input data for each covariate
    for i in 1:n_cov
        z[i, :, :] = chn[Symbol("z[$i]")]  # Extract z values for covariate i
        ϕ[i, :, :] = chn[Symbol("ϕ[$i]")]  # Extract ϕ values for covariate i
    end

    # Extract R2 and Ω_β values from chain
    R2 = chn[:R2].data
    Ω_β = chn[:Ω_β].data

    # Compute τ2
    τ2 = Ω_β .* R2 ./ (1 .- R2)
    # Compute the covariate scores (β) using the formula B = z * sqrt(ϕ * τ2)
    B = z .* sqrt.(ϕ .* repeat(reshape(τ2, 1, :, n_thr), n_cov))

    # Return the computed covariate scores
    return B
end

# Define the Turing model for estimating natural history
@model function model_natural_history(δy::Vector{Float64}, t::Vector{Float64}, id::Vector{Int64}, Z::Matrix{Float64})
    # Precompute parameters
    N = size(Z, 1)       # Number of individuals
    J = size(Z, 2)       # Number of covariates
    r_e = mean(δy ./ t)  # Broad estimate of the rate from data

    # Define the model priors
    σ   ~ truncated(Normal(0, 1), 0, Inf)  # Standard deviation, constrained to be positive
    Ω_β ~ truncated(Normal(0, 1), 0, Inf)  # Prior for Ω_β, constrained to be positive
    Ω_r ~ truncated(Normal(0, 1), 0, Inf)  # Prior for Ω_r, constrained to be positive

    # ------------------------- Horse Shoe prior (R2D2) -------------------------------
    # based on Jose Storopoli's implementation: https://discourse.julialang.org/t/regularized-horseshoe-prior/71599/2
    mean_R2 = 0.5
    prec_R2 = 2.0
    cons_D2 = 1.0
    z ~ filldist(Normal(), J)
    R2 ~ Beta(mean_R2 * prec_R2, (1 - mean_R2) * prec_R2) # R2 parameter
    ϕ ~ Dirichlet(J, cons_D2)
    τ2 = Ω_β * R2 / (1 - R2)
    β = z .* sqrt.(ϕ * τ2)
    # ---------------------------------------------------------------------------------

    β_0 ~ Normal(r_e, 2.0 * abs(r_e))  # Prior for the intercept β_0
    r_μ_i = (β_0 .+ Z * β)             # Compute the mean progression rate for each individual
    r_η_i ~ MvNormal(zeros(N), Ω_r)    # Random effects with multivariate normal distribution

    r_i = r_μ_i .+ r_η_i               # progression rate for each invidual
    δμ = r_i[id] .* t                  # Compute the mean change in score
    δy ~ MvNormal(δμ, σ)               # Likelihood of the observed data
end

## -----------------------------------------------------------------------------------------------------------------------------

# loading simulated data
data = CSV.read("./simulated_data/2_simulated_data.csv", DataFrame)

# model input
δy = data[:, :δy]
t  = data[:, :t ]
id = data[:, :id]

# loading covariate matrix
Z = Array(CSV.read("./simulated_data/2_simulated_covariatematrix.csv", DataFrame))

# Establish the Turing model
m = model_natural_history(δy, t, id, Z)                             

# Run sampling using the NUTS algorithm
chn = sample(m, NUTS(500, 0.65, adtype=AutoReverseDiff(false)), 200)  

# Estimate the weight of each covariate (β) in predicting progression rate. 
β_HSR2D2 = mean(compute_betas(chn, size(Z, 2)), dims=[2])[:]

## -----------------------------------------------------------------------------------------------------------------------------
#  -----------------------------------------------------------------------------------------------------------------------------
#  -----------------------------------------------------------------------------------------------------------------------------

# Comparing Horse-Shoe prior with Normal distribution
@model function model_natural_history_normal(δy::Vector{Float64}, t::Vector{Float64}, id::Vector{Int64}, Z::Matrix{Float64})
    # Precompute parameters
    N = size(Z, 1)       # Number of individuals
    J = size(Z, 2)       # Number of covariates
    r_e = mean(δy ./ t)  # Broad estimate of the rate from data

    # Define the model priors
    σ   ~ truncated(Normal(0, 1), 0, Inf)  # Standard deviation, constrained to be positive
    Ω_β ~ truncated(Normal(0, 1), 0, Inf)  # Prior for Ω_β, constrained to be positive
    Ω_r ~ truncated(Normal(0, 1), 0, Inf)  # Prior for Ω_r, constrained to be positive

    # ---------------------------------------------------------------------------------
    β ~ filldist(Normal(), J)
    # ---------------------------------------------------------------------------------

    β_0 ~ Normal(r_e, 2.0 * abs(r_e))  # Prior for the intercept β_0
    r_μ_i = (β_0 .+ Z * β)             # Compute the mean progression rate for each individual
    r_η_i ~ MvNormal(zeros(N), Ω_r)    # Random effects with multivariate normal distribution

    r_i = r_μ_i .+ r_η_i               # progression rate for each invidual
    δμ = r_i[id] .* t                  # Compute the mean change in score
    δy ~ MvNormal(δμ, σ)               # Likelihood of the observed data
end


# Establish the Turing model
m = model_natural_history_normal(δy, t, id, Z)                             

# Run sampling using the NUTS algorithm
chn = sample(m, NUTS(500, 0.65, adtype=AutoReverseDiff(false)), 200)  

# Estimate the weight of each covariate (β) in predicting progression rate. 
β_normal = [mean(chn[Symbol("β[$i]")]) for i in 1:length(β)]


## -----------------------------------------------------------------------------------------------------------------------------
# Comparing Horse-Shoe prior with Normal distribution

β = vcat(zeros(30), collect(0.01:0.01:0.2))  # 30 non relevant covariates, 20 relevant covariates

# estimating Residual Sum of Squares
RSS_HSR2D2 = sum((β_HSR2D2 .- β) .^ 2)
RSS_normal = sum((β_normal .- β) .^ 2)

# R2D2 HorseShoe prior leads to lower RSS as compared to normal distribution prior:
RSS_R2D2 / RSS_normal