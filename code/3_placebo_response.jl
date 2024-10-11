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


exp_func(t; τ=0.05) = (1.0 - exp(-t / τ))

# Define the Turing model for estimating placebo response
@model function model_placebo_response(δy::Vector{Float64}, t::Vector{Float64}, id::Vector{Int64}, r_μ_i::Vector{Float64}, Ω_r::Float64, Z::Matrix{Float64})

    # Precompute parameters
    t_func = exp_func.(t)  # temporal dependency of placebo response (quick improvement)
    N = size(Z, 1)         # Number of individuals
    J = size(Z, 2)         # Number of covariates
    
    # Define the model priors
    σ   ~ truncated(Normal(0, 1), 0, Inf)  # Standard deviation, constrained to be positive
    Ω_β ~ truncated(Normal(0, 1), 0, Inf)  # Prior for Ω_β, constrained to be positive
    Ω_α ~ truncated(Normal(0, 1), 0, Inf)  # Prior for Ω_r, constrained to be positive
    
    # ------------------------- Horse Shoe prior (R2D2) -------------------------------
    # based on Jose Storopoli's implementation: https://discourse.julialang.org/t/regularized-horseshoe-prior/71599/2
    mean_R2 = 0.5
    prec_R2 = 2.0
    cons_D2 = 1.0
    z ~ filldist(Normal(), J)
    R2 ~ Beta(mean_R2 * prec_R2, (1 - mean_R2) * prec_R2)
    ϕ ~ Dirichlet(J, cons_D2)
    τ2 = Ω_β * R2 / (1 - R2)
    β = z .* sqrt.(ϕ * τ2)
    # ---------------------------------------------------------------------------------

    β_0 ~ Normal(0, 5)                 # Prior for the intercept β_0
    α_μ_i = (β_0 .+ Z * β)             # Compute the mean placebo response for each individual
    α_η_i ~ MvNormal(zeros(N), Ω_α)    # Random effects for placebo response
    r_η_i ~ MvNormal(zeros(N), Ω_r)    # Random effects for progression rate

    α_i = α_η_i[id] .+ α_μ_i[id]       # placebo response for each invidual
    r_i = r_η_i[id] .+ r_μ_i           # progression rate for each invidual

    μ = r_i .* t .+ α_i .* t_func      # Compute the mean change in score
    δy ~ MvNormal(μ, σ)                # Likelihood of the observed data
end

## -----------------------------------------------------------------------------------------------------------------------------

# ----------------- fitting natural history model ------------------------------
# loading simulated data
data = CSV.read("./simulated_data/2_simulated_data.csv", DataFrame)

# model input
δy = data[:, :δy]
t  = data[:, :t ]
id = data[:, :id]

# loading covariate matrix
Z_r = Array(CSV.read("./simulated_data/2_simulated_covariatematrix.csv", DataFrame))

# Establish the Turing model
m = model_natural_history(δy, t, id, Z_r)                             

# Run sampling using the NUTS algorithm
chn = sample(m, NUTS(500, 0.65, adtype=AutoReverseDiff(false)), 200)  
# -----------------------------------------------------------------------------

# ---------------- fitting placebo reponse model ------------------------------

# loading simulated data
data = CSV.read("./simulated_data/3_simulated_data.csv", DataFrame)

# model input
δy = data[:, :δy]
t  = data[:, :t ]
id = data[:, :id]

# estimating parameters related to progression rate
β_r = mean(compute_betas(chn, size(Z_r,2)), dims=[2])[:]
Ω_r = mean(chn[:Ω_r])     # population variability in the progression rate
r_μ_i = (Z_r * β_r)[id]   # individual average progression rate                        
 
# loading covariate matrix
Z_α = Array(CSV.read("./simulated_data/3_simulated_covariatematrix.csv", DataFrame))

# Establish the Turing model
m = model_placebo_response(δy, t, id, r_μ_i, Ω_r, Z_α)

# Run sampling using the NUTS algorithm
chn = sample(m, NUTS(500, 0.65, adtype=AutoReverseDiff(false)), 200)

# Estimate the weight of each covariate (β) in predicting placebo response 
β_α = mean(compute_betas(chn, size(Z_α, 2)), dims=[2])[:]
# -----------------------------------------------------------------------------


