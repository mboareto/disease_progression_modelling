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

# save simulated data as a DataFrame for different types of models (linear, exponential, logistic)
data = DataFrame(x = Float64[], x0 = Float64[], t = Float64[], id = Int64[],  model = String[])
for model in ["linear", "exponential", "logistic"]

    # define the progression rate
    r = ones(n_datapoints)               # Linear progression in time
    if     model == "exponential"
        r = x0                          # Exponential progression in time
    elseif model == "logistic"
        r = x0 .* (maximum(x0) .- x0)   # Logistic progression in time
    end
    # Normalizing so that the average progression is 1 point per year
    r = r / mean(r)

    # Generate the observed scores
    x = x0 .+ r .* t .+ ϵ  

    data = vcat(data, DataFrame(x=x, x0=x0, t=t, id=id, model=model))
end

# Display the first few rows of the DataFrame
first(data, 5)
## -----------------------------------------------------------------------------------------------------------------------------

CSV.write("./simulated_data/1_simulated_data.csv", data)
## -----------------------------------------------------------------------------------------------------------------------------
