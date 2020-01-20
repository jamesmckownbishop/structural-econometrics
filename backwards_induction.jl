# This script computes expected values for an optimal stopping problem

using Random, Distributions

Random.seed!(9999)

δ = 0.96
μ = [-1, 0, -2, 0]
σ = [2, 1, 1, 1]
T = length(μ)

d = Normal()
v⁺ = zeros(T)
for i in T-1:-1:1
    z = - (μ[i+1] + v⁺[i+1]) / σ[i+1]
    v⁺[i] = δ * (σ[i+1] * (pdf(d, z) - z * (1 - cdf.(d, z))) + v⁺[i+1])
end

hazards = cdf.(Normal.(μ, σ), -v⁺)
