module Crank_Nicolson

export pin

using ..TDQMC
using ..TDQMC.Potential
using SparseArrays, LinearAlgebra


function Evolution_Matrix(P::Parameter, Dy::Dynamics)
    local λ = 2 * (P.Δx)^2 / P.Δt
    local Square_Δx = P.Δx^2
    local Former::SparseMatrixCSC = spzeros(typeof(λ), P.space_N, P.space_N)
    local later::SparseMatrixCSC = spzeros(typeof(λ), P.space_N, P.space_N)

    Former = spdiagm(-1 => ones(typeof(λ), P.space_N - 1), 1 => ones(typeof(λ), P.space_N - 1), 0 => V)
    later = spdiagm(-1 => -ones(typeof(λ), P.space_N - 1), 1 => -ones(typeof(λ), P.space_N - 1),  0 => V)

    return Former,later
end


function Crank_Evolution(Wave::Vector{T}) where T<:AbstractFloat
    local A





end







end