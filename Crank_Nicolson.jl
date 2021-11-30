module Crank_Nicolson

export pin

using ..TDQMC
using ..TDQMC.Potential_Matrix
using SparseArrays, LinearAlgebra


function Evolution_Matrix(P::Parameter, Dy::Dynamics)         #计算的是某一个时刻的一组系综粒子的演化矩阵, 相当于对应每一个系综粒子的波函数, 这里使用的数据结构为对应每一个系综粒子, 每一个系综粒子有空间划分格点这么多的n×n的稀疏矩阵, 内部元素的非零元素的个数大概为3n
    local λ = 2im * 2 * (P.Δx)^2 / P.Δt
    local Square_Δx = 2 * P.Δx^2
    local Former::Matrix{<:SparseMatrixCSC} = spzeros(typeof(λ), P.space_N, P.space_N)
    local later::Matrix{<:SparseMatrixCSC} = spzeros(typeof(λ), P.space_N, P.space_N)


    Former = spdiagm(-1 => ones(typeof(λ), P.space_N - 1), 1 => ones(typeof(λ), P.space_N - 1), 0 => (λ - 2) .- Square_Δx .* V_Operator())
    later = spdiagm(-1 => -ones(typeof(λ), P.space_N - 1), 1 => -ones(typeof(λ), P.space_N - 1), 0 => (λ + 2))

    return Former, later
end


function Crank_Evolution(Wave::Vector{T}) where {T<:AbstractFloat}         #求解本征值方程得到演化的本征矢量
    local A





end







end