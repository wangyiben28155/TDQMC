module Crank_Nicolson

export Crank_Evolution!

using ..TDQMC
using ..TDQMC.Potential_Matrix
using SparseArrays, LinearAlgebra


function ChangePart_Matrix!(P::Parameter, Dy::Dynamics, serial_num::Int64)             #这一部分要用到循环, 所以单独拿出来优化
    local Group_particle::Vector{<:AbstractFloat} = Dy.Trajectory[:, serial_num]
    local Matrix_Vtd::Vector{Vector{}} = Vtd_Operator(Group_particle, P, Dy.Time)      #提取第serial_num组系综粒子的位置坐标计算随时间变化的势能项,便于后面的线程并行化

    for i = 1:P.electron
        @inbounds Dy.Change_Matrix[i] = spdiagm.(complex(Matrix_Vtd[:, i]))       #把含时势能矩阵的每一列变成稀疏矩阵的对角元
    end

end

function Construct_matrix(P::Parameter, later_fix, T_later) where {T<:SparseMatrixCSC}

    for i = 1:P.electron
        @inbounds later[i] .-= Change_matrix[i]
        @inbounds former[i] .+= Change_matrix[i]
    end

end



function Crank_Evolution!(P::Parameter, Dy::Dynamics, serial_num::Int64)         #计算的是某一个时刻的一组系综粒子的演化矩阵, 相当于对应每一个系综粒子的波函数, 这里使用的数据结构为对应每一个系综粒子, 每一个系综粒子有空间划分格点这么多的n×n的稀疏矩阵, 内部元素的非零元素的个数大概为3n
    local Square_Δx = 2 * (P.Δx)^2
    local λ = Square_Δx * 1im / P.Δt
    local Constructure = ones(typeof(λ), P.space_N - 1)
    local later_fix::Vector{<:SparseMatrixCSC} = spdiagm(-1 => Constructure, 1 => Constructure, 0 => (λ - 2) .- Square_Δx .* V_ne.(P.sampling))
    local later_change::Vector{<:SparseMatrixCSC} = spzeros
    local former_fix::Vector{<:SparseMatrixCSC} = spdiagm(-1 => -Constructure, 1 => -Constructure, 0 => (λ + 2) .+ Square_Δx .* V_ne.(P.sampling))
    local former_change::Vector{<:SparseMatrixCSC} = 

    local Change_matrix::Vector{<:SparseMatrixCSC} = ChangePart_Matrix!(P, Dy, serial_num)
    local Fix_Matrix = FixPart

    for i = 1:P.electron
        @inbounds later[i] .-= Change_matrix[i]
        @inbounds former[i] .+= Change_matrix[i]
    end


end    

end