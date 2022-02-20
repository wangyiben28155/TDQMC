module Crank_Nicolson

export Construct_matrix!, Reset_matrix!

using ..TDQMC
using ..TDQMC.Potential_Matrix
using SparseArrays, LinearAlgebra


function ChangePart_Matrix!(pre_alloc::Vector{<:SparseMatrixCSC}, P::Parameter, Dy::Dynamics, serial_num::Int)             #=这一部分要用到循环,                                                                                                                        所以单独拿出来优化=#
    local Group_particle::Vector{<:AbstractFloat} = Dy.Trajectory[:, serial_num]
    local Matrix_Vtd::Matrix{<:AbstractFloat} = P.Twice_Δx² * Vtd_Operator(Group_particle, P, Dy.Time[serial_num])  #= 提取第serial_num组系综粒子的位置坐标
                                                                                            计算随时间变化的势能项, 便于后面的线程并行化 =#
    for i = 1:P.electron
        pre_alloc[i] = spdiagm(complex(Matrix_Vtd[:, i]))            #把含时势能矩阵的每一列变成稀疏矩阵的对角元
    end

end


function Reset_matrix!(P::Parameter, Dy::Dynamics, serial_num::Int,                 #这个函数主要用来赋值
    Change_matrix_former::Vector{T}, Change_matrix_later::Vector{T}) where {T<:SparseMatrixCSC}

    ChangePart_Matrix!(Change_matrix_former, P, Dy, serial_num)                     #将former个变量里的稀疏矩阵重置为每一时刻的差分稀疏矩阵
    Change_matrix_later[:] = -deepcopy(Change_matrix_former)
end


function Construct_matrix!(P::Parameter, later_fix::T, former_fix::T,
    Change_matrix_former::Vector{T}, Change_matrix_later::Vector{T}) where {T<:SparseMatrixCSC}

    for i = 1:P.electron
        Change_matrix_later[i] .+= later_fix
        Change_matrix_former[i] .+= former_fix
    end

end

end