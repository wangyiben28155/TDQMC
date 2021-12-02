module Crank_Nicolson

export CN_Evolution!

using ..TDQMC
using ..TDQMC.Potential_Matrix
using SparseArrays, LinearAlgebra


function ChangePart_Matrix(P::Parameter, Dy::Dynamics, serial_num::Int)             #这一部分要用到循环, 所以单独拿出来优化
    local Group_particle::Vector{<:AbstractFloat} = Dy.Trajectory[:, serial_num]
    local Matrix_Vtd::Matrix{T} = Vtd_Operator(Group_particle, P, Dy.Time)          #提取第serial_num组系综粒子的位置坐标计算随时间变化的势能项,便于后
    local Change_Matrix::Vector{<:SparseMatrixCSC} = fill(spezeros(eltype(Matrix_Vtd), P.space_N, P.space_N), P.electron)#面的线程并行化

    for i = 1:P.electron
        @inbounds Change_Matrix[i] = spdiagm.(complex(Matrix_Vtd[:, i]))            #把含时势能矩阵的每一列变成稀疏矩阵的对角元
    end
    return Change_Matrix
end


function Construct_matrix!(P::Parameter, later_fix::T, former_fix::T,
     Change_matrix_former::Vector{T}, Change_matrix_later::Vector{T}) where {T<:SparseMatrixCSC}

    for i = 1:P.electron
        @inbounds Change_matrix_later[i] .+= later_fix
        @inbounds Change_matrix_former[i] .+= former_fix
    end

end


function Reset_matrix!(P::Parameter, Dy::Dynamics, serial_num::Int,        #这个函数主要用来赋值
    Change_matrix_former::Vector{T}, Change_matrix_later::Vector{T}) where {T<:SparseMatrixCSC}

    Change_matrix_former = ChangePart_Matrix(P, Dy, serial_num)
    Change_matrix_later = -deepcopy(Change_matrix_former)
end


function CN_Evolution!(P::Parameter, Dy::Dynamics, serial_num::Int)         #计算的是某一个时刻的一组系综粒子的演化矩阵, 相当于对应每一个系综粒子的波函数, 这里使用的数据结构为对应每一个系综粒子, 每一个系综粒子有空间划分格点这么多的n×n的稀疏矩阵, 内部元素的非零元素的个数大概为3n
    local λ = P.Square_Δx * 1im / P.Δt
    local Constructure = ones(typeof(λ), P.space_N - 1)
    local later_fix::SparseMatrixCSC = spdiagm(-1 => Constructure, 1 => Constructure, 0 => (λ - 2) .- Square_Δx .* V_ne.(P.sampling))
    local former_fix::SparseMatrixCSC = spdiagm(-1 => -Constructure, 1 => -Constructure, 0 => (λ + 2) .+ Square_Δx .* V_ne.(P.sampling))

    local Change_matrix_former::Vector{<:SparseMatrixCSC} = ChangePart_Matrix(P, Dy, serial_num)     #这两个向量稀疏矩阵用来作为计算的两个矩阵
    local Change_matrix_later::Vector{<:SparseMatrixCSC} = -deepcopy(Change_matrix_former)           #所有的元素取负值
 


    for i = 1:P.step_t
        Reset_matrix!(P, Dy, serial_num, Change_matrix_former, Change_matrix_later)
        Construct_matrix!(P, later_fix, former_fix, Change_matrix_former, Change_matrix_later)
    
        Dy.Guide_Wave = Change_matrix_former .* Dy.Guide_Wave
        Dy.Guide_Wave = Change_matrix_later .\ Dy.Guide_Wave

        Dy.Time += P.Δt
    end


end    


end