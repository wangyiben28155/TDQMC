module Evolution

export CN_Evolution!

using ..TDQMC
using ..TDQMC.Crank_Nicolson

using SparseArrays

function CN_Evolution!(P::Parameter, Dy::Dynamics, serial_num::Int)         #计算的是某一个时刻的一组系综粒子的演化矩阵, 相当于对应每一个系综粒子的波函数, 这里使用的数据结构为对应每一个系综粒子, 每一个系综粒子有空间划分格点这么多的n×n的稀疏矩阵, 内部元素的非零元素的个数大概为3n
    local λ = P.Square_Δx * 1im / P.Δt
    local Constructure = ones(typeof(λ), P.space_N - 1)
    local later_fix::SparseMatrixCSC = spdiagm(-1 => Constructure, 1 => Constructure, 0 => (λ - 2.0) .- Square_Δx .* V_ne.(P.sampling))
    local former_fix::SparseMatrixCSC = spdiagm(-1 => -Constructure, 1 => -Constructure, 0 => (λ + 2.0) .+ Square_Δx .* V_ne.(P.sampling))

    local Change_matrix_former::Vector{<:SparseMatrixCSC} = fill(spzeros(eltype(later_fix), P.space_N, P.space_N), P.electron)     #这两个向量稀 
    #疏矩阵用来计算的两个矩阵
    local Change_matrix_later::Vector{<:SparseMatrixCSC} = -deepcopy(Change_matrix_former)           #所有的元素取负值



    for i = 1:P.step_t
        Reset_matrix!(P, Dy, serial_num, Change_matrix_former, Change_matrix_later)
        Construct_matrix!(P, later_fix, former_fix, Change_matrix_former, Change_matrix_later)

        Dy.Guide_Wave[:] = Change_matrix_former .* Dy.Guide_Wave
        Dy.Guide_Wave[:] = Change_matrix_later .\ Dy.Guide_Wave

        Dy.Time += P.Δt
    end


end



end