module Evolution

export CN_Evolution!

using ..TDQMC
using ..TDQMC.Potential_Matrix
using ..TDQMC.Crank_Nicolson
using ..TDQMC.Trajectory
using ..TDQMC.Quantity

using SparseArrays


function record(P::Parameter, Dy::Dynamics)                                      #这个函数虽然和function_1里的函数同名但是作用域是隔离的
    local df = DataFrame()
    

end


function CN_Evolution!(P::Parameter, Dy::Dynamics, serial_num::Int;
    later_fix::SparseMatrixCSC, former_fix::SparseMatrixCSC)         #计算的是某一个时刻的一组系综粒子的演化矩阵, 相当于对应每一个系综粒子的波函数, 这里使用的数据结构为对应每一个系综粒子, 每一个系综粒子有空间划分格点这么多的n×n的稀疏矩阵, 内部元素的非零元素的个数大概为3n

    local Change_matrix_former::Vector{<:SparseMatrixCSC} = [spzeros(eltype(later_fix), P.space_N, P.space_N) for i = 1:P.electron]     #这两个向量稀 
    #疏矩阵用来计算的两个矩阵
    local Change_matrix_later::Vector{<:SparseMatrixCSC} = -deepcopy(Change_matrix_former)           #所有的元素取负值
    local Vec_wave = view(Dy.Guide_Wave, :, serial_num)

    if imag(P.Δt) == 0.0
    
        for i = 1:P.step_t
            Movement!(P, Dy, serial_num, dt = P.Δt / ifelse(i == 1, 2.0, 1.0))
    
            Reset_matrix!(P, Dy, serial_num, Change_matrix_former, Change_matrix_later)
            Construct_matrix!(P, later_fix, former_fix, Change_matrix_former, Change_matrix_later)
    
            Vec_wave[:] = Change_matrix_former .* Vec_wave
            Vec_wave[:] = Change_matrix_later .\ Vec_wave
    
        end
    
        Dy.Time = P.step_t * P.Δt
    else
        return @error "the Time shold be a real number"
    end


end



end