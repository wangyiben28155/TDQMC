module Evolution

export CN_Evolution!

using ..TDQMC
using ..TDQMC.Potential_Matrix
using ..TDQMC.Crank_Nicolson
using ..TDQMC.Trajectory
using ..TDQMC.Quantity

using SparseArrays

function find_inbound!(P::Parameter, Dy::Dynamics, serial_num::Integer, Vec_Trajectory::SubArray)
    Dy.Index[serial_num] = findall(x -> abs(x) < P.scope, Vec_Trajectory)
    Dy.In_num[serial_num] = length(Dy.Index[serial_num])
end


function CN_Evolution!(P::Parameter, Dy::Dynamics, serial_num::Int;
    later_fix::SparseMatrixCSC, former_fix::SparseMatrixCSC)         #计算的是某一个时刻的一组系综粒子的演化矩阵, 相当于对应每一个系综粒子的波函数, 这里使用的数据结构为对应每一个系综粒子, 每一个系综粒子有空间划分格点这么多的n×n的稀疏矩阵, 内部元素的非零元素的个数大概为3n

    local Change_matrix_former::Vector{<:SparseMatrixCSC} = [spzeros(eltype(later_fix), P.space_N, P.space_N) for i = 1:P.electron]     #这两个向量稀 
    #疏矩阵用来计算的两个矩阵
    local Change_matrix_later::Vector{<:SparseMatrixCSC} = -deepcopy(Change_matrix_former)           #所有的元素取负值

    local Vec_wave = view(Dy.Guide_Wave, :, serial_num)
    local Vec_Trajectory = view(Dy.Trajectory, :, serial_num)

    
    Dy.Displace[1, serial_num, :] = Vec_Trajectory

    if imag(P.Δt) == 0.0
    
        for i = 1:P.step_t
            find_inbound!(P, Dy, serial_num, Vec_Trajectory)
        
            if Dy.In_num[serial_num] == 0
                break
            end
            
            Reset_matrix!(P, Dy, serial_num, Change_matrix_former, Change_matrix_later)
            Construct_matrix!(Dy, serial_num, later_fix, former_fix, Change_matrix_former, Change_matrix_later)
        
            Evolution!(Dy, serial_num, Vec_wave, Change_matrix_former, Change_matrix_later)
        
            Movement!(P, Dy, serial_num, Vec_Trajectory, dt = P.Δt)
        
            Dy.Time[serial_num] += P.Δt
            Dy.Displace[i+1, serial_num, :] = Dy.Trajectory[:, serial_num]
        end
    

    
    else
        return @error "the Time shold be a real number"
    end


end



end