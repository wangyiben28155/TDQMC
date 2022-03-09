module Evolution

export CN_Evolution!

using ..TDQMC
using ..TDQMC.Potential_Matrix
using ..TDQMC.Crank_Nicolson
using ..TDQMC.Trajectory
using ..TDQMC.Quantity

using SparseArrays, LinearAlgebra

@inline function find_inbound!(P::Parameter, Dy::Dynamics, serial_num::Integer, Vec_Trajectory::SubArray)
    Dy.Index[serial_num] = findall(x -> abs(x) < P.scope, Vec_Trajectory)
    Dy.In_num[serial_num] = length(Dy.Index[serial_num])
end

@inline function Spin_Effect(P::Parameter, Dy::Dynamics, serial_num::Integer)
    local Select_Spin::Vector = P.Spin[Dy.Index[serial_num]]
    local Spin_Effect::Matrix{<:Integer} = Select_Spin * Select_Spin'

    replace!(tril!(Spin_Effect), 0 => 1)

    return Spin_Effect
end

function Absorb(x::T; scope::T, ratio = 0.9) where {T<:AbstractFloat}
    local x₀ = ratio * scope

    if abs(x) < x₀
        return 1.0
    elseif x₀ <= abs(x) <= scope
        return cos((abs(x) - abs(x₀)) / (scope - x₀))^(1 / 6)
    end

end

function damping!(Vec_wave::SubArray{<:Vector}, Index::Vector{<:Integer}, Absorber::Vector)
    for i in Index
        Vec_wave[i] .*= Absorber
    end
end


function CN_Evolution!(P::Parameter, Dy::Dynamics, serial_num::Int;
    later_fix::SparseMatrixCSC, former_fix::SparseMatrixCSC)         #计算的是某一个时刻的一组系综粒子的演化矩阵, 相当于对应每一个系综粒子的波函数, 这里使用的数据结构为对应每一个系综粒子, 每一个系综粒子有空间划分格点这么多的n×n的稀疏矩阵, 内部元素的非零元素的个数大概为3n

    local Change_matrix_former::Vector{<:SparseMatrixCSC} = [spzeros(eltype(later_fix), P.space_N, P.space_N) for i = 1:P.electron]     #这两个向量稀 
    #疏矩阵用来计算的两个矩阵
    local Change_matrix_later::Vector{<:SparseMatrixCSC} = -deepcopy(Change_matrix_former)           #所有的元素取负值

    local Velocity_WaveFunc::Vector{Matrix{<:Complex}} = [zeros(eltype(later_fix), P.electron, P.electron) for i = 1:P.electron]
    local Vector_velocity::Vector{<:AbstractFloat} = zeros(eltype(Dy.Trajectory), P.electron)

    local Vec_wave = view(Dy.Guide_Wave, :, serial_num)
    local Vec_Trajectory = view(Dy.Trajectory, :, serial_num)

    local Spin_Matrix::Matrix{<:Integer} = zeros(Int, P.electron, P.electron)

    local record_Change::Integer = 0
    local Absorber = Absorb.(P.sampling, scope = P.scope)

    Dy.Displace[1, serial_num, :] = Vec_Trajectory

    find_inbound!(P, Dy, serial_num, Vec_Trajectory)
    Spin_Matrix = Spin_Effect(P, Dy, serial_num)

    if imag(P.Δt) == 0.0

        for i = 1:P.step_t
            record_Change = Dy.In_num[serial_num]
        
            find_inbound!(P, Dy, serial_num, Vec_Trajectory)
        
            if record_Change != Dy.In_num[serial_num]
                Spin_Matrix = Spin_Effect(P, Dy, serial_num)
        
                if Dy.In_num[serial_num] == 0
                    break
                end
            end
        
            Reset_matrix!(P, Dy, serial_num, Change_matrix_former, Change_matrix_later)
            Construct_matrix!(Dy, serial_num, later_fix, former_fix, Change_matrix_former, Change_matrix_later)
        
            Evolution!(Dy, serial_num, Vec_wave, Change_matrix_former, Change_matrix_later)

            damping!(Vec_wave, Dy.Index[serial_num], Absorber)
            
            Velocity!(P, Dy, serial_num, Velocity_WaveFunc, Vector_velocity; Spin_Matrix = Spin_Matrix)  #重置速度向量,这里只对Index索引的更新速度
            Movement!(P, Dy, serial_num, Vec_Trajectory, Vector_velocity = Vector_velocity)
        
            Dy.Time[serial_num] += P.Δt
            Dy.Displace[i+1, serial_num, :] = Dy.Trajectory[:, serial_num]
        end



    else
        return @error "the Time shold be a real number"
    end


end



end