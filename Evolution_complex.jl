module Ground_state

export CT_Evolution!                              #虚时演化

using ..TDQMC
using ..TDQMC.Potential_Matrix
using ..TDQMC.Crank_Nicolson
using ..TDQMC.Trajectory
using ..TDQMC.Quantity

using SparseArrays, NumericalIntegration, LinearAlgebra

@inline function find_inbound!(P::Parameter, Dy::Dynamics, serial_num::Integer, Vec_Trajectory::SubArray)
    Dy.Index[serial_num] = findall(x -> abs(x) < P.scope, Vec_Trajectory)
    Dy.In_num[serial_num] = length(Dy.Index[serial_num])
end


@inline function Normalizer(P::Parameter, Dy::Dynamics, serial_num::Integer)
    return sqrt.(integrate(P.sampling, abs2.(hcat(Dy.Guide_Wave[Dy.Index[serial_num], serial_num]...)), SimpsonEven()))
end


@inline function Normalization!(P::Parameter, Dy::Dynamics, serial_num::Integer, Vec_wave::SubArray{<:Vector})
    Vec_wave[Dy.Index[serial_num]] ./= Normalizer(P, Dy, serial_num)
end


@inline function Spin_Effect(P::Parameter, Dy::Dynamics, serial_num::Integer)
    local Select_Spin::Vector = P.Spin[Dy.Index[serial_num]]
    local Spin_Effect::Matrix{<:Integer} = Select_Spin * Select_Spin'

    replace!(tril!(Spin_Effect), 0 => 1)

    return Spin_Effect
end


function CT_Evolution!(P::Parameter, Dy::Dynamics, serial_num::Integer;
    later_fix::SparseMatrixCSC, former_fix::SparseMatrixCSC)

    local Change_matrix_former::Vector{<:SparseMatrixCSC} = [spzeros(eltype(later_fix), P.space_N, P.space_N) for i = 1:P.electron]
    local Change_matrix_later::Vector{<:SparseMatrixCSC} = -deepcopy(Change_matrix_former)

    local Velocity_WaveFunc::Vector{Matrix{<:Complex}} = [zeros(eltype(later_fix), P.electron, P.electron) for i = 1:P.electron]
    local Vector_velocity::Vector{<:AbstractFloat} = zeros(eltype(Dy.Trajectory), P.electron)

    local Vec_wave = view(Dy.Guide_Wave, :, serial_num)
    local Vec_Trajectory = view(Dy.Trajectory, :, serial_num)

    local Spin_Matrix::Matrix{<:Integer} = zeros(Int, P.electron, P.electron)

    local record_Change::Integer = 0

    Dy.Displace[1, serial_num, :] = Vec_Trajectory

    find_inbound!(P, Dy, serial_num, Vec_Trajectory)
    Spin_Matrix = Spin_Effect(P, Dy, serial_num)

    if imag(P.Δt) != 0.0                        #这里确保是复时演化

        for i in 1:P.step_t
        
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
        
            Normalization!(P, Dy, serial_num, Vec_wave)                                  #因为要重复使用,这里应该需要储存到一个变量里比较好.
        
            Velocity!(P, Dy, serial_num, Velocity_WaveFunc, Vector_velocity; Spin_Matrix = Spin_Matrix)  #重置速度向量,这里只对Index索引的更新速度
            Movement!(P, Dy, serial_num, Vec_Trajectory, Vector_velocity = Vector_velocity)
        
            Dy.Time[serial_num] += P.Δt
            Dy.Displace[i+1, serial_num, :] = Vec_Trajectory           #注意存储数据的时候不要数据竞争
        
        end

    else
        return @error "the Time shold be a Complex number"
    end
end

end