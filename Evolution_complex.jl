module Ground_state

export CT_Evolution!                              #虚时演化

using ..TDQMC
using ..TDQMC.Potential_Matrix
using ..TDQMC.Crank_Nicolson
using ..TDQMC.Trajectory
using ..TDQMC.Quantity

using SparseArrays, NumericalIntegration


@inline function Normalization(P::Parameter, Dy::Dynamics, serial_num::Integer)
    return sqrt.(integrate(P.sampling, abs2.(hcat(Dy.Guide_Wave[:, serial_num])...), SimpsonEven()))
end


function CT_Evolution!(P::Parameter, Dy::Dynamics, serial_num::Integer;
    later_fix::SparseMatrixCSC, former_fix::SparseMatrixCSC)


    local Change_matrix_former::Vector{<:SparseMatrixCSC} = [spzeros(eltype(later_fix), P.space_N, P.space_N) for i = 1:P.electron]
    local Change_matrix_later::Vector{<:SparseMatrixCSC} = -deepcopy(Change_matrix_former)

    #    local Vec_Trajectory = view(Dy.Trajectory, :, serial_num)
    #    local Trajectory_past = zeros(eltype(Dy.Trajectory), P.electron)                                #用来记录波函数的历史轨迹,用来决定迭代停止的判断条件
    local Vec_wave = view(Dy.Guide_Wave, :, serial_num)
    local Vec_Trajectory = view(Dy.Trajectory, :, serial_num)

    local Normalizer::Vector{<:AbstractFloat} = zeros(eltype(Dy.Trajectory), P.electron)


    Dy.Displace[1, serial_num, :] = Vec_Trajectory


    if imag(P.Δt) != 0.0                        #这里确保是复时演化

        for i in 1:P.step_t
        
            Reset_matrix!(P, Dy, serial_num, Change_matrix_former, Change_matrix_later)
            Construct_matrix!(P, later_fix, former_fix, Change_matrix_former, Change_matrix_later)
        
            Vec_wave[:] = Change_matrix_former .* Vec_wave
            Vec_wave[:] = Change_matrix_later .\ Vec_wave
        
            Normalizer .= Normalization(P, Dy, serial_num)                                  #因为要重复使用,这里应该需要储存到一个变量里比较好.
            Vec_wave ./= Normalizer
        
            Movement!(P, Dy, serial_num, dt = P.Δt)
        
            Dy.Time[serial_num] += P.Δt
            Dy.Displace[i+1, serial_num, :] = Vec_Trajectory           #注意存储数据的时候不要数据竞争
        
        end

        #record(P, Wave)                         #这里结束后说明迭代得到了稳定的波函数
    else
        return @error "the Time shold be a Complex number"
    end
end

end