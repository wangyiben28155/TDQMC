module Ground_state

export CTR!                              #虚时演化

using ..TDQMC
using ..TDQMC.Potential_Matrix
using ..TDQMC.Crank_Nicolson
using ..TDQMC.Trajectory
using ..TDQMC.Quantity

using SparseArrays, NumericalIntegration



@inline function Normalization(P::Parameter, Dy::Dynamics, serial_num::Integer)
    return sqrt.(integrate(P.sampling, abs2.(hcat(Dy.Guide_Wave[:, serial_num])...), SimpsonEven()))
end


function CTR!(P::Parameter, Dy::Dynamics, serial_num::Integer;
    later_fix::SparseMatrixCSC, former_fix::SparseMatrixCSC)

    local count::Int64 = 0

    local Change_matrix_former::Vector{<:SparseMatrixCSC} = [spzeros(eltype(later_fix), P.space_N, P.space_N) for i = 1:P.electron]     #这两个向量稀 
    #疏矩阵用来计算的两个矩阵
    local Change_matrix_later::Vector{<:SparseMatrixCSC} = -deepcopy(Change_matrix_former)           #所有的元素取负值
    local Vec_Trajectory = view(Dy.Trajectory, :, serial_num)
    local Trajectory_past = zeros(eltype(Dy.Trajectory), P.electron)                                #用来记录波函数的历史轨迹,用来决定迭代停止的判断条件
    local Vec_wave = view(Dy.Guide_Wave, :, serial_num)

    local Normalizer::Vector{<:AbstractFloat} = zeros(eltype(Dy.Trajectory), P.electron)


    if imag(P.Δt) != 0.0                        #这里确保是虚时演化

        while true

            if sum(@. abs(Trajectory_past - Vec_Trajectory)) >= 1e-15                     #用轨迹判断终止条件, 到最后应该有粒子的轨迹前后变化很小
                Trajectory_past[:] = Vec_Trajectory
                #正态分布函数本身就是归一化的,所以第一次运行的时候,wave_past是归一化的
            
                Reset_matrix!(P, Dy, serial_num, Change_matrix_former, Change_matrix_later)
                Construct_matrix!(P, later_fix, former_fix, Change_matrix_former, Change_matrix_later)
            
                Vec_wave[:] = Change_matrix_former .* Vec_wave
                Vec_wave[:] = Change_matrix_later .\ Vec_wave
            
                Normalizer .= Normalization(P, Dy, serial_num)                                  #因为要重复使用,这里应该需要储存到一个变量里比较好.
                Vec_wave ./= Normalizer
            
                Movement!(P, Dy, serial_num, dt = P.Δt)
            
                count += 1
            else
            
                Dy.Time[serial_num] = count * P.Δt
                break
            end

        end

        #record(P, Wave)                         #这里结束后说明迭代得到了稳定的波函数


    else
        return @error "the Time shold be a imaginary number"
    end
end

end