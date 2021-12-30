module Quantity

export Group_Energy

using ..TDQMC
using ..TDQMC.Discrete_Func_diff
using ..TDQMC.Find_nearest_k
using ..TDQMC.Potential_Matrix: V_ne

using Interpolations, LinearAlgebra

#在另一个模块两个类似的函数是用来计算矢量的
V_ee(x_difference::T; β::T = 1.0) where {T<:AbstractFloat} = 1.0 / sqrt(β + x_difference^2)


function kinetic_part(P::Parameter, Dy::Dynamics, serial_num::Integer)
    local Type = eltype(eltype(Dy.Guide_Wave))
    local location_vec::Vector = Dy.Trajectory[:, serial_num]
    local indexVec_k::Vector{<:UnitRange{<:Integer}} = find_k_index.(location_vec, x = P.sampling, k = 5)
    local k_itp_Wave_Vec::Vector{<:Vector{<:Complex}} = [zeros(Type, k) for i = 1:P.electron]    #选取的五个轨迹点最近邻点的波函数的值
    local derivative_Wave_Vec::Vector{<:Vector{<:Complex}} = deepcopy(k_itp_Wave_Vec)
    local Vector_Interp::Vector{<:Complex} = zeros(Type, P.electron)
    local Vector_Interp_derv::Vector{<:Complex} = zeros(Type, P.electron)

    for i = 1:P.electron
        k_itp_Wave_Vec[i] = Dy.Guide_Wave[i, serial_num][indexVec_k[i]]
    end

    derivative_Wave_Vec = Derivative_2.(k_itp_Wave_Vec, dL = P.Δx)

    for i = 1:P.electron
        interp_cubic_1 = CubicSplineInterpolation(P.sampling[indexVec_k[i]], k_itp_Wave_Vec[i])
        interp_cubic_2 = CubicSplineInterpolation(P.sampling[indexVec_k[i]], derivative_Wave_Vec[i])
        Vector_Interp[i] = interp_cubic_1(localtion_vec[i])
        Vector_Interp_derv[i] = interp_cubic_2(localtion_vec[i])
    end


    return sum(Vector_Interp_derv ./ Vector_Interp)

end




function potential_part(P::Parameter, Dy::Dynamics, serial_num::Integer)
    local location_vec::Vector = Dy.Trajectory[:, serial_num]
    local potential_vec::Vector = V_ne.(location_vec)       #核的势能项
    local location_Matrix::Matrix = repeat(location_vec', outer = (2, 1)) .- repeat(location_vec, outer = (1, 2))
    local interaction_energy::AbstractFloat = 0.0

    for i = 1:P.electron-1
        interaction_energy += sum(V_ee.(diag(location_Matrix, i)))    #这一部分是电子与电子之间相互作用的能量
    end

    interaction_energy += sum(potential_vec)

    return interaction_energy
end


function Group_Energy(P::Parameter, Dy::Dynamics, serial_num::Integer)
    return kinetic_part(P, Dy, serial_num) + potential_part(P, Dy, serial_num)
end


function TD_dipole!()


end

function acceleration()

end


end