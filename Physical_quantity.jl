module Quantity

export Group_Energy, HHG

using ..TDQMC
using ..TDQMC: b
using ..TDQMC.Discrete_Func_diff
using ..TDQMC.Find_nearest_k
using ..TDQMC.Potential_Matrix: V_ne

using Interpolations, LinearAlgebra, SparseArrays, FFTW

#在另一个模块两个类似的函数是用来计算矢量的
V_ee(x_difference::T; β::T=b) where {T<:AbstractFloat} = 1.0 / sqrt(β + x_difference^2)


@inline function sptril(A::Matrix, n::Int)
    return sparse(tril(A, n))
end

@inline function relative(A::AbstractVector)
    return repeat(A', outer=(length(A), 1)) - repeat(A, outer=(1, length(A)))
end


function kinetic_part(P::Parameter, Dy::Dynamics, serial_num::Integer; k::Integer=5)
    local Type = eltype(eltype(Dy.Guide_Wave))
    local location_vec::Vector = Dy.Trajectory[:, serial_num]
    local indexVec_k::Vector{<:UnitRange{<:Integer}} = find_k_index.(location_vec, x=P.sampling, k=k)
    local k_itp_Wave_Vec::Vector{<:Vector{<:Complex}} = [zeros(Type, k) for i = 1:P.electron]    #选取的五个轨迹点最近邻点的波函数的值
    local derivative_Wave_Vec::Vector{<:Vector{<:Complex}} = deepcopy(k_itp_Wave_Vec)
    local Vector_Interp::Vector{<:Complex} = zeros(Type, P.electron)
    local Vector_Interp_derv::Vector{<:Complex} = zeros(Type, P.electron)

    for i = 1:P.electron
        k_itp_Wave_Vec[i] = Dy.Guide_Wave[i, serial_num][indexVec_k[i]]
    end

    derivative_Wave_Vec = Derivative_2.(k_itp_Wave_Vec, dL=P.Δx)

    for i = 1:P.electron
        interp_cubic_1 = CubicSplineInterpolation(P.sampling[indexVec_k[i]], k_itp_Wave_Vec[i])
        interp_cubic_2 = CubicSplineInterpolation(P.sampling[indexVec_k[i]], derivative_Wave_Vec[i])
        Vector_Interp[i] = interp_cubic_1(location_vec[i])
        Vector_Interp_derv[i] = interp_cubic_2(location_vec[i])
    end


    return (-1 / 2) * real(sum(Vector_Interp_derv ./ Vector_Interp))               #这里要取实数,是因为最后能量计算得到的值为实数,有虚部是因为波函数最后还是带有一定的虚部, 影响最后的计算结果
end


function potential_part(P::Parameter, Dy::Dynamics, serial_num::Integer)
    local location_vec::Vector = Dy.Trajectory[:, serial_num]
    local potential_vec::Vector = V_ne.(location_vec)       #核的势能项
    local location_Matrix::SparseMatrixCSC = sptril(relative(location_vec), 0)
    local interaction_energy::AbstractFloat = 0.0

    for i = 1:P.electron-1
        interaction_energy += sum(V_ee.(diag(location_Matrix, -i)))    #这一部分是电子与电子之间相互作用的能量
    end

    interaction_energy += sum(potential_vec)

    return interaction_energy
end


function Group_Energy(P::Parameter, Dy::Dynamics, serial_num::Integer)
    return kinetic_part(P, Dy, serial_num) + potential_part(P, Dy, serial_num)
end


function HHG(P::Parameter, Dy::Dynamics)
    local a, b, c = size(Dy.Displace)                               #第一维为对时间的采样, 第二维为系综数,需要进行求和取平均,第三维为电子数
    local Type_0 = eltype(Dy.Displace)
    local floor_a = floor(Int, a / 2)
    local fₛ = a / (P.step_t * real(P.Δt))
    local Discrete_acc = zeros(Type_0, a, b, c)
    local Discrete_ft_dipole = zeros(Complex{Type_0}, floor_a + 1, b, c)
    local Discrete_ft_acc = zeros(Complex{Type_0}, floor_a + 1, b, c)       #预置元素为复数的数组
    local Total_ft_dipole = zeros(Complex{Type_0}, floor_a + 1, c)
    local Total_ft_acc = zeros(Complex{Type_0}, floor_a + 1, c)

    for j in 1:c, i in 1:b
        Discrete_acc[:, i, j] = Derivative_2(Dy.Displace[:, i, j], dL=P.Δt)
    end

    Discrete_ft_dipole[:, :, :] = rfft(Dy.Displace, 1)
    Discrete_ft_acc[:, :, :] = rfft(Discrete_acc, 1)
    Total_ft_dipole[:, :] = sum(Discrete_ft_dipole, dims=2) / b                 #对每个粒子的轨迹进行DFT之后求和之后取平均
    Total_ft_acc[:, :] = sum(Discrete_ft_acc, dims=2) / b

    return rfftfreq(a, fₛ), Total_ft_dipole, Total_ft_acc
end

end