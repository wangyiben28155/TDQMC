module Trajectory

export Movement

using ..TDQMC
using ..TDQMC.Discrete_Func_diff: Five_Point!               #使用此文件中定义的数值微分的函数
using ..TDQMC.Find_nearest_k
using Interpolations   #然后对不在格点上的轨迹进行插值求得其导数等等

function find_lattice(Wave_num::Integer, P::Parameter, Dy::Dynamics, serial_num::Integer; k::Integer = 5) where {T<:AbstractFloat}     #这个对应的是粒子的轨迹矢量,用来得到相应的波函数插值的部分
    local indexVec_k::Vector{<:Vector{<:Integer}} = find_k_index(Dy.Trajectory[:, serial_num], x = P.sampling, k = k)
    local Five_IP_Wave_Vector::Vector{<:Vector{<:Complex{T}}} = fill(zeros(Complex{T}, k), P.electron)     #选取五个点用来插值的波函数矢量

    for i = 1:P.electron
        @inbounds Five_IP_Wave_Vector[i] = Dy.Guide_Wave[Wave_num, serial_num][indexVec_k[i]]              
    end

    return Five_IP_Wave_Vector
end

function slater_determinant(P::Parameter, Dy::Dynamics, serial_num::Integer)           #通过此函数得到交叉关联的波函数, 按理来说有多少个电子就应该有多少个坐标,对应Electron_num维度的电子波函数
    local symmetric_Wave



    return symmertric_Wave
end

function derivative_Ontrajectory(P::Parameter, Dy::Dynamics, serial_num::Integer)
    local Wave_lattice = find_lattice(P, Dy, serial_num)            #得到的是波函数的取样点,用来插值的位置

end




function Movement(P::Parameter, Dy::Dynamics)                     # 这里我们使用欧拉法即可
    local 







end