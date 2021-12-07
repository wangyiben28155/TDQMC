module Trajectory

export Movement!

using ..TDQMC
using ..TDQMC.Find_nearest_k
using Interpolations, LinearAlgebra   #然后对不在格点上的轨迹进行插值求得其导数等等


function find_lattice_wave(Particle_num::Integer, serial_num::Integer, P::Parameter, Dy::Dynamics; k::Integer = 5) where {T<:AbstractFloat}     #这个对应的是粒子的轨迹矢量,用来得到相应的波函数插值的部分
    local localtion = Dy.Trajectory[Particle_num, serial_num]
    local indexVec_k::Vector{<:Integer} = find_k_index(localtion, x = P.sampling, k = k)
    local Type_0 = eltype(eltype(Dy.Guide_Wave))
    local k_itp_Wave_Vector::Vector{<:Vector{<:Complex{T}}} = [zeros(Type_0, k) for i in 1:P.electron]     #=选取五个点用来插值的波函数矢量,注意这种
                                                                                                           赋值的方式每一个重复的部分在后面赋值的时候不会产生关联.=#

    for i = 1:P.electron
        @inbounds k_itp_Wave_Vector[i] = Dy.Guide_Wave[i, serial_num][indexVec_k]
    end

    return P.sampling[indexVec_k], k_itp_Wave_Vector
end


function Interpolation_Wave(Particle_num::Integer, serial_num::Integer, P::Parameter, Dy::Dynamics)
    local localtion = Dy.Trajectory[Particle_num, serial_num]
    local xd, yd = find_lattice_wave(Particle_num, serial_num, P, Dy)
    local Type_1 = eltype(eltype(yd))
    local Vec_Wave::Vector{<:Complex} = zeros(Type_1, P.electron)
    local Vec_Derivate::Vector{<:Complex} = zeros(Type_1, P.electron)
    local interp_cubic::Function

    for i = 1:P.electron
        interp_cubic = CubicSplineInterpolation(xd, yd[i])                         #得到三次样条插值函数
        Vec_Wave[i] = interp_cubic(localtion)            #然后
        Vec_Derivate[i] = Interpolations.gradient(interp_cubic, localtion)
    end
    return Vec_Wave, Vec_Derivate
end



function Slater_determinant(P::Parameter, Dy::Dynamics, serial_num::Integer)           #通过此函数得到交叉关联的波函数, 按理来说有多少个电子就应该有多少个坐标,对应Electron_num维度的电子波函数
    local Type_2 = eltype(eltype(yd))
    local Vec_Wave::Vector{<:Complex}, Vec_Derivate::Vector{<:Complex} = (zeros(Type_2, P.electron), zeros(Type_1, P.electron))
    local symmetric_determinate::Matrix{<:Complex} = zeros(Type_2, (P.electron, P.electron))
    local Derivate_eachcoodinate::Matrix{<:Complex} = zeros(Type_2, (P.electron, P.electron))

    for i = 1:P.electron
        Vec_Wave, Vec_Derivate = Interpolation_Wave(i, serial_num, P, Dy)
        symmetric_determinate[:, i] = Vec_Wave
        Derivate_eachcoodinate[:, i] = Vec_Derivate
    end


    return symmetric_determinate, Derivate_eachcoodinate
end




function Accelaration(P::Parameter, Dy::Dynamics, serial_num::Integer)
    local symmetric_determinate, Derivate_eachcoodinate = Slater_determinant(P, Dy, serial_num)
    local Derivate_WaveFunc::Vector{<:Matrix{<:Complex}} = [symmetric_WaveFunc for i = 1:P.electron]
    local symmetric_WaveFunc::Vector{<:Matrix{<:Complex}} = fill(symmetric_WaveFunc, P.electron)
    local Vector_accelerate::Vector{<:AbstractFloat} = zeros(eltype(symmetric_determinate), P.electron)

    for i = 1:P.electron
        @inbounds Derivate_WaveFunc[i][:, i] = Derivate_eachcoodinate[:, i]
    end

    Vector_accelerate[:] = @. imag(det(Derivate_WaveFunc) / det(symmetric_WaveFunc))

    return Vector_accelerate
end




function Movement!(P::Parameter, Dy::Dynamics, serial_num::Integer; dt = P.Δt)                     # 这里我们使用欧拉法即可
    Dy.Trajectory[:, serial_num] .+= abs(dt) * Accelaration(P, Dy, serial_num)
end



end